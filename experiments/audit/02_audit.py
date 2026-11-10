"""Stage 1 audit: is per-residue anisotropy real, reproducible, and predictable?

Reads the dataset built by ``00_build_dataset.py``. Statistics are computed
per protein and aggregated across proteins with a bootstrap over proteins,
because residues within a protein are not independent.

Sections
--------
R1  Reproducibility (MD replica i -> replica j). The kill-criterion test.
    Win rate of Sigma_i over its isotropic surrogate, against an isotropic
    truth passed through the identical superposition (the artifact
    ceiling), and robustness to three superposition choices.
R2  How much of the anisotropy is real, and what it is worth.
    Optimal shape power alpha (0 = scalar optimal, 1 = take Sigma at face
    value), principal-axis agreement, and held-out Gaussian log-likelihood
    gain in nats per residue per frame.
R3  Predictability from structure alone. ANM (one parameter, no learning)
    against MD, placed between the isotropic floor and the replica ceiling.
    This is the bar a learned covariance model has to clear.
R4  External check against experiment: deposited crystallographic ANISOU.
R5  Where the signal lives: secondary structure, flexibility, termini.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from auf.audit import airm_distance, fractional_anisotropy, isotropize
from auf.audit.direction import (
    direction_wins,
    gaussian_loglik,
    optimal_shape_power,
    principal_axis_cos,
    shape_power,
    unit_trace,
)

ALPHAS = np.linspace(0.0, 2.0, 41)
_IU = (np.array([0, 1, 2, 0, 0, 1]), np.array([0, 1, 2, 1, 2, 2]))


def _floor_psd(u, eps=1e-4):
    w, v = np.linalg.eigh(u)
    return (v * np.maximum(w, eps)[..., None, :]) @ np.swapaxes(v, -1, -2)


def tls_plus(xyz, u, extra):
    """RSS of U ~ TLS(20) + c * extra (6 unique comps per residue)."""
    from auf.data.crystal import tls_design

    design = tls_design(xyz)
    if extra is not None:
        design = np.column_stack([design, extra[:, _IU[0], _IU[1]].reshape(-1)])
    target = u[:, _IU[0], _IU[1]].reshape(-1)
    coef, *_ = np.linalg.lstsq(design, target, rcond=None)
    return float(np.sum((design @ coef - target) ** 2))


# --------------------------------------------------------------------------
# Per-protein measurements
# --------------------------------------------------------------------------

def pairs(n):
    return list(itertools.permutations(range(n), 2))


def replica_wins(sig_rep):
    """(n_pairs, n_res) wins of Sigma_i over iso(Sigma_i) at predicting Sigma_j."""
    return np.stack([direction_wins(sig_rep[i], sig_rep[j]) for i, j in pairs(len(sig_rep))])


def artifact_control(z, rng, n_frames):
    """R1 on isotropic truth at each residue's real RMSF, same superposition."""
    from auf.data.atlas import covariances_from_trajectories

    ref = z["reference"]
    scale = np.sqrt(np.trace(z["sigma"], axis1=-2, axis2=-1) / 3.0)
    frames = ref + scale[:, None] * rng.normal(size=(3, n_frames, len(ref), 3))
    return replica_wins(covariances_from_trajectories("ctl", frames).sigma_rep).mean(0)


def loglik_gain(sig_rep, alpha=1.0):
    """Held-out log-lik gain of shaped Sigma_i over iso(Sigma_i) on replica j,
    nats per residue per frame, (n_res,) averaged over pairs. Trace kept."""
    out = []
    for i, j in pairs(len(sig_rep)):
        si, sj = sig_rep[i], sig_rep[j]
        tr = np.trace(si, axis1=-2, axis2=-1)[:, None, None]
        shaped = shape_power(si, alpha) * tr
        out.append(gaussian_loglik(shaped, sj) - gaussian_loglik(isotropize(si), sj))
    return np.mean(out, axis=0)


def pooled_alpha(preds, targets):
    """Best shape power when each pred_k must explain target_k (stacked)."""
    return optimal_shape_power(np.concatenate(preds), np.concatenate(targets), ALPHAS)


def measure(z, rng, control_frames):
    sig_rep = z["sigma_rep_mean"]
    n = sig_rep.shape[1]
    r = {"pdb_chain": str(z["pdb_chain"]), "n_res": n}
    rmsf = np.sqrt(np.trace(z["sigma"], axis1=-2, axis2=-1))

    # R1
    wins = replica_wins(sig_rep)
    ctl = artifact_control(z, rng, control_frames)
    r["R1_win"] = wins.mean()
    r["R1_control"] = ctl.mean()
    for mode in ("start", "core"):
        r[f"R1_win_{mode}"] = replica_wins(z[f"sigma_rep_{mode}"]).mean()

    # R2
    pp = pairs(len(sig_rep))
    alpha, _ = pooled_alpha([sig_rep[i] for i, _ in pp], [sig_rep[j] for _, j in pp])
    r["R2_alpha_rep"] = alpha
    r["R2_axis_cos_rep"] = np.mean([principal_axis_cos(sig_rep[i], sig_rep[j]).mean() for i, j in pp])
    ll = loglik_gain(sig_rep)
    r["R2_loglik_gain"] = ll.mean()
    r["R2_loglik_gain_shrunk"] = loglik_gain(sig_rep, min(alpha, 1.0)).mean()
    r["R2_fa"] = np.median(fractional_anisotropy(z["sigma"]))

    # R3: ANM vs MD. Compare each predictor's shape against replica j.
    anm = z["anm"]
    a_anm, _ = optimal_shape_power(anm, z["sigma"], ALPHAS)
    r["R3_anm_win"] = direction_wins(anm, z["sigma"]).mean()
    r["R3_anm_alpha"] = a_anm
    r["R3_anm_axis_cos"] = principal_axis_cos(anm, z["sigma"]).mean()
    d_iso, d_anm, d_anm_t, d_rep = [], [], [], []
    iso = np.broadcast_to(np.eye(3) / 3, anm.shape)
    for i, j in pp:
        t = unit_trace(sig_rep[j])
        d_iso.append(airm_distance(iso, t))
        d_anm.append(airm_distance(unit_trace(anm), t))
        d_anm_t.append(airm_distance(shape_power(anm, a_anm), t))
        d_rep.append(airm_distance(unit_trace(sig_rep[i]), t))
    r["R3_dist_iso"] = np.median(d_iso)
    r["R3_dist_anm"] = np.median(d_anm)
    r["R3_dist_anm_tuned"] = np.median(d_anm_t)
    r["R3_dist_replica"] = np.median(d_rep)

    # R4: crystal
    u = z["crystal_u"]
    ok = np.isfinite(u[:, 0, 0])
    if ok.sum() >= 10:
        ok[ok] &= np.linalg.eigvalsh(u[ok]).min(-1) > 0
    if ok.sum() >= 10:
        U, S, A = u[ok], z["sigma"][ok], anm[ok]
        r["R4_n"] = int(ok.sum())
        r["R4_rmsd"] = float(z["crystal_rmsd"])
        r["R4_fa_crystal"] = np.median(fractional_anisotropy(U))
        r["R4_fa_md"] = np.median(fractional_anisotropy(S))
        for name, P in (("md", S), ("anm", A)):
            r[f"R4_{name}_win"] = direction_wins(P, U).mean()
            r[f"R4_{name}_alpha"] = optimal_shape_power(P, U, ALPHAS)[0]
            r[f"R4_{name}_axis_cos"] = principal_axis_cos(P, U).mean()
        # Rigid body: how much of the crystal's anisotropy is one TLS group?
        from auf.data.crystal import fit_tls

        X = z["reference"][ok]
        T = _floor_psd(fit_tls(X, U))
        r["R4_tls_win"] = direction_wins(T, U).mean()
        r["R4_tls_axis_cos"] = principal_axis_cos(T, U).mean()
        # After TLS, does internal-motion direction help beyond a scalar of
        # the same per-residue size? Same parameter count on both sides.
        for name, P in (("md", S), ("anm", A)):
            rss_dir = tls_plus(X, U, P)
            rss_iso = tls_plus(X, U, isotropize(P))
            r[f"R4_tls_plus_{name}_better"] = float(rss_dir < rss_iso)
            r[f"R4_tls_plus_{name}_rss_ratio"] = rss_dir / rss_iso

    # R5 strata (per-residue arrays, aggregated later)
    ss = z["ss"].astype(str)
    term = np.zeros(n, bool)
    term[:5] = term[-5:] = True
    r["_res"] = {
        "win": wins.mean(0), "control": ctl, "loglik": ll, "rmsf": rmsf,
        "ss": ss, "terminal": term,
        "anm_win": direction_wins(anm, z["sigma"]).astype(float),
        "axis_cos": np.mean([principal_axis_cos(sig_rep[i], sig_rep[j]) for i, j in pp], 0),
    }
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in r.items()}


# --------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------

def boot_ci(x, rng, n_boot=4000, stat=np.mean):
    x = np.asarray([v for v in x if np.isfinite(v)])
    if len(x) == 0:
        return (np.nan, np.nan, np.nan, 0)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    b = stat(x[idx], axis=1)
    return (float(stat(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)), len(x))


def fmt(ci, pct=False, digits=2):
    m, lo, hi, n = ci
    if pct:
        return f"{100 * m:.0f}% [{100 * lo:.0f}, {100 * hi:.0f}] (n={n})"
    return f"{m:.{digits}f} [{lo:.{digits}f}, {hi:.{digits}f}] (n={n})"


def col(rows, key):
    return [r.get(key, np.nan) for r in rows]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--data", default="data/audit")
    p.add_argument("--crystal-extra", default="data/audit_highres",
                   help="extra entries used for R4 only (high resolution)")
    p.add_argument("--control-frames", type=int, default=3000)
    p.add_argument("--out", default="results/02_audit")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)

    main_files = sorted(Path(args.data).glob("*.npz"))
    rows = []
    for f in main_files:
        with np.load(f) as z:
            rows.append(measure(dict(z), rng, args.control_frames))
    seen = {r["pdb_chain"] for r in rows}
    extra = []
    if Path(args.crystal_extra).exists():
        for f in sorted(Path(args.crystal_extra).glob("*.npz")):
            if f.stem in seen:
                continue
            with np.load(f) as z:
                extra.append(measure(dict(z), rng, args.control_frames))
    cryst = [r for r in rows + extra if "R4_n" in r]

    L = []
    say = L.append
    say(f"# Stage 1 audit: {len(rows)} ATLAS proteins "
        f"({sum(r['n_res'] for r in rows)} residues), random draw")
    say("Mean over proteins, [95% bootstrap CI over proteins].\n")

    say("## R1  Reproducibility of direction (kill criterion)")
    say(f"Sigma_i beats iso(Sigma_i) at predicting Sigma_j:  {fmt(boot_ci(col(rows, 'R1_win'), rng), True)}")
    say(f"  same, isotropic truth through same superposition: {fmt(boot_ci(col(rows, 'R1_control'), rng), True)}")
    margin = np.array(col(rows, "R1_win")) - np.array(col(rows, "R1_control"))
    say(f"  margin over control:                              {fmt(boot_ci(margin, rng), True)}")
    say(f"  proteins clearing control by >10 pts:             "
        f"{(margin > 0.10).sum()}/{len(rows)}")
    say(f"  fit on starting structure:                        {fmt(boot_ci(col(rows, 'R1_win_start'), rng), True)}")
    say(f"  fit on rigid half only:                           {fmt(boot_ci(col(rows, 'R1_win_core'), rng), True)}")

    say("\n## R2  How much anisotropy is real, and what it is worth")
    say(f"fractional anisotropy (median per protein):        {fmt(boot_ci(col(rows, 'R2_fa'), rng))}")
    say(f"optimal shape power alpha, replica -> replica:     {fmt(boot_ci(col(rows, 'R2_alpha_rep'), rng))}")
    say(f"principal-axis |cos|, replica -> replica (rand .5): {fmt(boot_ci(col(rows, 'R2_axis_cos_rep'), rng))}")
    say(f"held-out log-lik gain over scalar, nats/res/frame: {fmt(boot_ci(col(rows, 'R2_loglik_gain'), rng), digits=3)}")
    say(f"  same with alpha-shrunk Sigma:                     {fmt(boot_ci(col(rows, 'R2_loglik_gain_shrunk'), rng), digits=3)}")

    say("\n## R3  Predictable from structure alone? (ANM, no learning)")
    say(f"ANM direction beats scalar vs MD:                  {fmt(boot_ci(col(rows, 'R3_anm_win'), rng), True)}")
    say(f"ANM optimal shape power:                           {fmt(boot_ci(col(rows, 'R3_anm_alpha'), rng))}")
    say(f"ANM principal-axis |cos| vs MD:                    {fmt(boot_ci(col(rows, 'R3_anm_axis_cos'), rng))}")
    say("median unit-trace AIRM to replica j (lower is better):")
    for k, name in (("R3_dist_iso", "scalar (floor)"), ("R3_dist_anm", "ANM raw"),
                    ("R3_dist_anm_tuned", "ANM, alpha tuned"), ("R3_dist_replica", "replica i (ceiling)")):
        say(f"  {name:<22} {fmt(boot_ci(col(rows, k), rng))}")
    gap = (np.array(col(rows, "R3_dist_iso")) - np.array(col(rows, "R3_dist_anm_tuned"))) / (
        np.array(col(rows, "R3_dist_iso")) - np.array(col(rows, "R3_dist_replica")))
    say(f"fraction of floor->ceiling gap closed by ANM:      {fmt(boot_ci(gap, rng, stat=np.median))} (median)")

    say(f"\n## R4  Against experiment: crystallographic ANISOU ({len(cryst)} proteins)")
    if cryst:
        say(f"FA crystal vs FA MD:                               {fmt(boot_ci(col(cryst, 'R4_fa_crystal'), rng))} vs "
            f"{fmt(boot_ci(col(cryst, 'R4_fa_md'), rng))}")
        for name in ("md", "anm"):
            say(f"{name.upper():<4} beats scalar vs crystal:                    {fmt(boot_ci(col(cryst, f'R4_{name}_win'), rng), True)}")
            say(f"     optimal shape power:                         {fmt(boot_ci(col(cryst, f'R4_{name}_alpha'), rng))}")
            say(f"     principal-axis |cos| (random 0.5):           {fmt(boot_ci(col(cryst, f'R4_{name}_axis_cos'), rng))}")
        say(f"TLS (1 rigid group, 20 params, in-sample) beats scalar: {fmt(boot_ci(col(cryst, 'R4_tls_win'), rng), True)}")
        say(f"     principal-axis |cos|:                        {fmt(boot_ci(col(cryst, 'R4_tls_axis_cos'), rng))}")
        for name in ("md", "anm"):
            better = col(cryst, f"R4_tls_plus_{name}_better")
            say(f"TLS + {name.upper()} shape fits better than TLS + scalar: "
                f"{int(np.nansum(better))}/{len(better)} proteins; median RSS ratio "
                f"{fmt(boot_ci(col(cryst, f'R4_tls_plus_{name}_rss_ratio'), rng, stat=np.median))}")

    say("\n## R5  Where the signal lives (residue level, pooled; descriptive)")
    res = {k: np.concatenate([r["_res"][k] for r in rows]) for k in rows[0]["_res"]}
    q = np.quantile(res["rmsf"], [1 / 3, 2 / 3])
    strata = {
        "helix": res["ss"] == "H", "strand": res["ss"] == "E", "coil": res["ss"] == "C",
        "termini (5)": res["terminal"], "interior": ~res["terminal"],
        "RMSF low third": res["rmsf"] <= q[0],
        "RMSF mid third": (res["rmsf"] > q[0]) & (res["rmsf"] <= q[1]),
        "RMSF high third": res["rmsf"] > q[1],
    }
    say(f"{'stratum':<17}{'n':>6}{'rep win':>9}{'control':>9}{'axis cos':>10}{'loglik':>9}{'ANM win':>9}")
    for name, m in strata.items():
        say(f"{name:<17}{m.sum():>6}{100 * res['win'][m].mean():>8.0f}%{100 * res['control'][m].mean():>8.0f}%"
            f"{res['axis_cos'][m].mean():>10.2f}{res['loglik'][m].mean():>9.3f}{100 * res['anm_win'][m].mean():>8.0f}%")

    report = "\n".join(L)
    print(report)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.md").write_text(report + "\n")
    flat = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows + extra]
    (out / "per_protein.json").write_text(json.dumps(flat, indent=1))
    print(f"\nwrote {out}/report.md, per_protein.json")


if __name__ == "__main__":
    main()
