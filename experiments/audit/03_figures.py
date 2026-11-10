"""Figures for the Stage 1 audit, from results/02_audit/per_protein.json."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
})


def strip(ax, cols, labels, colors, rng, ref=None, ref_label=None):
    for k, (vals, lab, c) in enumerate(zip(cols, labels, colors)):
        v = np.asarray([x for x in vals if np.isfinite(x)])
        jitter = rng.uniform(-0.18, 0.18, len(v))
        ax.scatter(k + jitter, v, s=14, color=c, alpha=0.55, linewidths=0)
        med = np.median(v)
        ax.plot([k - 0.3, k + 0.3], [med, med], color=INK, lw=2, solid_capstyle="round")
        ax.annotate(f"{med:.2f}", (k + 0.33, med), va="center", fontsize=8, color=INK)
    if ref is not None:
        ax.axhline(ref, color=INK2, lw=1, ls=(0, (3, 3)))
        ax.annotate(ref_label, (len(cols) - 0.5, ref), va="bottom", ha="right",
                    fontsize=8, color=INK2)
    ax.set_xticks(range(len(cols)), labels)
    ax.set_xlim(-0.6, len(cols) - 0.2)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--src", default="results/02_audit/per_protein.json")
    p.add_argument("--out", default="notes/figures")
    args = p.parse_args()
    rows = json.loads(Path(args.src).read_text())
    main_rows = [r for r in rows if "R1_win" in r]
    cryst = [r for r in rows if "R4_n" in r]
    rng = np.random.default_rng(0)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    def g(rs, k):
        return [r.get(k, np.nan) for r in rs]

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), dpi=200)

    ax = axes[0]
    strip(ax, [g(main_rows, "R1_control"), g(main_rows, "R1_win")],
          ["isotropic truth\n(artifact control)", "MD replicas"],
          [INK2, BLUE], rng)
    ax.set_ylim(0, 1)
    ax.set_title("A. Direction reproduces across MD replicas\n"
                 "share of residues where Σᵢ beats its scalar at predicting Σⱼ",
                 loc="left", fontsize=9)

    ax = axes[1]
    strip(ax, [g(main_rows, "R3_dist_iso"), g(main_rows, "R3_dist_anm_tuned"),
               g(main_rows, "R3_dist_replica")],
          ["scalar\n(floor)", "ANM\n(no learning)", "other replica\n(ceiling)"],
          [INK2, ORANGE, BLUE], rng)
    ax.set_title("B. Headroom for a learned model\n"
                 "median shape distance to MD replica j (lower is better)",
                 loc="left", fontsize=9)

    ax = axes[2]
    strip(ax, [g(cryst, "R4_md_axis_cos"), g(cryst, "R4_anm_axis_cos"),
               g(cryst, "R4_tls_axis_cos")],
          ["MD", "ANM", "rigid body\n(TLS, 20 params)"],
          [BLUE, ORANGE, AQUA], rng, ref=0.5, ref_label="random axes")
    ax.set_ylim(0.3, 1)
    ax.set_title(f"C. Against crystal ANISOU ({len(cryst)} proteins)\n"
                 "principal-axis agreement |cos|",
                 loc="left", fontsize=9)

    fig.tight_layout()
    fig.savefig(out / "stage1_audit.png", facecolor=SURFACE)
    print(f"wrote {out / 'stage1_audit.png'}")


if __name__ == "__main__":
    main()
