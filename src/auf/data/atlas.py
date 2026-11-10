"""ATLAS molecular-dynamics loader: per-residue Cα covariances.

ATLAS (https://www.dsimb.inserm.fr/ATLAS) provides three 100 ns replicas per
protein. The ``analysis`` archive carries each replica at 1000 frames
(100 ps spacing; 1001 including t=0), protein only, alongside ATLAS's own
per-replica RMSF, which is used here as a cross-check.

The pipeline is:

1. Read Cα coordinates for every replica, shape (n_frames, n_res, 3), in Å.
2. Superpose *all replicas jointly* onto one common reference (iterative
   Procrustes to the mean, by default). A shared frame is what makes
   replica-to-replica covariance comparisons meaningful: two covariances
   from separately fitted replicas differ by a rotation as well as by noise.
3. Per-residue 3x3 covariance about the per-replica mean, plus a pooled
   covariance over all frames, plus an effective sample size per residue
   from the integrated autocorrelation time.

Data licence: ATLAS is CC BY-NC 4.0. Nothing derived from it is committed.
"""
from __future__ import annotations

import io
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "ATLAS_API",
    "AtlasCovariances",
    "fetch_atlas",
    "read_ca_trajectories",
    "read_atlas_rmsf",
    "kabsch_rotation",
    "superpose",
    "per_residue_covariance",
    "integrated_autocorr_time",
    "effective_sample_size",
    "load_atlas_covariances",
]

ATLAS_API = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS"
DEFAULT_DATA_DIR = Path("data/atlas")


# --------------------------------------------------------------------------
# Download and read
# --------------------------------------------------------------------------

def fetch_atlas(
    pdb_chain: str,
    data_dir: str | Path = DEFAULT_DATA_DIR,
    kind: str = "analysis",
    timeout: float = 600.0,
) -> Path:
    """Download and unpack one ATLAS entry; return its directory.

    ``kind="analysis"`` is the 1000-frame archive (a few MB for a small
    protein). ``kind="protein"`` is the full-resolution trajectory archive
    and is typically hundreds of MB to GB. Already-unpacked entries are
    not downloaded again.
    """
    if kind not in {"analysis", "protein"}:
        raise ValueError(f"unknown ATLAS archive kind {kind!r}")
    out = Path(data_dir) / pdb_chain
    if (out / f"{pdb_chain}.pdb").exists() and any(out.glob(f"{pdb_chain}_R*.xtc")):
        return out
    out.mkdir(parents=True, exist_ok=True)
    url = f"{ATLAS_API}/{kind}/{pdb_chain}"
    tmp = out / f".{pdb_chain}_{kind}.zip.part"
    with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    with zipfile.ZipFile(tmp) as zf:
        for member in zf.namelist():
            # Archives are flat; refuse anything that would escape ``out``.
            name = Path(member).name
            if not name or name != member:
                continue
            with zf.open(member) as src, open(out / name, "wb") as dst:
                shutil.copyfileobj(src, dst)
    tmp.unlink()
    return out


def _replica_paths(entry_dir: Path, pdb_chain: str) -> list[Path]:
    paths = sorted(entry_dir.glob(f"{pdb_chain}_R*.xtc"))
    if not paths:
        raise FileNotFoundError(f"no {pdb_chain}_R*.xtc trajectories in {entry_dir}")
    return paths


def read_ca_trajectories(entry_dir: str | Path, pdb_chain: str) -> np.ndarray:
    """Cα coordinates for every replica, shape (n_rep, n_frames, n_res, 3), Å.

    Requires ``mdtraj`` (``pip install -e ".[atlas]"``) to decode XTC.
    """
    try:
        import mdtraj as md
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "reading ATLAS trajectories needs mdtraj: pip install -e '.[atlas]'"
        ) from exc

    entry_dir = Path(entry_dir)
    top = md.load_topology(str(entry_dir / f"{pdb_chain}.pdb"))
    ca = top.select("name CA")
    if ca.size == 0:
        raise ValueError(f"no CA atoms in {pdb_chain}.pdb")
    reps = []
    for path in _replica_paths(entry_dir, pdb_chain):
        traj = md.load(str(path), top=top, atom_indices=ca)
        reps.append(traj.xyz.astype(np.float64) * 10.0)  # nm -> Å
    n_frames = min(r.shape[0] for r in reps)
    return np.stack([r[:n_frames] for r in reps])


def read_atlas_rmsf(entry_dir: str | Path, pdb_chain: str) -> np.ndarray | None:
    """ATLAS's own per-replica RMSF, shape (n_rep, n_res), Å; None if absent."""
    path = Path(entry_dir) / f"{pdb_chain}_RMSF.tsv"
    if not path.exists():
        return None
    data = np.genfromtxt(path, delimiter="\t", names=True, dtype=None, encoding="utf-8")
    cols = [c for c in data.dtype.names if c.startswith("RMSF_")]
    return np.stack([np.asarray(data[c], dtype=float) for c in cols])


# --------------------------------------------------------------------------
# Superposition
# --------------------------------------------------------------------------

def kabsch_rotation(mobile: np.ndarray, target: np.ndarray,
                    weights: np.ndarray | None = None) -> np.ndarray:
    """Proper rotation(s) R minimising ``|| mobile @ R - target ||``.

    Both inputs are centred (..., n, 3). Returns (..., 3, 3) acting on row
    vectors, so the superposed coordinates are ``mobile @ R``.
    """
    if weights is not None:
        mobile = mobile * weights[:, None]
    h = np.swapaxes(mobile, -1, -2) @ target
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt))
    fix = np.ones(u.shape[:-2] + (3,))
    fix[..., -1] = d
    return (u * fix[..., None, :]) @ vt


def superpose(
    coords: np.ndarray,
    reference: np.ndarray | None = None,
    fit_mask: np.ndarray | None = None,
    max_iter: int = 50,
    tol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Superpose frames (n_frames, n_res, 3) onto a reference conformation.

    With ``reference=None`` this is iterative Procrustes: frames are fitted to
    their own running mean until the mean stops moving, which is the standard
    way to avoid biasing the covariance toward an arbitrary single frame.
    With an explicit reference (e.g. the ATLAS starting structure) a single
    fit is done.

    ``fit_mask`` (n_res,) bool restricts the residues used to determine the
    rotation, e.g. to exclude floppy termini; every residue is transformed.

    Returns the superposed coordinates and the final reference.
    """
    coords = np.asarray(coords, dtype=float)
    n_res = coords.shape[1]
    mask = np.ones(n_res, bool) if fit_mask is None else np.asarray(fit_mask, bool)
    if mask.sum() < 3:
        raise ValueError("need at least 3 residues to fit a rotation")

    centred = coords - coords[:, mask].mean(axis=1, keepdims=True)

    def fit(ref: np.ndarray) -> np.ndarray:
        ref_c = ref - ref[mask].mean(axis=0)
        rot = kabsch_rotation(centred[:, mask], np.broadcast_to(ref_c[mask], centred[:, mask].shape))
        return centred @ rot

    if reference is not None:
        ref = np.asarray(reference, dtype=float)
        out = fit(ref)
        return out, ref - ref[mask].mean(axis=0)

    ref = centred[0]
    out = fit(ref)
    for _ in range(max_iter):
        new_ref = out.mean(axis=0)
        new_ref = new_ref - new_ref[mask].mean(axis=0)
        shift = np.sqrt(np.mean(np.sum((new_ref - ref) ** 2, axis=-1)))
        ref = new_ref
        out = fit(ref)
        if shift < tol:
            break
    return out, ref


# --------------------------------------------------------------------------
# Covariance and effective sample size
# --------------------------------------------------------------------------

def per_residue_covariance(aligned: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Per-residue 3x3 positional covariance, (n_frames, n_res, 3) -> (n_res, 3, 3)."""
    aligned = np.asarray(aligned, dtype=float)
    dev = aligned - aligned.mean(axis=0, keepdims=True)
    return np.einsum("tri,trj->rij", dev, dev) / (aligned.shape[0] - ddof)


def integrated_autocorr_time(x: np.ndarray, c: float = 5.0) -> np.ndarray:
    """Integrated autocorrelation time along axis 0, Sokal's automatic window.

    ``tau = 1 + 2 sum_{k=1}^{M} rho(k)`` with the smallest M >= c * tau(M).
    Returns tau with the trailing shape of ``x``; tau >= 1.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    x = x - x.mean(axis=0, keepdims=True)
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, n=nfft, axis=0)
    acov = np.fft.irfft(f * np.conj(f), n=nfft, axis=0)[:n]
    var = acov[0]
    rho = np.divide(acov, var, out=np.zeros_like(acov), where=var > 0)
    taus = 2.0 * np.cumsum(rho, axis=0) - 1.0  # tau(M) for M = 0..n-1
    m = np.arange(n).reshape((n,) + (1,) * (x.ndim - 1))
    ok = m >= c * taus
    first = np.where(ok.any(axis=0), ok.argmax(axis=0), n - 1)
    tau = np.take_along_axis(taus, first[None], axis=0)[0]
    return np.maximum(tau, 1.0)


def effective_sample_size(aligned: np.ndarray, min_n: int = 5) -> np.ndarray:
    """Per-residue effective number of independent frames, shape (n_res,).

    Uses the *largest* autocorrelation time across the three Cartesian
    components, so the estimate is conservative: a smaller n_eff widens the
    sampling-noise null and makes anisotropy harder to claim.
    """
    aligned = np.asarray(aligned, dtype=float)
    tau = integrated_autocorr_time(aligned).max(axis=-1)
    n_eff = np.floor(aligned.shape[0] / tau).astype(int)
    return np.clip(n_eff, min_n, aligned.shape[0])


# --------------------------------------------------------------------------
# Top level
# --------------------------------------------------------------------------

@dataclass
class AtlasCovariances:
    """Per-residue Cα covariances for one ATLAS entry, in a common frame.

    Attributes
    ----------
    pdb_chain : ATLAS identifier, e.g. ``"2erl_A"``.
    sigma_rep : (n_rep, n_res, 3, 3) covariance of each replica about its own
        mean, Å². Use for replica-split tests.
    sigma : (n_res, 3, 3) pooled covariance over every frame of every
        replica about the pooled mean, Å². Includes between-replica spread.
    n_eff_rep : (n_rep, n_res) effective sample size per replica.
    n_frames : frames per replica.
    reference : (n_res, 3) reference conformation the frames were fitted to.
    atlas_rmsf : (n_rep, n_res) ATLAS's published RMSF, or None.
    """

    pdb_chain: str
    sigma_rep: np.ndarray
    sigma: np.ndarray
    n_eff_rep: np.ndarray
    n_frames: int
    reference: np.ndarray
    atlas_rmsf: np.ndarray | None = None

    @property
    def n_res(self) -> int:
        return self.sigma.shape[0]

    @property
    def n_eff(self) -> np.ndarray:
        """Pooled effective sample size per residue (sum over replicas)."""
        return self.n_eff_rep.sum(axis=0)

    @property
    def rmsf_rep(self) -> np.ndarray:
        """RMSF implied by the replica covariances, sqrt(tr Σ), Å."""
        return np.sqrt(np.trace(self.sigma_rep, axis1=-2, axis2=-1))

    def save(self, path: str | Path) -> None:
        extra = {} if self.atlas_rmsf is None else {"atlas_rmsf": self.atlas_rmsf}
        np.savez_compressed(
            path, pdb_chain=self.pdb_chain, sigma_rep=self.sigma_rep,
            sigma=self.sigma, n_eff_rep=self.n_eff_rep,
            n_frames=self.n_frames, reference=self.reference, **extra,
        )

    @classmethod
    def load(cls, path: str | Path) -> "AtlasCovariances":
        with np.load(path, allow_pickle=False) as z:
            return cls(
                pdb_chain=str(z["pdb_chain"]),
                sigma_rep=z["sigma_rep"], sigma=z["sigma"],
                n_eff_rep=z["n_eff_rep"], n_frames=int(z["n_frames"]),
                reference=z["reference"],
                atlas_rmsf=z["atlas_rmsf"] if "atlas_rmsf" in z.files else None,
            )


def covariances_from_trajectories(
    pdb_chain: str,
    ca: np.ndarray,
    reference: np.ndarray | None = None,
    fit_mask: np.ndarray | None = None,
    atlas_rmsf: np.ndarray | None = None,
) -> AtlasCovariances:
    """Core computation, separated from I/O so it can be tested directly.

    ``ca`` is (n_rep, n_frames, n_res, 3). All replicas are fitted jointly to
    one reference so their covariances share a coordinate frame.
    """
    ca = np.asarray(ca, dtype=float)
    n_rep, n_frames, n_res, _ = ca.shape
    flat = ca.reshape(n_rep * n_frames, n_res, 3)
    aligned, ref = superpose(flat, reference=reference, fit_mask=fit_mask)
    aligned = aligned.reshape(n_rep, n_frames, n_res, 3)
    return AtlasCovariances(
        pdb_chain=pdb_chain,
        sigma_rep=np.stack([per_residue_covariance(a) for a in aligned]),
        sigma=per_residue_covariance(aligned.reshape(-1, n_res, 3)),
        n_eff_rep=np.stack([effective_sample_size(a) for a in aligned]),
        n_frames=n_frames,
        reference=ref,
        atlas_rmsf=atlas_rmsf,
    )


def load_atlas_covariances(
    pdb_chain: str,
    data_dir: str | Path = DEFAULT_DATA_DIR,
    reference: str = "mean",
    download: bool = True,
    cache: bool = True,
) -> AtlasCovariances:
    """Per-residue Cα covariances for an ATLAS entry, downloading if needed.

    Parameters
    ----------
    reference : ``"mean"`` for iterative Procrustes to the average structure
        (default), or ``"start"`` to fit once onto the ATLAS starting
        structure (closer to how ATLAS computes its RMSF).
    download : fetch the ``analysis`` archive if it is not on disk.
    cache : read/write ``<data_dir>/<pdb_chain>/covariances_<reference>.npz``.
    """
    if reference not in {"mean", "start"}:
        raise ValueError("reference must be 'mean' or 'start'")
    entry = Path(data_dir) / pdb_chain
    cache_path = entry / f"covariances_{reference}.npz"
    if cache and cache_path.exists():
        return AtlasCovariances.load(cache_path)
    if download:
        entry = fetch_atlas(pdb_chain, data_dir)

    ca = read_ca_trajectories(entry, pdb_chain)
    ref = None
    if reference == "start":
        import mdtraj as md

        start = md.load(str(entry / f"{pdb_chain}.pdb"))
        ref = start.xyz[0, start.top.select("name CA")].astype(float) * 10.0
    result = covariances_from_trajectories(
        pdb_chain, ca, reference=ref, atlas_rmsf=read_atlas_rmsf(entry, pdb_chain)
    )
    if cache:
        result.save(cache_path)
    return result
