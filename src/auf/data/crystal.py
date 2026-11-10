"""Experimental anisotropic displacement parameters (ANISOU) for ATLAS entries.

High-resolution crystal structures deposit a 3x3 anisotropic displacement
tensor U per atom. It is the only *experimental* per-residue covariance
available at scale, and it is independent of MD. It is not the same object
(it also carries lattice disorder and rigid-body/TLS motion), so agreement in
*direction* beyond the isotropic surrogate is the test, not equality.

Residue mapping uses ATLAS's ``<pdb_chain>_corresp.tsv``: row k is ATLAS
residue k, ``PDB_num`` is the author numbering, and ``gap == 1`` marks
residues modelled in by ATLAS and absent from the crystal.
"""
from __future__ import annotations

import csv
import urllib.request
from pathlib import Path

import numpy as np

__all__ = ["fetch_pdb", "parse_ca_anisou", "crystal_for_atlas", "fit_tls", "tls_design"]

RCSB = "https://files.rcsb.org/download"


def fetch_pdb(pdb_id: str, data_dir: str | Path, timeout: float = 120.0) -> Path | None:
    """Download the deposited PDB-format file; None if RCSB has none."""
    out = Path(data_dir) / f"{pdb_id.upper()}.pdb"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(f"{RCSB}/{pdb_id.upper()}.pdb", timeout=timeout) as r:
            out.write_bytes(r.read())
    except Exception:  # large structures have no PDB-format file
        return None
    return out


def parse_ca_anisou(path: str | Path, chain: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Map author residue id (resSeq + iCode) -> (CA xyz, U 3x3) for one chain.

    First model and first altloc only. U is in Å^2.
    """
    xyz: dict[str, np.ndarray] = {}
    u: dict[str, np.ndarray] = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec == "ENDMDL":
                break
            if rec not in ("ATOM  ", "ANISOU"):
                continue
            if line[12:16].strip() != "CA" or line[21] != chain:
                continue
            if line[16] not in (" ", "A"):
                continue
            key = (line[22:26].strip() + line[26].strip())
            if rec == "ATOM  " and key not in xyz:
                xyz[key] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elif rec == "ANISOU" and key not in u:
                v = [int(line[28 + 7 * k: 35 + 7 * k]) for k in range(6)]
                u11, u22, u33, u12, u13, u23 = (x * 1e-4 for x in v)
                u[key] = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
    return {k: (xyz[k], u[k]) for k in u if k in xyz}


def crystal_for_atlas(
    pdb_chain: str, entry_dir: str | Path, pdb_dir: str | Path
) -> tuple[np.ndarray, np.ndarray] | None:
    """Crystal CA coordinates and U per ATLAS residue, NaN where missing.

    Returns (xyz (n_res, 3), U (n_res, 3, 3)) or None if the entry has no
    usable ANISOU records.
    """
    pdb_id, chain = pdb_chain.split("_")
    corresp = Path(entry_dir) / f"{pdb_chain}_corresp.tsv"
    if not corresp.exists():
        return None
    with open(corresp) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    path = fetch_pdb(pdb_id, pdb_dir)
    if path is None:
        return None
    records = parse_ca_anisou(path, chain)
    if not records:
        return None
    n = len(rows)
    xyz = np.full((n, 3), np.nan)
    u = np.full((n, 3, 3), np.nan)
    for k, row in enumerate(rows):
        if row.get("gap", "0") != "0":
            continue
        rec = records.get(str(row["PDB_num"]).strip())
        if rec is not None:
            xyz[k], u[k] = rec
    if np.isfinite(xyz[:, 0]).sum() < 10:
        return None
    return xyz, u


def _skew(r: np.ndarray) -> np.ndarray:
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    o = np.zeros_like(x)
    return np.stack([np.stack([o, z, -y], -1),
                     np.stack([-z, o, x], -1),
                     np.stack([y, -x, o], -1)], -2)


_IU = (np.array([0, 1, 2, 0, 0, 1]), np.array([0, 1, 2, 1, 2, 2]))


def tls_design(xyz: np.ndarray) -> np.ndarray:
    """Linear map from 21 TLS parameters (T 6, L 6, S 9) to the 6 unique
    components of U at each position: (n * 6, 21)."""
    r = np.asarray(xyz, float) - np.mean(xyz, axis=0)
    a = _skew(r)  # rigid rotation by angle vector w moves r by a @ w
    n = len(r)
    cols = []
    basis6 = []
    for i, j in zip(*_IU):
        e = np.zeros((3, 3))
        e[i, j] = e[j, i] = 1.0
        basis6.append(e)
    for e in basis6:  # T
        cols.append(np.broadcast_to(e, (n, 3, 3)))
    for e in basis6:  # L
        cols.append(a @ e @ np.swapaxes(a, -1, -2))
    for i in range(3):  # S (full 3x3; its trace is not identifiable)
        for j in range(3):
            e = np.zeros((3, 3))
            e[i, j] = 1.0
            m = a @ e
            cols.append(m + np.swapaxes(m, -1, -2))
    return np.stack([c[:, _IU[0], _IU[1]].reshape(-1) for c in cols], axis=-1)


def fit_tls(xyz: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Least-squares whole-chain TLS fit to per-residue U; returns U_TLS.

    One rigid group, 20 identifiable parameters: the crystallographers'
    parsimonious model of anisotropy (Schomaker & Trueblood 1968).
    """
    design = tls_design(xyz)
    target = np.asarray(u)[:, _IU[0], _IU[1]].reshape(-1)
    coef, *_ = np.linalg.lstsq(design, target, rcond=None)
    flat = (design @ coef).reshape(len(xyz), 6)
    out = np.zeros((len(xyz), 3, 3))
    out[:, _IU[0], _IU[1]] = flat
    out[:, _IU[1], _IU[0]] = flat
    return out
