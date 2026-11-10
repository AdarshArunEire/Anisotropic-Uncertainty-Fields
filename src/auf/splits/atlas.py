"""Homology-aware train / validation / test splits for ATLAS entries.

Proteins are grouped by ECOD homology group (H-group, "X.H" of the ECOD_ID),
which is stricter than a 30% sequence-identity cut: members of an H-group
are believed homologous even when sequence identity is undetectable. A
multi-domain protein links every H-group it contains, so groups are the
connected components of the protein-H-group graph. Whole components are
assigned to one split, so no homologue crosses a split boundary.

Assignment is a deterministic function of the component's sorted member
ids and a seed (hash-based), so it does not depend on any result.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

__all__ = ["atlas_homology_components", "assign_splits", "load_splits"]


def _h_groups(ecod_id: str) -> set[str]:
    out = set()
    for part in ecod_id.replace('"', "").split(","):
        f = part.strip().split(".")
        if len(f) >= 2 and f[0] not in ("", "-"):
            out.add(f"{f[0]}.{f[1]}")
    return out


def atlas_homology_components(info_tsv: str | Path) -> list[list[str]]:
    """Connected components of ATLAS entries linked by shared ECOD H-groups."""
    with open(info_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        parent[find(a)] = find(b)

    for r in rows:
        node = "P:" + r["PDB"]
        find(node)
        for h in _h_groups(r.get("ECOD_ID", "")):
            union(node, "H:" + h)
    comps: dict[str, list[str]] = {}
    for r in rows:
        comps.setdefault(find("P:" + r["PDB"]), []).append(r["PDB"])
    return sorted((sorted(c) for c in comps.values()), key=lambda c: c[0])


def assign_splits(
    components: list[list[str]],
    fractions: tuple[float, float, float] = (0.6, 0.2, 0.2),
    seed: str = "auf-splits-v1",
    giant: float = 0.05,
) -> dict[str, str]:
    """Map entry -> 'train' / 'val' / 'test' by hashing each component.

    Components holding more than ``giant`` of all entries (superfamilies
    chained together by multi-domain proteins) go to train whole, since
    splitting them would leak homologues; the remaining components are
    hashed with fractions adjusted so the overall split hits ``fractions``.
    """
    total = sum(len(c) for c in components)
    big = [c for c in components if len(c) > giant * total]
    n_big = sum(len(c) for c in big)
    rest = total - n_big
    f_train = max(fractions[0] * total - n_big, 0.0) / rest
    f_val = fractions[1] * total / rest
    cut_train, cut_val = f_train, f_train + f_val
    out = {e: "train" for c in big for e in c}
    for comp in components:
        if len(comp) > giant * total:
            continue
        h = hashlib.sha256((seed + "|" + ",".join(comp)).encode()).hexdigest()
        u = int(h[:12], 16) / float(1 << 48)
        split = "train" if u < cut_train else "val" if u < cut_val else "test"
        for e in comp:
            out[e] = split
    return out


def load_splits(path: str | Path = "splits/atlas_splits.tsv") -> dict[str, str]:
    with open(path) as fh:
        return {r["pdb_chain"]: r["split"] for r in csv.DictReader(fh, delimiter="\t")}
