"""Phase A2 splits: MMseqs2 clusters (30% identity, 80% coverage) -> train/val/test.

Reads the accepted ``<afset>/*.npz`` (sequence of mapped residues), clusters
with ``mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0``, and assigns
whole clusters to splits with ``auf.splits.assign_splits`` (hash of sorted
member names, seed ``auf-afset-splits-v1``, 60/20/20), as pre-registered.
Writes ``splits/afset_splits.tsv`` (name, uniprot, cluster representative,
split), which is committed so later steps never re-cluster.

Usage: python3 experiments/phaseA/a2_splits.py --mmseqs /path/to/mmseqs
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from auf.splits.atlas import assign_splits


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--afset", default="data/afset")
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--out", default="splits/afset_splits.tsv")
    args = p.parse_args()

    with open(Path(args.afset) / "manifest.tsv") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["status"] == "accepted"]
    seqs = {}
    for r in rows:
        d = np.load(Path(args.afset) / f"{r['name']}.npz")
        seqs[r["name"]] = str(d["sequence"])
    acc = {r["name"]: r["uniprot"] for r in rows}

    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "seqs.fasta"
        fa.write_text("".join(f">{n}\n{s}\n" for n, s in sorted(seqs.items())))
        subprocess.run([args.mmseqs, "easy-cluster", str(fa), str(Path(tmp) / "clu"),
                        str(Path(tmp) / "work"), "--min-seq-id", "0.3", "-c", "0.8",
                        "--cov-mode", "0", "--threads", "4", "-v", "1"], check=True)
        rep = {}
        with open(Path(tmp) / "clu_cluster.tsv") as fh:
            for line in fh:
                r, m = line.split()
                rep[m] = r

    comps: dict[str, list[str]] = {}
    for m, r in rep.items():
        comps.setdefault(r, []).append(m)
    components = sorted((sorted(c) for c in comps.values()), key=lambda c: c[0])
    split = assign_splits(components, seed="auf-afset-splits-v1")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("name\tuniprot\tcluster\tsplit\n")
        for n in sorted(seqs):
            fh.write(f"{n}\t{acc[n]}\t{rep[n]}\t{split[n]}\n")
    counts = {s: sum(v == s for v in split.values()) for s in ("train", "val", "test")}
    print(f"{len(seqs)} proteins, {len(components)} clusters: {counts}")


if __name__ == "__main__":
    main()
