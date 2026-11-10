"""Second Stage-1 draw for the crystal comparison: ATLAS entries at <= 1.2 Å.

Reconstructed 2026-11-10: the original list (data/highres_entries.txt) was in
an ignored folder and lost with its session. The draw follows the description
in notes/stage1_audit.md (seeded draw of 60 entries at <= 1.2 Å resolution,
selected on resolution only): entries of 60-300 residues at <= 1.2 Å that are
not in the random 150, 60 drawn with seed 0. It may differ from the original
list. Both lists are now committed under manifests/.

Usage: python3 experiments/audit/00b_highres_entries.py
then   python3 experiments/audit/00_build_dataset.py \\
           --entries $(cat manifests/atlas_highres60.txt) --out data/audit_highres
"""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import numpy as np

INFO = Path("data/2024_11_18_ATLAS_info.tsv")
spec = importlib.util.spec_from_file_location("b", Path(__file__).with_name("00_build_dataset.py"))
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)

random150 = b.sample_entries(INFO, 150, 60, 300, 0)
with open(INFO) as fh:
    rows = [r for r in csv.DictReader(fh, delimiter="\t")
            if 60 <= int(r["length"]) <= 300 and r["PDB_resolution"] not in ("", "NA")
            and float(r["PDB_resolution"]) <= 1.2 and r["PDB"] not in set(random150)]
rng = np.random.default_rng(0)
hi = [rows[i]["PDB"] for i in rng.choice(len(rows), size=min(60, len(rows)), replace=False)]
out = Path("manifests")
out.mkdir(exist_ok=True)
(out / "atlas_random150.txt").write_text("\n".join(random150) + "\n")
(out / "atlas_highres60.txt").write_text("\n".join(hi) + "\n")
print(f"{len(rows)} eligible at <= 1.2 Å; wrote 150 + {len(hi)} entries to manifests/")
