"""Phase A2: AlphaFold error set (criteria in notes/phaseA_preregistration.md).

1. RCSB search: protein entities released >= 2023-01-01, one protein entity
   and one instance in the entry, 60-500 residues, X-ray <= 2.5 Å or
   cryo-EM <= 3.0 Å.
2. Seeded shuffle; process entries until ``--target`` proteins are accepted,
   keeping one entry per UniProt accession (first in shuffled order).
3. Per entry: SIFTS (PDBe) residue mapping to UniProt; AlphaFold DB model,
   pLDDT and PAE; experimental Cα from mmCIF. Accepted if
   - single UniProt mapping for the chain;
   - AlphaFold model created *before* the entry's release date;
   - <= 10% of the modelled range missing in the experimental chain;
   - AlphaFold covers >= 90% of the observed residues;
   - mean pLDDT over mapped residues >= 50.
4. Writes ``<out>/<PDB>_<chain>.npz`` with matched Cα arrays, pLDDT, PAE
   (mapped residues only), sequence and metadata, plus ``manifest.tsv``
   (accepted and rejected, with reason).

Splits (MMseqs2, 30% identity / 80% coverage) are made by a2_splits.py.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import random
import urllib.request
from multiprocessing import Pool
from pathlib import Path

import numpy as np

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_CIF = "https://files.rcsb.org/download/{}.cif.gz"
SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{}"
AFDB = "https://alphafold.ebi.ac.uk/api/prediction/{}"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "MSE": "M"}


def _get(url: str, timeout: float = 60.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "auf-research/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def search_entities() -> list[str]:
    def term(attr, op, value):
        return {"type": "terminal", "service": "text",
                "parameters": {"attribute": attr, "operator": op, "value": value}}

    method = lambda m, res: {"type": "group", "logical_operator": "and", "nodes": [
        term("exptl.method", "exact_match", m),
        term("rcsb_entry_info.resolution_combined", "less_or_equal", res)]}
    query = {"type": "group", "logical_operator": "and", "nodes": [
        term("rcsb_accession_info.initial_release_date", "greater_or_equal", "2023-01-01"),
        term("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
        term("rcsb_entry_info.polymer_entity_count", "equals", 1),
        term("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
        term("entity_poly.rcsb_sample_sequence_length", "range",
             {"from": 60, "to": 500, "include_lower": True, "include_upper": True}),
        {"type": "group", "logical_operator": "or", "nodes": [
            method("X-RAY DIFFRACTION", 2.5), method("ELECTRON MICROSCOPY", 3.0)]},
    ]}
    ids, start = [], 0
    while True:
        body = {"query": query, "return_type": "polymer_entity",
                "request_options": {"paginate": {"start": start, "rows": 10000},
                                    "results_content_type": ["experimental"]}}
        req = urllib.request.Request(RCSB_SEARCH, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.load(r)
        ids += [x["identifier"] for x in d.get("result_set", [])]
        start += 10000
        if start >= d["total_count"]:
            return sorted(ids)


def ca_from_cif(cif_bytes: bytes, chain: str) -> dict[int, np.ndarray]:
    """Cα by label_seq_id (1-based index into the entity sequence)."""
    import gemmi

    doc = gemmi.cif.read_string(gzip.decompress(cif_bytes).decode())
    st = gemmi.make_structure_from_block(doc.sole_block())
    out = {}
    for ch in st[0]:
        if ch.name != chain:
            continue
        for res in ch:
            ca = res.find_atom("CA", "*")
            if ca is not None and res.het_flag != "H" and res.label_seq is not None:
                out.setdefault(int(res.label_seq), np.array(ca.pos.tolist()))
    return out


def ca_from_af(pdb_text: str) -> tuple[dict[int, np.ndarray], dict[int, float], dict[int, str]]:
    xyz, plddt, aa = {}, {}, {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            i = int(line[22:26])
            xyz[i] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            plddt[i] = float(line[60:66])
            aa[i] = THREE.get(line[17:20], "X")
    return xyz, plddt, aa


def process(args) -> dict:
    entity_id, out_dir = args
    pdb_id = entity_id.split("_")[0].lower()
    rec = {"entity": entity_id, "status": "rejected"}
    try:
        entry = json.loads(_get(RCSB_ENTRY.format(pdb_id)))
        release = entry["rcsb_accession_info"]["initial_release_date"][:10]
        rec["release"] = release
        rec["method"] = entry["exptl"][0]["method"]
        rec["resolution"] = (entry.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]

        sifts = json.loads(_get(SIFTS.format(pdb_id))).get(pdb_id, {}).get("UniProt", {})
        if len(sifts) != 1:
            rec["reason"] = f"uniprot mappings: {len(sifts)}"
            return rec
        acc, info = next(iter(sifts.items()))
        maps = info["mappings"]
        chains = {m["chain_id"] for m in maps}
        if len(chains) != 1:
            rec["reason"] = "multiple chains mapped"
            return rec
        chain = chains.pop()
        rec.update(uniprot=acc, chain=chain)

        pred = json.loads(_get(AFDB.format(acc)))
        if not pred:
            rec["reason"] = "no AFDB model"
            return rec
        pred = pred[0]
        created = pred.get("modelCreatedDate", "")[:10]
        rec["af_created"] = created
        if not created or created >= release:
            rec["reason"] = "AF model not older than release"
            return rec

        # experimental residue (label_seq) -> uniprot position via SIFTS segments.
        # "Modelled range" = first..last observed residue, so disordered
        # termini do not count as missing (pre-registration wording).
        exp_ca = ca_from_cif(_get(RCSB_CIF.format(pdb_id)), chain)
        mapped, labels = {}, []
        for m in maps:
            l0, l1 = m["start"]["residue_number"], m["end"]["residue_number"]
            u0 = m["unp_start"]
            for lab in range(l0, l1 + 1):
                if lab in exp_ca:
                    mapped[u0 + (lab - l0)] = exp_ca[lab]
                    labels.append(lab)
        modelled = (max(labels) - min(labels) + 1) if labels else 0
        if modelled == 0 or len(mapped) < 0.9 * modelled:
            rec["reason"] = f"missing residues {1 - len(mapped) / max(modelled, 1):.2f}"
            return rec

        af_xyz, af_plddt, af_aa = ca_from_af(_get(pred["pdbUrl"]).decode())
        shift = pred.get("uniprotStart", 1) - 1
        common = sorted(u for u in mapped if (u - shift) in af_xyz)
        if len(common) < 0.9 * len(mapped):
            rec["reason"] = "AF coverage < 90%"
            return rec
        af_idx = [u - shift for u in common]
        plddt = np.array([af_plddt[i] for i in af_idx])
        if plddt.mean() < 50:
            rec["reason"] = f"mean pLDDT {plddt.mean():.0f}"
            return rec

        pae_doc = json.loads(_get(pred["paeDocUrl"]))
        pae_doc = pae_doc[0] if isinstance(pae_doc, list) else pae_doc
        pae = np.asarray(pae_doc.get("predicted_aligned_error") or pae_doc.get("distance"), float)
        if pae.ndim == 1:  # old flat format
            nn = int(round(np.sqrt(pae.size)))
            pae = pae.reshape(nn, nn)
        sel = np.array(af_idx) - 1
        name = f"{pdb_id.upper()}_{chain}"
        np.savez_compressed(
            Path(out_dir) / f"{name}.npz",
            name=name, uniprot=acc, release=release, af_created=created,
            method=rec["method"], resolution=rec["resolution"] or np.nan,
            uniprot_pos=np.array(common),
            exp_xyz=np.stack([mapped[u] for u in common]),
            af_xyz=np.stack([af_xyz[i] for i in af_idx]),
            plddt=plddt, pae=pae[np.ix_(sel, sel)],
            sequence="".join(af_aa[i] for i in af_idx),
        )
        rec.update(status="accepted", name=name, n_res=len(common),
                   mean_plddt=round(float(plddt.mean()), 1))
        return rec
    except Exception as exc:
        rec["reason"] = f"error: {type(exc).__name__}: {str(exc)[:80]}"
        return rec


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--target", type=int, default=1500)
    p.add_argument("--out", default="data/afset")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ids = search_entities()
    random.Random(args.seed).shuffle(ids)
    (out / "candidates.txt").write_text("\n".join(ids) + "\n")
    print(f"{len(ids)} candidate entities", flush=True)

    seen_acc, rows, accepted = set(), [], 0
    batch = args.workers * 4
    with Pool(args.workers) as pool:
        for b in range(0, len(ids), batch):
            for rec in pool.map(process, [(e, out) for e in ids[b:b + batch]]):
                if rec["status"] == "accepted":
                    if rec["uniprot"] in seen_acc:
                        (out / f"{rec['name']}.npz").unlink(missing_ok=True)
                        rec.update(status="rejected", reason="duplicate uniprot")
                    else:
                        seen_acc.add(rec["uniprot"])
                        accepted += 1
                rows.append(rec)
            print(f"processed {len(rows)}, accepted {accepted}", flush=True)
            if accepted >= args.target:
                break

    cols = ["entity", "status", "reason", "name", "uniprot", "chain", "release",
            "af_created", "method", "resolution", "n_res", "mean_plddt"]
    with open(out / "manifest.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"done: {accepted} accepted of {len(rows)} processed")


if __name__ == "__main__":
    main()
