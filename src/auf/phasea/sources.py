"""Direction sources and controls for one protein (roadmap A3/A4).

``load_protein`` reads an A2 ``.npz`` and computes the error vector;
``build_sources`` returns {name: 3N x 3N covariance (Å², before the global
scale)} using parameters frozen on train (``params/phasea_a3.json``).
"""
from __future__ import annotations

import numpy as np

from auf.phasea import pae as P

# Sources entering A4. Gate candidates: the three PAE constructions (one is
# chosen on validation) and ANM. Others are controls or secondary.
PAE_SOURCES = ["pae_reldisp", "pae_domains", "pae_lowrank"]
SOURCES = PAE_SOURCES + ["anm", "anm_recal", "plddt_scalar"]
CONTROLS = ["oracle_scalar", "pae_rigid_control"]
# Broad prior of the practical-standard control (fixed, not fitted): each PAE
# domain may move 10 Å and rotate ~20 degrees.
CONTROL_SIG_T, CONTROL_SIG_R = 10.0, 0.35


def load_protein(path, mode: str = "core") -> dict:
    d = np.load(path)
    exp, af, plddt = d["exp_xyz"], d["af_xyz"], d["plddt"]
    af_sup, e, used_core = P.error_vector(exp, af, plddt, mode=mode)
    return {"name": str(d["name"]), "exp": exp, "af": af_sup, "plddt": plddt,
            "pae": d["pae"], "e": e, "used_core": used_core,
            "method": str(d["method"]), "n": len(exp)}


def build_sources(p: dict, params: dict, names=None, oracle: bool = True) -> dict:
    names = list(SOURCES + CONTROLS) if names is None else list(names)
    x, plddt, pae = p["af"], p["plddt"], p["pae"]
    noise = P.diag_cov(P.plddt_variance(plddt, params["plddt_table"]))
    out = {}
    need_a = {"pae_reldisp", "pae_lowrank"} & set(names)
    a = P.pae_reldisp_cov(pae) if need_a else None
    need_dom = {"pae_domains", "pae_rigid_control"} & set(names)
    labels = P.pae_domains(pae) if need_dom else None
    for n in names:
        if n == "pae_reldisp":
            out[n] = a
        elif n == "pae_domains":
            out[n] = P.rigid_domain_cov(x, labels, params["dom_sig_t"],
                                        params["dom_sig_r"]) + noise
        elif n == "pae_lowrank":
            out[n] = params["lowrank_w"] * P.lowrank(a, params["lowrank_r"]) + noise
        elif n == "anm":
            out[n] = P.anm_spectrum_cov(x, 1.0)
        elif n == "anm_recal":
            out[n] = P.anm_spectrum_cov(x, params["anm_power"])
        elif n == "plddt_scalar":
            out[n] = noise
        elif n == "oracle_scalar":
            if oracle:
                v = np.maximum(np.sum(p["e"].reshape(-1, 3) ** 2, 1) / 3.0, 1e-4)
                out[n] = P.diag_cov(v)
        elif n == "pae_rigid_control":
            out[n] = P.rigid_domain_cov(x, labels, CONTROL_SIG_T, CONTROL_SIG_R) + noise
        else:
            raise KeyError(n)
    return out
