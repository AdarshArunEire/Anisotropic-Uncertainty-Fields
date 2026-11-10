import numpy as np

from auf.audit.subspace import project_out_rigid
from auf.phasea.pae import (error_vector, kabsch_superpose, lddt_ca, pae_reldisp_cov,
                            random_rotation_cov, rigid_domain_cov, whitened_quad)


def _coords(n=30, seed=0):
    rng = np.random.default_rng(seed)
    return np.cumsum(rng.normal(size=(n, 3)) * 2.2, axis=0)


def test_kabsch_recovers_rigid_motion():
    x = _coords()
    th = 0.7
    r = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    y = x @ r.T + 5.0
    assert np.allclose(kabsch_superpose(y, x), x, atol=1e-8)


def test_reldisp_reproduces_variogram():
    # Euclidean-embeddable PAE: PAE_ij^2 / 3 = squared distance of points p.
    rng = np.random.default_rng(1)
    p = rng.normal(size=(20, 4))
    d2 = np.sum((p[:, None] - p[None]) ** 2, -1)
    pae = np.sqrt(3 * d2)
    c = pae_reldisp_cov(pae, return_residue=True)
    var_diff = np.diag(c)[:, None] + np.diag(c)[None] - 2 * c
    assert np.allclose(var_diff, d2, atol=1e-8)


def test_single_rigid_domain_is_pure_rigid_motion():
    x = _coords()
    c = rigid_domain_cov(x, np.zeros(len(x), int), 2.0, 0.1)
    assert np.allclose(project_out_rigid(x, c), 0, atol=1e-8)
    c2 = rigid_domain_cov(x, (np.arange(len(x)) >= 15).astype(int), 2.0, 0.1)
    assert np.linalg.norm(project_out_rigid(x, c2)) > 1.0


def test_error_vector_has_no_rigid_part_and_whitening():
    x = _coords()
    rng = np.random.default_rng(2)
    af = x + rng.normal(size=x.shape)
    af_sup, e, _ = error_vector(x, af, np.full(len(x), 90.0))
    assert np.allclose(project_out_rigid(af_sup, e), e, atol=1e-8)
    q, logdet, d = whitened_quad(np.eye(3 * len(x)), e, af_sup)
    assert d == 3 * len(x) - 6 and abs(logdet) < 1e-3
    assert np.isclose(q, e @ e, rtol=1e-4)


def test_random_rotation_keeps_spectrum():
    x = _coords()
    c = rigid_domain_cov(x, (np.arange(len(x)) >= 15).astype(int), 2.0, 0.1) + np.eye(90) * 0.1
    r = random_rotation_cov(c, x, np.random.default_rng(0))
    pc = project_out_rigid(x, c)
    assert np.allclose(np.sort(np.linalg.eigvalsh(pc))[6:], np.sort(np.linalg.eigvalsh(r))[6:])


def test_lddt_identity_and_drop():
    x = _coords()
    assert lddt_ca(x, x) == 1.0
    assert lddt_ca(x + np.random.default_rng(3).normal(size=x.shape) * 2, x) < 0.8
