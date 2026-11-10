"""Loader tests on synthetic trajectories: no network, no mdtraj."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from auf.audit import isotropize, wishart_isotropic_null_airm
from auf.data.atlas import (
    AtlasCovariances,
    covariances_from_trajectories,
    effective_sample_size,
    integrated_autocorr_time,
    kabsch_rotation,
    per_residue_covariance,
    superpose,
)


def rigid_scramble(frames, rng):
    """Apply an independent random rotation and translation to every frame."""
    rots = Rotation.random(len(frames), random_state=rng).as_matrix()
    shifts = rng.normal(scale=20.0, size=(len(frames), 1, 3))
    return np.einsum("tri,tji->trj", frames, rots) + shifts


def synthetic_ensemble(rng, n_frames=4000, n_res=30):
    """Rigid scaffold plus per-residue anisotropic Gaussian fluctuations."""
    ref = rng.normal(scale=8.0, size=(n_res, 3))
    axes = Rotation.random(n_res, random_state=rng).as_matrix()
    evals = np.stack([np.full(n_res, 1.0), np.full(n_res, 0.2), np.full(n_res, 0.05)], -1)
    sigma = (axes * evals[:, None, :]) @ np.swapaxes(axes, -1, -2)
    chol = np.linalg.cholesky(sigma)
    noise = np.einsum("rij,trj->tri", chol, rng.normal(size=(n_frames, n_res, 3)))
    return ref + noise, sigma


def test_kabsch_recovers_known_rotation():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(20, 3))
    x -= x.mean(0)
    r = Rotation.random(random_state=rng).as_matrix()
    rot = kabsch_rotation(x, x @ r)
    np.testing.assert_allclose(rot, r, atol=1e-10)
    assert np.isclose(np.linalg.det(rot), 1.0)


def test_kabsch_never_returns_reflection():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(20, 3))
    x -= x.mean(0)
    mirror = x * np.array([1.0, 1.0, -1.0])
    assert np.isclose(np.linalg.det(kabsch_rotation(x, mirror)), 1.0)


def test_superpose_removes_rigid_motion():
    rng = np.random.default_rng(2)
    base = rng.normal(scale=8.0, size=(25, 3))
    frames = rigid_scramble(np.broadcast_to(base, (50, 25, 3)), rng)
    aligned, _ = superpose(frames)
    np.testing.assert_allclose(aligned.std(axis=0), 0.0, atol=1e-8)


def test_superposition_is_invariant_to_rigid_scramble():
    """Per-frame rigid motion must not change the result at all."""
    rng = np.random.default_rng(3)
    frames, _ = synthetic_ensemble(rng)
    clean = per_residue_covariance(superpose(frames)[0])
    scrambled = per_residue_covariance(superpose(rigid_scramble(frames, rng))[0])
    np.testing.assert_allclose(
        np.linalg.eigvalsh(scrambled), np.linalg.eigvalsh(clean), atol=1e-10
    )


def test_covariance_recovered_when_superposition_bias_is_small():
    """With many residues the fitted rotation absorbs little of any one
    residue's fluctuation, and per-residue eigenvalues come back intact."""
    rng = np.random.default_rng(3)
    frames, sigma = synthetic_ensemble(rng, n_res=300)
    est = per_residue_covariance(superpose(rigid_scramble(frames, rng))[0])
    np.testing.assert_allclose(
        np.linalg.eigvalsh(est), np.linalg.eigvalsh(sigma), rtol=0.15
    )


def test_superposition_bias_is_directional_on_small_proteins():
    """Documents a known artifact rather than asserting correctness.

    On a small protein the fitted rotation absorbs part of each residue's
    fluctuation and leaks others' into it. The bias is in *eigenvalues*, so
    it is a shape distortion, not just a scale one: anything measured about
    anisotropy on small proteins carries it. If this test starts failing
    because the bias vanished, the superposition changed.
    """
    rng = np.random.default_rng(3)
    frames, sigma = synthetic_ensemble(rng, n_res=30)
    est = per_residue_covariance(superpose(rigid_scramble(frames, rng))[0])
    ratio = np.linalg.eigvalsh(est) / np.linalg.eigvalsh(sigma)
    assert ratio.max() > 1.3 or ratio.min() < 0.8


def test_joint_fit_puts_replicas_in_a_common_frame():
    """Replicas scrambled independently must come back with matching
    covariances, not covariances that differ by a rotation."""
    rng = np.random.default_rng(4)
    frames, sigma = synthetic_ensemble(rng, n_frames=6000)
    reps = np.stack([rigid_scramble(f, rng) for f in np.split(frames, 2)])
    out = covariances_from_trajectories("syn_A", reps)
    a, b = out.sigma_rep
    assert np.median(np.linalg.norm(a - b, axis=(-2, -1))) < 0.1
    assert out.sigma_rep.shape == (2, 30, 3, 3)
    assert out.n_eff_rep.shape == (2, 30)


def test_autocorr_time_of_white_noise_is_one():
    rng = np.random.default_rng(5)
    tau = integrated_autocorr_time(rng.normal(size=(20000, 4)))
    np.testing.assert_allclose(tau, 1.0, atol=0.15)


def test_autocorr_time_of_ar1_matches_theory():
    rng = np.random.default_rng(6)
    phi, n = 0.9, 100_000
    x = np.empty(n)
    x[0] = 0.0
    eps = rng.normal(size=n)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + eps[t]
    expected = (1 + phi) / (1 - phi)  # 19
    assert abs(integrated_autocorr_time(x) - expected) / expected < 0.15


def test_effective_sample_size_shrinks_with_correlation():
    rng = np.random.default_rng(7)
    white = rng.normal(size=(2000, 5, 3))
    smooth = np.cumsum(white, axis=0) * 0.0
    kernel = np.ones(20) / 20
    for r in range(5):
        for i in range(3):
            smooth[:, r, i] = np.convolve(white[:, r, i], kernel, mode="same")
    assert np.all(effective_sample_size(smooth) < effective_sample_size(white) / 5)


def test_save_load_roundtrip(tmp_path):
    rng = np.random.default_rng(8)
    frames, _ = synthetic_ensemble(rng, n_frames=600, n_res=10)
    out = covariances_from_trajectories(
        "syn_A", frames.reshape(2, 300, 10, 3), atlas_rmsf=np.ones((2, 10))
    )
    out.save(tmp_path / "c.npz")
    back = AtlasCovariances.load(tmp_path / "c.npz")
    assert back.pdb_chain == "syn_A"
    np.testing.assert_array_equal(back.sigma_rep, out.sigma_rep)
    np.testing.assert_array_equal(back.atlas_rmsf, out.atlas_rmsf)


def test_isotropic_null_shrinks_with_frames_and_is_scale_free():
    rng = np.random.default_rng(9)
    few = wishart_isotropic_null_airm(10, n_draws=500, rng=rng)
    many = wishart_isotropic_null_airm(1000, n_draws=500, rng=rng)
    assert np.median(many) < np.median(few) / 5


def test_isotropic_null_does_not_flag_isotropic_truth():
    """A sphere sampled at n frames must look no more anisotropic than the
    null says it should, i.e. p-values are roughly uniform."""
    from auf.audit import airm_distance

    rng = np.random.default_rng(10)
    n = 40
    x = rng.normal(size=(300, n, 3))
    s = np.einsum("rti,rtj->rij", x, x) / n
    obs = airm_distance(s, isotropize(s))
    null = wishart_isotropic_null_airm(n, n_draws=4000, rng=rng)
    p = (null[None, :] >= obs[:, None]).mean(axis=1)
    assert 0.03 < (p < 0.05).mean() < 0.08


def test_tls_fit_recovers_rigid_body_u():
    from auf.data.crystal import fit_tls, tls_design

    rng = np.random.default_rng(30)
    xyz = rng.normal(scale=10.0, size=(60, 3))
    coef = rng.normal(size=21) * 0.01
    coef[:6] = [0.2, 0.2, 0.2, 0.01, 0.0, 0.02]
    flat = (tls_design(xyz) @ coef).reshape(60, 6)
    u = np.zeros((60, 3, 3))
    iu = (np.array([0, 1, 2, 0, 0, 1]), np.array([0, 1, 2, 1, 2, 2]))
    u[:, iu[0], iu[1]] = flat
    u[:, iu[1], iu[0]] = flat
    np.testing.assert_allclose(fit_tls(xyz, u), u, atol=1e-10)


def test_anm_rigid_modes_removed_and_psd():
    from auf.data.anm import anm_covariance, anm_hessian

    rng = np.random.default_rng(31)
    ca = np.cumsum(rng.normal(scale=2.2, size=(40, 3)), axis=0)
    w = np.linalg.eigvalsh(anm_hessian(ca, cutoff=20.0))
    assert np.sum(np.abs(w) < 1e-8) == 6
    assert np.all(np.linalg.eigvalsh(anm_covariance(ca, cutoff=20.0)) > 0)


def test_parse_ca_anisou(tmp_path):
    from auf.data.crystal import parse_ca_anisou

    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      2  CA  ALA A  10       1.000   2.000   3.000  1.00 10.00           C\n"
        "ANISOU    2  CA  ALA A  10     1000   2000   3000    100    200    300       C\n"
        "ATOM      9  CA  GLY B  10       9.000   9.000   9.000  1.00 10.00           C\n"
    )
    rec = parse_ca_anisou(pdb, "A")
    xyz, u = rec["10"]
    np.testing.assert_allclose(xyz, [1, 2, 3])
    np.testing.assert_allclose(u, [[0.1, 0.01, 0.02], [0.01, 0.2, 0.03], [0.02, 0.03, 0.3]])
    assert "10" not in parse_ca_anisou(pdb, "B")


def test_splits_keep_homologues_together(tmp_path):
    from auf.splits import assign_splits, atlas_homology_components

    tsv = tmp_path / "info.tsv"
    tsv.write_text(
        "PDB\tECOD_ID\n"
        "a_A\t1.1.1.1\n"
        "b_A\t1.1.2.1, 2.5.1.1\n"   # links H-groups 1.1 and 2.5
        "c_A\t2.5.3.3\n"
        "d_A\t7.7.1.1\n"
    )
    comps = atlas_homology_components(tsv)
    assert ["a_A", "b_A", "c_A"] in comps and ["d_A"] in comps
    s = assign_splits(comps, giant=1.0)
    assert s["a_A"] == s["b_A"] == s["c_A"]
    assert assign_splits(comps, giant=1.0) == s  # deterministic
