import numpy as np
import pytest

from auf.audit import (
    airm_distance,
    chi2_calibration_test,
    fractional_anisotropy,
    isotropize,
    logfrob_distance,
    mahalanobis_sq,
    wishart_null_airm,
)


def random_spd(rng, d=3, n=None, scale=1.0):
    shape = (d, d) if n is None else (n, d, d)
    a = rng.normal(size=shape)
    sigma = np.swapaxes(a, -1, -2) @ a + d * np.eye(d)
    return scale * sigma


def test_isotropize_preserves_trace():
    rng = np.random.default_rng(0)
    s = random_spd(rng, n=16)
    iso = isotropize(s)
    np.testing.assert_allclose(
        np.trace(s, axis1=-2, axis2=-1), np.trace(iso, axis1=-2, axis2=-1)
    )


def test_isotropize_is_isotropic():
    rng = np.random.default_rng(1)
    iso = isotropize(random_spd(rng, n=8))
    assert np.allclose(fractional_anisotropy(iso), 0.0, atol=1e-10)


def test_fractional_anisotropy_bounds():
    rng = np.random.default_rng(2)
    fa = fractional_anisotropy(random_spd(rng, n=64))
    assert np.all(fa >= 0.0)
    assert np.all(fa < 1.0)


def test_fractional_anisotropy_degenerate_tends_to_one():
    sigma = np.diag([1.0, 1e-8, 1e-8])
    assert fractional_anisotropy(sigma) > 0.99


@pytest.mark.parametrize("dist", [airm_distance, logfrob_distance])
def test_distance_identity_is_zero(dist):
    rng = np.random.default_rng(3)
    s = random_spd(rng, n=8)
    np.testing.assert_allclose(dist(s, s), 0.0, atol=1e-8)


@pytest.mark.parametrize("dist", [airm_distance, logfrob_distance])
def test_distance_symmetric(dist):
    rng = np.random.default_rng(4)
    a, b = random_spd(rng, n=8), random_spd(rng, n=8)
    np.testing.assert_allclose(dist(a, b), dist(b, a), rtol=1e-8)


def test_airm_is_scale_invariant():
    """AIRM is affine invariant, so a shared rescaling changes nothing."""
    rng = np.random.default_rng(5)
    a, b = random_spd(rng, n=8), random_spd(rng, n=8)
    np.testing.assert_allclose(
        airm_distance(a, b), airm_distance(7.5 * a, 7.5 * b), rtol=1e-8
    )


def test_mahalanobis_matches_explicit_form():
    rng = np.random.default_rng(6)
    sigma = random_spd(rng)
    e = rng.normal(size=3)
    expected = e @ np.linalg.inv(sigma) @ e
    np.testing.assert_allclose(mahalanobis_sq(e, sigma), expected, rtol=1e-10)


def test_calibration_passes_for_correct_covariance():
    rng = np.random.default_rng(7)
    sigma = random_spd(rng)
    chol = np.linalg.cholesky(sigma)
    e = (chol @ rng.normal(size=(3, 20000))).T
    res = chi2_calibration_test(e, np.broadcast_to(sigma, (20000, 3, 3)))
    assert res.ks_pvalue > 0.01
    assert 0.95 < res.scale_factor < 1.05


def test_calibration_detects_overconfidence():
    """Errors twice as large as advertised must register as overconfident."""
    rng = np.random.default_rng(8)
    sigma = random_spd(rng)
    chol = np.linalg.cholesky(4.0 * sigma)
    e = (chol @ rng.normal(size=(3, 20000))).T
    res = chi2_calibration_test(e, np.broadcast_to(sigma, (20000, 3, 3)))
    assert res.scale_factor > 3.5
    assert res.ks_pvalue < 1e-6


def test_calibration_detects_wrong_shape_at_right_magnitude():
    """The isotropic surrogate has the right trace but the wrong shape.

    This is the test that distinguishes a magnitude failure from an
    anisotropy failure, and it is the one the published models never run.
    """
    rng = np.random.default_rng(9)
    sigma = np.diag([4.0, 1.0, 0.25])
    chol = np.linalg.cholesky(sigma)
    e = (chol @ rng.normal(size=(3, 40000))).T
    iso = np.broadcast_to(isotropize(sigma), (40000, 3, 3))
    res = chi2_calibration_test(e, iso)
    assert res.ks_pvalue < 1e-6


def test_wishart_null_is_positive_and_shrinks_with_more_frames():
    rng = np.random.default_rng(10)
    sigma = random_spd(rng)
    few = wishart_null_airm(sigma, n_eff=20, n_draws=200, rng=rng)
    many = wishart_null_airm(sigma, n_eff=2000, n_draws=200, rng=rng)
    assert np.all(few > 0)
    assert np.median(many) < np.median(few)


def test_wishart_null_rejects_degenerate_n_eff():
    rng = np.random.default_rng(11)
    with pytest.raises(ValueError):
        wishart_null_airm(random_spd(rng), n_eff=2)


# --- directional metrics -------------------------------------------------

from auf.audit.direction import (  # noqa: E402
    direction_wins,
    gaussian_loglik,
    optimal_shape_power,
    principal_axis_cos,
    shape_power,
    unit_trace,
)


def test_shape_power_endpoints():
    rng = np.random.default_rng(20)
    s = random_spd(rng, n=6)
    np.testing.assert_allclose(shape_power(s, 1.0), unit_trace(s), atol=1e-10)
    np.testing.assert_allclose(shape_power(s, 0.0), np.broadcast_to(np.eye(3) / 3, s.shape), atol=1e-10)


def test_optimal_shape_power_recovers_overstated_anisotropy():
    """A predictor with the right axes but squared elongation should get
    alpha ~ 0.5: half of its (log) anisotropy is real."""
    rng = np.random.default_rng(21)
    truth = random_spd(rng, n=50)
    pred = shape_power(truth, 2.0)
    alpha, _ = optimal_shape_power(pred, truth)
    assert abs(alpha - 0.5) < 0.06


def test_direction_wins_ignores_magnitude():
    rng = np.random.default_rng(22)
    truth = random_spd(rng, n=50) * np.array([1.0, 0.5, 0.1])[None, None, :]
    truth = truth @ np.swapaxes(truth, -1, -2)
    assert direction_wins(100.0 * truth, truth).all()


def test_principal_axis_cos_random_is_half():
    rng = np.random.default_rng(23)
    a, b = random_spd(rng, n=20000), random_spd(rng, n=20000)
    assert abs(principal_axis_cos(a, b).mean() - 0.5) < 0.01


def test_gaussian_loglik_is_maximised_by_true_covariance():
    rng = np.random.default_rng(24)
    s = random_spd(rng, n=10)
    best = gaussian_loglik(s, s)
    assert np.all(best > gaussian_loglik(isotropize(s), s))
    assert np.all(best > gaussian_loglik(1.3 * s, s))


# --- downstream fusion ----------------------------------------------------

from auf.audit.fusion import block_diag_3, fuse_distances, posterior_blocks  # noqa: E402


def _toy_structure(rng, n=12):
    return np.cumsum(rng.normal(scale=3.0, size=(n, 3)), axis=0)


def test_fusion_ignores_useless_data():
    rng = np.random.default_rng(40)
    mu = _toy_structure(rng)
    prior = block_diag_3(random_spd(rng, n=len(mu)) * 0.1)
    pairs = np.array([[0, 5], [2, 9]])
    d = np.linalg.norm(mu[pairs[:, 0]] - mu[pairs[:, 1]], axis=-1) + 3.0
    x, cov = fuse_distances(mu, prior, pairs, d, noise_sd=1e4)
    np.testing.assert_allclose(x, mu, atol=1e-5)
    np.testing.assert_allclose(cov, prior, atol=1e-6)


def test_fusion_matches_exact_distances_with_tiny_noise():
    rng = np.random.default_rng(41)
    mu = _toy_structure(rng)
    truth = mu + rng.normal(scale=0.3, size=mu.shape)
    pairs = np.array([[0, 5], [2, 9], [3, 11], [1, 7]])
    d = np.linalg.norm(truth[pairs[:, 0]] - truth[pairs[:, 1]], axis=-1)
    x, _ = fuse_distances(mu, block_diag_3(np.broadcast_to(0.1 * np.eye(3), (12, 3, 3))), pairs, d, 1e-3)
    fit = np.linalg.norm(x[pairs[:, 0]] - x[pairs[:, 1]], axis=-1)
    np.testing.assert_allclose(fit, d, atol=1e-3)


def test_fusion_batches_and_blocks():
    rng = np.random.default_rng(42)
    mu = _toy_structure(rng)
    prior = block_diag_3(random_spd(rng, n=len(mu)) * 0.1)
    pairs = np.array([[0, 5], [2, 9]])
    d = np.linalg.norm(mu[pairs[:, 0]] - mu[pairs[:, 1]], axis=-1)
    x, cov = fuse_distances(mu, prior, pairs, np.stack([d, d + 0.5, d - 0.5]), 0.5)
    assert x.shape == (3, 12, 3) and cov.shape == (3, 36, 36)
    blocks = posterior_blocks(cov)
    assert blocks.shape == (3, 12, 3, 3)
    # Unmeasured residues keep their prior block exactly.
    np.testing.assert_allclose(blocks[0, 4], prior[12:15, 12:15], atol=1e-10)


# --- subspace capture -----------------------------------------------------

from auf.audit.subspace import capture_ratio, project_out_rigid, rigid_basis  # noqa: E402


def test_rigid_basis_spans_rigid_motions():
    rng = np.random.default_rng(50)
    x = rng.normal(scale=5.0, size=(20, 3))
    q = rigid_basis(x)
    assert q.shape == (60, 6)
    w = np.array([0.01, -0.02, 0.03])
    moved = (x - x.mean(0)) + np.cross(w, x - x.mean(0)) + np.array([1.0, 2.0, 3.0])
    delta = (moved - (x - x.mean(0))).ravel()
    np.testing.assert_allclose(project_out_rigid(x, delta), 0.0, atol=1e-10)


def test_capture_ratio_random_is_one_and_aligned_is_large():
    rng = np.random.default_rng(51)
    n = 30
    x = rng.normal(scale=5.0, size=(n, 3))
    a = rng.normal(size=(3 * n, 3 * n))
    cov = a @ a.T
    rand = capture_ratio(cov, rng.normal(size=(4000, 3 * n)), ks=(10,), rigid_ref=x)[10]
    assert abs(rand.mean() - 1.0) < 0.05
    top = np.linalg.eigh(project_out_rigid(x, cov))[1][:, -1]
    assert capture_ratio(cov, top, ks=(10,), rigid_ref=x)[10] > 0.99 * (3 * n - 6) / 10
