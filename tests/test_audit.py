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
