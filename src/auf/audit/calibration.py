"""Calibration tests for predicted covariances.

The published covariance models are scored with correlation and distribution
distances. None of them ask whether a predicted Sigma is the *right size and
shape* for the errors actually observed. These are the missing tests.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

__all__ = ["mahalanobis_sq", "chi2_calibration_test", "wishart_null_airm",
           "CalibrationResult"]


def mahalanobis_sq(residual: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """``e^T Sigma^-1 e`` for batched residuals and covariances.

    Parameters
    ----------
    residual : (..., d) observed error vectors.
    sigma : (..., d, d) predicted covariances.
    """
    residual = np.asarray(residual, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    solved = np.linalg.solve(sigma, residual[..., None])[..., 0]
    return np.einsum("...i,...i->...", residual, solved)


@dataclass
class CalibrationResult:
    """Outcome of a chi-squared calibration test."""

    dim: int
    n: int
    mean_mahalanobis: float
    expected_mean: float
    ks_statistic: float
    ks_pvalue: float

    @property
    def scale_factor(self) -> float:
        """Multiplicative error in predicted covariance magnitude.

        1.0 is correct; >1 means the model is overconfident.
        """
        return self.mean_mahalanobis / self.expected_mean

    def __str__(self) -> str:
        return (
            f"n={self.n} d={self.dim} "
            f"mean_mahalanobis={self.mean_mahalanobis:.3f} "
            f"(expected {self.expected_mean:.1f}, "
            f"scale={self.scale_factor:.3f}) "
            f"KS={self.ks_statistic:.4f} p={self.ks_pvalue:.3g}"
        )


def chi2_calibration_test(
    residual: np.ndarray, sigma: np.ndarray
) -> CalibrationResult:
    """Test whether standardised residuals follow chi-squared with d dof.

    If the predicted covariances are calibrated, ``e^T Sigma^-1 e`` is
    chi2_d distributed. A scale factor above 1 means overconfidence. A
    significant KS statistic with a scale factor near 1 means the magnitude is
    right but the *shape* is wrong, which is the anisotropy failing rather
    than the magnitude.
    """
    residual = np.asarray(residual, dtype=float).reshape(-1, np.shape(residual)[-1])
    sigma = np.asarray(sigma, dtype=float).reshape(-1, *np.shape(sigma)[-2:])
    d = residual.shape[-1]
    m2 = mahalanobis_sq(residual, sigma)
    ks = stats.kstest(m2, "chi2", args=(d,))
    return CalibrationResult(
        dim=d,
        n=int(m2.size),
        mean_mahalanobis=float(m2.mean()),
        expected_mean=float(d),
        ks_statistic=float(ks.statistic),
        ks_pvalue=float(ks.pvalue),
    )


def wishart_null_airm(
    sigma: np.ndarray,
    n_eff: int,
    n_draws: int = 200,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Sampling-noise null for covariance distances.

    A covariance estimated from a finite trajectory of ``n_eff`` effectively
    independent frames differs from the truth even when the truth is known
    exactly. This draws Wishart resamples of ``sigma`` and returns their AIRM
    distance to it, giving the distance scale that means *nothing*.

    Any claimed improvement in covariance prediction has to clear this.
    """
    from auf.audit.spd import airm_distance

    rng = np.random.default_rng() if rng is None else rng
    sigma = np.asarray(sigma, dtype=float)
    d = sigma.shape[-1]
    if n_eff <= d:
        raise ValueError("n_eff must exceed the dimension for a proper Wishart")
    draws = stats.wishart.rvs(
        df=n_eff, scale=sigma / n_eff, size=n_draws, random_state=rng
    )
    draws = np.asarray(draws, dtype=float).reshape(n_draws, d, d)
    return airm_distance(np.broadcast_to(sigma, draws.shape), draws)
