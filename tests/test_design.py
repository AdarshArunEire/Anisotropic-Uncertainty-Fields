import numpy as np

from auf.design import candidate_pairs, greedy_variance_reduction, plddt_pairs


def _x(n=25, seed=0):
    return np.cumsum(np.random.default_rng(seed).normal(size=(n, 3)) * 2.2, axis=0)


def test_greedy_first_pick_matches_brute_force():
    x = _x()
    rng = np.random.default_rng(1)
    g = rng.normal(size=(75, 75))
    prior = g @ g.T / 75
    cands = candidate_pairs(x)
    best = greedy_variance_reduction(x, prior, cands, 1, 1.0)[0]
    gains = []
    for a, b in cands:
        u = (x[a] - x[b]) / np.linalg.norm(x[a] - x[b])
        h = np.zeros(75)
        h[3 * a:3 * a + 3], h[3 * b:3 * b + 3] = u, -u
        ph = prior @ h
        gains.append(ph @ ph / (h @ ph + 1.0))
    assert tuple(best) == tuple(cands[int(np.argmax(gains))])


def test_greedy_picks_distinct_pairs():
    x = _x()
    cands = candidate_pairs(x)
    out = greedy_variance_reduction(x, np.eye(75), cands, 10, 1.0)
    assert len({tuple(p) for p in out}) == 10


def test_plddt_pairs_prefers_low_confidence_and_no_reuse():
    x = _x()
    pl = np.full(25, 95.0)
    pl[[3, 20]] = 40.0
    out = plddt_pairs(candidate_pairs(x), pl, 4)
    assert 3 in out[0] or 20 in out[0]
    assert len(set(out.ravel())) == out.size


def _helix(n=20):
    t = np.arange(n) * np.deg2rad(100)
    return np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)], 1)


def test_ca_secondary_structure_detects_helix_and_strand():
    from auf.design import ca_secondary_structure
    assert (ca_secondary_structure(_helix())[1:-1] == "H").all()
    strand = np.stack([3.3 * np.arange(12), np.tile([0.0, 0.9], 6), np.zeros(12)], 1)
    assert (ca_secondary_structure(strand)[1:-1] == "E").all()


def test_ensemble_mrmr_avoids_redundant_duplicate():
    from auf.design import ensemble_mrmr
    x = _x()
    cands = candidate_pairs(x)
    out = ensemble_mrmr(x, np.eye(75), cands, 5, 1.0, np.random.default_rng(0))
    assert len({tuple(p) for p in out}) == 5


def test_linear_posterior_trace_decreases_with_data():
    from auf.design import linear_posterior_trace
    x = _x()
    cands = candidate_pairs(x)
    t0 = linear_posterior_trace(x, np.eye(75), cands[:0], 1.0)
    t5 = linear_posterior_trace(x, np.eye(75), cands[:5], 1.0)
    t20 = linear_posterior_trace(x, np.eye(75), cands[:20], 1.0)
    assert t0 == 75 and t0 > t5 > t20 > 0
