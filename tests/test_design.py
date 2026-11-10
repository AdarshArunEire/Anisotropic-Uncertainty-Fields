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
