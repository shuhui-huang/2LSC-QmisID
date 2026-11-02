import awkward as ak
import numpy as np
import pytest
import uproot

from qmisid import expr
from qmisid.hist import Hist, write_root
from qmisid import likelihood as LH
from qmisid.closure import flip_weight, flip_weight_derivs


# --------------------------------------------------------------------------- expr
def _arrays():
    return {
        "a": np.array([1.0, 2.0, 3.0, -1.5], dtype=np.float32),
        "b": np.array([0, 1, 1, 0], dtype=np.int32),
        "c": np.array([2.0, 2.0, 5.0, 0.5]),
        "v": ak.Array([[1, 2], [3], [], [4, 5, 6]]),
    }


@pytest.mark.parametrize("text,expected", [
    ("a<2 && b==1", [0, 0, 0, 0]),
    ("a<2 || b==1", [1, 1, 1, 1]),
    ("!(a>1.5)", [1, 0, 0, 1]),
    ("!b", [1, 0, 0, 1]),
    ("a*0.001", [0.001, 0.002, 0.003, -0.0015]),
    ("fabs(a)", [1, 2, 3, 1.5]),
    ("TMath::Abs(a)+1", [2, 3, 4, 2.5]),
    ("a>0 ? c : -c", [2, 2, 5, -0.5]),
    ("pow(c,2)-1", [3, 3, 24, -0.75]),
    ("1+2*3==7", [1, 1, 1, 1]),
    ("a-1-1", [-1, 0, 1, -3.5]),
    ("(l>0)", None),
])
def test_expr_values(text, expected):
    arr = _arrays()
    if expected is None:
        with pytest.raises(expr.ExprError):
            expr.evaluate(text, arr)
        return
    np.testing.assert_allclose(expr.evaluate(text, arr), expected, rtol=1e-12, atol=1e-15)


def test_expr_precedence_and_or():
    arr = _arrays()
    # C++: && binds tighter than ||
    np.testing.assert_array_equal(expr.evaluate("b==1 || a>2 && c>4", arr), [0, 1, 1, 0])
    np.testing.assert_array_equal(expr.evaluate("(b==1 || a>2) && c>4", arr), [0, 0, 1, 0])


def test_expr_jagged_index_out_of_range_is_zero():
    arr = _arrays()
    # entry 2 has no element [0] -> whole expression 0 (ROOT skips the entry)
    np.testing.assert_array_equal(expr.evaluate("v[0]>0 || 1", arr), [1, 1, 0, 1])
    np.testing.assert_array_equal(expr.evaluate("v[1]", arr), [2, 0, 0, 5])


def test_expr_double_precision():
    # TTreeFormula computes in double: float32 branch * 0.001 must not be float32 arithmetic
    arr = {"m": np.array([91187.6], dtype=np.float32)}
    v = expr.evaluate("m*0.001", arr)
    assert v.dtype == np.float64
    assert v[0] == float(np.float32(91187.6)) * 0.001


def test_expr_variables():
    assert expr.variables("fabs(l1_eta)<2.47 && el_pt[0]*0.001>20 && TMath::Abs(x)") == {"l1_eta", "el_pt", "x"}


# --------------------------------------------------------------------------- hist
def test_hist_root_conventions(tmp_path):
    h = Hist([0, 1, 2], [0, 10], name="h")
    h.fill([-1, 0, 0.999, 1, 2, np.nan], [5, 5, 5, 5, 5, 5], weight=[1, 2, 3, 4, 5, 6])
    # underflow, [0,1), [1,2), overflow (2 is overflow in ROOT, NaN too)
    np.testing.assert_allclose(h.sumw[:, 1], [1, 5, 4, 11])
    np.testing.assert_allclose(h.sumw2[:, 1], [1, 13, 16, 61])
    f = tmp_path / "h.root"
    write_root(str(f), {"h": h})
    with uproot.open(f) as fh:
        r = fh["h"]
        assert r.classname == "TH2F"
        np.testing.assert_allclose(r.values(flow=True), h.sumw)
        np.testing.assert_allclose(r.variances(flow=True), h.sumw2)
    back = Hist.read(str(f), "h")
    np.testing.assert_allclose(back.sumw, h.sumw)


def test_hist_divide_like_root():
    n = Hist([0, 1, 2], name="n").fill([0.5, 0.5, 1.5], weight=[1, 2, 3])
    d = Hist([0, 1, 2], name="d").fill([0.5, 1.5, 1.5], weight=[4, 5, 0])
    r = n.divide(d)
    c1, c2, e1, e2 = 3.0, 4.0, 5.0, 16.0
    assert r.values()[0] == pytest.approx(c1 / c2)
    assert r.variances()[0] == pytest.approx((e1 * c2 ** 2 + e2 * c1 ** 2) / c2 ** 4)
    assert r.values()[1] == pytest.approx(3 / 5)


def test_merge_overflow():
    h = Hist([0, 1, 2]).fill([0.5, 1.5, 7, 9])
    h.merge_overflow()
    np.testing.assert_allclose(h.values(), [1, 3])
    assert h.sumw[-1] == 0


# --------------------------------------------------------------------------- likelihood
def test_find_bins_conventions():
    edges = [20, 50, 100]
    x = np.array([20.0, 35, 50, 100, 150, 10])
    np.testing.assert_array_equal(LH.find_bins(edges, x, "right_closed", True), [-1, 0, 0, 1, 1, -1])
    np.testing.assert_array_equal(LH.find_bins(edges, x, "left_closed", True), [0, 0, 1, 1, 1, -1])


def _asimov_cells(eps, n_eta, n_pt, n_per_cell=1e6, model="approx"):
    ss = np.zeros((n_eta, n_eta, n_pt, n_pt))
    os_ = np.zeros_like(ss)
    for e1 in range(n_eta):
        for e2 in range(e1, n_eta):
            for p1 in range(n_pt):
                for p2 in range(n_pt):
                    a, b = eps[e1, p1], eps[e2, p2]
                    mu = a + b - 2 * a * b if model == "exact" else a + b
                    ss[e1, e2, p1, p2] = n_per_cell * mu
                    os_[e1, e2, p1, p2] = n_per_cell * (1 - mu)
    return ss, os_


@pytest.mark.parametrize("model", ["approx", "exact"])
def test_fit_recovers_asimov(model):
    eps = np.array([[1e-4, 3e-4, 1e-3], [0, 0, 0], [5e-4, 2e-3, 2e-2]])
    ss, os_ = _asimov_cells(eps, 3, 3, model=model)
    res = LH.fit_rates(ss, os_, [0, 1, 2, 3], [0, 1, 2, 3], skip_eta_bins=[1], model=model)
    assert res.valid
    np.testing.assert_allclose(res.rates, eps, rtol=1e-3, atol=1e-12)
    assert np.all(res.errors[[0, 2]] > 0)
    assert np.all(res.errors[1] == 0)


def test_gradient_matches_numeric():
    rng = np.random.default_rng(1)
    ss = rng.poisson(5, (3, 3, 2, 2)).astype(float)
    os_ = rng.poisson(5000, (3, 3, 2, 2)).astype(float)
    for model in ("approx", "exact"):
        nll, grad = LH.make_nll(ss, os_, 3, 2, skip=[], model=model)
        x = rng.uniform(1e-4, 1e-2, 6)
        g = grad(x)
        h = 1e-9
        num = np.array([(nll(x + h * np.eye(6)[i]) - nll(x - h * np.eye(6)[i])) / (2 * h) for i in range(6)])
        np.testing.assert_allclose(g, num, rtol=1e-4, atol=1e-3)


def test_flip_weight_derivatives():
    e1, e2 = np.array([1e-3, 0.02]), np.array([5e-4, 0.1])
    d1, d2 = flip_weight_derivs(e1, e2)
    h = 1e-8
    np.testing.assert_allclose(d1, (flip_weight(e1 + h, e2) - flip_weight(e1 - h, e2)) / (2 * h), rtol=1e-5)
    np.testing.assert_allclose(d2, (flip_weight(e1, e2 + h) - flip_weight(e1, e2 - h)) / (2 * h), rtol=1e-5)
    # P/(1-P) with P = e1 + e2 - 2 e1 e2
    p = e1 + e2 - 2 * e1 * e2
    np.testing.assert_allclose(flip_weight(e1, e2), p / (1 - p))
