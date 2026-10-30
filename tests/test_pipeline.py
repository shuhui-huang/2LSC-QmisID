"""End-to-end tests on small toy ntuples with known (binned) charge-flip rates."""

import json
import os

import numpy as np
import pytest
import uproot

from qmisid import closure, config, expr, likelihood, systematics, toys, truth
from qmisid.toys import BINNED_RATES

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="module")
def setup(tmp_path_factory):
    d = tmp_path_factory.mktemp("toys")
    lists = toys.make_all(str(d / "nt"), n_data=400_000, n_z=300_000, n_tt=60_000, seed=11)
    out = d / "out"
    ov = {
        "output_dir": str(out),
        "samples": {k: {"path": str(d / "nt"), "filelist": v} for k, v in lists.items()},
        "truth": {"output_dir": str(out / "truth")},
        "likelihood": {"output_dir_data": str(out / "lh_data"), "output_dir_mc": str(out / "lh_mc"),
                       "reference": {"file": str(out / "truth" / "Zjets_Default_nJets_le1.root"),
                                     "hist": "Zjets_Default_nJets_le1_all_num"}},
        "closure": {"output_dir": str(out / "closure")},
        "systematics": {"output_dir": str(out / "syst"), "plots": False},
    }
    cfg = config.load_config(os.path.join(REPO, "configs", "default.yaml"), ov)
    cfg["samples"]["Data"]["is_data"] = True
    for k in ("Zjets", "ttbar"):
        cfg["samples"][k]["weight"] = "weight"
    return cfg, d


def test_truth_recovers_generated_rates(setup):
    cfg, _ = setup
    acc, outfile = truth.run(cfg, "Zjets")
    rate = acc.rate.values()
    # expected (binomial) uncertainty from the generated rate, robust for empty bins
    den, den_w2 = acc.den.values(), acc.den.variances()
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.sqrt(BINNED_RATES * den_w2) / den
    m = np.ones_like(rate, bool)
    m[1] = False
    pulls = (rate[m] - BINNED_RATES[m]) / err[m]
    assert np.all(np.abs(pulls) < 4.5), pulls
    with uproot.open(outfile) as f:
        keys = {k.split(";")[0] for k in f.keys()}
    # object names of ChargeFlipTruthv5New.C
    for k in ("Zjets_Default_nJets_le1_all_num", "Zjets_Default_nJets_le1_wrongTrack_num1", "m_histo_all_num",
              "Zjets_Default_nJets_le1_histo_RF_SS_num", "raw/all_den"):
        assert k in keys
    truth.run(cfg, "ttbar")


def test_likelihood_recovers_generated_rates(setup):
    cfg, _ = setup
    acc, res, outfile = likelihood.run(cfg, "Data")
    assert res.valid
    m = res.free.reshape(res.rates.shape)
    pulls = (res.rates[m] - BINNED_RATES[m]) / res.errors[m]
    assert np.all(np.abs(pulls) < 4.5), pulls
    assert np.mean(pulls ** 2) < 4
    s = json.load(open(outfile.replace(".root", ".json")))
    assert s["fit"]["status"] == 3
    with uproot.open(outfile) as f:
        assert "Data_Default_nJets_le1_all_num" in f and "SS" in f and "OS" in f
        np.testing.assert_allclose(f["Data_Default_nJets_le1_all_num"].values(), res.rates, rtol=1e-6)


def test_closure_matches_macro_event_loop(setup):
    """Columnar closure == literal per-event re-implementation of ClosureTestv5New.C."""
    cfg, d = setup
    lh_file = os.path.join(cfg["likelihood"]["output_dir_data"], "Data_Default_nJets_le1.root")
    sec = config.section(cfg, "closure", "Data")
    sec["variables"] = [{"name": "PtLep1", "expr": "l1_pt*0.001", "bins": [20, 50, 100, 200, 2600]}]
    sec["merge_overflow"] = "nominal"
    sec["subtract_background"] = False
    sec["float32_kinematics"] = True
    rates = closure.RateMap.from_file(lh_file)
    acc = closure.ClosureAccumulator(cfg, config.get_sample(cfg, "Data"), rates, sec=sec)
    fname = os.path.join(cfg["samples"]["Data"]["path"], "toy_2lSC_QMisID_data.root")
    with uproot.open(fname) as f:
        arr = f["qmisid_cr"].arrays(library="np", entry_stop=50_000)
    acc.process(arr, len(arr["l1_pt"]))
    res = acc.finalize()["PtLep1"]["hists"]

    # literal macro logic
    h = rates.h
    edges_x, edges_y = h.edges
    ss = np.zeros(6)
    osn = np.zeros(6)
    osu = np.zeros(6)
    vedges = np.array([20, 50, 100, 200, 2600])
    sel = expr.evaluate_mask(sec["selection"], arr)
    for i in range(len(arr["l1_pt"])):
        f = lambda x: float(np.float32(x))  # noqa: E731  (C++ float variables)
        mll = f(float(arr["m_l1l2"][i]) * 0.001)
        if mll < 81 or mll > 101 or not sel[i]:
            continue
        pt1, pt2 = f(float(arr["l1_pt"][i]) * 0.001), f(float(arr["l2_pt"][i]) * 0.001)
        eta1, eta2 = abs(float(arr["l1_eta"][i])), abs(float(arr["l2_eta"][i]))
        b = np.searchsorted(vedges, pt1, side="right")
        if arr["same_charge"][i]:
            ss[b] += 1
        elif arr["opposite_charge"][i]:
            def eps(eta, pt):
                ix = np.searchsorted(edges_x, eta, side="right")
                iy = np.searchsorted(edges_y, pt, side="right")
                return h.sumw[ix, iy], np.sqrt(h.sumw2[ix, iy])
            (e1, s1), (e2, s2) = eps(eta1, pt1), eps(eta2, pt2)
            W = lambda a, c: (a + c - 2 * a * c) / (1 - a - c + 2 * a * c)  # noqa: E731
            osn[b] += W(e1, e2)
            osu[b] += W(min(e1 + s1, 1), min(e2 + s2, 1))
    ss[-2] += ss[-1]
    osn[-2] += osn[-1]
    np.testing.assert_allclose(res["SS"].values(), ss[1:-1])
    np.testing.assert_allclose(res["OS"].values(), osn[1:-1], rtol=1e-10)
    np.testing.assert_allclose(res["OSUp"].values(flow=True)[1:-1], osu[1:-1], rtol=1e-10)


def test_closure_closes_on_toys(setup):
    cfg, _ = setup
    lh_file = os.path.join(cfg["likelihood"]["output_dir_data"], "Data_Default_nJets_le1.root")
    acc, res, yields, _ = closure.run(cfg, "Data", lh_file, plot=False)
    y = yields["PtLep1"]
    # SS observed vs OS-weighted: compatible within the SS statistical uncertainty (+ rate stat)
    tot = np.hypot(y["SS_stat"], y["rate_stat_cov"])
    assert abs(y["SS"] - y["OS_weighted"]) < 3 * tot
    # covariance-based rate uncertainty is smaller than the coherent +-1 sigma shift
    assert y["rate_stat_cov"] < y["rate_up"]


def test_systematics_runs(setup):
    cfg, _ = setup
    cfg = config.deep_merge(cfg, {"systematics": {"stat": {"n_bootstrap": 8}}})
    cfg["systematics"]["sources"]["mll_window"]["mode"] = "beyond_stat"
    runner, res = systematics.run(cfg)
    comps = res["components"]
    for k in ("mll_window", "bkg_norm", "lh_model", "non_closure", "method_bias", "process_dependence"):
        assert k in comps and comps[k].shape == (4, 4)
        assert np.all(comps[k][1] == 0)  # crack not fitted
    assert np.all(res["total"] >= res["stat"] - 1e-15)
    # ttbar toys are generated with 1.3x the Z rates -> process dependence ~ +30 %
    r = np.array(runner.result["details"]["process_dependence"]["relative"])
    assert 0.0 < r[0, 0] < 0.8
    base = os.path.join(cfg["systematics"]["output_dir"], "Data_Default_nJets_le1_syst")
    with uproot.open(base + ".root") as f:
        assert "Data_Default_nJets_le1_all_num" in f
        assert "Data_Default_nJets_le1_rates_bkg_norm_up" in f
        np.testing.assert_allclose(f["Data_Default_nJets_le1_all_num"].errors(), res["total"], rtol=1e-5)
    # the syst file is directly usable by the closure step (alternative maps)
    alt = {"bkg_norm_up": (base + ".root", "Data_Default_nJets_le1_rates_bkg_norm_up")}
    _, _, yields, _ = closure.run(cfg, "Data", base + ".root", "Data_Default_nJets_le1_all_num", alt=alt, plot=False,
                                  outdir=str(os.path.join(cfg["systematics"]["output_dir"], "apply")))
    assert yields["PtLep1"]["OS_bkg_norm_up"] >= yields["PtLep1"]["OS_weighted"] * 0.999


def test_cli_smoke(setup, tmp_path):
    from qmisid.cli import main

    cfg, _ = setup
    lh_file = os.path.join(cfg["likelihood"]["output_dir_data"], "Data_Default_nJets_le1.root")
    assert main(["plot", lh_file, "-o", str(tmp_path / "rates")]) == 0
    assert os.path.exists(tmp_path / "rates_1D.pdf")
