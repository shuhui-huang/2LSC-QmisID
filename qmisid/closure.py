"""Step 3: closure test / application of the rates (port of ClosureTestv5New.C).

Opposite-sign events are weighted with the probability ratio that the pair is
reconstructed same-sign because of a charge flip,

    w = P / (1 - P),   P = e1 + e2 - 2 e1 e2,

with ``e1, e2`` the rates looked up at the two electrons' (|eta|, pT)
(``TH2::FindBin``; flow bins give 0). The weighted OS distributions are compared
to the observed SS ones. The same machinery applies the rates in any region
(``qmisid closure`` on a signal-region ntuple gives the QmisID prediction).

Uncertainty bands:
* ``up/down``: all rates moved coherently by +-1 sigma (clipped to [0, 1]), as
  in the macro;
* ``stat_cov``: linear propagation of the fit covariance, using per-parameter
  derivative histograms ``J_k = sum w dW/d eps_k`` (sigma^2 = J C J^T).
"""

from __future__ import annotations

import json
import logging
import os

import numpy as np

from . import config as C
from . import expr, io
from .hist import Hist, write_root

log = logging.getLogger(__name__)


def flip_weight(e1, e2):
    p = e1 + e2 - 2.0 * e1 * e2
    with np.errstate(divide="ignore", invalid="ignore"):
        return p / (1.0 - p)


def flip_weight_derivs(e1, e2):
    p = e1 + e2 - 2.0 * e1 * e2
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 1.0 / (1.0 - p) ** 2
    return d * (1.0 - 2.0 * e2), d * (1.0 - 2.0 * e1)


def _var_edges(v):
    if "bins" in v:
        return np.asarray(v["bins"], float)
    return np.linspace(float(v["low"]), float(v["high"]), int(v["nbins"]) + 1)


class RateMap:
    """A rate histogram with fast per-electron lookup."""

    def __init__(self, hist, cov=None):
        self.h = hist
        self.nx, self.ny = hist.nbins
        self.cov = cov

    @classmethod
    def from_file(cls, filename, hist_name=None, with_cov=True):
        import uproot

        if hist_name is None:
            hist_name = os.path.basename(filename).replace(".root", "") + "_all_num"
        h = Hist.read(filename, hist_name)
        cov = None
        if with_cov:
            with uproot.open(filename) as fh:
                key = hist_name.split("_all_num")[0] + "_covariance"
                if key in fh:
                    cov = fh[key].values()
        return cls(h, cov)

    def lookup(self, eta, pt):
        ix, iy = self.h.bin_indices(eta, pt)
        val = self.h.sumw[ix, iy]
        err = np.sqrt(self.h.sumw2[ix, iy])
        inside = (ix >= 1) & (ix <= self.nx) & (iy >= 1) & (iy <= self.ny)
        k = np.where(inside, (ix - 1) * self.ny + (iy - 1), -1)
        return val, err, k


class ClosureAccumulator:
    def __init__(self, cfg, sample, rates, alt_rates=None, sec=None, name=None):
        self.cfg = cfg
        self.sample = sample
        self.sec = sec if sec is not None else C.section(cfg, "closure", sample["name"])
        s = self.sec
        self.name = name or sample["name"]
        self.rates = rates
        self.alt = alt_rates or {}
        self.selection = s.get("selection") or "1"
        self.ss_expr = s.get("same_sign", "same_charge")
        self.os_expr = s.get("opposite_sign", "opposite_charge")
        self.mll_expr = s.get("mll", "m_l1l2*0.001")
        self.window = s.get("mll_window", [81, 101])
        self.sidebands = s.get("sidebands")  # [low, high] outer edges, or None
        self.subtract = bool(s.get("subtract_background", False)) and self.sidebands is not None
        if self.subtract:
            lo, hi = self.window
            slo, shi = self.sidebands
            self.scale = (hi - lo) / ((lo - slo) + (shi - hi))
        lep = s["leptons"]
        self.eta_expr, self.pt_expr = lep["eta"], lep["pt"]
        self.weight_expr = None if sample.get("is_data") else sample.get("weight")
        self.rate_clip = bool(s.get("clip_variations", True))
        self.f32 = bool(s.get("float32_kinematics", False))  # macro stores kinematics in floats
        self.variables = s["variables"]
        self.merge_overflow = s.get("merge_overflow", "all")  # all | nominal (macro) | none
        npar = rates.nx * rates.ny
        self.npar = npar
        self.h = {}
        self.jac = {}
        for v in self.variables:
            edges = _var_edges(v)
            t = [v.get("xtitle", v["name"])]
            mk = lambda suf: Hist(edges, name=f"{v['name']}_{suf}", title=f"{v['name']}_{suf}", axis_titles=t)  # noqa: E731
            d = {"SS": mk("SS"), "OS": mk("OS"), "OSUp": mk("OSUp"), "OSDown": mk("OSDown")}
            for a in self.alt:
                d[f"OS_{a}"] = mk(f"OS_{a}")
            self.h[v["name"]] = d
            self.jac[v["name"]] = np.zeros((npar, len(edges) + 1))
        self.n_ss = 0.0
        self.n_os_w = 0.0

    def expressions(self):
        out = [self.selection, self.ss_expr, self.os_expr, self.mll_expr, self.weight_expr,
               *self.eta_expr, *self.pt_expr]
        out += [v["expr"] for v in self.variables]
        return out

    def process(self, arrays, n, cache=None):
        cache = {} if cache is None else cache

        def ev(e):
            if e not in cache:
                cache[e] = expr.evaluate(e, arrays, n)
            return cache[e]

        def kin(e):
            v = ev(e)
            return v.astype(np.float32).astype(np.float64) if self.f32 else v

        w = ev(self.weight_expr) if self.weight_expr else np.ones(n)
        sel = ev(self.selection) != 0
        if self.window is not None:
            mll = kin(self.mll_expr)
            lo, hi = self.window
            peak = ~((mll < lo) | (mll > hi))
            if self.subtract:
                slo, shi = self.sidebands
                side = ~((mll < slo) | (mll > shi)) & ~peak
                w = np.where(peak, w, np.where(side, -self.scale * w, 0.0))
                sel &= peak | side
            else:
                sel &= peak
        is_ss = (ev(self.ss_expr) != 0) & sel
        is_os = (ev(self.os_expr) != 0) & ~(ev(self.ss_expr) != 0) & sel
        self.n_ss += float(w[is_ss].sum())

        eta = [kin(self.eta_expr[i]) for i in (0, 1)]
        pt = [kin(self.pt_expr[i]) for i in (0, 1)]
        o = np.nonzero(is_os)[0]
        e1, s1, k1 = self.rates.lookup(eta[0][o], pt[0][o])
        e2, s2, k2 = self.rates.lookup(eta[1][o], pt[1][o])
        wo = w[o]
        qw = flip_weight(e1, e2)
        clip = (lambda x: np.clip(x, 0.0, 1.0)) if self.rate_clip else (lambda x: x)  # noqa: E731
        qw_up = flip_weight(clip(e1 + s1), clip(e2 + s2))
        qw_dn = flip_weight(clip(e1 - s1), clip(e2 - s2))
        d1, d2 = flip_weight_derivs(e1, e2)
        alt_w = {}
        for a, rm in self.alt.items():
            a1, _, _ = rm.lookup(eta[0][o], pt[0][o])
            a2, _, _ = rm.lookup(eta[1][o], pt[1][o])
            alt_w[a] = flip_weight(a1, a2)
        self.n_os_w += float((wo * qw).sum())

        for v in self.variables:
            x = kin(v["expr"])
            d = self.h[v["name"]]
            d["SS"].fill(x[is_ss], weight=w[is_ss])
            xo = x[o]
            d["OS"].fill(xo, weight=wo * qw)
            d["OSUp"].fill(xo, weight=wo * qw_up)
            d["OSDown"].fill(xo, weight=wo * qw_dn)
            for a, aw in alt_w.items():
                d[f"OS_{a}"].fill(xo, weight=wo * aw)
            # derivative histograms for covariance propagation
            xb = d["OS"].bins_of(xo)[0]
            J = self.jac[v["name"]]
            for k, dd in ((k1, d1), (k2, d2)):
                ok = k >= 0
                if ok.any():
                    flat = np.ravel_multi_index((k[ok], xb[ok]), J.shape)
                    J += np.bincount(flat, weights=wo[ok] * dd[ok], minlength=J.size).reshape(J.shape)

    # ------------------------------------------------------------------------
    def finalize(self):
        res = {}
        cov = self.rates.cov
        for v in self.variables:
            d = self.h[v["name"]]
            J = self.jac[v["name"]]
            if self.merge_overflow in ("all", "nominal"):
                keys = ["SS", "OS"] if self.merge_overflow == "nominal" else list(d)
                for key in keys:
                    d[key].merge_overflow()
                if self.merge_overflow == "all":
                    J[:, -2] += J[:, -1]
                    J[:, -1] = 0
            stat = None
            if cov is not None and cov.shape[0] == J.shape[0]:
                var = np.einsum("ib,ij,jb->b", J, cov, J)
                stat = np.sqrt(np.clip(var, 0, None))
            res[v["name"]] = {"hists": d, "stat_cov": stat, "jacobian": J}
        return res

    def yields(self, res, var_name=None):
        v = var_name or self.variables[0]["name"]
        d = res[v]["hists"]
        nom = d["OS"].integral()
        out = {
            "SS": d["SS"].integral(), "SS_stat": d["SS"].integral_error(),
            "OS_weighted": nom, "OS_weighted_stat": d["OS"].integral_error(),
            "rate_up": d["OSUp"].integral() - nom, "rate_down": nom - d["OSDown"].integral(),
        }
        J = res[v]["jacobian"]
        if self.rates.cov is not None and self.rates.cov.shape[0] == J.shape[0]:
            jt = J[:, 1:-1].sum(axis=1)
            out["rate_stat_cov"] = float(np.sqrt(max(jt @ self.rates.cov @ jt, 0.0)))
        for a in self.alt:
            out[f"OS_{a}"] = d[f"OS_{a}"].integral()
        out["SS/OS"] = out["SS"] / nom if nom else float("nan")
        return out


def ratio_with_errors(ss, os_nom, os_up, os_down):
    """SS/OS ratio with the macro's errors: rel. SS stat (+) rel. rate band."""
    num, ne = ss.values(), ss.errors()
    den = os_nom.values()
    eu = os_up.values() - den
    ed = den - os_down.values()
    with np.errstate(divide="ignore", invalid="ignore"):
        r = num / den
        up = np.abs(r) * np.sqrt((ne / num) ** 2 + (eu / den) ** 2)
        dn = np.abs(r) * np.sqrt((ne / num) ** 2 + (ed / den) ** 2)
    return r, dn, up


def output_dir(cfg, sample_name):
    sec = C.section(cfg, "closure", sample_name)
    return sec.get("output_dir") or C.output_path(cfg, "Closure")


def run(cfg, sample_name, rates_file, rates_hist=None, max_events=None, outdir=None, alt=None, plot=True, label=None):
    sample = C.get_sample(cfg, sample_name)
    sec = C.section(cfg, "closure", sample_name)
    rates = RateMap.from_file(rates_file, rates_hist)
    alt_maps = {k: RateMap.from_file(f, h, with_cov=False) for k, (f, h) in (alt or {}).items()}
    acc = ClosureAccumulator(cfg, sample, rates, alt_maps, sec=sec)
    files = io.resolve_inputs(sample)
    tree = sec.get("tree", cfg.get("tree"))
    log.info("Closure/apply for %s with rates %s:%s", sample_name, rates_file, rates.h.name)
    log.info("Selection: %s | Mll window %s | bkg subtraction %s", acc.selection, acc.window, acc.subtract)
    for arrays, n in io.iterate(files, tree, acc.expressions(), step_size=sec.get("step_size", "50 MB"), max_events=max_events):
        acc.process(arrays, n)
    res = acc.finalize()
    outdir = outdir or output_dir(cfg, sample_name)
    os.makedirs(outdir, exist_ok=True)
    objs = {}
    for vname, r in res.items():
        d = r["hists"]
        for key, h in d.items():
            objs[f"h_{sample_name}_{vname}_{key}"] = h.clone(f"h_{sample_name}_{vname}_{key}")
        ratio = d["SS"].divide(d["OS"], name=f"Ratio_{vname}_OS_Standard")
        rr, rdn, rup = ratio_with_errors(d["SS"], d["OS"], d["OSUp"], d["OSDown"])
        ratio.set_values(np.nan_to_num(rr), np.nan_to_num(0.5 * (rdn + rup)))
        objs[f"Ratio_{vname}_OS_Standard"] = ratio
        if r["stat_cov"] is not None:
            hs = d["OS"].clone(f"h_{sample_name}_{vname}_OS_statcov")
            hs.set_values(d["OS"].values(), r["stat_cov"][1:-1])
            objs[hs.name] = hs
    outfile = os.path.join(outdir, f"closure_{sample_name}.root")
    write_root(outfile, objs, kind="D")
    yields = {vname: acc.yields(res, vname) for vname in res}
    with open(os.path.join(outdir, f"closure_{sample_name}_yields.json"), "w") as fh:
        json.dump({"rates": f"{rates_file}:{rates.h.name}", "selection": acc.selection,
                   "mll_window": acc.window, "yields": yields}, fh, indent=1)
    y0 = yields[acc.variables[0]["name"]]
    log.info("SS events: %.6g ; OS weighted: %.6g +%.4g/-%.4g (coherent +-1 sigma)%s ; SS/OS = %.4f",
             y0["SS"], y0["OS_weighted"], y0["rate_up"], y0["rate_down"],
             f" ; +-{y0['rate_stat_cov']:.4g} (covariance)" if "rate_stat_cov" in y0 else "", y0["SS/OS"])
    if plot:
        from . import plotting

        for v in acc.variables:
            r = res[v["name"]]
            plotting.closure_plot(r, v, os.path.join(outdir, f"h_{sample_name}_{v['name']}"),
                                  label=label or sec.get("label", sample_name), lumi=sec.get("lumi_label"),
                                  logy=v.get("logy", False))
    log.info("Wrote %s", outfile)
    return acc, res, yields, outfile
