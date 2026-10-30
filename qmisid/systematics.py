"""Step 4: systematic uncertainties on the likelihood charge-flip rates.

Strategy (see docs/systematics.md for the full discussion). Every source gives,
per (|eta|, pT) bin, an absolute shift ``delta_s`` of the nominal data rate:

* ``stat``      - statistical uncertainty of the nominal fit (HESSE covariance;
                  optional Poisson-bootstrap cross-check that includes the
                  sideband-subtraction weights);
* ``refit``     - the data fit is repeated with modified settings (Z-mass
                  window / sidebands, background normalisation, likelihood
                  model, any selection change ...); delta = eps_var - eps_nom;
* ``closure``   - data non-closure in the fit region: for each rate bin k the
                  observed SS yield of pairs with an electron in bin k is
                  compared with the OS-weighted prediction (beyond stat.);
* ``mc_method`` - method bias: likelihood vs truth-matched rates in the same
                  MC sample (Z+jets), applied as a relative uncertainty;
* ``process``   - extrapolation from Z->ee to the signal-region composition:
                  truth rates in ttbar vs Z+jets MC, as a relative uncertainty.

Two-sided variations are symmetrised (``envelope`` = max |shift|, or
``half_diff``); sources are added in quadrature. All data variations are
filled in a single pass over the ntuples (accumulators with identical
event-level settings are shared), MC jobs in one pass per sample.
"""

from __future__ import annotations

import copy
import json
import logging
import os

import numpy as np

from . import config as C
from . import io
from . import likelihood as LH
from . import truth as TR
from .closure import flip_weight
from .hist import Hist, write_root

log = logging.getLogger(__name__)

_ACC_KEYS = ("selection", "same_sign", "opposite_sign", "mll", "mll_window", "bin_convention", "binning", "leptons")


def _acc_key(sec):
    return json.dumps({k: sec.get(k) for k in _ACC_KEYS}, sort_keys=True, default=str)


def _scale(sec):
    low, low_z, up_z, up = [float(v) for v in sec.get("mll_window", [61, 81, 101, 121])]
    s = sec.get("sideband_scale", "auto")
    return (up_z - low_z) / ((low_z - low) + (up - up_z)) if s in (None, "auto") else float(s)


def _symmetrize(deltas, how):
    arr = np.asarray(list(deltas.values()))
    if len(arr) == 1 or how == "envelope":
        return np.max(np.abs(arr), axis=0)
    if how == "half_diff":
        return 0.5 * np.abs(arr[0] - arr[-1])
    if how == "rms":
        return np.sqrt(np.mean(arr ** 2, axis=0))
    raise ValueError(f"Unknown symmetrisation {how!r}")


def _merge_granularity(num, den, how):
    """Sum numerator/denominator over eta and/or pT before taking the ratio."""
    if how == "bin":
        return num, den
    if how == "pt":  # one value per pT bin
        return np.broadcast_to(num.sum(0, keepdims=True), num.shape), np.broadcast_to(den.sum(0, keepdims=True), den.shape)
    if how == "eta":
        return np.broadcast_to(num.sum(1, keepdims=True), num.shape), np.broadcast_to(den.sum(1, keepdims=True), den.shape)
    if how == "inclusive":
        return np.full(num.shape, num.sum()), np.full(den.shape, den.sum())
    raise ValueError(f"Unknown granularity {how!r}")


# ---------------------------------------------------------------------------
# non-closure in the fit region, from the pair-cell counts
# ---------------------------------------------------------------------------
def cell_closure(acc, rates, subtract=None, scale=None, skip=()):
    """Per rate-bin closure: SS observed vs OS x w(eps_a, eps_b) for all pairs
    with at least one electron in the bin. Returns dict of arrays (n_eta, n_pt)."""
    ss, os_ = acc.pairs(subtract, scale)
    ss_w2, os_w2 = acc.pairs_w2(subtract, scale)
    ne, npt = acc.n_eta, acc.n_pt
    E1, E2, P1, P2 = np.meshgrid(np.arange(ne), np.arange(ne), np.arange(npt), np.arange(npt), indexing="ij")
    eps = rates.ravel()
    a = E1 * npt + P1
    b = E2 * npt + P2
    w = flip_weight(eps[a], eps[b])
    pred = os_ * w
    pred_w2 = os_w2 * w * w
    obs = np.zeros(ne * npt)
    exp = np.zeros(ne * npt)
    obs_v = np.zeros(ne * npt)
    exp_v = np.zeros(ne * npt)
    for k in range(ne * npt):
        m = (a == k) | (b == k)
        obs[k] = ss[m].sum()
        exp[k] = pred[m].sum()
        obs_v[k] = ss_w2[m].sum()
        exp_v[k] = pred_w2[m].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        r = obs / exp - 1.0
        sig = np.sqrt(obs_v / exp ** 2 + obs ** 2 * exp_v / exp ** 4)
    shape = (ne, npt)
    # goodness of fit over all cells
    keep = (E2 >= E1)
    for s in skip:
        keep &= (E1 != s) & (E2 != s)
    var = ss_w2 + pred_w2
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_cells = np.where(keep & (var > 0), (ss - pred) ** 2 / var, 0.0)
    return {"obs": obs.reshape(shape), "pred": exp.reshape(shape), "nonclosure": r.reshape(shape),
            "nonclosure_stat": sig.reshape(shape), "chi2": float(chi2_cells.sum()),
            "ncells": int((keep & (var > 0)).sum())}


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
class SystematicsRunner:
    def __init__(self, cfg, sample_name=None, max_events=None):
        self.cfg = cfg
        self.ss = copy.deepcopy(cfg.get("systematics", {}))
        self.sample_name = sample_name or self.ss.get("sample", "Data")
        self.sample = C.get_sample(cfg, self.sample_name)
        self.nom_sec = C.section(cfg, "likelihood", self.sample_name)
        self.max_events = max_events
        self.name = C.run_name(self.sample_name, cfg)
        self.outdir = self.ss.get("output_dir") or C.output_path(cfg, "Systematics")
        self.sources = self.ss.get("sources", {}) or {}
        self.skip = self.nom_sec.get("skip_eta_bins", [])
        self.eta_edges = np.asarray(self.nom_sec["binning"]["eta"], float)
        self.pt_edges = np.asarray(self.nom_sec["binning"]["pt"], float)

    # -- data pass ---------------------------------------------------------------
    def _data_pass(self):
        stat = self.ss.get("stat", {}) or {}
        nboot = int(stat.get("n_bootstrap", 0))
        accs = {}
        nominal = LH.LikelihoodAccumulator(self.cfg, self.sample, name=self.name, sec=self.nom_sec, n_bootstrap=nboot)
        accs[_acc_key(self.nom_sec)] = nominal
        self.var_secs = {}
        for src, spec in self.sources.items():
            if spec.get("type") != "refit":
                continue
            for vname, ov in (spec.get("variations") or {}).items():
                sec = C.deep_merge(self.nom_sec, ov or {})
                self.var_secs[(src, vname)] = sec
                key = _acc_key(sec)
                if key not in accs:
                    accs[key] = LH.LikelihoodAccumulator(self.cfg, self.sample, name=f"{self.name}_{src}_{vname}", sec=sec, n_bootstrap=nboot)
        self.accs = accs
        self.nominal_acc = nominal
        exprs = []
        for a in accs.values():
            exprs += a.expressions()
        files = io.resolve_inputs(self.sample)
        tree = self.nom_sec.get("tree", self.cfg.get("tree"))
        log.info("Data pass over %d files with %d distinct accumulators, %d bootstrap replicas", len(files), len(accs), nboot)
        rng = np.random.default_rng(int(stat.get("seed", 12345)))
        for arrays, n in io.iterate(files, tree, exprs, step_size=self.nom_sec.get("step_size", "50 MB"), max_events=self.max_events):
            cache = {}
            for a in accs.values():
                a.process(arrays, n, cache, external_boot=True)
            # the same Poisson weights for every accumulator -> correlated replicas
            for r in range(nboot):
                pw = rng.poisson(1.0, size=n).astype(np.float64)
                for a in accs.values():
                    a.fill_boot(r, pw)

    def _fit(self, sec, start, boot=False):
        acc = self.accs[_acc_key(sec)]
        sub, scale = bool(sec.get("subtract_background", True)), _scale(sec)
        ss, os_ = acc.pairs(sub, scale)
        res = LH.fit_rates(ss, os_, self.eta_edges, self.pt_edges, self.skip, start, sec.get("minuit", {}), sec.get("model", "approx"))
        b = LH.bootstrap_rates(acc, res, sec, sub, scale) if (boot and acc.n_boot) else None
        return res, b

    # -- MC passes -----------------------------------------------------------------
    def _truth_cfg(self, extra_cut):
        tsec = copy.deepcopy(self.cfg.get("truth", {}))
        tsec["binning"] = copy.deepcopy(self.nom_sec["binning"])
        if extra_cut:
            ec = tsec.get("extra_cut")
            tsec["extra_cut"] = f"({ec}) && ({extra_cut})" if ec else extra_cut
        return C.deep_merge(self.cfg, {"truth": tsec})

    def _mc_pass(self, sample_name, truth_cuts=(None,), want_lh=False):
        """One pass over an MC sample: truth rates for each requested extra
        cut, plus (optionally) the likelihood pair counts."""
        sample = C.get_sample(self.cfg, sample_name)
        out = {"truth": {}}
        accs = []
        for cut in truth_cuts:
            t = TR.TruthAccumulator(self._truth_cfg(cut), sample, name=f"{sample_name}_truth")
            out["truth"][cut] = t
            accs.append(t)
        if want_lh:
            lsec = C.section(self.cfg, "likelihood", sample_name)
            lh = LH.LikelihoodAccumulator(self.cfg, sample, name=f"{sample_name}_lh", sec=lsec)
            out["lh"] = (lh, lsec)
            accs.append(lh)
        exprs = []
        for a in accs:
            exprs += a.expressions()
        files = io.resolve_inputs(sample)
        tree = self.nom_sec.get("tree", self.cfg.get("tree"))
        log.info("MC pass over %s (%d files)", sample_name, len(files))
        for arrays, n in io.iterate(files, tree, exprs, max_events=self.max_events):
            cache = {}
            for a in accs:
                if isinstance(a, LH.LikelihoodAccumulator):
                    a.process(arrays, n, cache)
                else:
                    a.process(arrays, n)
        for t in out["truth"].values():
            t.finalize()
        return out

    @staticmethod
    def _relative(r, rs, mode, cap=None):
        """Relative uncertainty from a measured relative difference ``r`` with
        statistical uncertainty ``rs``: ``abs`` |r|; ``beyond_stat``
        sqrt(max(r^2 - rs^2, 0)); ``max_stat`` max(|r|, rs)."""
        r = np.asarray(r, float)
        rs = np.asarray(rs, float)
        bad = ~np.isfinite(r) | ~np.isfinite(rs)
        if bad.any():
            log.warning("Undefined relative difference in %d bin(s) (empty reference?): set to 0", int(bad.sum()))
            r = np.where(bad, 0.0, r)
            rs = np.where(bad, 0.0, rs)
        if mode == "abs":
            rel = np.abs(r)
        elif mode == "beyond_stat":
            rel = np.sqrt(np.clip(r ** 2 - rs ** 2, 0, None))
        elif mode == "max_stat":
            rel = np.maximum(np.abs(r), rs)
        else:
            raise ValueError(f"Unknown mode {mode!r}")
        if cap is not None:
            rel = np.minimum(rel, float(cap))
        return rel

    # -- main --------------------------------------------------------------------------
    def run(self):
        self._data_pass()
        nom_acc = self.nominal_acc
        ref = None
        if self.nom_sec.get("reference"):
            try:
                ref = LH.read_reference(self.nom_sec["reference"])
            except (FileNotFoundError, KeyError, OSError) as exc:
                log.warning("Reference not readable (%s): default starting values", exc)
        nominal = LH.fit_rates(*nom_acc.pairs(), self.eta_edges, self.pt_edges, self.skip, ref,
                               self.nom_sec.get("minuit", {}), self.nom_sec.get("model", "approx"))
        if not nominal.valid:
            log.error("Nominal fit did not converge!\n%s", nominal.message)
        eps = nominal.rates
        log.info("Nominal fit: fval=%.4f status=%d\n%s", nominal.fval, nominal.status,
                 LH.format_rates(nominal, self.eta_edges, self.pt_edges))
        free2d = nominal.free.reshape(eps.shape)

        components = {}      # name -> delta (n_eta, n_pt), absolute
        variations = {}      # name -> {vname: varied rates}
        details = {}

        # statistical
        stat_cfg = self.ss.get("stat", {}) or {}
        stat_hesse = nominal.errors.copy()
        boot = None
        if nom_acc.n_boot:
            boot = LH.bootstrap_rates(nom_acc, nominal, self.nom_sec)
            stat_boot = np.nanstd(boot, axis=0, ddof=1)
            details["stat_bootstrap"] = stat_boot.tolist()
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(free2d, stat_boot / stat_hesse, np.nan)
            log.info("Bootstrap/Hesse stat ratio per bin:\n%s", np.array2string(ratio, precision=3))
        stat = stat_boot if (stat_cfg.get("method") == "bootstrap" and boot is not None) else stat_hesse
        details["stat_hesse"] = stat_hesse.tolist()

        # MC jobs: one pass per (sample, truth cut), with the likelihood if any source needs it
        mc_cuts, mc_lh = {}, {}
        for spec in self.sources.values():
            if spec.get("enabled", True) is False:
                continue
            cut = spec.get("truth_extra_cut")
            if spec.get("type") == "mc_method":
                smp = spec.get("mc_sample", "Zjets")
                mc_cuts.setdefault(smp, []).append(cut)
                mc_lh[smp] = True
            elif spec.get("type") == "process":
                for smp in (spec.get("reference_sample", "Zjets"), spec.get("target_sample", "ttbar")):
                    mc_cuts.setdefault(smp, []).append(cut)
        mc_cache = {}

        def mc(sample_name, truth_extra_cut=None):
            if sample_name not in mc_cache:
                cuts = list(dict.fromkeys(mc_cuts.get(sample_name, [truth_extra_cut])))
                mc_cache[sample_name] = self._mc_pass(sample_name, cuts, mc_lh.get(sample_name, False))
            job = mc_cache[sample_name]
            return {"truth": job["truth"][truth_extra_cut], "lh": job.get("lh")}

        for src, spec in self.sources.items():
            typ = spec.get("type")
            if spec.get("enabled", True) is False:
                continue
            if typ == "refit":
                vr, deltas, dstat = {}, {}, {}
                mode = spec.get("mode", "abs")
                use_boot = boot is not None and mode != "abs"
                if mode != "abs" and boot is None:
                    log.warning("%s: mode %s needs systematics.stat.n_bootstrap > 0; using 'abs'", src, mode)
                    mode = "abs"
                for vname in (spec.get("variations") or {}):
                    sec = self.var_secs[(src, vname)]
                    res, vboot = self._fit(sec, eps, boot=use_boot)
                    if not res.valid:
                        log.warning("Variation %s/%s: fit not valid", src, vname)
                    vr[vname] = res.rates
                    d = res.rates - eps
                    if vboot is not None:
                        # stat. uncertainty of the shift from correlated bootstrap replicas
                        sd = np.nanstd(vboot - boot, axis=0, ddof=1)
                        dstat[vname] = sd
                        # (the formulas of _relative are homogeneous: use absolute shifts)
                        d = np.sign(d) * self._relative(np.where(free2d, d, 0), np.where(free2d, sd, 0), mode)
                    deltas[vname] = d
                    log.info("Variation %s/%s: fval=%.4f status=%d", src, vname, res.fval, res.status)
                variations[src] = vr
                details[src] = {"shift": {k: (vr[k] - eps).tolist() for k in vr},
                                "shift_stat": {k: v.tolist() for k, v in dstat.items()}, "mode": mode}
                components[src] = _symmetrize(deltas, spec.get("symmetrize", "envelope"))
            elif typ == "closure":
                cc = cell_closure(nom_acc, eps, skip=self.skip)
                r = np.where(free2d, cc["nonclosure"], 0)
                s = np.where(free2d, cc["nonclosure_stat"], 0)
                rel = self._relative(r, s, spec.get("mode", "beyond_stat"), spec.get("cap"))
                components[src] = rel * eps
                details[src] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in cc.items()}
                log.info("Data closure chi2/ncells = %.1f / %d", cc["chi2"], cc["ncells"])
            elif typ == "mc_method":
                smp = spec.get("mc_sample", "Zjets")
                job = mc(smp, spec.get("truth_extra_cut"))
                lh_acc, lsec = job["lh"]
                lres = LH.fit_rates(*lh_acc.pairs(), self.eta_edges, self.pt_edges, self.skip, eps,
                                    lsec.get("minuit", {}), lsec.get("model", "approx"))
                t = job["truth"]
                num, den = _merge_granularity(t.allnum.values(), t.den.values(), "bin")
                with np.errstate(divide="ignore", invalid="ignore"):
                    tr = num / den
                    r = lres.rates / tr - 1.0
                    # stat of the ratio (truth binomial-ish + LH Hesse)
                    tr_err = t.rate.errors()
                    rs = np.abs(lres.rates / tr) * np.sqrt((lres.errors / lres.rates) ** 2 + (tr_err / tr) ** 2)
                r, rs = np.where(free2d, r, 0), np.where(free2d, rs, 0)
                rel = self._relative(r, rs, spec.get("mode", "beyond_stat"), spec.get("cap"))
                components[src] = rel * eps
                variations[src] = {"mc_likelihood": lres.rates, "mc_truth": np.nan_to_num(tr)}
                details[src] = {"relative": np.nan_to_num(r).tolist(), "relative_stat": np.nan_to_num(rs).tolist()}
            elif typ == "process":
                # Truth flips observed in the target sample (e.g. ttbar) vs the flips
                # predicted for the same electrons with the reference (Z) truth rates:
                # r = sum(N_flip^tgt) / sum(N_e^tgt * eps_ref) - 1, summed over the
                # merged bins -> independent of the kinematic composition.
                refj = mc(spec.get("reference_sample", "Zjets"), spec.get("truth_extra_cut"))["truth"]
                tgtj = mc(spec.get("target_sample", "ttbar"), spec.get("truth_extra_cut"))["truth"]
                g = spec.get("granularity", "bin")
                eps_ref = refj.rate.values()
                eps_ref_v = refj.rate.variances()
                n_t, d_t = tgtj.allnum.values(), tgtj.den.values()
                n_t_v = tgtj.allnum.variances()
                pred = d_t * eps_ref
                pred_v = d_t ** 2 * eps_ref_v
                obs_g, pred_g = _merge_granularity(n_t, pred, g)
                obs_gv, pred_gv = _merge_granularity(n_t_v, pred_v, g)
                with np.errstate(divide="ignore", invalid="ignore"):
                    r = np.where(free2d, obs_g / pred_g - 1.0, 0)
                    rs = np.where(free2d, np.sqrt(obs_gv / pred_g ** 2 + obs_g ** 2 * pred_gv / pred_g ** 4), 0)
                components[src] = self._relative(r, rs, spec.get("mode", "abs"), spec.get("cap")) * eps
                variations[src] = {"reference_truth": np.nan_to_num(eps_ref), "target_truth": np.nan_to_num(tgtj.rate.values())}
                details[src] = {"relative": np.nan_to_num(r).tolist(), "relative_stat": np.nan_to_num(rs).tolist(), "granularity": g}
            elif typ == "external":
                # relative uncertainty map given in a ROOT file (hist content = relative unc.)
                h = Hist.read(spec["file"], spec["hist"])
                components[src] = np.abs(h.values()) * eps
            else:
                raise ValueError(f"Unknown systematic type {typ!r} for source {src}")
            components[src] = np.where(free2d, np.nan_to_num(components[src]), 0.0)

        syst = np.sqrt(sum(d ** 2 for d in components.values())) if components else np.zeros_like(eps)
        total = np.sqrt(stat ** 2 + syst ** 2)
        self.result = dict(nominal=nominal, stat=stat, components=components, syst=syst, total=total,
                           variations=variations, details=details, boot=boot)
        self._write()
        return self.result

    # -- output ---------------------------------------------------------------------
    def _write(self):
        R = self.result
        os.makedirs(self.outdir, exist_ok=True)
        base = os.path.join(self.outdir, f"{self.name}_syst")
        eta, pt = self.eta_edges, self.pt_edges
        eps = R["nominal"].rates
        mk = LH.rates_hist
        objs = {
            f"{self.name}_all_num": mk(f"{self.name}_all_num", eta, pt, eps, R["total"]),
            f"{self.name}_all_num_stat": mk(f"{self.name}_all_num_stat", eta, pt, eps, R["stat"]),
            f"{self.name}_all_num_syst": mk(f"{self.name}_all_num_syst", eta, pt, eps, R["syst"]),
        }
        for src, d in R["components"].items():
            objs[f"{self.name}_syst_{src}"] = mk(f"{self.name}_syst_{src}", eta, pt, eps, d)
            up = np.clip(eps + d, 0, 1)
            dn = np.clip(eps - d, 0, 1)
            objs[f"{self.name}_rates_{src}_up"] = mk(f"{self.name}_rates_{src}_up", eta, pt, up, np.zeros_like(eps))
            objs[f"{self.name}_rates_{src}_down"] = mk(f"{self.name}_rates_{src}_down", eta, pt, dn, np.zeros_like(eps))
        for src, vr in R["variations"].items():
            for vname, v in vr.items():
                key = f"{self.name}_var_{src}_{vname}"
                objs[key] = mk(key, eta, pt, v, np.zeros_like(eps))
        # statistical uncertainty as uncorrelated eigen-variations of the fit covariance
        # (one nuisance parameter each; scaled so that their sum in quadrature gives the
        # quoted per-bin stat. uncertainty, i.e. HESSE or bootstrap)
        for k, vec in enumerate(self.stat_eigen()):
            up = np.clip(eps + vec, 0, 1)
            dn = np.clip(eps - vec, 0, 1)
            objs[f"{self.name}_rates_stat_eig{k}_up"] = mk(f"{self.name}_rates_stat_eig{k}_up", eta, pt, up, np.zeros_like(eps))
            objs[f"{self.name}_rates_stat_eig{k}_down"] = mk(f"{self.name}_rates_stat_eig{k}_down", eta, pt, dn, np.zeros_like(eps))
        write_root(base + ".root", objs)
        import uproot

        cov = R["nominal"].cov
        hc = Hist(np.arange(cov.shape[0] + 1), np.arange(cov.shape[0] + 1), name=f"{self.name}_covariance")
        hc.set_values(cov)
        with uproot.update(base + ".root") as f:
            f[hc.name] = hc.to_root("D")

        summary = {
            "name": self.name, "eta_edges": eta.tolist(), "pt_edges": pt.tolist(),
            "nominal": eps.tolist(), "stat": R["stat"].tolist(), "syst": R["syst"].tolist(), "total": R["total"].tolist(),
            "components": {k: v.tolist() for k, v in R["components"].items()},
            "details": R["details"],
            "fit": {"fval": R["nominal"].fval, "status": R["nominal"].status, "valid": R["nominal"].valid},
            "config_sources": self.sources,
        }
        with open(base + ".json", "w") as fh:
            json.dump(summary, fh, indent=1, default=str)
        with open(base + ".md", "w") as fh:
            fh.write(self.table())
        if self.ss.get("plots", True):
            from . import plotting

            comps = {"stat": R["stat"], **R["components"]}
            plotting.syst_breakdown(eta, pt, eps, comps, R["total"], base + "_breakdown", skip_eta_bins=self.skip)
            hs = [objs[f"{self.name}_all_num"]]
            labels = ["nominal"]
            for src, vr in R["variations"].items():
                for vname, v in vr.items():
                    hs.append(LH.rates_hist("x", eta, pt, v, np.zeros_like(v)))
                    labels.append(f"{src}:{vname}")
            if len(hs) > 1:
                plotting.compare_rates(hs, labels, base + "_variations", ref=0, skip_eta_bins=self.skip, ratio_label="var. / nominal")
        log.info("Wrote %s.{root,json,md}", base)

    def stat_eigen(self):
        """Eigen-decomposition of the stat. covariance -> list of shift maps."""
        R = self.result
        cov = np.array(R["nominal"].cov)
        eps = R["nominal"].rates
        free = R["nominal"].free
        boot = R.get("boot")
        if boot is not None and self.ss.get("stat", {}).get("method") == "bootstrap":
            b = boot.reshape(len(boot), -1)[:, free]
            b = b[np.all(np.isfinite(b), axis=1)]
            c = np.cov(b, rowvar=False)
            cov = np.zeros((free.size, free.size))
            cov[np.ix_(free, free)] = c
        sub = cov[np.ix_(free, free)]
        if sub.size == 0 or not np.any(sub):
            return []
        w, v = np.linalg.eigh(sub)
        out = []
        for k in np.argsort(w)[::-1]:
            if w[k] <= 0:
                continue
            full = np.zeros(free.size)
            full[free] = np.sqrt(w[k]) * v[:, k]
            out.append(full.reshape(eps.shape))
        return out

    def table(self):
        R = self.result
        eps = R["nominal"].rates
        eta, pt = self.eta_edges, self.pt_edges
        comps = {"stat": R["stat"], **R["components"]}
        head = "| |eta| bin | pT bin [GeV] | rate | " + " | ".join(comps) + " | syst | total |\n"
        head += "|" + "---|" * (len(comps) + 5) + "\n"
        rows = []
        for e in range(len(eta) - 1):
            if e in self.skip:
                continue
            for p in range(len(pt) - 1):
                rel = lambda d: f"{100 * d[e, p] / eps[e, p]:.1f}%" if eps[e, p] > 0 else "-"  # noqa: E731
                rows.append(f"| [{eta[e]:g}, {eta[e + 1]:g}] | [{pt[p]:g}, {pt[p + 1]:g}] | {eps[e, p]:.3e} | "
                            + " | ".join(rel(d) for d in comps.values()) + f" | {rel(R['syst'])} | {rel(R['total'])} |")
        return f"# Charge-flip rate uncertainties: {self.name}\n\nRelative uncertainties per bin.\n\n" + head + "\n".join(rows) + "\n"


def run(cfg, sample_name=None, max_events=None):
    r = SystematicsRunner(cfg, sample_name, max_events)
    res = r.run()
    log.info("\n%s", r.table())
    return r, res
