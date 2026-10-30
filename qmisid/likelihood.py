"""Step 2: likelihood charge-flip rates from Z->ee (port of Likelihoodv5New.C).

Events with two selected electrons in a window around the Z mass are counted
in "pair cells" ``(eta_lo, eta_hi, pt_lo, pt_hi)`` (the lepton with the lower
|eta| bin first), separately for same-sign (SS) and opposite-sign (OS) pairs.
The rates ``eps_k`` (one per |eta| x pT bin) maximise the Poisson likelihood

    ln L = sum_cells  N_SS ln(N mu) - N mu,    N = N_SS + N_OS,
    mu = eps_a + eps_b                    (model "approx", as in the macro)
    mu = eps_a + eps_b - 2 eps_a eps_b    (model "exact")

The non-Z background is removed with a sideband subtraction: events in the
sidebands enter with weight ``-scale`` (macro: -0.5 = 20 GeV peak / 40 GeV of
sidebands). Peak and sideband counts are stored separately, so subtraction
choices can be changed after the (expensive) event processing.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field

import numpy as np

from . import config as C
from . import expr, io
from .hist import Hist, write_root

log = logging.getLogger(__name__)

TINY_N = 1e-18  # value used by the macro for cells with N <= 0


# ---------------------------------------------------------------------------
# binning helpers
# ---------------------------------------------------------------------------
def find_bins(edges, x, convention="right_closed", overflow_to_last=False):
    """Bin index (0..n-1) or -1 if outside.

    ``right_closed``: x in (lo, hi]  (FillVectorsEtaPt in the macro)
    ``left_closed`` : x in [lo, hi)  (ROOT TH1 convention)
    """
    edges = np.asarray(edges, float)
    x = np.asarray(x, float)
    nb = len(edges) - 1
    idx = np.full(x.shape, -1, np.int64)
    if convention == "right_closed":
        for e in edges:
            idx += x > e
    else:
        for e in edges:
            idx += x >= e
    if overflow_to_last:
        over = x > edges[-1] if convention == "right_closed" else x >= edges[-1]
        idx = np.where(over, nb - 1, idx)
    idx = np.where((idx < 0) | (idx >= nb) | np.isnan(x), -1, idx)
    return idx


# ---------------------------------------------------------------------------
# accumulation
# ---------------------------------------------------------------------------
class LikelihoodAccumulator:
    """Collects the SS/OS pair-cell counts for one configuration.

    Several accumulators (e.g. systematic variations) can be fed from the same
    chunks, so the ntuples are read only once.
    """

    def __init__(self, cfg, sample, name=None, sec=None, n_bootstrap=None, seed=12345):
        self.cfg = cfg
        self.sample = sample
        self.sec = sec if sec is not None else C.section(cfg, "likelihood", sample["name"])
        self.name = name or C.run_name(sample["name"], cfg)
        s = self.sec
        self.eta_edges = np.asarray(s["binning"]["eta"], float)
        self.pt_edges = np.asarray(s["binning"]["pt"], float)
        self.n_eta, self.n_pt = len(self.eta_edges) - 1, len(self.pt_edges) - 1
        self.eta_expr = s["leptons"]["eta"]
        self.pt_expr = s["leptons"]["pt"]
        self.selection = s.get("selection") or "1"
        self.ss_expr = s.get("same_sign", "same_charge")
        self.os_expr = s.get("opposite_sign", "opposite_charge")
        self.mll_expr = s.get("mll", "m_l1l2*0.001")
        self.window = [float(v) for v in s.get("mll_window", [61, 81, 101, 121])]
        low, low_z, up_z, up = self.window
        self.subtract = bool(s.get("subtract_background", True))
        scale = s.get("sideband_scale", "auto")
        self.scale = (up_z - low_z) / ((low_z - low) + (up - up_z)) if scale in (None, "auto") else float(scale)
        self.convention = s.get("bin_convention", "right_closed")
        # the macro stores pT, |eta| and Mll in C++ floats: emulate the rounding (legacy)
        self.f32 = bool(s.get("float32_kinematics", False))
        self.weight_expr = None if sample.get("is_data") else sample.get("weight")
        self.n_boot = int(s.get("n_bootstrap", 0) if n_bootstrap is None else n_bootstrap)
        self.rng = np.random.default_rng(seed)

        shape = (2, 2, self.n_eta, self.n_eta, self.n_pt, self.n_pt)  # (SS/OS, peak/side, ...)
        self.counts = np.zeros(shape)
        self.counts_w2 = np.zeros(shape)
        self.boot = np.zeros((self.n_boot,) + shape) if self.n_boot else None
        self._boot_idx = None
        nb = int(round(up - low)) or 1
        mll_edges = np.linspace(low, up, nb + 1)
        self.mll_ss = Hist(mll_edges, name="SS", title="SS", axis_titles=["M_{ee} [GeV]"])
        self.mll_os = Hist(mll_edges, name="OS", title="OS", axis_titles=["M_{ee} [GeV]"])
        self.n_used = 0.0
        self.n_outside = 0.0
        self.n_bad_bins = 0
        self.n_no_charge = 0

    def expressions(self):
        return [self.selection, self.ss_expr, self.os_expr, self.mll_expr, self.weight_expr,
                *self.eta_expr, *self.pt_expr]

    def process(self, arrays, n, cache=None, external_boot=False):
        """Fill from one chunk. ``cache`` (dict) shares evaluated expressions
        between accumulators fed with the same chunk. With ``external_boot``
        the caller provides the bootstrap weights through :meth:`fill_boot`
        (so that several accumulators use the same replicas)."""
        cache = {} if cache is None else cache
        self._boot_idx = None

        def ev(e):
            if e not in cache:
                cache[e] = expr.evaluate(e, arrays, n)
            return cache[e]

        def kin(e):
            v = ev(e)
            return v.astype(np.float32).astype(np.float64) if self.f32 else v

        w = ev(self.weight_expr) if self.weight_expr else np.ones(n)
        sel = ev(self.selection) != 0
        mll = kin(self.mll_expr)
        is_ss = ev(self.ss_expr) != 0
        is_os = (ev(self.os_expr) != 0) & ~is_ss
        low, low_z, up_z, up = self.window
        inwin = ~((mll < low) | (mll > up))
        self.n_outside += float(w[sel & ~inwin].sum())
        m = sel & inwin
        nocharge = m & ~is_ss & ~is_os
        if nocharge.any():
            self.n_no_charge += int(nocharge.sum())
            m &= ~nocharge

        idx = np.nonzero(m)[0]
        bins = []
        for i in (0, 1):
            be = find_bins(self.eta_edges, kin(self.eta_expr[i])[idx], self.convention)
            bp = find_bins(self.pt_edges, kin(self.pt_expr[i])[idx], self.convention, overflow_to_last=True)
            bins.append((be, bp))
        good = (bins[0][0] >= 0) & (bins[1][0] >= 0) & (bins[0][1] >= 0) & (bins[1][1] >= 0)
        if not good.all():
            self.n_bad_bins += int((~good).sum())
            idx = idx[good]
            bins = [(be[good], bp[good]) for be, bp in bins]
        if len(idx) == 0:
            return
        e1, p1 = bins[0]
        e2, p2 = bins[1]
        swap = e1 > e2
        elo, ehi = np.where(swap, e2, e1), np.where(swap, e1, e2)
        plo, phi = np.where(swap, p2, p1), np.where(swap, p1, p2)
        ww = w[idx]
        mm = mll[idx]
        charge = np.where(is_ss[idx], 0, 1)
        region = np.where((mm > low_z) & (mm < up_z), 0, 1)
        flat = np.ravel_multi_index((charge, region, elo, ehi, plo, phi), self.counts.shape)
        size = self.counts.size
        self.counts += np.bincount(flat, weights=ww, minlength=size).reshape(self.counts.shape)
        self.counts_w2 += np.bincount(flat, weights=ww * ww, minlength=size).reshape(self.counts.shape)
        self.n_used += float(ww.sum())
        ss_i = charge == 0
        self.mll_ss.fill(mm[ss_i], weight=ww[ss_i])
        self.mll_os.fill(mm[~ss_i], weight=ww[~ss_i])
        if self.n_boot:
            self._boot_idx, self._boot_flat, self._boot_w = idx, flat, ww
            if not external_boot:
                for r in range(self.n_boot):
                    self.fill_boot(r, self.rng.poisson(1.0, size=n))

    def fill_boot(self, r, pw):
        """Add bootstrap replica ``r`` for the current chunk; ``pw`` are the
        Poisson(1) weights of all events in the chunk."""
        if self._boot_idx is None:
            return
        self.boot[r] += np.bincount(self._boot_flat, weights=self._boot_w * pw[self._boot_idx],
                                    minlength=self.counts.size).reshape(self.counts.shape)

    # -- derived counts ---------------------------------------------------------
    def pairs(self, subtract=None, scale=None):
        """Return background-subtracted ``(ss, os)`` cell counts."""
        subtract = self.subtract if subtract is None else subtract
        scale = self.scale if scale is None else scale
        f = -scale if subtract else 1.0
        ss = self.counts[0, 0] + f * self.counts[0, 1]
        os_ = self.counts[1, 0] + f * self.counts[1, 1]
        return ss, os_

    def pairs_w2(self, subtract=None, scale=None):
        subtract = self.subtract if subtract is None else subtract
        scale = self.scale if scale is None else scale
        f2 = scale * scale if subtract else 1.0
        return self.counts_w2[0, 0] + f2 * self.counts_w2[0, 1], self.counts_w2[1, 0] + f2 * self.counts_w2[1, 1]

    def boot_pairs(self, r, subtract=None, scale=None):
        subtract = self.subtract if subtract is None else subtract
        scale = self.scale if scale is None else scale
        f = -scale if subtract else 1.0
        b = self.boot[r]
        return b[0, 0] + f * b[0, 1], b[1, 0] + f * b[1, 1]

    def save_pairs(self, filename):
        np.savez_compressed(
            filename, counts=self.counts, counts_w2=self.counts_w2,
            boot=self.boot if self.boot is not None else np.zeros((0,)),
            eta_edges=self.eta_edges, pt_edges=self.pt_edges, window=np.asarray(self.window),
            subtract=self.subtract, scale=self.scale,
        )

    def merge(self, other):
        self.counts += other.counts
        self.counts_w2 += other.counts_w2
        self.n_used += other.n_used
        self.n_outside += other.n_outside
        if self.boot is not None and other.boot is not None:
            self.boot += other.boot
        self.mll_ss.add(other.mll_ss)
        self.mll_os.add(other.mll_os)


# ---------------------------------------------------------------------------
# fit
# ---------------------------------------------------------------------------
@dataclass
class FitResult:
    rates: np.ndarray            # (n_eta, n_pt)
    errors: np.ndarray           # (n_eta, n_pt)
    cov: np.ndarray              # (n_eta*n_pt, n_eta*n_pt), eta-major
    fval: float
    edm: float
    valid: bool
    accurate_covar: bool
    nfcn: int
    start: np.ndarray
    free: np.ndarray             # bool mask of fitted parameters (eta-major)
    message: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def status(self):
        """Mimic TMinuit's icstat: 3 = full accurate covariance."""
        if not self.valid:
            return 1
        return 3 if self.accurate_covar else 2


def default_start(n_eta, n_pt):
    """Default starting values of the macro (used when no reference is given).

    Note: the macro fills this pattern pT-major into an eta-major parameter
    vector; here it is filled consistently eta-major.
    """
    init = 0.00005
    out = np.zeros((n_eta, n_pt))
    for p in range(n_pt):
        v = (p + 1) * init
        if p == 2:
            v = 0.0005
        if p >= 3:
            v = 0.001
        for e in range(n_eta):
            out[e, p] = v * 2 ** e
    return out


def read_reference(ref):
    """Starting values from a (truth) rate histogram: ``{file:, hist:}``."""
    h = Hist.read(ref["file"], ref["hist"])
    return h.values()


def _cells(n_eta, n_pt, skip):
    """Cell index arrays over the macro's loop (eta_hi >= eta_lo)."""
    E1, E2, P1, P2 = np.meshgrid(np.arange(n_eta), np.arange(n_eta), np.arange(n_pt), np.arange(n_pt), indexing="ij")
    keep = E2 >= E1
    for s in skip:
        keep &= (E1 != s) & (E2 != s)
    a = (E1 * n_pt + P1)[keep]
    b = (E2 * n_pt + P2)[keep]
    return keep, a, b


def make_nll(ss, os_, n_eta, n_pt, skip=(), model="approx"):
    """Return ``nll(eps_full)`` for the given cell counts (vectorised)."""
    keep, a, b = _cells(n_eta, n_pt, skip)
    ss_c = ss[keep]
    all_c = ss_c + os_[keep]
    all_c = np.where(all_c <= 0, TINY_N, all_c)
    has_ss = ss_c != 0
    exact = model == "exact"

    def nll(eps):
        ea, eb = eps[a], eps[b]
        mu = ea + eb - 2.0 * ea * eb if exact else ea + eb
        lam = all_c * mu
        with np.errstate(divide="ignore", invalid="ignore"):
            logterm = np.where(has_ss, ss_c * np.log(lam), 0.0)
        return -float(np.sum(logterm - lam))

    def grad(eps):
        ea, eb = eps[a], eps[b]
        mu = ea + eb - 2.0 * ea * eb if exact else ea + eb
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.where(has_ss, ss_c / mu, 0.0) - all_c  # d(lnL)/d(mu)
        da = d * ((1 - 2 * eb) if exact else 1.0)
        db = d * ((1 - 2 * ea) if exact else 1.0)
        g = np.zeros(len(eps))
        np.add.at(g, a, da)
        np.add.at(g, b, db)
        return -g

    return nll, grad


def fit_rates(ss, os_, eta_edges, pt_edges, skip_eta_bins=(), start=None, minuit_cfg=None, model="approx"):
    """Maximum-likelihood fit of the rates with iminuit (MIGRAD [+ HESSE])."""
    from iminuit import Minuit

    mc = {"errordef": 0.5, "strategy": 1, "step": 1e-10, "limits": [1e-10, 0.99], "hesse": True,
          "ncall": 10_000_000, "tol": 0.1, "use_gradient": False}
    mc.update(minuit_cfg or {})
    n_eta, n_pt = len(eta_edges) - 1, len(pt_edges) - 1
    npar = n_eta * n_pt
    skip = list(skip_eta_bins or [])
    free = np.ones((n_eta, n_pt), bool)
    for s in skip:
        free[s, :] = False
    free = free.ravel()
    lo, hi = mc["limits"]
    s0 = default_start(n_eta, n_pt).ravel()
    if start is not None:
        st = np.asarray(start, float).ravel()
        if st.shape != s0.shape:
            raise ValueError(f"Reference has {st.size} bins, expected {npar}")
        ok = (st > lo) & (st < hi) & np.isfinite(st)
        s0 = np.where(ok, st, s0)
    nll_full, grad_full = make_nll(ss, os_, n_eta, n_pt, skip, model)
    fidx = np.nonzero(free)[0]

    def expand(x):
        e = np.zeros(npar)
        e[fidx] = x
        return e

    def fcn(x):
        return nll_full(expand(np.asarray(x)))

    def gcn(x):
        return grad_full(expand(np.asarray(x)))[fidx]

    names = [f"epsilon_{k + 1}" for k in fidx]
    m = Minuit(fcn, s0[fidx], grad=gcn if mc["use_gradient"] else None, name=names)
    m.errordef = mc["errordef"]
    m.strategy = mc["strategy"]
    m.tol = mc["tol"]
    m.limits = [(lo, hi)] * len(fidx)
    m.errors = [mc["step"]] * len(fidx)
    m.print_level = int(mc.get("print_level", 0))
    m.migrad(ncall=int(mc["ncall"]))
    if mc["hesse"] and m.valid:
        try:
            m.hesse()
        except Exception as exc:  # pragma: no cover
            log.warning("HESSE failed: %s", exc)
    vals = expand(np.asarray(m.values))
    errs = expand(np.asarray(m.errors))
    cov = np.zeros((npar, npar))
    if m.covariance is not None:
        cov[np.ix_(fidx, fidx)] = np.asarray(m.covariance)
    fm = m.fmin
    return FitResult(
        rates=vals.reshape(n_eta, n_pt), errors=errs.reshape(n_eta, n_pt), cov=cov,
        fval=float(fm.fval), edm=float(fm.edm), valid=bool(m.valid),
        accurate_covar=bool(fm.has_accurate_covar), nfcn=int(fm.nfcn), start=s0.reshape(n_eta, n_pt),
        free=free, message=str(fm) if not m.valid else "",
    )


def bootstrap_rates(acc, nominal, sec, subtract=None, scale=None):
    """Refit the Poisson-bootstrap replicas; return (n_boot, n_eta, n_pt)."""
    out = []
    for r in range(acc.n_boot):
        ss, os_ = acc.boot_pairs(r, subtract, scale)
        res = fit_rates(ss, os_, acc.eta_edges, acc.pt_edges, sec.get("skip_eta_bins", []), nominal.rates,
                        {**sec.get("minuit", {}), "hesse": False}, sec.get("model", "approx"))
        out.append(res.rates if res.valid else np.full_like(res.rates, np.nan))
    return np.asarray(out)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def rates_hist(name, eta_edges, pt_edges, rates, errors):
    h = Hist(eta_edges, pt_edges, name=name, title=name, axis_titles=["Lepton #eta", "Lepton p_{T} [GeV]"])
    h.set_values(rates, errors)
    h.entries = float(rates.size)
    return h


def output_file(cfg, sample_name, name=None):
    sample = C.get_sample(cfg, sample_name)
    sec = C.section(cfg, "likelihood", sample_name)
    key = "output_dir_data" if sample.get("is_data") else "output_dir_mc"
    outdir = sec.get(key) or C.output_path(cfg, "Likelihood_Data" if sample.get("is_data") else "Likelihood_MC")
    return os.path.join(outdir, f"{name or C.run_name(sample_name, cfg)}.root")


def fit_accumulator(acc, sec, reference=None):
    ss, os_ = acc.pairs()
    if reference is None and sec.get("reference"):
        try:
            reference = read_reference(sec["reference"])
            log.info("Starting values from %s:%s", sec["reference"]["file"], sec["reference"]["hist"])
        except (FileNotFoundError, KeyError, OSError) as exc:
            log.warning("Could not read reference rates (%s); using default starting values", exc)
    res = fit_rates(ss, os_, acc.eta_edges, acc.pt_edges, sec.get("skip_eta_bins", []), reference,
                    sec.get("minuit", {}), sec.get("model", "approx"))
    return res


def write_output(acc, res, outfile, boot=None):
    h = rates_hist(f"{acc.name}_all_num", acc.eta_edges, acc.pt_edges, res.rates, res.errors)
    npar = res.cov.shape[0]
    cov = Hist(np.arange(npar + 1), np.arange(npar + 1), name=f"{acc.name}_covariance", title="rate covariance (eta-major)")
    cov.set_values(res.cov)
    objs = {h.name: h, "OS": acc.mll_os, "SS": acc.mll_ss}
    write_root(outfile, {h.name: h, "OS": acc.mll_os, "SS": acc.mll_ss})
    import uproot

    with uproot.update(outfile) as fout:
        fout[cov.name] = cov.to_root("D")
        if boot is not None and len(boot):
            bstd = np.nanstd(boot, axis=0, ddof=1)
            hb = rates_hist(f"{acc.name}_bootstrap_std", acc.eta_edges, acc.pt_edges, np.nanmean(boot, axis=0), bstd)
            fout[hb.name] = hb.to_root("F")
    acc.save_pairs(outfile.replace(".root", "_pairs.npz"))
    if boot is not None and len(boot):
        np.save(outfile.replace(".root", "_bootstrap.npy"), boot)
    ss, os_ = acc.pairs()
    summary = {
        "name": acc.name, "selection": acc.selection, "mll_window": acc.window, "subtract_background": acc.subtract,
        "sideband_scale": acc.scale, "bin_convention": acc.convention,
        "eta_edges": acc.eta_edges.tolist(), "pt_edges": acc.pt_edges.tolist(),
        "rates": res.rates.tolist(), "errors": res.errors.tolist(), "start": res.start.tolist(),
        "fit": {"fval": res.fval, "edm": res.edm, "valid": res.valid, "status": res.status,
                "accurate_covar": res.accurate_covar, "nfcn": res.nfcn},
        "events_used_sumw": acc.n_used, "events_outside_window_sumw": acc.n_outside,
        "events_dropped_bad_bins": acc.n_bad_bins, "events_dropped_no_charge": acc.n_no_charge,
        "total_ss_after_subtraction": float(ss.sum()), "total_os_after_subtraction": float(os_.sum()),
    }
    with open(outfile.replace(".root", ".json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    return objs


def format_rates(res, eta_edges, pt_edges):
    lines = []
    for e in range(len(eta_edges) - 1):
        for p in range(len(pt_edges) - 1):
            k = e * (len(pt_edges) - 1) + p
            if not res.free[k]:
                continue
            lines.append(f"  epsilon_{k + 1:<3d} |eta| [{eta_edges[e]:.2f},{eta_edges[e + 1]:.2f}] "
                         f"pT [{pt_edges[p]:g},{pt_edges[p + 1]:g}]: {res.rates[e, p]:.6e} +- {res.errors[e, p]:.3e}")
    return "\n".join(lines)


def run(cfg, sample_name, max_events=None, outfile=None):
    sample = C.get_sample(cfg, sample_name)
    sec = C.section(cfg, "likelihood", sample_name)
    acc = LikelihoodAccumulator(cfg, sample, sec=sec)
    files = io.resolve_inputs(sample)
    tree = sec.get("tree", cfg.get("tree"))
    log.info("Likelihood rates for %s (%d files, tree %s)", acc.name, len(files), tree)
    log.info("Selection: %s | window %s | subtract bkg %s (scale %.3f) | model %s",
             acc.selection, acc.window, acc.subtract, acc.scale, sec.get("model", "approx"))
    for arrays, n in io.iterate(files, tree, acc.expressions(), step_size=sec.get("step_size", "50 MB"), max_events=max_events):
        acc.process(arrays, n)
    log.info("Events used (sum w): %.6g ; outside window: %.6g", acc.n_used, acc.n_outside)
    if acc.n_bad_bins:
        log.warning("%d events dropped: lepton outside the |eta|/pT binning", acc.n_bad_bins)
    if acc.n_no_charge:
        log.warning("%d events dropped: neither SS nor OS", acc.n_no_charge)
    res = fit_accumulator(acc, sec)
    if not res.valid:
        log.error("Minuit did not converge:\n%s", res.message)
    log.info("Fit: fval=%.6f edm=%.3g status=%d nfcn=%d\n%s", res.fval, res.edm, res.status, res.nfcn,
             format_rates(res, acc.eta_edges, acc.pt_edges))
    boot = bootstrap_rates(acc, res, sec) if acc.n_boot else None
    outfile = outfile or output_file(cfg, sample_name)
    write_output(acc, res, outfile, boot)
    log.info("Wrote %s", outfile)
    return acc, res, outfile
