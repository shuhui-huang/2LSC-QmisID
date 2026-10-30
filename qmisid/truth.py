"""Step 1: truth-matched charge-flip rates in MC (port of ChargeFlipTruthv5New.C).

For each lepton ``i`` (leading/sub-leading) the macro fills, in bins of
(|eta|, pT):

* denominator ``den_i``: all leptons of events passing the baseline selection;
* numerators: leptons with a wrong-track charge flip (``IFF==3``), a photon
  conversion flip (dummy in v5), or either (``total``);
* the event-type bookkeeping histograms (RR, RQ, QQ, RF, QF, FF, ... with
  R = real/prompt, Q = charge flipped, F = fake).

The rate is ``(num_1 + num_2) / (den_1 + den_2)``. The output ROOT file keeps
the object names of the macro (``<name>_all_num`` holds the total rate, a ROOT
``Clone()``-naming quirk), so files are drop-in replacements. Raw counts are
additionally stored under ``raw/``.
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

_LEPS = (0, 1)


class TruthAccumulator:
    def __init__(self, cfg, sample, name=None):
        self.cfg = cfg
        self.sample = sample
        self.sec = C.section(cfg, "truth", sample["name"])
        self.name = name or C.run_name(sample["name"], cfg)
        b = self.sec["binning"]
        self.eta_edges = np.asarray(b["eta"], float)
        self.pt_edges = np.asarray(b["pt"], float)
        lep = self.sec["leptons"]
        self.eta_expr = lep["eta"]
        self.pt_expr = lep["pt"]
        base = self.sec.get("baseline") or "1"
        extra = self.sec.get("extra_cut") or ""
        self.base_expr = f"({base}) && ({extra})" if extra else base
        self.wt_expr = self.sec["wrong_track"]
        self.pc_expr = self.sec.get("photon_conversion") or ["0", "0"]
        self.real_expr = self.sec["real"]
        self.ss_expr = self.sec.get("same_charge", "same_charge==1")
        self.split = bool(self.sec.get("split_charge", False))
        # as in the macro, the e-/e+ selections are combined with the extra cut only
        andx = (lambda s: f"({s}) && ({extra})") if extra else (lambda s: s)  # noqa: E731
        self.ele_expr = [andx(s) for s in (self.sec.get("electron_selection") or ["0", "0"])]
        self.pos_expr = [andx(s) for s in (self.sec.get("positron_selection") or ["0", "0"])]
        use_w = self.sec.get("use_weight", True)
        self.weight_expr = sample.get("weight") if (use_w and not sample.get("is_data")) else None
        self.h = {}
        for key in _HIST_KEYS:
            self.h[key] = self._new(key)

    def _new(self, label):
        return Hist(
            self.eta_edges, self.pt_edges, name=f"{self.name}_{label}", title=f"{self.name}_{label}",
            axis_titles=["Lepton #eta", "Lepton p_{T} [GeV]"],
        )

    def expressions(self):
        out = [self.base_expr, self.ss_expr, self.weight_expr]
        for lst in (self.eta_expr, self.pt_expr, self.wt_expr, self.pc_expr, self.real_expr):
            out.extend(lst)
        if self.split:
            out.extend(self.ele_expr)
            out.extend(self.pos_expr)
        return out

    def process(self, arrays, n):
        ev = lambda e: expr.evaluate_mask(e, arrays, n)  # noqa: E731
        w = expr.evaluate(self.weight_expr, arrays, n) if self.weight_expr else np.ones(n)
        base = ev(self.base_expr)
        ss = ev(self.ss_expr)
        eta = [expr.evaluate(self.eta_expr[i], arrays, n) for i in _LEPS]
        pt = [expr.evaluate(self.pt_expr[i], arrays, n) for i in _LEPS]
        wt = [ev(self.wt_expr[i]) for i in _LEPS]
        pc = [ev(self.pc_expr[i]) for i in _LEPS]
        real = [ev(self.real_expr[i]) for i in _LEPS]
        flip = [wt[i] | pc[i] for i in _LEPS]
        fake = [~wt[i] & ~pc[i] & ~real[i] for i in _LEPS]

        any_h = self.h["all_den1"]
        bins = [any_h.bins_of(eta[i], pt[i]) for i in _LEPS]  # computed once per lepton

        def fill(key, lep, mask):
            sel = np.nonzero(mask)[0]
            if len(sel):
                self.h[key].fill_bins([b[sel] for b in bins[lep]], w[sel], [eta[lep][sel], pt[lep][sel]])

        for i in _LEPS:
            j = 1 - i
            fill(f"wrongTrack_num{i + 1}", i, base & wt[i])
            fill(f"photonConversion_num{i + 1}", i, base & pc[i])
            fill(f"total_num{i + 1}", i, base & flip[i])
            fill(f"all_den{i + 1}", i, base)
            # event categories (lepton i kinematics)
            fill("histo_RR_num", i, base & real[0] & real[1])
            fill("histo_QQ_num", i, base & flip[0] & flip[1])
            fill("histo_FF_num", i, base & fake[0] & fake[1])
            fill("histo_RQ_num", i, base & flip[i] & real[j])
            fill("histo_RnoR_num", i, base & real[i] & ~real[j] & ss)
            fill("histo_QF_forQflip_num", i, base & flip[i] & fake[j])
            fill("histo_QF_forFakes_num", i, base & fake[i] & flip[j])
            fill("histo_RF_num", i, base & fake[i] & real[j])
            fill("histo_FF_SS_num", i, base & fake[0] & fake[1] & ss)
            fill("histo_QF_forFakes_SS_num", i, base & fake[i] & flip[j] & ss)
            fill("histo_RF_SS_num", i, base & fake[i] & real[j] & ss)
            if self.split:
                ele = ev(self.ele_expr[i])
                pos = ev(self.pos_expr[i])
                fill("electron_num", i, ele & flip[i])
                fill("positron_num", i, pos & flip[i])
                fill("electron_den", i, ele)
                fill("positron_den", i, pos)

    # ------------------------------------------------------------------------
    def finalize(self):
        h = self.h
        n = self.name
        den = h["all_den1"].clone(f"{n}_all_den").add(h["all_den2"])
        wt = h["wrongTrack_num1"].clone(f"{n}_wrongTrack_num").add(h["wrongTrack_num2"])
        pc = h["photonConversion_num1"].clone(f"{n}_photon_conversion_num").add(h["photonConversion_num2"])
        allnum = wt.clone(f"{n}_all_num").add(pc)
        self.den, self.allnum = den, allnum

        ratios = [
            h["wrongTrack_num1"].divide(h["all_den1"]),
            h["wrongTrack_num2"].divide(h["all_den2"]),
            h["photonConversion_num1"].divide(h["all_den1"]),
            h["photonConversion_num2"].divide(h["all_den2"]),
            h["total_num1"].divide(h["all_den1"]),
            h["total_num2"].divide(h["all_den2"]),
            wt.divide(den),
            pc.divide(den),
            allnum.divide(den),
            h["electron_num"].divide(h["electron_den"]),
            h["positron_num"].divide(h["positron_den"]),
            h["histo_RnoR_num"].divide(den),
            h["histo_RQ_num"].divide(den),
            h["histo_RF_num"].divide(den),
            h["histo_RF_SS_num"].divide(den),
        ]
        stats = [wt.clone("m_histo_wrongTrack_num"), pc.clone("m_histo_photonConversion_num"), allnum.clone("m_histo_all_num")]
        out = {}
        for r in ratios + stats:
            out[r.name] = r
        raw = dict(h)
        raw["all_den"] = den
        raw["wrongTrack_num"] = wt
        raw["photonConversion_num"] = pc
        raw["all_num"] = allnum
        for k, v in raw.items():
            out[f"raw/{k}"] = v.clone(k)
        self.rate = allnum.divide(den, name=f"{n}_all_num")
        return out

    def summary(self):
        h = self.h
        lines = ["### Event categories (sum of weights, both leptons):"]
        for k in ("histo_RR_num", "histo_QQ_num", "histo_RQ_num", "histo_QF_forQflip_num", "histo_QF_forFakes_num",
                  "histo_RF_num", "histo_FF_num", "histo_QF_forFakes_SS_num", "histo_RF_SS_num", "histo_FF_SS_num", "histo_RnoR_num"):
            lines.append(f"  {k[6:-4]:<16s}: {h[k].integral():.6g}")
        lines.append("### Charge flips:")
        for k in ("wrongTrack_num1", "wrongTrack_num2", "photonConversion_num1", "photonConversion_num2",
                  "total_num1", "total_num2", "all_den1", "all_den2"):
            lines.append(f"  {k:<22s}: {h[k].integral():.6g}")
        return "\n".join(lines)


_HIST_KEYS = [
    "wrongTrack_num1", "wrongTrack_num2", "photonConversion_num1", "photonConversion_num2",
    "total_num1", "total_num2", "all_den1", "all_den2",
    "electron_num", "positron_num", "electron_den", "positron_den",
    "histo_RR_num", "histo_QQ_num", "histo_FF_num", "histo_RQ_num", "histo_RnoR_num",
    "histo_QF_forQflip_num", "histo_QF_forFakes_num", "histo_RF_num",
    "histo_FF_SS_num", "histo_RF_SS_num", "histo_QF_forFakes_SS_num",
]


def output_file(cfg, sample_name, name=None):
    sec = C.section(cfg, "truth", sample_name)
    outdir = sec.get("output_dir") or C.output_path(cfg, "Truth_MC")
    return os.path.join(outdir, f"{name or C.run_name(sample_name, cfg)}.root")


def run(cfg, sample_name, max_events=None, outfile=None):
    sample = C.get_sample(cfg, sample_name)
    acc = TruthAccumulator(cfg, sample)
    files = io.resolve_inputs(sample)
    sec = acc.sec
    tree = sec.get("tree", cfg.get("tree"))
    log.info("Truth rates for %s (%d files, tree %s)", acc.name, len(files), tree)
    log.info("Baseline: %s", acc.base_expr)
    for arrays, n in io.iterate(files, tree, acc.expressions(), step_size=sec.get("step_size", "50 MB"), max_events=max_events):
        acc.process(arrays, n)
    objs = acc.finalize()
    log.info("\n%s", acc.summary())
    outfile = outfile or output_file(cfg, sample_name)
    write_root(outfile, objs)
    with open(outfile.replace(".root", ".json"), "w") as fh:
        json.dump({
            "name": acc.name, "files": files, "baseline": acc.base_expr, "weight": acc.weight_expr,
            "eta_edges": acc.eta_edges.tolist(), "pt_edges": acc.pt_edges.tolist(),
            "rate": acc.rate.values().tolist(), "rate_error": acc.rate.errors().tolist(),
        }, fh, indent=1)
    log.info("Wrote %s", outfile)
    return acc, outfile
