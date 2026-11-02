"""Plots (matplotlib + mplhep, ATLAS style): rates, rate comparisons, closure,
systematic breakdowns. Ports of PlotRates.C / PlotRatesCompare.C and the
plotting part of ClosureTestv5New.C."""

from __future__ import annotations

import os

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

try:
    import mplhep as hep

    hep.style.use("ATLAS")
except Exception:  # pragma: no cover
    hep = None

COLORS = ["#3f3f3f", "#d62728", "#2ca02c", "#1f77b4", "#ff7f0e", "#9467bd", "#8c564b", "#17becf"]


def _label(ax, text="Internal", lumi=None, extra=None):
    if hep is not None:
        hep.atlas.text(text, ax=ax, loc=0)
    y = 0.88
    for line in [lumi, extra]:
        if line:
            ax.text(0.04, y, line, transform=ax.transAxes, fontsize=12, va="top")
            y -= 0.06


def _save(fig, outbase, formats=("pdf", "png")):
    os.makedirs(os.path.dirname(os.path.abspath(outbase)), exist_ok=True)
    for f in formats:
        fig.savefig(f"{outbase}.{f}", bbox_inches="tight")
    plt.close(fig)


def _pt_label(lo, hi, last):
    return rf"$p_T > {lo:g}$ GeV" if last else rf"$p_T \in [{lo:g}, {hi:g}]$ GeV"


def rates_1d(hist, outbase, title="Likelihood", lumi=r"$\sqrt{s}$ = 13 TeV, 140 fb$^{-1}$", skip_eta_bins=(1,),
             ylim=(1e-5, 1.0)):
    """Rate vs |eta| for each pT bin (PlotRate1D)."""
    eta, pt = hist.edges
    vals, errs = hist.values(), hist.errors()
    fig, ax = plt.subplots(figsize=(7, 5.5))
    centers = 0.5 * (eta[1:] + eta[:-1])
    half = 0.5 * np.diff(eta)
    for p in range(len(pt) - 1):
        m = np.ones(len(centers), bool)
        m[list(skip_eta_bins)] = False
        ax.errorbar(centers[m], vals[m, p], xerr=half[m], yerr=errs[m, p], fmt="^", color=COLORS[p % 8],
                    label=_pt_label(pt[p], pt[p + 1], p == len(pt) - 2), ms=5, lw=1)
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.set_xlim(eta[0], eta[-1])
    ax.set_xlabel(r"$|\eta|$")
    ax.set_ylabel(r"$\epsilon_{\mathrm{QmisID}}$")
    ax.legend(title=title, loc="lower right", fontsize=11)
    _label(ax, lumi=lumi)
    _save(fig, outbase)


def rates_2d(hist, outbase, title="Charge-flip rates"):
    """2D map with values printed (PlotRate2D)."""
    eta, pt = hist.edges
    vals = hist.values()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ptp = pt.copy()
    ptp[-1] = min(pt[-1], 3 * pt[-2]) if len(pt) > 2 else pt[-1]  # compress the last bin visually
    with np.errstate(divide="ignore"):
        z = np.where(vals > 0, vals, np.nan)
    mesh = ax.pcolormesh(eta, ptp, z.T, norm=matplotlib.colors.LogNorm(), cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="Charge-flip rate")
    for i in range(len(eta) - 1):
        for j in range(len(pt) - 1):
            if vals[i, j] > 0:
                ax.text(0.5 * (eta[i] + eta[i + 1]), np.sqrt(ptp[j] * ptp[j + 1]), f"{vals[i, j]:.2e}",
                        ha="center", va="center", fontsize=8, rotation=30, color="w")
    ax.set_yscale("log")
    ax.set_yticks(ptp)
    ax.set_yticklabels([f"{v:g}" for v in pt])
    ax.set_xlabel(r"$|\eta|$")
    ax.set_ylabel(r"$p_T$ [GeV]")
    ax.set_title(title, fontsize=12)
    _save(fig, outbase)


def compare_rates(hists, labels, outbase, ref=-1, skip_eta_bins=(1,), ratio_label=None, ylim=None):
    """Overlay several rate maps per pT bin with a ratio to ``hists[ref]``
    (PlotRatesCompare)."""
    eta, pt = hists[0].edges
    centers = 0.5 * (eta[1:] + eta[:-1])
    half = 0.5 * np.diff(eta)
    m = np.ones(len(centers), bool)
    m[list(skip_eta_bins)] = False
    refi = ref % len(hists)
    refv = hists[refi].values()
    cmap = plt.get_cmap("tab20" if len(hists) > 8 else "tab10")
    nh = len(hists)
    for p in range(len(pt) - 1):
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        allv = []
        for k, (h, lab) in enumerate(zip(hists, labels)):
            v, e = h.values(), h.errors()
            j = k if k < refi else k - 1  # index among the non-reference maps
            off = 0 * half if k == refi else ((j + 1) / nh - 0.5) * 0.8 * half  # spread markers inside the bin
            col = "k" if k == refi else cmap(k % cmap.N)
            xe = half[m] if k == refi else None
            ax.errorbar(centers[m] + off[m], v[m, p], xerr=xe, yerr=e[m, p], fmt="o" if k == refi else ".", color=col, label=lab, ms=5)
            allv.append(v[m, p])
            if k != refi:
                with np.errstate(divide="ignore", invalid="ignore"):
                    rr = v[:, p] / refv[:, p]
                    re = e[:, p] / refv[:, p]
                rx.errorbar(centers[m] + off[m], rr[m], yerr=re[m], fmt=".", color=col, ms=5)
        for b in eta[1:-1]:
            ax.axvline(b, color="0.85", lw=0.8, zorder=0)
            rx.axvline(b, color="0.85", lw=0.8, zorder=0)
        ax.set_yscale("log")
        if ylim is None:
            pos = np.concatenate(allv)
            pos = pos[pos > 0]
            if len(pos):
                ax.set_ylim(pos.min() / 5, pos.max() * 50)
        else:
            ax.set_ylim(*ylim)
        ax.set_ylabel(r"$\epsilon_{\mathrm{QmisID}}$")
        ax.legend(title=_pt_label(pt[p], pt[p + 1], p == len(pt) - 2), loc="upper left", fontsize=9, ncol=2 if nh > 5 else 1,
                  bbox_to_anchor=(0.0, 0.88), title_fontsize=10)
        _label(ax)
        rx.axhline(1, color="k", lw=1)
        rx.set_ylim(0, 2.5)
        rx.set_ylabel(ratio_label or f"X / {labels[refi]}", fontsize=11)
        rx.set_xlabel(r"$|\eta|$")
        rx.set_xlim(eta[0], eta[-1])
        _save(fig, f"{outbase}_pt{p + 1}")


def closure_plot(res, var, outbase, label="", lumi=None, logy=False):
    d = res["hists"]
    ss, os_, up, dn = d["SS"], d["OS"], d["OSUp"], d["OSDown"]
    edges = ss.edges[0]
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, gridspec_kw={"height_ratios": [7, 3], "hspace": 0.05})
    ax.stairs(ss.values(), edges, color="k", lw=2, label="SS")
    ax.stairs(os_.values(), edges, color="tab:blue", lw=2, label="OS weighted")
    lo_band = dn.values()
    hi_band = up.values()
    stat = res.get("stat_cov")
    ax.stairs(hi_band, edges, baseline=lo_band, fill=True, alpha=0.3, color="tab:blue", hatch="///", label=r"rates $\pm1\sigma$")
    ymax = max(ss.values().max(), os_.values().max(), 1e-9)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(0.1, ymax * 1e4)
    else:
        ax.set_ylim(0, 2.0 * ymax)
    ax.set_ylabel("Events")
    ax.legend(loc="upper right", fontsize=11)
    _label(ax, lumi=lumi, extra=label)
    from .closure import ratio_with_errors

    r, rdn, rup = ratio_with_errors(ss, os_, up, dn)
    centers = 0.5 * (edges[1:] + edges[:-1])
    rx.errorbar(centers, r, yerr=[np.nan_to_num(rdn), np.nan_to_num(rup)], xerr=0.5 * np.diff(edges), fmt="o", color="tab:blue", ms=4)
    if stat is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = stat[1:-1] / os_.values()
        rx.stairs(1 + np.nan_to_num(rel), edges, baseline=1 - np.nan_to_num(rel), fill=True, alpha=0.25, color="gray",
                  label="rate stat. (cov.)")
        rx.legend(fontsize=9, loc="upper right")
    rx.axhline(1, color="k", lw=1)
    rx.set_ylim(0, 2.0)
    rx.set_ylabel("SS/OS")
    rx.set_xlabel(var.get("xtitle", var["name"]))
    rx.set_xlim(edges[0], edges[-1])
    _save(fig, outbase)


def syst_breakdown(eta, pt, nominal, components, total, outbase, skip_eta_bins=(1,), title="Relative uncertainty"):
    """Relative uncertainty per source vs |eta|, one panel per pT bin."""
    npt = len(pt) - 1
    ncol = min(npt, 2)
    nrow = int(np.ceil(npt / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 4.8 * nrow), squeeze=False)
    m = np.ones(len(eta) - 1, bool)
    m[list(skip_eta_bins)] = False
    with np.errstate(divide="ignore", invalid="ignore"):
        for p in range(npt):
            ax = axes[p // ncol][p % ncol]
            for k, (name, delta) in enumerate(components.items()):
                rel = np.where(m, delta[:, p] / nominal[:, p], np.nan)
                ax.stairs(np.nan_to_num(rel), eta, color=COLORS[(k + 1) % 8], lw=1.8, label=name)
            relt = np.where(m, total[:, p] / nominal[:, p], np.nan)
            ax.stairs(np.nan_to_num(relt), eta, color="k", lw=2.5, ls="--", label="total")
            ax.set_title(_pt_label(pt[p], pt[p + 1], p == npt - 1), fontsize=12)
            ax.set_xlabel(r"$|\eta|$")
            ax.set_ylabel(title, fontsize=11)
            ax.set_xlim(eta[0], eta[-1])
            ax.set_ylim(bottom=0)
            if p == 0:
                ax.legend(fontsize=8, ncol=2)
    for q in range(npt, nrow * ncol):
        axes[q // ncol][q % ncol].axis("off")
    fig.tight_layout()
    _save(fig, outbase)
