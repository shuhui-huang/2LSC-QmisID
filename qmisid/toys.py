"""Toy ntuples in the v5 ``qmisid_cr`` format, with known charge-flip rates.

Used by the tests and to cross-check the Python port against the original
ROOT macros. Branches (as read by the macros):

    l1_pt, l2_pt, l1_eta, l2_eta (float, MeV)   m_l1l2 (float, MeV)
    l1_IDTight, l2_IDTight (int)  l1_Tight, l2_Tight (char)
    same_charge, opposite_charge (int)
    MC only: weight (double), l1_IFF, l2_IFF (int; 2 prompt, 3 charge flip, 4 fake)
"""

from __future__ import annotations

import os

import numpy as np
import uproot

MZ, GZ = 91.1876, 2.4952


ETA_EDGES = np.array([0, 1.37, 1.52, 2.0, 2.6])
PT_EDGES = np.array([20, 50, 100, 200, 2600])
# piecewise-constant generating rates (similar to the v5 Z+jets truth rates)
BINNED_RATES = np.array([
    [4.0e-5, 6.0e-5, 1.2e-4, 8.0e-4],
    [2.0e-4, 3.0e-4, 6.0e-4, 2.0e-3],
    [4.5e-4, 6.0e-4, 1.4e-3, 1.0e-2],
    [8.0e-4, 1.1e-3, 2.5e-3, 2.0e-2],
])


def true_rate(abseta, pt, scale=1.0, model="binned"):
    """Generating charge-flip probability (|eta|, pT [GeV]).

    ``binned``: constant inside the analysis bins (the likelihood is then
    unbiased); ``smooth``: strong variation inside the bins (exhibits the
    binning/method bias that the ``mc_method`` systematic measures).
    """
    if model == "smooth":
        return np.clip(scale * 2.5e-5 * (1 + 60 * (abseta / 2.5) ** 3) * (1 + (pt / 60.0) ** 1.5), 0, 0.2)
    ie = np.clip(np.searchsorted(ETA_EDGES, abseta, side="right") - 1, 0, len(ETA_EDGES) - 2)
    ip = np.clip(np.searchsorted(PT_EDGES, pt, side="right") - 1, 0, len(PT_EDGES) - 2)
    return np.clip(scale * BINNED_RATES[ie, ip], 0, 0.5)


def _kinematics(rng, n, mll):
    pt_a = np.abs(rng.normal(0.5, 0.12, n)) * mll + 5.0
    tail = rng.random(n) < 0.04
    pt_a[tail] += rng.exponential(150.0, tail.sum())
    pt_b = pt_a * rng.uniform(0.45, 1.0, n)
    pt_b = np.maximum(pt_b, 20.5 + rng.exponential(5.0, n))
    pt_a = np.maximum(pt_a, pt_b + 0.1)
    eta = []
    for _ in range(2):
        e = rng.normal(0, 1.3, n)
        bad = np.abs(e) >= 2.55
        while bad.any():
            e[bad] = rng.normal(0, 1.3, bad.sum())
            bad = np.abs(e) >= 2.55
        e[e == 0] = 1e-4
        eta.append(e)
    return pt_a, pt_b, eta[0], eta[1]


def generate(n_z, n_bkg, seed=1, rate_scale=1.0, mc=True, process="Zjets", fake_frac=1e-5, rate_model="binned"):
    rng = np.random.default_rng(seed)
    if process == "ttbar":
        mll_z = np.clip(rng.gamma(4.0, 30.0, n_z), 15, 1000)
    else:
        mll_z = MZ + GZ / 2 * np.tan(np.pi * (rng.random(n_z) - 0.5))
        mll_z = mll_z + rng.normal(0, 1.5, n_z)
        keep = (mll_z > 40) & (mll_z < 200)
        mll_z = mll_z[keep]
    nz = len(mll_z)
    mll_b = 40 + rng.exponential(150.0, n_bkg)
    mll_b = mll_b[mll_b < 200]
    nb = len(mll_b)
    mll = np.concatenate([mll_z, mll_b])
    n = nz + nb
    pt1, pt2, eta1, eta2 = _kinematics(rng, n, mll)

    # truth classification: prompt (2) / flip (3) / fake (4)
    iff = []
    flipped = []
    for pt, eta in ((pt1, eta1), (pt2, eta2)):
        p = true_rate(np.abs(eta), pt, rate_scale, rate_model)
        f = rng.random(n) < p
        cls = np.where(f, 3, 2)
        flipped.append(f)
        iff.append(cls)
    is_bkg = np.arange(n) >= nz
    fake = (rng.random(n) < fake_frac) | is_bkg & (rng.random(n) < 0.7)
    which = rng.integers(0, 2, n)
    iff[0] = np.where(fake & (which == 0), 4, iff[0])
    iff[1] = np.where(fake & (which == 1), 4, iff[1])
    # observed charge relation: signal pairs are true OS; one flip -> SS
    ss = flipped[0] ^ flipped[1]
    # fakes: random charge relation; continuum background: mostly OS
    rand_ss = rng.random(n) < np.where(is_bkg, 0.05, 0.45)
    ss = np.where(fake | is_bkg, rand_ss, ss)

    tight1 = rng.random(n) < 0.92
    tight2 = rng.random(n) < 0.92
    out = {
        "l1_pt": (pt1 * 1000).astype(np.float32), "l2_pt": (pt2 * 1000).astype(np.float32),
        "l1_eta": eta1.astype(np.float32), "l2_eta": eta2.astype(np.float32),
        "m_l1l2": (mll * 1000).astype(np.float32),
        "l1_IDTight": tight1.astype(np.int32), "l2_IDTight": tight2.astype(np.int32),
        "l1_Tight": tight1.astype(np.int8), "l2_Tight": tight2.astype(np.int8),
        "same_charge": ss.astype(np.int32), "opposite_charge": (~ss).astype(np.int32),
    }
    if mc:
        w = rng.normal(1.0, 0.2, n)
        w[rng.random(n) < 0.02] *= -1
        out["weight"] = w.astype(np.float64)
        out["l1_IFF"] = iff[0].astype(np.int32)
        out["l2_IFF"] = iff[1].astype(np.int32)
    perm = rng.permutation(n)
    return {k: v[perm] for k, v in out.items()}


def write(filename, branches, tree="qmisid_cr", basket=200_000):
    """Write a TTree (explicitly: recent uproot writes RNTuple for ``f[k] = dict``)."""
    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    n = len(next(iter(branches.values())))
    with uproot.recreate(filename) as f:
        t = f.mktree(tree, {k: v.dtype for k, v in branches.items()})
        for start in range(0, n, basket):
            t.extend({k: v[start:start + basket] for k, v in branches.items()})


def make_all(outdir, n_data=2_000_000, n_z=1_500_000, n_tt=300_000, seed=7, rate_model="binned"):
    """Write data / Z+jets / ttbar toy ntuples and the matching ``.list`` files."""
    os.makedirs(outdir, exist_ok=True)
    files = {
        "Data": ("toy_2lSC_QMisID_data.root", generate(n_data, int(0.01 * n_data), seed, 1.0, mc=False, rate_model=rate_model)),
        "Zjets": ("toy_2lSC_QMisID_Vjets_Zee.root", generate(n_z, int(0.002 * n_z), seed + 1, 1.0, mc=True, rate_model=rate_model)),
        "ttbar": ("toy_2lSC_QMisID_ttbar.root", generate(n_tt, 0, seed + 2, 1.3, mc=True, process="ttbar", fake_frac=0.03, rate_model=rate_model)),
    }
    lists = {}
    for key, (fname, br) in files.items():
        write(os.path.join(outdir, fname), br)
        lst = os.path.join(outdir, f"{key}{'' if key == 'Data' else 'MC'}.list")
        with open(lst, "w") as fh:
            fh.write(fname + "\n")
        lists[key] = lst
    return lists
