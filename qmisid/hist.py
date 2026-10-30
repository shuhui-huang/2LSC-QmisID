"""Minimal ROOT-compatible weighted histograms (1D / 2D) on numpy.

The bin convention is ROOT's: bin ``i`` covers ``[edge[i-1], edge[i])``, index 0
is the underflow and ``n+1`` the overflow. Contents and sum of squared weights
are stored for all bins including the flow bins, so that ratios, overflow
merging and the written ``TH1F``/``TH2F`` objects behave exactly like the
originals produced by the ROOT macros.
"""

from __future__ import annotations

import numpy as np
import uproot
from uproot.writing import identify as _ident


def _is_uniform(edges):
    d = np.diff(edges)
    return len(d) > 0 and np.allclose(d, d[0], rtol=1e-12, atol=0)


def _find_bins(edges, x):
    """ROOT ``TAxis::FindBin`` for an array: 0 underflow, n+1 overflow.

    Fixed-width axes use ROOT's formula ``1 + int(n (x - lo) / (hi - lo))``;
    variable-width axes a binary search (``TMath::BinarySearch``)."""
    x = np.asarray(x, dtype=np.float64)
    nb = len(edges) - 1
    if _is_uniform(edges):
        lo, hi = float(edges[0]), float(edges[-1])
        with np.errstate(invalid="ignore"):
            idx = np.floor(nb * (x - lo) / (hi - lo))
            idx = np.where(x < lo, -1, np.where(x >= hi, nb, idx))
        idx = np.nan_to_num(idx, nan=nb).astype(np.int64) + 1
        return idx
    if nb <= 16:  # few edges: comparisons are faster than a binary search
        idx = np.zeros(x.shape, np.int64)
        for e in edges:
            idx += x >= e
    else:
        idx = np.searchsorted(edges, x, side="right")
    idx[np.isnan(x)] = nb + 1  # ROOT puts NaN in the overflow
    return idx


class Hist:
    """Weighted N-dimensional (N=1,2) histogram with ROOT bin conventions."""

    def __init__(self, *edges, name="h", title="", axis_titles=None):
        self.edges = [np.asarray(e, dtype=np.float64) for e in edges]
        shape = tuple(len(e) + 1 for e in self.edges)  # n bins + 2 flow
        self.sumw = np.zeros(shape)
        self.sumw2 = np.zeros(shape)
        self.name = name
        self.title = title
        self.axis_titles = list(axis_titles or [""] * len(edges))
        self.entries = 0.0
        # in-range statistics, as kept by ROOT (fTsumw, fTsumw2, fTsumwx, ...)
        self.stats = np.zeros(7)

    # -- construction helpers -------------------------------------------------
    @property
    def ndim(self):
        return len(self.edges)

    @property
    def nbins(self):
        return tuple(len(e) - 1 for e in self.edges)

    def clone(self, name=None):
        h = Hist(*self.edges, name=name or self.name, title=self.title, axis_titles=self.axis_titles)
        h.sumw = self.sumw.copy()
        h.sumw2 = self.sumw2.copy()
        h.entries = self.entries
        h.stats = self.stats.copy()
        return h

    def empty_like(self, name=None):
        return Hist(*self.edges, name=name or self.name, title=self.title, axis_titles=self.axis_titles)

    # -- filling ----------------------------------------------------------------
    def fill(self, *coords, weight=None):
        """Fill with arrays of coordinates. Entries with weight 0 are skipped
        (as ``TTree::Draw`` does for a selection evaluating to 0)."""
        coords = [np.asarray(c, dtype=np.float64) for c in coords]
        n = len(coords[0])
        w = np.ones(n) if weight is None else np.broadcast_to(np.asarray(weight, dtype=np.float64), (n,))
        keep = w != 0
        if not keep.all():
            coords = [c[keep] for c in coords]
            w = w[keep]
        if len(w) == 0:
            return self
        idx = [_find_bins(e, c) for e, c in zip(self.edges, coords)]
        return self.fill_bins(idx, w, coords)

    def bins_of(self, *coords):
        """Bin indices (with flow) of coordinates, to be reused with :meth:`fill_bins`."""
        return [_find_bins(e, np.asarray(c, dtype=np.float64)) for e, c in zip(self.edges, coords)]

    def fill_bins(self, idx, weight, coords=None):
        """Fill with precomputed bin indices (see :meth:`bins_of`)."""
        w = np.asarray(weight, dtype=np.float64)
        keep = w != 0
        if not keep.all():
            idx = [i[keep] for i in idx]
            w = w[keep]
            if coords is not None:
                coords = [c[keep] for c in coords]
        if len(w) == 0:
            return self
        flat = idx[0] if self.ndim == 1 else np.ravel_multi_index(idx, self.sumw.shape)
        size = self.sumw.size
        self.sumw += np.bincount(flat, weights=w, minlength=size).reshape(self.sumw.shape)
        self.sumw2 += np.bincount(flat, weights=w * w, minlength=size).reshape(self.sumw.shape)
        self.entries += len(w)
        if coords is None:
            return self
        inr = np.ones(len(w), bool)
        for e, i in zip(self.edges, idx):
            inr &= (i >= 1) & (i <= len(e) - 1)
        wi = w[inr]
        x = coords[0][inr]
        self.stats[0] += wi.sum()
        self.stats[1] += (wi * wi).sum()
        self.stats[2] += (wi * x).sum()
        self.stats[3] += (wi * x * x).sum()
        if self.ndim == 2:
            y = coords[1][inr]
            self.stats[4] += (wi * y).sum()
            self.stats[5] += (wi * y * y).sum()
            self.stats[6] += (wi * x * y).sum()
        return self

    # -- arithmetic ---------------------------------------------------------------
    def add(self, other, c=1.0):
        """``TH1::Add(other, c)``."""
        self.sumw = self.sumw + c * other.sumw
        self.sumw2 = self.sumw2 + c * c * other.sumw2
        self.entries += other.entries
        self.stats = self.stats + c * other.stats
        return self

    def __iadd__(self, other):
        return self.add(other)

    def divide(self, den, name=None):
        """Return ``self / den`` with ``TH1::Divide(den)`` error propagation
        (uncorrelated): e^2 = (e1^2 c2^2 + e2^2 c1^2) / c2^4; empty den -> 0."""
        out = self.clone(name)
        c1, c2 = self.sumw, den.sumw
        e1, e2 = self.sumw2, den.sumw2
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(c2 != 0, c1 / c2, 0.0)
            var = np.where(c2 != 0, (e1 * c2 * c2 + e2 * c1 * c1) / (c2 ** 4), 0.0)
        out.sumw, out.sumw2 = ratio, var
        return out

    def merge_overflow(self):
        """Add the overflow into the last visible bin (1D), as the closure macro."""
        assert self.ndim == 1
        self.sumw[-2] += self.sumw[-1]
        self.sumw2[-2] += self.sumw2[-1]
        self.sumw[-1] = 0.0
        self.sumw2[-1] = 0.0
        return self

    # -- access -------------------------------------------------------------------
    def _inner(self, a):
        sl = tuple(slice(1, -1) for _ in self.edges)
        return a[sl]

    def values(self, flow=False):
        return self.sumw.copy() if flow else self._inner(self.sumw).copy()

    def variances(self, flow=False):
        return self.sumw2.copy() if flow else self._inner(self.sumw2).copy()

    def errors(self, flow=False):
        return np.sqrt(self.variances(flow))

    def integral(self):
        """``TH1::Integral()`` (visible bins only)."""
        return float(self._inner(self.sumw).sum())

    def integral_error(self):
        return float(np.sqrt(self._inner(self.sumw2).sum()))

    def set_values(self, values, errors=None):
        """Set visible-bin contents (and errors) from arrays of shape ``nbins``."""
        sl = tuple(slice(1, -1) for _ in self.edges)
        self.sumw[sl] = values
        if errors is not None:
            self.sumw2[sl] = np.asarray(errors) ** 2
        return self

    def lookup(self, *coords, flow=True):
        """Content and error at ``FindBin(coords)`` for arrays of coordinates
        (``GetBinContent(FindBin(x, y))``; flow bins included)."""
        idx = tuple(_find_bins(e, c) for e, c in zip(self.edges, coords))
        return self.sumw[idx], np.sqrt(self.sumw2[idx])

    def bin_indices(self, *coords):
        return tuple(_find_bins(e, c) for e, c in zip(self.edges, coords))

    # -- ROOT I/O -----------------------------------------------------------------
    def to_root(self, kind="F"):
        """Return an uproot-writable ``TH1F/TH2F`` (or ``D``) model."""
        dtype = ">f4" if kind == "F" else ">f8"
        axes = []
        for i, e in enumerate(self.edges):
            uniform = np.allclose(np.diff(e), (e[-1] - e[0]) / (len(e) - 1))
            axes.append(
                _ident.to_TAxis(
                    fName=["xaxis", "yaxis", "zaxis"][i],
                    fTitle=self.axis_titles[i] if i < len(self.axis_titles) else "",
                    fNbins=len(e) - 1,
                    fXmin=float(e[0]),
                    fXmax=float(e[-1]),
                    fXbins=None if uniform else np.asarray(e, dtype=">f8"),
                )
            )
        data = np.asarray(self.sumw.ravel(order="F"), dtype=dtype)
        sumw2 = np.asarray(self.sumw2.ravel(order="F"), dtype=">f8")
        s = self.stats
        if self.ndim == 1:
            return _ident.to_TH1x(
                fName=self.name, fTitle=self.title, data=data, fEntries=float(self.entries),
                fTsumw=s[0], fTsumw2=s[1], fTsumwx=s[2], fTsumwx2=s[3], fSumw2=sumw2, fXaxis=axes[0],
            )
        return _ident.to_TH2x(
            fName=self.name, fTitle=self.title, data=data, fEntries=float(self.entries),
            fTsumw=s[0], fTsumw2=s[1], fTsumwx=s[2], fTsumwx2=s[3], fTsumwy=s[4], fTsumwy2=s[5],
            fTsumwxy=s[6], fSumw2=sumw2, fXaxis=axes[0], fYaxis=axes[1],
        )

    @classmethod
    def from_uproot(cls, h, name=None):
        edges = [ax.edges() for ax in h.axes]
        out = cls(*edges, name=name or h.member("fName"), title=h.member("fTitle"))
        out.sumw = np.asarray(h.values(flow=True), dtype=np.float64)
        var = h.variances(flow=True)
        out.sumw2 = np.asarray(var if var is not None else np.abs(out.sumw), dtype=np.float64)
        try:
            out.entries = float(h.member("fEntries"))
        except Exception:  # pragma: no cover
            pass
        return out

    @classmethod
    def read(cls, filename, key):
        with uproot.open(filename) as fh:
            if key not in fh:
                raise KeyError(f"Histogram '{key}' not found in {filename}. Available: {fh.keys()}")
            return cls.from_uproot(fh[key], name=key)


def write_root(filename, objects, kind="F"):
    """Write ``{key: Hist}`` into a new ROOT file (RECREATE)."""
    import os

    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
    with uproot.recreate(filename) as fout:
        for key, obj in objects.items():
            if isinstance(obj, Hist):
                fout[key] = obj.to_root(kind)
            else:
                fout[key] = obj
