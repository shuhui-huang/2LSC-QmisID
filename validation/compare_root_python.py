#!/usr/bin/env python3
"""Compare histograms of two ROOT files (e.g. ROOT-macro vs Python outputs).

    python validation/compare_root_python.py ROOT_FILE PY_FILE [--rtol 1e-6] [--keys k1 k2]

For every histogram present in both files (TH1/TH2), the bin contents and
errors (visible + flow bins) are compared. Prints one line per histogram and
exits non-zero if any histogram differs beyond the tolerance.
"""

import argparse
import sys

import numpy as np
import uproot


def compare(fa, fb, rtol, atol, etol, keys=None, pull=None):
    ok = True
    rows = []
    with uproot.open(fa) as A, uproot.open(fb) as B:
        ka = {k.split(";")[0] for k, c in A.classnames().items() if c.startswith("TH")}
        kb = {k.split(";")[0] for k, c in B.classnames().items() if c.startswith("TH")}
        common = sorted(ka & kb) if keys is None else keys
        only_a, only_b = sorted(ka - kb), sorted(kb - ka)
        for k in common:
            ha, hb = A[k], B[k]
            va, vb = ha.values(flow=True), hb.values(flow=True)
            if va.shape != vb.shape:
                rows.append((k, "SHAPE", str(va.shape), str(vb.shape)))
                ok = False
                continue
            ea = np.sqrt(np.abs(ha.variances(flow=True))) if ha.variances(flow=True) is not None else np.zeros_like(va)
            eb = np.sqrt(np.abs(hb.variances(flow=True))) if hb.variances(flow=True) is not None else np.zeros_like(vb)
            dv = np.abs(va - vb)
            de = np.abs(ea - eb)
            sv = np.max(dv / np.maximum(np.abs(va), atol)) if dv.size else 0.0
            se = np.max(de / np.maximum(np.abs(ea), atol)) if de.size else 0.0
            if pull is not None:
                # fitted quantities: compare the difference in units of the uncertainty
                with np.errstate(divide="ignore", invalid="ignore"):
                    pulls = np.where(ea > 0, dv / ea, 0.0)
                sv = float(np.max(pulls)) if pulls.size else 0.0
                good = sv < pull and np.allclose(ea, eb, rtol=etol, atol=atol)
            else:
                good = np.allclose(va, vb, rtol=rtol, atol=atol) and np.allclose(ea, eb, rtol=etol, atol=atol)
            ok &= good
            rows.append((k, "OK" if good else "DIFF", f"{sv:.2e}", f"{se:.2e}"))
    w = max([len(r[0]) for r in rows] + [10])
    vlab = "max |diff|/error" if pull is not None else "max_rel_diff(values)"
    print(f"{'histogram':<{w}}  status  {vlab:>20}  max_rel_diff(errors)")
    for r in rows:
        print(f"{r[0]:<{w}}  {r[1]:<6}  {r[2]:>20}  {r[3]:>20}")
    if only_a:
        print("only in first :", only_a)
    if only_b:
        print("only in second:", only_b)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--rtol", type=float, default=1e-5, help="relative tolerance on values")
    ap.add_argument("--etol", type=float, default=None, help="relative tolerance on errors (default: rtol)")
    ap.add_argument("--atol", type=float, default=1e-12)
    ap.add_argument("--keys", nargs="*", default=None)
    ap.add_argument("--pull", type=float, default=None,
                    help="for fitted maps: max allowed |difference| in units of the bin error")
    a = ap.parse_args()
    ok = compare(a.a, a.b, a.rtol, a.atol, a.etol if a.etol is not None else a.rtol, a.keys, a.pull)
    print("RESULT:", "AGREE" if ok else "DIFFER")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
