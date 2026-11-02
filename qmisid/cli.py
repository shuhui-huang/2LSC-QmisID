"""Command line interface: ``qmisid <step> ...`` (or ``python -m qmisid``)."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from . import config as C


def _common(p, sample_required=True):
    p.add_argument("-c", "--config", required=True, help="YAML configuration")
    p.add_argument("-s", "--sample", action="append", required=sample_required,
                   help="sample name(s) from the config (repeatable)")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override a config value, e.g. --set likelihood.model=exact")
    p.add_argument("-n", "--max-events", type=int, default=None, help="process at most N events (debugging)")
    p.add_argument("-o", "--output", default=None, help="output file/directory (default from config)")


def _load(args):
    return C.load_config(args.config, C.parse_set_options(args.set))


def _parse_ref(text):
    if ":" in text and not text.endswith(".root"):
        f, h = text.rsplit(":", 1)
        return f, h
    return text, None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qmisid", description="Electron charge-flip (QmisID) rate framework")
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("truth", help="1) truth-matched rates in MC")
    _common(p)
    p = sub.add_parser("likelihood", help="2) likelihood rates from Z->ee (data or MC)")
    _common(p)
    p = sub.add_parser("closure", help="3) closure test / apply rates to OS events")
    _common(p)
    p.add_argument("--rates", required=True, help="rates ROOT file[:histogram]")
    p.add_argument("--alt", action="append", default=[], metavar="LABEL=FILE[:HIST]",
                   help="additional rate maps to propagate (e.g. systematic variations)")
    p.add_argument("--syst", default=None, help="<name>_syst.root from 'qmisid syst': propagate every source")
    p.add_argument("--no-plots", action="store_true")
    p = sub.add_parser("syst", help="4) systematic uncertainties on the data rates")
    _common(p, sample_required=False)
    p = sub.add_parser("plot", help="plot rate maps and comparisons")
    p.add_argument("inputs", nargs="+", help="FILE[:HIST] (default hist <stem>_all_num)")
    p.add_argument("--labels", nargs="*", default=None)
    p.add_argument("-o", "--output", required=True, help="output base name (no extension)")
    p.add_argument("--ref", type=int, default=-1, help="index of the reference map for ratios")
    p = sub.add_parser("make-toys", help="write toy ntuples with known rates")
    p.add_argument("outdir")
    p.add_argument("--n-data", type=int, default=2_000_000)
    p.add_argument("--n-z", type=int, default=1_500_000)
    p.add_argument("--n-tt", type=int, default=300_000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--rate-model", choices=["binned", "smooth"], default="binned")

    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s", datefmt="%H:%M:%S")
    for noisy in ("matplotlib", "fontTools", "PIL"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    if args.cmd == "make-toys":
        from . import toys

        lists = toys.make_all(args.outdir, args.n_data, args.n_z, args.n_tt, args.seed, args.rate_model)
        print("\n".join(f"{k}: {v}" for k, v in lists.items()))
        return 0

    if args.cmd == "plot":
        from . import plotting
        from .hist import Hist

        hists, labels = [], []
        for i, inp in enumerate(args.inputs):
            f, h = _parse_ref(inp)
            h = h or os.path.basename(f).replace(".root", "") + "_all_num"
            hists.append(Hist.read(f, h))
            labels.append(args.labels[i] if args.labels and i < len(args.labels) else h)
        if len(hists) == 1:
            plotting.rates_1d(hists[0], args.output + "_1D", title=labels[0])
            plotting.rates_2d(hists[0], args.output + "_2D", title=labels[0])
        else:
            plotting.compare_rates(hists, labels, args.output, ref=args.ref)
        return 0

    cfg = _load(args)

    if args.cmd == "truth":
        from . import truth

        for s in args.sample:
            truth.run(cfg, s, args.max_events, args.output if len(args.sample) == 1 else None)
    elif args.cmd == "likelihood":
        from . import likelihood

        for s in args.sample:
            likelihood.run(cfg, s, args.max_events, args.output if len(args.sample) == 1 else None)
    elif args.cmd == "closure":
        from . import closure

        f, h = _parse_ref(args.rates)
        alt = {}
        for a in args.alt:
            lab, ref = a.split("=", 1)
            alt[lab] = _parse_ref(ref)
        if args.syst:
            import uproot

            with uproot.open(args.syst) as fh:
                for k in fh.keys(cycle=False):
                    if "_rates_" in k and (k.endswith("_up") or k.endswith("_down")):
                        alt[k.split("_rates_", 1)[1]] = (args.syst, k)
        for s in args.sample:
            closure.run(cfg, s, f, h, args.max_events, args.output, alt=alt, plot=not args.no_plots)
    elif args.cmd == "syst":
        from . import systematics

        samples = args.sample or [None]
        for s in samples:
            systematics.run(cfg, s, args.max_events)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
