"""Chunked, columnar reading of the input TTrees with uproot."""

from __future__ import annotations

import logging
import os
import time

import numpy as np
import uproot

from . import expr

log = logging.getLogger(__name__)


def read_filelist(listfile, path=""):
    """Read a ``.list`` file (one ROOT file per line, ``#`` comments allowed).

    Relative entries are prefixed with ``path`` exactly like the ROOT macros
    (``m_path + file``). A list entry may also be a glob pattern.
    """
    import glob

    files = []
    with open(listfile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            full = line if (os.path.isabs(line) or "://" in line or not path) else os.path.join(path, line)
            matches = sorted(glob.glob(full)) if any(c in full for c in "*?[") else [full]
            files.extend(matches)
    if not files:
        raise FileNotFoundError(f"No input files found in list {listfile}")
    return files


def resolve_inputs(sample_cfg):
    """Return the list of ROOT files for a sample config block.

    ``sample_cfg`` may define ``files`` (list) and/or ``filelist`` (path to a
    ``.list`` file), plus an optional ``path`` prefix.
    """
    path = sample_cfg.get("path", "") or ""
    files = []
    for f in sample_cfg.get("files", []) or []:
        files.append(f if (os.path.isabs(f) or "://" in f or not path) else os.path.join(path, f))
    if sample_cfg.get("filelist"):
        files.extend(read_filelist(sample_cfg["filelist"], path))
    if not files:
        raise ValueError("Sample has neither 'files' nor 'filelist'")
    return files


def needed_branches(expressions):
    out = set()
    for e in expressions:
        if e is None:
            continue
        out |= expr.variables(e)
    return out


def iterate(files, tree, expressions, step_size="50 MB", max_events=None, report_every=10.0):
    """Yield ``(arrays, n)`` chunks with every branch referenced by ``expressions``.

    Only the branches actually used are read. ``arrays`` maps branch name to a
    numpy array (flat branches) or an awkward array (jagged branches).
    """
    branches = sorted(needed_branches(expressions))
    seen = 0
    t0 = time.time()
    last = t0
    for fname in files:
        with uproot.open(fname) as fh:
            if tree not in fh:
                raise KeyError(f"Tree '{tree}' not found in {fname}; available: {fh.keys()}")
            t = fh[tree]
            missing = [b for b in branches if b not in t.keys()]
            if missing:
                raise KeyError(f"Branches {missing} not found in {fname}:{tree}")
            nent = t.num_entries
            log.info("Opening %s:%s with %d entries", fname, tree, nent)
            stop = None
            if max_events is not None:
                if seen >= max_events:
                    break
                stop = min(nent, max_events - seen)
            for chunk in t.iterate(branches, step_size=step_size, entry_stop=stop, library="ak", how=dict):
                n = len(chunk[branches[0]]) if branches else 0
                arrays = {}
                for k, v in chunk.items():
                    a = v.to_numpy() if v.ndim == 1 else v
                    if isinstance(a, np.ndarray):
                        if a.dtype == np.bool_:
                            a = a.astype(np.int8)
                        elif a.dtype.kind == "f" and a.dtype != np.float64:
                            a = a.astype(np.float64)  # TTreeFormula works in double
                    arrays[k] = a
                seen += n
                now = time.time()
                if now - last > report_every:
                    log.info("  processed %d events (%.0f ev/s)", seen, seen / (now - t0))
                    last = now
                yield arrays, n
    dt = time.time() - t0
    log.info("Read %d events in %.1f s (%.0f ev/s)", seen, dt, seen / dt if dt > 0 else 0.0)
