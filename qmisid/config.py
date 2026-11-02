"""YAML configuration handling.

A config may inherit from another one with ``base: other.yaml`` (path relative
to the file). Dictionaries are deep-merged; lists and scalars are replaced.
Command-line ``--set a.b.c=value`` overrides are applied last.
"""

from __future__ import annotations

import copy
import os

import yaml


def deep_merge(base, override):
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(path, overrides=None):
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}
    if "base" in cfg:
        base_path = cfg.pop("base")
        if not os.path.isabs(base_path):
            base_path = os.path.join(os.path.dirname(os.path.abspath(path)), base_path)
        cfg = deep_merge(load_config(base_path), cfg)
    if overrides:
        cfg = deep_merge(cfg, overrides)
    cfg.setdefault("_source", os.path.abspath(path))
    return cfg


def parse_set_options(items):
    """Turn ``["likelihood.minuit.errordef=0.5", ...]`` into a nested dict."""
    out = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"--set expects key=value, got {item!r}")
        key, val = item.split("=", 1)
        node = out
        parts = key.split(".")
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        node[parts[-1]] = yaml.safe_load(val)
    return out


def get_sample(cfg, name):
    samples = cfg.get("samples", {})
    if name not in samples:
        raise KeyError(f"Sample '{name}' not defined in config; available: {list(samples)}")
    s = dict(samples[name])
    s.setdefault("name", name)
    s.setdefault("is_data", False)
    return s


def section(cfg, name, sample=None):
    """Config section ``name`` with optional per-sample overrides
    (``<section>.per_sample.<sample>``) merged in, plus global defaults."""
    sec = copy.deepcopy(cfg.get(name, {}))
    per = sec.pop("per_sample", {}) or {}
    if sample is not None and sample in per:
        sec = deep_merge(sec, per[sample])
    for key in ("binning", "leptons", "tree", "step_size"):
        if key not in sec and key in cfg:
            sec[key] = copy.deepcopy(cfg[key])
        elif key in sec and isinstance(sec[key], dict) and isinstance(cfg.get(key), dict):
            sec[key] = deep_merge(cfg[key], sec[key])
    return sec


def run_name(sample_name, cfg):
    """Output name, e.g. ``Data_Default_nJets_le1`` like the ROOT macros."""
    return f"{sample_name}{cfg.get('tag', '')}"


def output_path(cfg, *parts):
    base = cfg.get("output_dir", "results_py")
    return os.path.join(base, *parts)
