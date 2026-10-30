# 2LSC-QmisID: electron charge-flip (QmisID) background estimation

## Table of contents
1. [Introduction](#introduction)
2. [Python framework (`qmisid`)](#python)
   1. [Setup](#py-setup)
   2. [Quick start on toy ntuples](#py-quickstart)
   3. [Running on the v5 ntuples](#py-run)
   4. [Configuration](#py-config)
   5. [Outputs](#py-outputs)
   6. [Validation against the ROOT macros](#py-validation)
   7. [Differences with respect to the ROOT macros](#py-differences)
3. [Systematic uncertainties](#systematics)
4. [Legacy ROOT macros](#legacy)

# Introduction <a name="introduction"></a>

The charge-flip ("QmisID") background in same-sign dielectron regions is
estimated from opposite-sign events weighted with per-electron charge-flip
rates. The rates are binned in |eta| x pT and measured

1. in MC with **truth matching** (IFF class of the electrons),
2. in data (and MC) with the **likelihood method** on Z->ee events,

and validated with a **closure test** (OS events weighted with the rates vs.
SS events). Step 4 assigns **systematic uncertainties** to the rates.

The original ROOT macros are in `Scripts/` (see [Legacy ROOT macros](#legacy)).
The `py-dev` branch adds a columnar Python port with all their functionality,
plus the systematics step.

# Python framework (`qmisid`) <a name="python"></a>

* Columnar: [uproot](https://github.com/scikit-hep/uproot5) reads only the
  needed TTree branches in chunks; selections and histograms are computed on
  whole numpy/awkward arrays. There is no event loop.
* The selection strings of the macros (`TTreeFormula` syntax, e.g.
  `(fabs(l1_eta)<1.37 || fabs(l1_eta)>1.52) && l1_IFF==3`) are used unchanged:
  `qmisid/expr.py` parses and evaluates them on arrays.
* The likelihood fit uses [iminuit](https://github.com/scikit-hep/iminuit)
  (Minuit2) with the same settings as the macro (parameter limits, steps,
  strategy).
* One YAML file (`configs/default.yaml`) replaces the hard-coded settings, so
  new ntuple versions, binnings or selections need no code changes.
* Output ROOT files keep the object names of the macros (for example
  `<name>_all_num` holds the rate map). Old and new files can be swapped,
  including as input of the ROOT closure macro.

| step | command | ROOT macro |
|---|---|---|
| 1. truth rates (MC) | `qmisid truth -c CFG -s Zjets -s ttbar` | `ChargeFlipTruthv5New.C` |
| 2. likelihood rates | `qmisid likelihood -c CFG -s Data` (or `-s Zjets`) | `Likelihoodv5New.C` |
| 3. closure / apply | `qmisid closure -c CFG -s Data --rates FILE[:HIST]` | `ClosureTestv5New.C` |
| 4. systematics | `qmisid syst -c CFG` | (new) |
| plots | `qmisid plot FILE[:HIST] [FILE2 ...] -o OUT` | `PlotRates.C`, `PlotRatesCompare.C` |
| toy ntuples | `qmisid make-toys DIR` | |

Package layout:

```
qmisid/
  expr.py         ROOT cut-string parser/evaluator (numpy/awkward)
  io.py           file lists, chunked uproot reading
  hist.py         ROOT-convention weighted TH1/TH2 (sumw2, flow bins, TH1F/TH2F writing)
  truth.py        step 1
  likelihood.py   step 2 (pair-cell counting, NLL, iminuit fit, bootstrap)
  closure.py      step 3 (OS weighting, rate bands, covariance propagation)
  systematics.py  step 4
  plotting.py     matplotlib/mplhep plots
  toys.py         toy ntuples with known rates
  cli.py          `qmisid` command
configs/          default.yaml (recommended), legacy_v5New.yaml (macro-identical), toys.yaml
validation/       ROOT <-> Python cross-check scripts
docs/             systematics.md (strategy for step 4)
tests/            pytest suite
```

## Setup <a name="py-setup"></a>

Python >= 3.9 with numpy, awkward, uproot, iminuit, pyyaml, matplotlib and
mplhep. ROOT is **not** needed.

```bash
git clone https://github.com/shuhui-huang/2LSC-QmisID.git
cd 2LSC-QmisID && git checkout py-dev
python3 -m venv venv && source venv/bin/activate     # or an LCG view with the packages above
pip install -e .          # installs the `qmisid` command
pytest                    # optional: 30 tests, about 5 s
```

## Quick start on toy ntuples <a name="py-quickstart"></a>

The toy ntuples have the v5 `qmisid_cr` format and known, piecewise-constant
charge-flip rates. They are useful to try the chain and for tests.

```bash
qmisid make-toys toys                                   # Data, Z+jets, ttbar toys + .list files
qmisid truth      -c configs/toys.yaml -s Zjets -s ttbar
qmisid likelihood -c configs/toys.yaml -s Data
qmisid closure    -c configs/toys.yaml -s Data \
       --rates results_toys/Likelihood_Data/Data_Default_nJets_le1.root
qmisid syst       -c configs/toys.yaml --set systematics.stat.n_bootstrap=50
qmisid plot results_toys/Likelihood_Data/Data_Default_nJets_le1.root \
            results_toys/Truth_MC_Weighted_Full/Zjets_Default_nJets_le1.root \
            --labels "Data likelihood" "Z+jets truth" -o results_toys/plots/compare
```

## Running on the v5 ntuples <a name="py-run"></a>

`configs/default.yaml` points to `/eos/user/s/shuhui/public/Run3_HHML/v5_ntuples_wIFF/`
and to `Inputlists/{Data,ZjetsMC,ttbarMC}.list` (one file name per line,
relative to `path`). Update the lists for the v5 files (the macros used
`../Inputlists_v5Run3/`), or override them from the command line:

```bash
CFG=configs/default.yaml
# 1) truth-matched rates (also used as starting values of the likelihood fit)
qmisid truth -c $CFG -s Zjets -s ttbar
# 2) likelihood rates in data and in Z+jets MC
qmisid likelihood -c $CFG -s Data -s Zjets
# 3) closure in data (and MC with data or MC rates)
qmisid closure -c $CFG -s Data  --rates results_py/Likelihood_Data/Data_Default_nJets_le1.root
qmisid closure -c $CFG -s Zjets --rates results_py/Likelihood_MC/Zjets_Default_nJets_le1.root
# 4) systematics (data pass + Z+jets and ttbar MC passes), then propagate to the yields
qmisid syst -c $CFG --set systematics.stat.n_bootstrap=100
qmisid closure -c $CFG -s Data \
    --rates results_py/Systematics/Data_Default_nJets_le1_syst.root:Data_Default_nJets_le1_all_num_stat \
    --syst  results_py/Systematics/Data_Default_nJets_le1_syst.root -o results_py/Closure_syst
```

Useful options: `--set key.sub=value` overrides any config entry, for example
`--set likelihood.model=exact` or `--set samples.Data.filelist=my.list`.
`-n N` processes only N events (for debugging), and `-v` gives verbose logs.

## Configuration <a name="py-config"></a>

`configs/default.yaml` is commented. Its main blocks are:

* `binning`, `leptons`: |eta| and pT edges, and the per-lepton expressions
  `[leading, sub-leading]`.
* `samples`: `path`, `filelist` and/or `files`, `is_data`, `weight` (MC weight
  expression).
* `truth`: `baseline`, `extra_cut`, `wrong_track`, `photon_conversion`,
  `real` (per lepton), `same_charge`, `split_charge` (e-/e+).
* `likelihood`: `selection`, `same_sign`, `opposite_sign`, `mll`,
  `mll_window: [low, low_Z, up_Z, up]`, `subtract_background`,
  `sideband_scale`, `model` (`approx`/`exact`), `bin_convention`,
  `skip_eta_bins`, `reference` (starting values), `minuit`, `n_bootstrap`.
* `closure`: `selection`, `mll_window`, `sidebands`/`subtract_background`,
  `variables` (name, expression, binning), `merge_overflow`.
* `systematics`: see [docs/systematics.md](docs/systematics.md).

Any block accepts `per_sample: {<sample>: {...}}` overrides (for example, no Mll
window for ttbar in the closure). A config can inherit from another one with
`base: other.yaml`.

## Outputs <a name="py-outputs"></a>

| step | file | content |
|---|---|---|
| truth | `Truth_MC_Weighted_Full/<name>.root` | the 18 objects of the macro (same names: `<name>_all_num` = total rate, `<name>_wrongTrack_num1`, ..., `m_histo_all_num`), plus raw numerators and denominators in `raw/`; `<name>.json` |
| likelihood | `Likelihood_{Data,MC}/<name>.root` | `<name>_all_num` (TH2F rates), `SS`/`OS` (Mll), `<name>_covariance`; `<name>.json` (fit summary), `<name>_pairs.npz` (pair-cell counts, for refits without rereading) |
| closure | `Closure/closure_<sample>.root` | per variable: `h_<s>_<var>_SS`, `_OS`, `_OSUp`, `_OSDown`, `_OS_statcov`, `Ratio_<var>_OS_Standard`; plots `h_<s>_<var>.{pdf,png}`; yields in `closure_<s>_yields.json` |
| systematics | `Systematics/<name>_syst.{root,json,md}` | see [docs/systematics.md](docs/systematics.md) |

## Validation against the ROOT macros <a name="py-validation"></a>

`validation/run_root_crosscheck.sh NTUPLE_DIR LIST_DIR` runs the original
`ChargeFlipTruthv5New.C` and `Likelihoodv5New.C` (with only the input paths
patched) and the Python port with `configs/legacy_v5New.yaml` on the same
ntuples, then compares every histogram (`validation/compare_root_python.py`).
It needs ROOT, so run it on lxplus with the real ntuples. It was run here with
ROOT 6.40/04 on toy ntuples (2.0M data, 1.5M Z+jets, 0.3M ttbar events):

| output | agreement |
|---|---|
| truth, all 18 histograms (ttbar, Z+jets) | relative difference <= 1.6e-5 (ROOT stores TH2F contents in single precision) |
| likelihood `SS`/`OS` Mll histograms, data | identical |
| likelihood `SS`/`OS` Mll histograms, weighted MC | <= 1e-5 relative (TH1F single precision) |
| likelihood rates, data and Z+jets | <= 0.004 sigma (TMinuit vs Minuit2 convergence tolerance); -ln L at both minima equal to 1e-7 |
| likelihood rate errors (MIGRAD, no HESSE) | <= 0.2 % (data), <= 5 % in one poorly constrained MC bin |

The closure macro depends on private `/afs` tools (`MiniTreeAnalyzer`, ...),
so it cannot run outside lxplus. Its logic is instead checked by a test that
reimplements the macro's event loop line by line
(`tests/test_pipeline.py::test_closure_matches_macro_event_loop`).

Timing on toys with 10M data and 5M MC events (2-CPU cloud machine; Python
times include about 1 s of start-up):

| step | ROOT macro | Python |
|---|---|---|
| truth, 5M MC events | 58 s | 2-3 s |
| likelihood, 10M data events | 5 s | 3 s |
| likelihood, 5M MC events | 3 s | 2 s |
| systematics: all 6 sources (6 data refits from one data pass, Z+jets and ttbar MC passes) | n/a | 10 s |
| same, with 20 correlated bootstrap replicas | n/a | 19 s |

The truth step gains the most: the macro runs 30 `TTree::Draw` passes over the
tree, while the port reads it once.

## Differences with respect to the ROOT macros <a name="py-differences"></a>

`configs/legacy_v5New.yaml` reproduces the macros, including the issues below.
`configs/default.yaml` fixes them and is recommended for new results.

| # | macro behaviour | effect | default.yaml |
|---|---|---|---|
| 1 | truth crack veto on **signed** eta: `(l1_eta<1.37 \|\| l1_eta>1.52)` | electrons with -1.52 < eta < -1.37 pass: non-zero crack bin, inconsistent with the likelihood selection; also no \|eta\| < 2.47 cut | `fabs(eta)`, plus \|eta\| < 2.47 |
| 2 | truth pT binning ends at 2500 GeV, likelihood at 2600 GeV | electrons with 2500 < pT < 2600 GeV are missing from the truth rates | 2600 everywhere |
| 3 | `SetErrorDef(1)` for a -ln L | likelihood rate errors too large by sqrt(2) | `errordef: 0.5`, plus HESSE |
| 4 | pair-cell counts stored as `float` | above 2^23 (8.4M) entries in a cell, sideband weights of -0.5 are rounded; above 2^24, increments of 1 are lost | always double precision |
| 5 | pT <= 20 GeV or \|eta\| = 0 give bin index -1 | out-of-bounds vector access (undefined behaviour) | events dropped and counted in the log |
| 6 | default starting values filled pT-major into an eta-major parameter vector | wrong start values (only without a reference file) | filled eta-major |
| 7 | closure: no tight-ID, \|eta\| or crack selection (`l1_Tight` is read but not used); no background subtraction of the SS events in the Z window | data closure is not like-for-like with the rate measurement, which likely contributes to SS/OS = 1.9 in `run/val_v5_Data_DataLike_all.log` | same selection as the likelihood, plus sideband subtraction |
| 8 | closure: `EtaLep1new`/`EtaLep2new` filled with \|eta\| (title matching on `Contains("EtaLep1")`) | the signed-eta plots show \|eta\| | signed eta |
| 9 | closure: overflow merged into SS and OS but not into OSUp/OSDown | rate band and yield +- wrong in the last bin | merged for all |
| 10 | closure band: all rates shifted coherently by +-1 sigma | ignores the fit correlations; conservative (3.6x the covariance result on the toys) | also covariance propagation (`_OS_statcov`, `rate_stat_cov`) |

Numerical conventions are kept in both configs: the likelihood bins are
`(lo, hi]` with pT overflow in the last bin, the closure uses ROOT `FindBin`,
and the MC weights are applied the same way. In legacy mode, kinematics are
rounded to `float` as in the macro.

# Systematic uncertainties <a name="systematics"></a>

The strategy is in [docs/systematics.md](docs/systematics.md). In short:

* **stat**: HESSE covariance of the data fit, cross-checked with a Poisson
  bootstrap that includes the sideband weights; also provided as eigen-variations;
* **mll_window**, **bkg_norm**, **lh_model**: refits of the data with a
  different Z window/sidebands, +-50 % sideband normalisation, and the exact
  pair probability. All are filled in one pass over the data, and a correlated
  bootstrap removes their statistical noise;
* **non_closure**: SS observed vs. OS-weighted prediction per rate bin in the fit region;
* **method_bias**: likelihood vs. truth rates in Z+jets MC;
* **process_dependence**: ttbar vs. Z+jets truth rates, independent of the
  kinematic composition.

Sources are symmetrised and added in quadrature. The output contains the total,
per-source varied rate maps (for one nuisance parameter per source) and stat.
eigen-variations. `qmisid closure --syst` propagates all of them to the
predicted yields.

# Legacy ROOT macros <a name="legacy"></a>

The scripts used for QmisID rate estimation with ROOT are in `Scripts/`.
- To get the MC truth-matched rate, run `Scripts/ChargeFlipTruthv5New.C`
- To get the likelihood rate, run `Scripts/Likelihoodv5New.C`
- Closure test: `Scripts/ClosureTestv5New.C`
- Plots: `Scripts/PlotRates.C`, `Scripts/PlotRatesCompare.C`

Example commands are in `run/`:
- `run/run_commands_ChargeFlipTruth.sh`
- `run/run_commands_Likelihood.sh`
- `run/run_commands_plots.sh`

Setup for the ROOT macros (first time):
```
setupATLAS
lsetup git
git clone ssh://git@gitlab.cern.ch:7999/shuhui/QmisID.git
cd QmisID/
source setup.sh
```
In a fresh terminal, run `source setup.sh` again.
