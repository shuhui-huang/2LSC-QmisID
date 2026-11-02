# Systematic uncertainties on the charge-flip rates: strategy

This note describes what `qmisid syst` (`qmisid/systematics.py`) computes and why.
The configuration lives in the `systematics:` block of `configs/default.yaml`.

## What is being estimated

The QmisID background in a same-sign (SS) region is predicted by weighting
opposite-sign (OS) events with `w = P/(1-P)`, `P = e1 + e2 - 2 e1 e2`, where
`e1, e2` are the charge-flip rates of the two electrons, binned in |eta| x pT and
measured in data with the Z->ee likelihood fit. The rates carry two kinds of
uncertainty:

1. **Statistical**: from the finite Z->ee sample (and the sideband subtraction).
2. **Systematic**: from how the rates are measured (background, mass window,
   likelihood model, method bias), and from transferring them from the Z->ee
   control region to the analysis region (process dependence, closure).

For every source `s`, the script computes an absolute shift `delta_s(eta, pT)`
of the nominal data rate. Two-sided variations are symmetrised, and sources are
added in quadrature:

    total(eta, pT) = sqrt( stat^2 + sum_s delta_s^2 )

## Sources

| Source (config key) | Type | What is varied / compared | Default |
|---|---|---|---|
| `stat` | fit | HESSE covariance of the nominal fit; optional Poisson bootstrap | HESSE |
| `mll_window` | `refit` | Z peak and sidebands: wide `[56,76,106,126]`, narrow `[66,86,96,116]` (nominal `[61,81,101,121]`) | envelope |
| `bkg_norm` | `refit` | sideband-subtraction normalisation +-50 % (scale 0.75 / 0.25, nominal 0.5) | envelope |
| `lh_model` | `refit` | likelihood `N(e1+e2)` (macro) vs. exact `N(e1+e2-2e1e2)` | one-sided |
| `non_closure` | `closure` | data, fit region: SS observed vs. OS-weighted prediction, per rate bin | beyond stat. |
| `method_bias` | `mc_method` | Z+jets MC: likelihood rates vs. truth-matched rates of the same events | beyond stat. |
| `process_dependence` | `process` | truth flips in ttbar MC vs. flips predicted for the same electrons with the Z+jets truth rates | inclusive, abs |

### Statistical (`stat`)

* **HESSE** (default): the inverse Hessian of the -ln L at the minimum, with the
  correct `errordef = 0.5` (the macro used 1, which inflated the errors by
  sqrt(2)).
* **Bootstrap** (`stat.n_bootstrap: N`): every event gets a Poisson(1) weight per
  replica, all pair-cell counts are refilled and the fit is repeated. Unlike
  HESSE, this includes the effect of the negative sideband weights (the variance
  of a subtracted count is `N_peak + scale^2 N_side`, not `N_peak - scale N_side`)
  and of MC event weights. With `stat.method: bootstrap` the bootstrap spread is
  used as the stat. uncertainty. On the toys, bootstrap/HESSE is 0.8-1.3 per
  bin (40 replicas).
* The stat. covariance is also written as **eigen-variations**
  (`<name>_rates_stat_eig<k>_up/down`): uncorrelated nuisance parameters that
  reproduce the full covariance (including the bin-to-bin anti-correlations from
  the fit), ready for a profile-likelihood fit.

### Refit variations (`type: refit`)

Each variation is an override of the `likelihood:` block (any key: `mll_window`,
`sideband_scale`, `subtract_background`, `model`, `selection`, `binning` ...), so
new variations are added in the config without code changes (e.g. a tighter
ID, a pile-up split, a data-period split).

* The ntuples are read **once**: variations that only change post-processing
  quantities (`model`, `sideband_scale`, `subtract_background`, `minuit`) reuse
  the nominal pair counts, because peak and sideband counts are stored
  separately. Variations that change the event selection or the mass window
  get their own accumulator, filled from the same chunks.
* A narrower window removes events, so part of the shift is statistical noise.
  With bootstrap replicas enabled, the **same Poisson weights** are applied to
  the nominal and to every variation, so the statistical uncertainty of the
  *shift* is known (correlated bootstrap). With `mode: beyond_stat` only the part
  of the shift beyond its stat. uncertainty is kept:
  `delta = sqrt(max(shift^2 - sigma_shift^2, 0))` (a Barlow-type test). Without
  bootstrap the full shift is used (`mode: abs`).
* `symmetrize: envelope` takes the largest |shift| of the variations;
  `half_diff` takes |up - down|/2.

### Non-closure in the fit region (`type: closure`)

The likelihood assumes that the SS probability of a pair factorises into
per-electron rates in the chosen bins. Its validity is checked on the data used
for the fit: for each rate bin k, all pairs with at least one electron in k are
selected and the observed (sideband-subtracted) SS yield is compared with the
OS-weighted prediction. This uses the pair-cell counts, so no second pass over
the data is needed. The relative non-closure `r_k` beyond its statistical
uncertainty becomes the uncertainty of bin k. A chi2 over all cells is also
reported as a goodness-of-fit measure.

### Method bias (`type: mc_method`)

On Z+jets MC the likelihood fit is run exactly as on data and compared with the
truth-matched rates (IFF class) of the same events (same selection and binning,
truth restricted to the fit window). The relative difference
`r = e_LH / e_truth - 1` measures the bias of the method itself (e.g. rate
variations inside a bin that the pair model cannot absorb, residual
background, fakes); it is applied to the data rate as `|r| e_data`. The MC
statistics of both estimates enter `sigma_r`, and with `mode: beyond_stat` only
the significant part is kept. With the binned toy rates this source is zero;
with rates that vary strongly inside the bins (`qmisid make-toys --rate-model
smooth`, 1.5M Z+jets MC events) it picks up 8-12 % in two central bins. Its
sensitivity is limited by the MC statistics.

### Process dependence (`type: process`)

The rates are measured in Z->ee but applied to SS regions whose electrons may
come from ttbar, ttV, ... (more hadronic activity, different isolation and
material interactions). The truth-matched flips in the target MC (ttbar) are
compared with the flips expected for the same electrons from the Z+jets truth
rates:

    r = sum N_flip(ttbar) / sum_bins [ N_e(ttbar) x e_truth(Z) ] - 1

Because the prediction uses the target sample's own electrons, differences in
the pT/eta composition cancel and only genuine rate differences remain.
`granularity` sets where the sums run: per bin (`bin`), per pT bin (`pt`), per
|eta| bin (`eta`) or globally (`inclusive`). The ttbar truth sample in the
`qmisid_cr` tree has very few flips (~10 weighted events in the v5 log), so the
default is `inclusive` with `mode: abs`. When more ttbar MC is available
(e.g. the inclusive tree, `truth.per_sample.ttbar.tree`), switch to `pt`.

### Not included by default (easy to add as `refit` / `external` sources)

* electron vs. positron rates (`truth.split_charge`, needs charge branches);
* dependence on N_jets / HT / pile-up (a `refit` with an extra `selection`);
* ECIDS / BDT working point;
* an `external` relative uncertainty map from a ROOT file.

## Correlations and how to use the output

`results_py/Systematics/<name>_syst.root` contains:

| object | content |
|---|---|
| `<name>_all_num` | nominal rates, error = total uncertainty (drop-in for the old files) |
| `<name>_all_num_stat`, `<name>_all_num_syst` | nominal rates with stat-only / syst-only errors |
| `<name>_syst_<src>` | nominal rates, error = delta of source `src` |
| `<name>_rates_<src>_up/down` | varied rate maps (nominal +- delta, clipped to [0, 1]) |
| `<name>_rates_stat_eig<k>_up/down` | stat. eigen-variations |
| `<name>_var_<src>_<variation>` | the raw rate map of each variation (diagnostics) |
| `<name>_covariance` | stat. covariance (eta-major parameter order) |

plus `<name>_syst.json` (all numbers), `<name>_syst.md` (table) and plots
(`_breakdown`: relative uncertainty per source; `_variations_pt<i>`: rate maps).

Recommended treatment in a fit: one nuisance parameter per systematic source,
fully correlated across (eta, pT) bins, and one per stat. eigen-variation. The
event-level propagation to any region is done by the closure/apply step:

    qmisid closure -c configs/default.yaml -s <sample> \
        --rates results_py/Systematics/Data_Default_nJets_le1_syst.root:Data_Default_nJets_le1_all_num_stat \
        --syst  results_py/Systematics/Data_Default_nJets_le1_syst.root

which fills the OS-weighted prediction for every `_rates_*_up/down` map; the
yields per source are in `closure_<sample>_yields.json`. On the toys, the
quadrature sum of the stat. eigen-variations on the yield equals the
covariance propagation (4.19 % vs 4.19 %).

## Caveats

* Shifts are evaluated with fits started from the nominal result; a variation
  whose fit is not valid is reported in the log.
* `beyond_stat` assumes Gaussian statistics; for bins with very few SS events
  prefer the bootstrap.
* Rates in the crack bin are not fitted and get no uncertainty.
* The weighted-MC likelihood fit uses HESSE errors that are only approximate for
  weighted events; the `method_bias` stat. uncertainty is therefore approximate
  (use the bootstrap for a precise estimate).
