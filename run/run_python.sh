#!/usr/bin/env bash
# Python (qmisid) equivalent of run_commands_*.sh. See README.md.
set -e
CFG=${CFG:-configs/default.yaml}     # configs/legacy_v5New.yaml reproduces the ROOT macros
cd "$(dirname "$0")/.."

# 1) truth-matched rates in MC (ChargeFlipTruthv5New.C("ttbar"/"Zjets","_Default",true,false,"_nJets_le1",...))
qmisid truth -c $CFG -s Zjets -s ttbar
# 2) likelihood rates (Likelihoodv5New.C("Data") / ("Zjets"))
qmisid likelihood -c $CFG -s Data -s Zjets
# 3) closure tests (ClosureTestv5New.C("Data","results/Likelihood_Data_v5/Data","Data","all",...))
qmisid closure -c $CFG -s Data  --rates results_py/Likelihood_Data/Data_Default_nJets_le1.root
qmisid closure -c $CFG -s Zjets --rates results_py/Likelihood_MC/Zjets_Default_nJets_le1.root -o results_py/Closure_Zjets_MCrates
qmisid closure -c $CFG -s Zjets --rates results_py/Likelihood_Data/Data_Default_nJets_le1.root -o results_py/Closure_Zjets_Datarates
# 4) systematics on the data rates, then propagation to the closure yields
qmisid syst -c $CFG --set systematics.stat.n_bootstrap=100
qmisid closure -c $CFG -s Data \
  --rates results_py/Systematics/Data_Default_nJets_le1_syst.root:Data_Default_nJets_le1_all_num_stat \
  --syst results_py/Systematics/Data_Default_nJets_le1_syst.root -o results_py/Closure_syst
# plots (PlotRates.C / PlotRatesCompare.C)
qmisid plot results_py/Likelihood_Data/Data_Default_nJets_le1.root -o results_py/Plots/Data_Likelihood
qmisid plot results_py/Likelihood_Data/Data_Default_nJets_le1.root \
            results_py/Truth_MC_Weighted_Full/Zjets_Default_nJets_le1.root \
            --labels "Data likelihood" "Z+jets truth" -o results_py/Plots/DatavsMC
