#!/usr/bin/env bash
# Run the original ROOT macros and the Python port (legacy config) on the same
# ntuples and compare every output histogram.
#
#   bash validation/run_root_crosscheck.sh NTUPLE_DIR LIST_DIR [WORKDIR]
#
# NTUPLE_DIR : directory with the ROOT ntuples (the macros' hard-coded EOS path)
# LIST_DIR   : directory with Data.list, ZjetsMC.list, ttbarMC.list
# WORKDIR    : scratch area (default: results_crosscheck)
#
# Needs `root` and the qmisid package in the environment. Example on toys:
#   qmisid make-toys toys && PYTHON=python3 bash validation/run_root_crosscheck.sh toys toys
set -euo pipefail
PYTHON=${PYTHON:-python3}   # interpreter with qmisid + uproot installed

NT=$(cd "$1" && pwd)/
LISTS=$(cd "$2" && pwd)
WORK=${3:-results_crosscheck}
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$WORK"/root/{Scripts,Inputlists_v5Run3,run}
WORK=$(cd "$WORK" && pwd)

# --- ROOT macros (paths patched, logic untouched) ---------------------------
cp "$REPO"/Scripts/ChargeFlipTruthv5New.C "$REPO"/Scripts/Likelihoodv5New.C "$WORK"/root/Scripts/
cp "$LISTS"/Data.list "$LISTS"/ZjetsMC.list "$LISTS"/ttbarMC.list "$WORK"/root/Inputlists_v5Run3/
sed -i "s#/eos/user/s/shuhui/public/Run3_HHML/v5_ntuples_wIFF/#${NT}#g" "$WORK"/root/Scripts/*.C
sed -i 's#"../run_v5New/Truth_MC_Weighted_Full/"#"Truth_MCRun3_Weighted_Full/"#' "$WORK"/root/Scripts/Likelihoodv5New.C
(
  cd "$WORK"/root/run
  for p in ttbar Zjets; do
    s=$(date +%s)
    root -l -b -q "../Scripts/ChargeFlipTruthv5New.C(\"$p\",\"_Default\",true,false,\"_nJets_le1\",\"l1_IDTight==1 && l2_IDTight==1\",\"all\",false)" > log_truth_$p.log 2>&1
    echo "ROOT truth $p: $(( $(date +%s) - s )) s"
  done
  for p in Data Zjets; do
    s=$(date +%s)
    root -l -b -q "../Scripts/Likelihoodv5New.C(\"$p\")" > log_lh_$p.log 2>&1
    echo "ROOT likelihood $p: $(( $(date +%s) - s )) s"
  done
)

# --- Python, legacy configuration --------------------------------------------
PY="$WORK"/python
SET=(--set "samples.Data.path=${NT}" --set "samples.Data.filelist=${LISTS}/Data.list"
     --set "samples.Zjets.path=${NT}" --set "samples.Zjets.filelist=${LISTS}/ZjetsMC.list"
     --set "samples.ttbar.path=${NT}" --set "samples.ttbar.filelist=${LISTS}/ttbarMC.list"
     --set "truth.output_dir=${PY}/Truth_MCRun3_Weighted_Full"
     --set "likelihood.output_dir_data=${PY}/results_Run3/Likelihood_Data_v5/Data"
     --set "likelihood.output_dir_mc=${PY}/results_Run3/Likelihood_MC_v5"
     --set "likelihood.reference.file=${PY}/Truth_MCRun3_Weighted_Full/ttbar_Default_nJets_le1.root")
for p in ttbar Zjets; do
  s=$(date +%s); "$PYTHON" -m qmisid truth -c "$REPO"/configs/legacy_v5New.yaml -s $p "${SET[@]}" > "$WORK"/py_truth_$p.log 2>&1
  echo "Python truth $p: $(( $(date +%s) - s )) s"
done
for p in Data Zjets; do
  s=$(date +%s); "$PYTHON" -m qmisid likelihood -c "$REPO"/configs/legacy_v5New.yaml -s $p "${SET[@]}" > "$WORK"/py_lh_$p.log 2>&1
  echo "Python likelihood $p: $(( $(date +%s) - s )) s"
done

# --- comparison ----------------------------------------------------------------
CMP="$PYTHON $REPO/validation/compare_root_python.py"
R="$WORK"/root/run
status=0
for p in ttbar Zjets; do
  echo "=== truth $p (TH2F single-precision storage in ROOT -> rtol 5e-5)"
  $CMP "$R"/Truth_MCRun3_Weighted_Full/${p}_Default_nJets_le1.root "$PY"/Truth_MCRun3_Weighted_Full/${p}_Default_nJets_le1.root --rtol 5e-5 | grep -v "^only" || status=1
done
for p in Data Zjets; do
  if [ $p = Data ]; then rf="$R"/results_Run3/Likelihood_Data_v5/Data/Data_Default_nJets_le1.root; pf="$PY"/results_Run3/Likelihood_Data_v5/Data/Data_Default_nJets_le1.root
  else rf="$R"/results_Run3/Likelihood_MC_v5/Zjets_Default_nJets_le1.root; pf="$PY"/results_Run3/Likelihood_MC_v5/Zjets_Default_nJets_le1.root; fi
  echo "=== likelihood $p: Mll histograms"
  $CMP "$rf" "$pf" --keys SS OS --rtol 5e-5 || status=1
  echo "=== likelihood $p: rates (TMinuit vs Minuit2: |diff| < 0.01 sigma, MIGRAD errors within 10%)"
  $CMP "$rf" "$pf" --keys ${p}_Default_nJets_le1_all_num --pull 0.01 --etol 0.1 || status=1
done
exit $status
