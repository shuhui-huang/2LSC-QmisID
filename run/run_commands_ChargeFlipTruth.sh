
# for events in inclusive region: njets>=0
#root -l -b -q '../Scripts/ChargeFlipTruth.C("ttbar","_Default",false,false,"_nJets_ge0","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_ttbar_mc20_njets_ge0.log
# root -l -b -q '../Scripts/ChargeFlipTruth.C("Zjets","_Default",false,false,"_nJets_ge0","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_Zjets_mc20_njets_ge0.log

# for events in QmisID CR
#root -l -b -q '../Scripts/ChargeFlipTruthv5.C("ttbar","_Default",false,false,"_nJets_le1","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_ttbar_mc20_njets_le1.log
#root -l -b -q '../Scripts/ChargeFlipTruth.C("Zjets","_Default",false,false,"_nJets_le1","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_Zjets_mc20_njets_le1.log


# root -l -b -q '../Scripts/ChargeFlipTruthv5New.C("ttbar","_Default",true,false,"_nJets_le1","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_ttbar_mc20_njets_le1_weighted.log
root -l -b -q '../Scripts/ChargeFlipTruthv5New.C("Zjets","_Default",true,false,"_nJets_le1","l1_IDTight==1 && l2_IDTight==1","all",false)' | tee log_truth_Zjets_mc20_njets_le1_weighted.log
