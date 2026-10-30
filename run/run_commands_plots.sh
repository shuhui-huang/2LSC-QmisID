#root -l -b -q ../Scripts/PlotRates.C
#root -l -b -q ../Scripts/PlotRatesCompare.C


# ClosureTest(TString process,TString inputRates,TString MCversion,TString var,TString outputDir,TString TestHt, bool testRealSS, bool testplus, TString TestDep="")

# root -l -q -b '../Scripts/ClosureTestv5.C("Data","results/Likelihood_Data_v5/Data","Data","all","ValidationPlots_wDataLikelihodd",false,false,"")' | tee val_v5_Data_DataLike_all.log
# root -l -q -b '../Scripts/ClosureTestv5.C("Zjets","results/Likelihood_Data_v5/Data","Data","all","ValidationPlots_wDataLikelihodd",false,false,"")' | tee val_v5_Zjets_DataLike_all.log


# root -l -q -b '../Scripts/ClosureTest.C("Data","Likelihood_Data_test/Data","Data","all","ValidationPlots_wDataLikelihodd",false,false,"")' | tee val_Data_DataLike_all.log
#root -l -q -b '../Scripts/ClosureTest.C("Zjets","Likelihood_MC_test","Zjets","all","ValidationPlots_wMCLikelihood_test",false,false,"")' | tee val_Zjets_mcaLike_all_test1.log
#root -l -q -b '../Scripts/ClosureTest.C("Zjets","Truth_MC_Weighted_Full","Zjets","all","ValidationPlots_wMCLikelihood_test",false,false,"")' | tee val_Zjets_mcaLike_all_test1.log

root -l -q -b '../Scripts/ClosureTestv5New.C("Data","results/Likelihood_Data_v5/Data","Data","all","ValidationPlots_wDataLikelihodd",false,false,"")' | tee val_v5_Data_DataLike_all.log
