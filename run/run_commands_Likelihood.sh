#root -l -b -q Likelihood.C | tee testLikelihood.log 
#root -l -b -q '../Scripts/Likelihoodv5.C("Data")' | tee log_Likelihood_Data.log 
#root -l -b -q '../Scripts/Likelihoodv5.C("Zjets")' | tee log_Likelihood_ZjetsMC.log 
# root -l -b -q '../Scripts/Likelihoodv5.C("ttbar")' | tee log_Likelihood_ttbarMC.log 
root -l -b -q '../Scripts/Likelihoodv5New.C("Zjets")' | tee log_Likelihoodv5New_ZjetsMC.log 
