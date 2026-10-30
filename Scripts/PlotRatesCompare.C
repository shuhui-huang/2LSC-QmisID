#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
//#include <iomanip>
//#include <vector>

#include "TPie.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/Yields.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/Channel.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/EventCut.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/VariableDistr2D.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/VariableDistr.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/PhysicsSample.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/PhysicsProcess.C"

#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/MiniTreeAnalyzer.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/MyStyle.C"

void PlotRate1D(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TString outName);
void PlotRate2D(TH2F *hRates,TString outName);
void GetHistPlot(TString FileName, TString TreeName, TString Version);
std::vector<TH1F*> GetRateHist(TString FileName, TString TreeName, TString Version);
void PlotHistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, TString name1, TString name2);
void Plot4HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, vector<TH1F*> v_hrate4, TString name1, TString name2, TString name3, TString name4);
void Plot2HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, TString name1, TString name2, TString meth1, TString meth2);
void Plot5HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, vector<TH1F*> v_hrate4,vector<TH1F*> v_hrate5, TString name1, TString name2, TString name3, TString name4, TString name5);
void Plot3HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, TString name1, TString name2, TString name3);

void PlotRatesCompare(){
    SetMyStyle();

    //std::vector<TString> data={"Data2015","Data2016","Data2017","Data2018"};
    std::vector<TString> data={"Data2015","Data2016","Data2017"};
    //std::vector<TString> mcversion={"mc16a","mc16d","mc16e"}; 
    //std::vector<TString> mcversion={"mc16e"}; 
    std::vector<TString> mcversion={"mc"}; 
  
    TString pathData="../run/results/Likelihood_Data_v5/Data/";
    TString pathmc="../run/results/Likelihood_MC_v5/";
    TString pathmctruth = "../run/Truth_MC_Weighted_Full/";

    // std::vector<TH1F*> v_hRate_mcLike  = GetRateHist(pathmc+"Zjets_Default_nJets_le1.root","Zjets_Default_nJets_le1_all_num","mc Likelihood");
    // std::vector<TH1F*> v_hRate_ttbarLike  = GetRateHist(pathmc+"ttbar_Default_nJets_le1.root","ttbar_Default_nJets_le1_all_num","mc Likelihood");
    //std::vector<TH1F*> v_hRate_mcTruth = GetRateHist(mctruthpath+"Zjets_Default_nJets_ge0.root","Zjets_Default_nJets_ge0_all_num","");
    // std::vector<TH1F*> v_hRate_dtLike = GetRateHist(pathData+"Data_Default_nJets_le1.root","Data_Default_nJets_le1_all_num","Data");
    // std::vector<TH1F*> v_hRate_dtLikeyulei = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/Data_V5_DNN_nJets_le1.root","Data_Default_nJets_le1_all_num","Data");
    // std::vector<TH1F*> v_hRate_dtLikev5New = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/run_v5New/results/Likelihood_Data_v5/Data/Data_Default_nJets_le1.root","Data_Default_nJets_le1_all_num","Data");
    std::vector<TH1F*> v_hRate_dtLikev5NewRun3 = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/run_v5New/results_Run3/Likelihood_Data_v5/Data/Data_Default_nJets_le1.root","Data_Default_nJets_le1_all_num","Data");
    // std::vector<TH1F*> v_hRate_mcLikev5New = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/run_v5New/results/Likelihood_MC_v5/Zjets_Default_nJets_le1.root","Zjets_Default_nJets_le1_all_num","Data");

    // std::vector<TH1F*> v_hRate_ZjetsTruth = GetRateHist(pathmctruth+"Zjets_Default_nJets_le1.root","Zjets_Default_nJets_le1_all_num","Data2015");
    // std::vector<TH1F*> v_hRate_ZjetsTruthv5New = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/run_v5New/Truth_MC_Weighted_Full/Zjets_Default_nJets_le1.root","Zjets_Default_nJets_le1_all_num","Data2015");
    std::vector<TH1F*> v_hRate_ZjetsTruthv5NewRun3 = GetRateHist("/afs/cern.ch/work/s/shuhui/Run3/HHML/QmisID/run_v5New/Truth_MCRun3_Weighted_Full/Zjets_Default_nJets_le1.root","Zjets_Default_nJets_le1_all_num","Data2015");
    //  std::vector<TH1F*> v_hRate_ttbarTruth = GetRateHist(pathmctruth+"ttbar_Default_nJets_le1.root","ttbar_Default_nJets_le1_all_num","Data2015");
    


    // mc vs data Likelihood
    // Plot2HistCompare (v_hRate_mcLike,v_hRate_dtLike,"mc16_all","Data_all","Likelihood","Likelihood");

    // different mc periods.
    // Plot4HistCompare (v_hRate_ZjetsTruth,v_hRate_ttbarTruth,v_hRate_mcLike,v_hRate_dtLike,"Zjets Truth-matched","ttbar Truth-matched","Zjets Likelihood","Data Likelihood");


    // PlotHistCompare(v_hRate_dtLike,v_hRate_mcLike, "Data", "MC");
    // PlotHistCompare(v_hRate_dtLikev5New, v_hRate_dtLike, "Data", "Data");
    // PlotHistCompare(v_hRate_dtLike, v_hRate_dtLikev5New, "Data", "Data");
    PlotHistCompare(v_hRate_dtLikev5NewRun3, v_hRate_ZjetsTruthv5NewRun3, "Data", "MC");
    // PlotHistCompare(v_hRate_mcLike, v_hRate_mcLikev5New, "MC", "MC");
    // PlotHistCompare(v_hRate_mcLike,v_hRate_ZjetsTruth , "MC", "MC");
    // PlotHistCompare(v_hRate_mcLikev5New,v_hRate_ZjetsTruthv5New , "MC", "MC");
    // PlotHistCompare(v_hRate_dtLikev5New,v_hRate_mcLikev5New , "Data", "MC");
    // PlotHistCompare(v_hRate_ttbarLike,v_hRate_ttbarTruth , "ttbar", "MC");


return;
}


void Plot5HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, vector<TH1F*> v_hrate4,vector<TH1F*> v_hrate5, TString name1, TString name2, TString name3, TString name4, TString name5){

    std::cout << "hrate1: " << name1 << " has size" << v_hrate1.size() << std::endl;
    std::cout << "hrate2: " << name2 << " has size" << v_hrate2.size() << std::endl;
    std::cout << "hrate3: " << name3 << " has size" << v_hrate3.size() << std::endl;
    std::cout << "hrate4: " << name4 << " has size" << v_hrate4.size() << std::endl;
    std::cout << "hrate5: " << name5 << " has size" << v_hrate5.size() << std::endl;
    
    TString pt1 = "#in [20, 50] GeV"; 
    TString pt2 = "#in [50, 100] GeV"; 
    TString pt3 = "#in [100, 200] GeV"; 
    TString pt4 = "> 200 GeV";
    std::vector<TString> v_ptbin = {pt1,pt2,pt3,pt4};

    for (unsigned int i=0; i<v_ptbin.size(); i++){
 
        TString ptrange = v_ptbin.at(i); 
        TH1F *h1 = v_hrate1.at(i);    
        TH1F *h2 = v_hrate2.at(i);    
        TH1F *h3 = v_hrate3.at(i);    
        TH1F *h4 = v_hrate4.at(i);    
        TH1F *h5 = v_hrate5.at(i);    

        TCanvas *c1 = new TCanvas("c", "canvas", 850, 810);
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.3);
        gStyle->SetErrorX(0.5);
        MiniTreeAnalyzer newanalyzer;
        h1->SetAxisRange(1e-5, 10, "Y");
        h1->GetYaxis()->SetTitleOffset(1.5);
        h1->GetYaxis()->SetLabelSize(0.02);
        h1->GetYaxis()->SetTitleSize(0.02);
        h1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
        h1->GetXaxis()->SetTitle("|#eta|");
        
        h1->SetMarkerColor(kRed+1);
        h1->SetMarkerStyle(kCircle);
        h1->SetMarkerSize(1);
        h1->SetLineColor(kRed+1);
        h1->SetLineStyle(2);
        h1->SetLineWidth(1);
        h1->Draw("HIST");
        
        h2->SetMarkerColor(kGreen+3);
        h2->SetMarkerStyle(kCircle);
        h2->SetMarkerSize(1);
        h2->SetLineColor(kGreen+3);
        h2->SetLineStyle(3);
        h2->SetLineWidth(1);
        
        h3->SetMarkerColor(kBlue);
        h3->SetMarkerStyle(kCircle);
        h3->SetMarkerSize(1);
        h3->SetLineColor(kBlue);
        h3->SetLineStyle(4);
        h3->SetLineWidth(1);
        
        h4->SetMarkerColor(kOrange-3);
        h4->SetMarkerStyle(kCircle);
        h4->SetMarkerSize(1);
        h4->SetLineColor(kOrange-3);
        h4->SetLineStyle(5);
        h4->SetLineWidth(1);
        
        h5->SetMarkerColor(kGray+3);
        h5->SetMarkerStyle(kCircle);
        h5->SetMarkerSize(1);
        h5->SetLineColor(kGray+3);
        h5->SetLineStyle(1);
        h5->SetLineWidth(1);
        
        c1->SetLogy(1);
       
         
        TH1F *f1 = (TH1F*)h1->Clone();
        TH1F *f2 = (TH1F*)h2->Clone();
        TH1F *f3 = (TH1F*)h3->Clone();
        TH1F *f4 = (TH1F*)h4->Clone();
        f1->Divide(h1,h5);
        f2->Divide(h2,h5);
        f3->Divide(h3,h5);
        f4->Divide(h4,h5);
        
        TH1F *f0_1 = (TH1F*)h4->Clone("f0_1");
        f0_1->SetLineColor(kWhite);
        f0_1->SetMarkerColor(kWhite);
        f0_1->SetMarkerStyle(kDot);
        TRatioPlot *rp = new TRatioPlot(f0_1,f0_1);
        
        //rp->SetSeparationMargin(0.02);
        rp->Draw();
        rp->GetLowerPad()->cd();
        for (int i=0; i<6; i++){
            rp->GetLowerRefGraph()->RemovePoint(0);
        } 
        cout<< "N points: "<< rp->GetLowerRefGraph()->GetN() << endl;
        rp->GetLowerRefGraph()->SetPoint(0,2.5,1);
        rp->GetLowerRefGraph()->Draw("HIST E0");
        f1->Draw("SAME E1p");
        f2->Draw("SAME E1p");
        f3->Draw("SAME E1p");
        f4->Draw("SAME E1p");
        rp->GetUpperPad()->cd();
        h1->Draw("SAME E1p");
        h2->Draw("SAME E1p");
        h3->Draw("SAME E1p");
        h4->Draw("SAME E1p");
        h5->Draw("SAME E1p");
        
        rp->GetUpperRefYaxis()->SetTitleOffset(1.5);
        rp->GetUpperRefYaxis()->SetLabelSize(0.03);
        rp->GetUpperRefYaxis()->SetTitleSize(0.03);
        rp->GetUpperRefYaxis()->SetTitle("#epsilon_{QmisID}");
        
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(2.0);
        rp->GetLowerRefYaxis()->SetLabelSize(0.02);
        rp->GetLowerRefYaxis()->SetTitleSize(0.02);
        //rp->GetLowerRefYaxis()->SetTitle(name1+" / "+name2);
        rp->GetLowerRefYaxis()->SetTitle("#epsilon_{QmisID} / #epsilon_{"+name5+"}");
        rp->GetLowerRefXaxis()->SetLimits(0.,2.5);
        rp->GetLowerRefXaxis()->SetTitleSize(0.03);
        rp->GetLowerRefXaxis()->SetLabelSize(0.03);
        rp->GetLowerRefXaxis()->SetTitle("|#eta|");
        
        newanalyzer.GetATLAS("Internal",0.13,0.81,false,0.05);
        newanalyzer.GetLabel(0.132,0.75,"13 TeV, 140 fb^{-1}",0.04);
        newanalyzer.GetLabel(0.132,0.69,"p_{T}^{l} "+ptrange,0.04);
        newanalyzer.GetLabel(0.132,0.63,"Z+jets",0.03);
        newanalyzer.GetLabel(0.132,0.58,"w/ECIDS",0.03);
        TLegend *legend = new TLegend(0.54,0.64,0.84,0.79, "");
        //legend->SetNColumns(2);
        legend->AddEntry(h1,name1+" (Likelihood)","lp");
        legend->AddEntry(h2,name2+" (Likelihood)","lp");
        legend->AddEntry(h3,name3+" (Likelihood)","lp");
        legend->AddEntry(h4,name4+" (Likelihood)","lp");
        legend->AddEntry(h5,name5+" (Likelihood)","lp");
        legend->SetTextSize(0.03);
        legend->Draw("same");
        TString ipt = std::to_string(i);
        c1->SaveAs("Plots/compare5_"+name1+"_"+name5+"_pt"+ipt+".pdf");
    } 


return;
}

void Plot2HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, TString name1, TString name2, TString meth1, TString meth2){
    std::cout << "hrate1: " << name1 << " has size" << v_hrate1.size() << std::endl;
    std::cout << "hrate2: " << name2 << " has size" << v_hrate2.size() << std::endl;
    //if (name1.Contains("mc")) name1="mc16_all";
    //if (name2.Contains("mc")) name2="mc16_all";
    
    TString pt1 = "[0,60] GeV"; 
    TString pt2 = "[60,90] GeV"; 
    TString pt3 = "[90,130] GeV"; 
    TString pt4 = "> 130 GeV";
    std::vector<TString> v_ptbin = {pt1,pt2,pt3,pt4};

    for (unsigned int i=0; i<v_ptbin.size(); i++){
 
        TString ptrange = v_ptbin.at(i); 
        TH1F *h1 = v_hrate1.at(i);    
        TH1F *h2 = v_hrate2.at(i);    

        TCanvas *c1 = new TCanvas("c", "canvas", 850, 810);
        c1->SetLeftMargin(0.15);
        gStyle->SetErrorX(0.5);
        MiniTreeAnalyzer newanalyzer;
        h1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
        h1->GetXaxis()->SetTitle("|#eta|");
        h1->GetYaxis()->SetTitleOffset(1.5);
        h1->GetYaxis()->SetLabelSize(0.03);
        h1->GetYaxis()->SetTitleSize(0.03);
        h1->SetAxisRange(1e-5, 10, "Y");
        h1->SetMarkerColor(kRed+1);
        h1->SetMarkerStyle(kCircle);
        h1->SetMarkerSize(1);
        h1->SetLineColor(kRed+1);
        h1->SetLineStyle(2);
        h1->SetLineWidth(1);
        h1->Draw("E1p");
        h2->SetMarkerColor(kGray+3);
        h2->SetMarkerStyle(kCircle);
        h2->SetMarkerSize(1);
        h2->SetLineColor(kGray+3);
        h2->SetLineStyle(1);
        h2->SetLineWidth(1);
        h2->Draw("SAME E1p");
        c1->SetLogy(1);
        
        TH1F *f1 = (TH1F*)h1->Clone("h1");
        TH1F *f2 = (TH1F*)h2->Clone("h2");
        f1->Divide(h1,h2);
        
        TH1F *f0_1 = (TH1F*)h1->Clone("f0_1");
        f0_1->SetLineColor(kWhite);
        f0_1->SetMarkerColor(kWhite);
        f0_1->SetMarkerStyle(kDot);
        TRatioPlot *rp = new TRatioPlot(f0_1,f0_1);
       
        
        //rp->SetSeparationMargin(0.02);
        rp->Draw();
        rp->GetLowerPad()->cd();
        for (int i=0; i<6; i++){
            rp->GetLowerRefGraph()->RemovePoint(0);
        } 
        cout<< "N points: "<< rp->GetLowerRefGraph()->GetN() << endl;
        rp->GetLowerRefGraph()->SetPoint(0,2.5,1);
        rp->GetLowerRefGraph()->Draw("HIST E0");
        f1->Draw("SAME E1p");
        rp->GetUpperPad()->cd();
        h1->Draw("SAME E1p");
        h2->Draw("SAME E1p");
        
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(2.0);
        rp->GetLowerRefYaxis()->SetLabelSize(0.02);
        rp->GetLowerRefYaxis()->SetTitleSize(0.02);
        //rp->GetLowerRefYaxis()->SetTitle(name1+" / "+name2/);
        rp->GetLowerRefYaxis()->SetTitle("#epsilon_{QmisID} / #epsilon_{"+name2+"("+meth2+")}");
        rp->GetLowerRefXaxis()->SetLimits(0.,2.5);
        rp->GetLowerRefXaxis()->SetTitleSize(0.03);
        rp->GetLowerRefXaxis()->SetLabelSize(0.03);
        rp->GetLowerRefXaxis()->SetTitle("|#eta|");
        
        newanalyzer.GetATLAS("Internal",0.13,0.81,false,0.05);
        newanalyzer.GetLabel(0.132,0.75,"13 TeV, 140 fb^{-1}",0.04);
        newanalyzer.GetLabel(0.132,0.69,"p_{T}^{l} "+ptrange,0.04);
        newanalyzer.GetLabel(0.132,0.63,"Z+jets",0.03);
        newanalyzer.GetLabel(0.132,0.58,"w/ECIDS",0.03);
        TLegend *legend = new TLegend(0.54,0.64,0.84,0.79, "");
        //legend->SetNColumns(2);
        legend->AddEntry(h1,name1+" ("+meth1+")","lp");
        legend->AddEntry(h2,name2+" ("+meth2+")","lp");
        legend->SetTextSize(0.03);
        legend->Draw("same");
        TString ipt = std::to_string(i);
        c1->SaveAs("Plots/compare2_"+name1+"_"+name2+"_pt"+ipt+".pdf");
    } 
    return;
}

void Plot3HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, TString name1, TString name2, TString name3){

    std::cout << "hrate1: " << name1 << " has size" << v_hrate1.size() << std::endl;
    std::cout << "hrate2: " << name2 << " has size" << v_hrate2.size() << std::endl;
    std::cout << "hrate3: " << name3 << " has size" << v_hrate3.size() << std::endl;
    
    TString pt1 = "[0,60] GeV"; 
    TString pt2 = "[60,90] GeV"; 
    TString pt3 = "[90,130] GeV"; 
    TString pt4 = "> 130 GeV";
    std::vector<TString> v_ptbin = {pt1,pt2,pt3,pt4};

    for (unsigned int i=0; i<v_ptbin.size(); i++){
 
        TString ptrange = v_ptbin.at(i); 
        TH1F *h1 = v_hrate1.at(i);    
        TH1F *h2 = v_hrate2.at(i);    
        TH1F *h3 = v_hrate3.at(i);    

        TCanvas *c1 = new TCanvas("c", "canvas", 850, 810);
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.3);
        gStyle->SetErrorX(0.5);
        MiniTreeAnalyzer newanalyzer;
        h1->SetAxisRange(1e-5, 10, "Y");
        h1->GetYaxis()->SetTitleOffset(1.5);
        h1->GetYaxis()->SetLabelSize(0.02);
        h1->GetYaxis()->SetTitleSize(0.02);
        h1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
        h1->GetXaxis()->SetTitle("|#eta|"); 
        h1->SetMarkerColor(kRed+1);
        h1->SetMarkerStyle(kCircle);
        h1->SetMarkerSize(1);
        h1->SetLineColor(kRed+1);
        h1->SetLineStyle(2);
        h1->SetLineWidth(1);
        h1->Draw("HIST");
        
        h2->SetMarkerColor(kGreen+3);
        h2->SetMarkerStyle(kCircle);
        h2->SetMarkerSize(1);
        h2->SetLineColor(kGreen+3);
        h2->SetLineStyle(3);
        h2->SetLineWidth(1);
        
        h3->SetMarkerColor(kGray+3);
        h3->SetMarkerStyle(kCircle);
        h3->SetMarkerSize(1);
        h3->SetLineColor(kGray+3);
        h3->SetLineStyle(1);
        h3->SetLineWidth(1);
        
        c1->SetLogy(1);
       
         
        TH1F *f1 = (TH1F*)h1->Clone();
        TH1F *f2 = (TH1F*)h2->Clone();
        TH1F *f3 = (TH1F*)h3->Clone();
        f1->Divide(h1,h3);
        f2->Divide(h2,h3);
        
        TH1F *f0_1 = (TH1F*)h1->Clone("f0_1");
        f0_1->SetLineColor(kWhite);
        f0_1->SetMarkerColor(kWhite);
        f0_1->SetMarkerStyle(kDot);
        TRatioPlot *rp = new TRatioPlot(f0_1,f0_1);
        
        //rp->SetSeparationMargin(0.02);
        rp->Draw();
        rp->GetLowerPad()->cd();
        for (int i=0; i<6; i++){
            rp->GetLowerRefGraph()->RemovePoint(0);
        } 
        cout<< "N points: "<< rp->GetLowerRefGraph()->GetN() << endl;
        rp->GetLowerRefGraph()->SetPoint(0,2.5,1);
        rp->GetLowerRefGraph()->Draw("HIST E0");
        f1->Draw("SAME E1p");
        f2->Draw("SAME E1p");
        rp->GetUpperPad()->cd();
        h1->Draw("SAME E1p");
        h2->Draw("SAME E1p");
        h3->Draw("SAME E1p");
        
        rp->GetUpperRefYaxis()->SetTitleOffset(1.5);
        rp->GetUpperRefYaxis()->SetLabelSize(0.03);
        rp->GetUpperRefYaxis()->SetTitleSize(0.03);
        rp->GetUpperRefYaxis()->SetTitle("#epsilon_{QmisID}");
        
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(2.0);
        rp->GetLowerRefYaxis()->SetLabelSize(0.02);
        rp->GetLowerRefYaxis()->SetTitleSize(0.02);
        //rp->GetLowerRefYaxis()->SetTitle(name1+" / "+name2);
        rp->GetLowerRefYaxis()->SetTitle("#epsilon_{QmisID} / #epsilon_{"+name3+"}");
        rp->GetLowerRefXaxis()->SetLimits(0.,2.5);
        rp->GetLowerRefXaxis()->SetTitleSize(0.03);
        rp->GetLowerRefXaxis()->SetLabelSize(0.03);
        rp->GetLowerRefXaxis()->SetTitle("|#eta|");
        
        newanalyzer.GetATLAS("Internal",0.13,0.81,false,0.05);
        newanalyzer.GetLabel(0.132,0.75,"13 TeV, 140 fb^{-1}",0.04);
        newanalyzer.GetLabel(0.132,0.69,"p_{T}^{l} "+ptrange,0.04);
        newanalyzer.GetLabel(0.132,0.63,"Z+jets",0.03);
        newanalyzer.GetLabel(0.132,0.58,"w/ECIDS",0.03);
        TLegend *legend = new TLegend(0.54,0.64,0.84,0.79, "");
        //legend->SetNColumns(2);
        legend->AddEntry(h1,name1+" (Likelihood)","lp");
        legend->AddEntry(h2,name2+" (Likelihood)","lp");
        legend->AddEntry(h3,name3+" (Likelihood)","lp");
        legend->SetTextSize(0.03);
        legend->Draw("same");
        TString ipt = std::to_string(i);
        c1->SaveAs("Plots/compare4_"+name1+"_"+name3+"_pt"+ipt+".pdf");
    } 


    return;
}
void Plot4HistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, vector<TH1F*> v_hrate3, vector<TH1F*> v_hrate4, TString name1, TString name2, TString name3, TString name4){

    std::cout << "hrate1: " << name1 << " has size" << v_hrate1.size() << std::endl;
    std::cout << "hrate2: " << name2 << " has size" << v_hrate2.size() << std::endl;
    std::cout << "hrate3: " << name3 << " has size" << v_hrate3.size() << std::endl;
    std::cout << "hrate4: " << name4 << " has size" << v_hrate4.size() << std::endl;
    
    TString pt1 = "#in [20, 50] GeV"; 
    TString pt2 = "#in [50, 100] GeV"; 
    TString pt3 = "#in [100, 200] GeV"; 
    TString pt4 = "> 200 GeV";
    std::vector<TString> v_ptbin = {pt1,pt2,pt3,pt4};

    for (unsigned int i=0; i<v_ptbin.size(); i++){
 
        TString ptrange = v_ptbin.at(i); 
        TH1F *h1 = v_hrate1.at(i);    
        TH1F *h2 = v_hrate2.at(i);    
        TH1F *h3 = v_hrate3.at(i);    
        TH1F *h4 = v_hrate4.at(i);    

        TCanvas *c1 = new TCanvas("c", "canvas", 850, 810);
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.3);
        gStyle->SetErrorX(0.5);
        MiniTreeAnalyzer newanalyzer;
        h1->SetAxisRange(1e-5, 10, "Y");
        h1->GetYaxis()->SetTitleOffset(1.5);
        h1->GetYaxis()->SetLabelSize(0.02);
        h1->GetYaxis()->SetTitleSize(0.02);
        h1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
        h1->GetXaxis()->SetTitle("|#eta|"); 
        std::cout << "--debug third eta lower edge=" << h1->GetXaxis()->GetBinLowEdge(4) << std::endl;

        h1->SetMarkerColor(kRed+1);
        h1->SetMarkerStyle(kOpenDiamond);
        h1->SetMarkerSize(1);
        h1->SetLineColor(kRed+1);
        h1->SetLineStyle(1);
        h1->SetLineWidth(2);
        h1->Draw("HIST");
        
        h2->SetMarkerColor(kGreen+3);
        h2->SetMarkerStyle(kOpenTriangleUp);
        h2->SetMarkerSize(1);
        h2->SetLineColor(kGreen+3);
        h2->SetLineStyle(1);
        h2->SetLineWidth(2);
        
        h3->SetMarkerColor(kBlue);
        h3->SetMarkerStyle(kOpenSquare);
        h3->SetMarkerSize(1);
        h3->SetLineColor(kBlue);
        h3->SetLineStyle(1);
        h3->SetLineWidth(2);
        
        h4->SetMarkerColor(kGray+3);
        h4->SetMarkerStyle(kOpenCircle);
        h4->SetMarkerSize(1);
        h4->SetLineColor(kGray+3);
        h4->SetLineStyle(1);
        h4->SetLineWidth(2);
        
        c1->SetLogy(1);
       
         
        TH1F *f1 = (TH1F*)h1->Clone();
        TH1F *f2 = (TH1F*)h2->Clone();
        TH1F *f3 = (TH1F*)h3->Clone();
        TH1F *f4 = (TH1F*)h4->Clone();
        f1->Divide(h1,h4);
        f2->Divide(h2,h4);
        f3->Divide(h3,h4);
        
        TH1F *f0_1 = (TH1F*)h1->Clone("f0_1");
        f0_1->SetLineColor(kWhite);
        f0_1->SetMarkerColor(kWhite);
        f0_1->SetMarkerStyle(kDot);
        TRatioPlot *rp = new TRatioPlot(f0_1,f0_1);
        
        //rp->SetSeparationMargin(0.02);
        rp->Draw();
        rp->GetLowerPad()->cd();
        for (int i=0; i<4; i++){
            rp->GetLowerRefGraph()->RemovePoint(0);
        } 
        cout<< "N points: "<< rp->GetLowerRefGraph()->GetN() << endl;
        rp->GetLowerRefGraph()->SetPoint(0,2.6,1);
        rp->GetLowerRefGraph()->Draw("HIST E0");
        f1->Draw("SAME E1p");
        f2->Draw("SAME E1p");
        f3->Draw("SAME E1p");
        rp->GetUpperPad()->cd();
        h1->Draw("SAME E1p");
        h2->Draw("SAME E1p");
        h3->Draw("SAME E1p");
        h4->Draw("SAME E1p");
        
        rp->GetUpperRefYaxis()->SetTitleOffset(1.5);
        rp->GetUpperRefYaxis()->SetLabelSize(0.03);
        rp->GetUpperRefYaxis()->SetTitleSize(0.03);
        rp->GetUpperRefYaxis()->SetTitle("#epsilon_{QmisID}");
        
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(2.0);
        rp->GetLowerRefYaxis()->SetLabelSize(0.02);
        rp->GetLowerRefYaxis()->SetTitleSize(0.02);
        //rp->GetLowerRefYaxis()->SetTitle(name1+" / "+name2);
        rp->GetLowerRefYaxis()->SetTitle("#epsilon_{QmisID} / #epsilon_{"+name4+"}");
        rp->GetLowerRefXaxis()->SetLimits(0.,2.5);
        rp->GetLowerRefXaxis()->SetTitleSize(0.03);
        rp->GetLowerRefXaxis()->SetLabelSize(0.03);
        rp->GetLowerRefXaxis()->SetTitle("|#eta|");
        
        newanalyzer.GetATLAS("Internal",0.18,0.81,false,0.05);
        newanalyzer.GetLabel(0.182,0.75,"13 TeV, 140 fb^{-1}",0.04);
        newanalyzer.GetLabel(0.182,0.69,"p_{T}^{l} "+ptrange,0.04);
        // newanalyzer.GetLabel(0.132,0.63,"Z+jets",0.03);
        // newanalyzer.GetLabel(0.132,0.58,"w/ECIDS",0.03);
        TLegend *legend = new TLegend(0.54,0.7,0.84,0.85, "");
        //legend->SetNColumns(2);
        legend->AddEntry(h1,name1,"lp");
        legend->AddEntry(h2,name2,"lp");
        legend->AddEntry(h3,name3,"lp");
        legend->AddEntry(h4,name4,"lp");
        legend->SetTextSize(0.03);
        legend->Draw("same");
        TString ipt = std::to_string(i);
        c1->SaveAs("../Plots/compare4Hist_pt"+ipt+".pdf");
    } 


    return;
}

void PlotHistCompare(std::vector<TH1F*> v_hrate1, vector<TH1F*> v_hrate2, TString name1, TString name2){
    if (name1.Contains("mc")) name1="mc16";
    TCanvas *c1 = new TCanvas("c", "canvas", 1200, 1000);
    //c1->Divide(1,2);
  //  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    //c1->cd(1);
    //c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.1);
    gStyle->SetErrorX(0.5);
    MiniTreeAnalyzer newanalyzer;
    
    std::cout << "hrate1 has size" << v_hrate1.size() << endl;
    TH1F *h1_1 = v_hrate1.at(0);    
    TH1F *h1_2 = v_hrate1.at(1);    
    TH1F *h1_3 = v_hrate1.at(2);    
    TH1F *h1_4 = v_hrate1.at(3);    

    std::cout << "hrate2 has size" << v_hrate2.size() << endl;
    TH1F *h2_1 = v_hrate2.at(0);    
    TH1F *h2_2 = v_hrate2.at(1);    
    TH1F *h2_3 = v_hrate2.at(2);    
    TH1F *h2_4 = v_hrate2.at(3);    

    h1_1->GetXaxis()->SetTitle("|#eta|");
    h1_1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
    h1_1->GetYaxis()->SetTitleOffset(1.2);
    h1_1->GetYaxis()->SetLabelSize(0.03);
    h1_1->GetYaxis()->SetTitleSize(0.03);
    if (name1.Contains("Data2015")) h1_1->SetAxisRange(1e-5, 1, "Y");
    else h1_1->SetAxisRange(1e-5, 1, "Y");
    h1_1->SetMarkerColor(kGray+3);
    h1_1->SetMarkerStyle(kOpenTriangleUp);
    h1_1->SetMarkerSize(1);
    h1_1->SetLineColor(kGray+3);
    h1_1->SetLineStyle(2);
    h1_1->SetLineWidth(1);
    h1_1->Draw("E1p");
    h1_2->SetMarkerColor(kRed+1);
    h1_2->SetMarkerStyle(kOpenTriangleUp);
    h1_2->SetMarkerSize(1);
    h1_2->SetLineColor(kRed+1);
    h1_2->SetLineStyle(2);
    h1_2->SetLineWidth(1);
    // h1_2->Draw("SAMEE1p");
    h1_3->SetMarkerColor(kGreen+3);
    h1_3->SetMarkerStyle(kOpenTriangleUp);
    h1_3->SetMarkerSize(1);
    h1_3->SetLineColor(kGreen+3);
    h1_3->SetLineStyle(2);
    h1_3->SetLineWidth(1);
    // h1_3->Draw("SAMEE1p");
    h1_4->SetMarkerColor(kBlue);
    h1_4->SetMarkerStyle(kOpenTriangleUp);
    h1_4->SetMarkerSize(1);
    h1_4->SetLineColor(kBlue);
    h1_4->SetLineStyle(2);
    h1_4->SetLineWidth(1);
    // h1_4->Draw("SAMEE1p");

    h2_1->SetMarkerColor(kGray+3);
    h2_1->SetMarkerStyle(kFullCircle);
    h2_1->SetMarkerSize(1);
    h2_1->SetLineColor(kGray+3);
    h2_1->SetLineStyle(1);
    h2_1->SetLineWidth(1);
    // h2_1->Draw("SAMEE1p");
    h2_2->SetMarkerColor(kRed+1);
    h2_2->SetMarkerStyle(kFullCircle);
    h2_2->SetMarkerSize(1);
    h2_2->SetLineColor(kRed+1);
    h2_2->SetLineStyle(1);
    h2_2->SetLineWidth(1);
    // h2_2->Draw("SAMEE1p");
    h2_3->SetMarkerColor(kGreen+3);
    h2_3->SetMarkerStyle(kFullCircle);
    h2_3->SetMarkerSize(1);
    h2_3->SetLineColor(kGreen+3);
    h2_3->SetLineStyle(1);
    h2_3->SetLineWidth(1);
    // h2_3->Draw("SAMEE1p");
    h2_4->SetMarkerColor(kBlue);
    h2_4->SetMarkerStyle(kFullCircle);
    h2_4->SetMarkerSize(1);
    h2_4->SetLineColor(kBlue);
    h2_4->SetLineStyle(1);
    h2_4->SetLineWidth(1);
    // h2_4->Draw("SAMEE1p");

    c1->SetLogy(1);
    // c1->SetTicks(1, 0); 

    TH1F *f0_1 = (TH1F*)h1_1->Clone("h1_1"); 
    TH1F *f0_2 = (TH1F*)h2_1->Clone("h2_1");
    f0_1->SetLineColor(kWhite); 
    f0_2->SetLineColor(kWhite); 
    //c1->cd(2);
    // TRatioPlot *rp = new TRatioPlot(f0_1,f0_2);

    f0_1->SetLineColor(kWhite);
    f0_1->SetMarkerColor(kWhite);
    f0_1->SetMarkerStyle(kDot);
    TRatioPlot *rp = new TRatioPlot(f0_1,f0_1);
    
    //rp->SetSeparationMargin(0.02);
    rp->Draw();
    rp->GetLowerPad()->cd();
    for (int i=0; i<4; i++){
        rp->GetLowerRefGraph()->RemovePoint(0);
    } 
    cout<< "N points: "<< rp->GetLowerRefGraph()->GetN() << endl;
    rp->GetLowerRefGraph()->SetPoint(0,2.6,1);
    rp->GetLowerRefGraph()->Draw("HIST E0");
    rp->GetLowerPad()->cd();

    TH1F *f1 = (TH1F*)h2_1->Clone("h2_1");
    TH1F *f2 = (TH1F*)h2_2->Clone("h2_2");
    TH1F *f3 = (TH1F*)h2_3->Clone("h2_3");
    TH1F *f4 = (TH1F*)h2_4->Clone("h2_4");;
    
    f1->Divide(h2_1,h1_1);
    f2->Divide(h2_2,h1_2);
    f3->Divide(h2_3,h1_3);
    f4->Divide(h2_4,h1_4);
    
    //rp->GetLowerRefYaxis()->SetRangeUser(0.0,2.0);
    f1->Draw("SAME E1p");
    f2->Draw("SAME E1p");
    f3->Draw("SAME E1p");
    f4->Draw("SAME E1p");
    //c1->Update();
    rp->GetUpperPad()->cd();
    h1_1->Draw("SAMEE1p");
    h1_2->Draw("SAMEE1p");
    h1_3->Draw("SAMEE1p");
    h1_4->Draw("SAMEE1p");
    h2_1->Draw("SAMEE1p");
    h2_2->Draw("SAMEE1p");
    h2_3->Draw("SAMEE1p");
    h2_4->Draw("SAMEE1p");
    
    // rp->GetLowerRefGraph()->SetMinimum(0.0);
    // rp->GetLowerRefGraph()->SetMaximum(2.0);

    rp->GetLowerRefGraph()->SetMinimum(0.5);
    rp->GetLowerRefGraph()->SetMaximum(1.5);
                        
    rp->GetLowerRefGraph()->SetPoint(0,2.6,1);
    rp->GetLowerRefGraph()->Draw("HIST E0");

    rp->GetLowerRefYaxis()->SetLabelSize(0.02);
    rp->GetLowerRefYaxis()->SetTitleSize(0.03);
    rp->GetLowerRefYaxis()->SetNdivisions(101); //404
    rp->GetLowerRefYaxis()->SetTitle("Ratio");
    rp->GetLowerRefXaxis()->SetLimits(0.,2.5);
    rp->GetLowerRefXaxis()->SetTitle("|#eta|");
    rp->GetLowerRefXaxis()->SetTitleSize(0.03);
    rp->GetLowerRefXaxis()->SetLabelSize(0.03);
    
    newanalyzer.GetATLAS("Internal",0.185,0.84,false,0.04);
    newanalyzer.GetLabel(0.186,0.79,"13 TeV, 140 fb^{-1}",0.03);

    TLegend *legend = new TLegend(0.186,0.60,0.59,0.75);
    legend->SetNColumns(2);
    if (name1.Contains("Data")) {
        legend->AddEntry(f0_1, "Data Likelihood", "l");
        legend->AddEntry(f0_1, "MC Likelihood", "l");
        // legend->AddEntry(f0_1, "Data Likelihood", "l");
    }
    else if (name1.Contains("MC")) {
        legend->AddEntry(f0_1, "MC Likelihood", "l");
        legend->AddEntry(f0_1, "MC Truth-matched", "l");
        // legend->AddEntry(f0_1, "MC Likelihood", "l");

    }

    legend->AddEntry(h1_1, "p_{T} #in [20, 50] GeV","lp");
    legend->AddEntry(h2_1, "p_{T} #in [20, 50] GeV","lp");

    legend->AddEntry(h1_2, "p_{T} #in [50, 100] GeV","lp");
    legend->AddEntry(h2_2, "p_{T} #in [50, 100] GeV","lp");

    legend->AddEntry(h1_3, "p_{T} #in [100, 200] GeV","lp");
    legend->AddEntry(h2_3, "p_{T} #in [100, 200] GeV","lp");

    legend->AddEntry(h1_4, "p_{T} > 200 GeV","lp");
    legend->AddEntry(h2_4, "p_{T} > 200 GeV","lp");

    legend->SetTextSize(0.03);
    legend->Draw("same");
    c1->SaveAs("Plots/"+name1+"vs"+name2+".pdf");
    c1->SaveAs("Plots/"+name1+"vs"+name2+".png");
    return;
}

std::vector<TH1F*> GetRateHist(TString FileName, TString TreeName, TString Version){
    TFile *f0 = new TFile(FileName);
    std::cout << "Openfile: " << FileName << std::endl;
    TH2F  *h0 = (TH2F*)f0->Get(TreeName);
    std::cout << "Opentree: " << TreeName << std::endl;

    Int_t   NEtaBin = 4;
    Int_t   NPtBin  = 4;
  Float_t EtaBinning[]  = {0, 1.37, 1.52, 2.0, 2.6};
  Float_t PtBinning[]   = {20, 50, 100, 200, 2600};
    //TH1F *hEta = new TH1F("Rates in Eta binning", "",NEtaBin,EtaBinning);
    TH1F *hPt1  = new TH1F("h1", "h1",NEtaBin,EtaBinning);
    TH1F *hPt2  = new TH1F("h2", "h2",NEtaBin,EtaBinning);
    TH1F *hPt3  = new TH1F("h3", "h3",NEtaBin,EtaBinning);
    TH1F *hPt4  = new TH1F("h4", "h4",NEtaBin,EtaBinning);
    TH2F *hRates = new TH2F("Rates", "Rates",NEtaBin,EtaBinning,NPtBin,PtBinning);
    Float_t Rates = 0;
    Float_t eRates = 0;
    for(int i = 1;i <= h0->GetNbinsX();i++) { 
        for(int j = 1;j <= h0->GetNbinsY();j++) {
            Rates = h0->GetBinContent(i,j);
            eRates = h0->GetBinError(i,j);
            if (i==2) {
                Rates=0;
                eRates=0;
                //shuhui: set dummy values for crack-veto region
            }
            hRates->SetBinContent(i,j,Rates);
            std::cout << "EtaBin: " << i << ", PtBin: " << j << ", rates: " << Rates << ", errors: " << eRates << std::endl;        
            if(j==1) {
                hPt1->SetBinContent(i,Rates);
                hPt1->SetBinError(i,eRates);
            }
            else if(j==2) {
                 hPt2->SetBinContent(i,Rates);
                 hPt2->SetBinError(i,eRates);
            }
            else if(j==3) {
                 hPt3->SetBinContent(i,Rates);
                 hPt3->SetBinError(i,eRates);
            }
            else if(j==4) {
                hPt4->SetBinContent(i,Rates);
                hPt4->SetBinError(i,eRates);
            }
        }
    }
    vector<TH1F*> v_hRate;
    v_hRate.push_back(hPt1);
    v_hRate.push_back(hPt2);
    v_hRate.push_back(hPt3);
    v_hRate.push_back(hPt4);

    return v_hRate;
}

void GetHistPlot(TString FileName, TString TreeName, TString Version){
    TFile *f0 = new TFile(FileName);
    std::cout << "Openfile: " << FileName << std::endl;
    TH2F  *h0 = (TH2F*)f0->Get(TreeName);
    std::cout << "Opentree: " << TreeName << std::endl;

    Int_t   NEtaBin = 6;
    Int_t   NPtBin  = 4;
    Float_t EtaBinning[]= {0,0.6,1.1,1.52,1.7,2.3,2.5};
    Float_t PtBinning[] = {0,60,90,130,2500};
    //TH1F *hEta = new TH1F("Rates in Eta binning", "",NEtaBin,EtaBinning);
    TH1F *hPt1  = new TH1F("h1", "h1",NEtaBin,EtaBinning);
    TH1F *hPt2  = new TH1F("h2", "h2",NEtaBin,EtaBinning);
    TH1F *hPt3  = new TH1F("h3", "h3",NEtaBin,EtaBinning);
    TH1F *hPt4  = new TH1F("h4", "h4",NEtaBin,EtaBinning);
    TH2F *hRates = new TH2F("Rates", "Rates",NEtaBin,EtaBinning,NPtBin,PtBinning);
    Float_t Rates = 0;
    Float_t eRates = 0;
    for(int i = 1;i <= h0->GetNbinsX();i++) { 
        for(int j = 1;j <= h0->GetNbinsY();j++) {
            Rates = h0->GetBinContent(i,j);
            eRates = h0->GetBinError(i,j);
            hRates->SetBinContent(i,j,Rates);
            std::cout << "EtaBin: " << i << ", PtBin: " << j << ", rates: " << Rates << ", errors: " << eRates << std::endl;        
            if(j==1) {
                hPt1->SetBinContent(i,Rates);
                hPt1->SetBinError(i,eRates);
            }
            else if(j==2) {
                 hPt2->SetBinContent(i,Rates);
                 hPt2->SetBinError(i,eRates);
            }
            else if(j==3) {
                 hPt3->SetBinContent(i,Rates);
                 hPt3->SetBinError(i,eRates);
            }
            else if(j==4) {
                hPt4->SetBinContent(i,Rates);
                hPt4->SetBinError(i,eRates);
            }
        }
    }


    PlotRate1D(hPt1,hPt2,hPt3,hPt4,Version);
    return ;
}

void PlotRate1D(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TString outName){

    TCanvas *c1 = new TCanvas;
    c1->cd();
    c1->SetRightMargin(0.15);
    //gStyle->SetPaintTextFormat("4.12f");
    //gStyle->SetTextFont(62);i
    gStyle->SetErrorX(0.5);
    MiniTreeAnalyzer newanalyzer;
    h1->GetXaxis()->SetTitle("|#eta|");
    h1->GetYaxis()->SetTitle("Chargeflip Rates");
    h1->GetYaxis()->SetTitleOffset(1.2);
    if (outName=="Data15") h1->SetAxisRange(1e-12, 1, "Y");
    else if (outName.Contains("mc")) h1->SetAxisRange(1e-5, 1, "Y");
    else h1->SetAxisRange(1e-5, 1, "Y");    
    h1->SetMarkerColor(kRed+1);
    h1->SetMarkerStyle(kFullCircle);
    h1->SetMarkerSize(1);
    h1->SetLineColor(kRed+1);
    h2->SetLineStyle(1);
    h2->SetLineWidth(1);
    h1->Draw("E1p");
    h2->SetMarkerColor(kBlue);
    h2->SetMarkerStyle(kFullCircle);
    h2->SetMarkerSize(1);
    h2->SetLineColor(kBlue);
    h2->SetLineStyle(1);
    h2->SetLineWidth(1);
    h2->Draw("SAMEE1p");
    h3->SetMarkerColor(kGreen+3);
    h3->SetMarkerStyle(kFullCircle);
    h3->SetMarkerSize(1);
    h3->SetLineColor(kGreen+3);
    h3->SetLineStyle(1);
    h3->SetLineWidth(1);
    h3->Draw("SAMEE1p");
    h4->SetMarkerColor(kOrange-3);
    h4->SetMarkerStyle(kFullCircle);
    h4->SetMarkerSize(1);
    h4->SetLineColor(kOrange-3);
    h4->SetLineStyle(1);
    h4->SetLineWidth(1);
    h4->Draw("SAMEE1p");
    c1->SetLogy(1);
    newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.04);
    newanalyzer.GetLabel(0.186,0.81,"",0.03);
    
    TLegend *legend = new TLegend(0.65,0.18,0.83,0.38, outName+" likelihood");
    legend->AddEntry(h1,"p_{T}[0,60]GeV","lp");
    legend->AddEntry(h2,"p_{T}[60,90]GeV","lp");
    legend->AddEntry(h3,"p_{T}[90,130]GeV","lp");
    legend->AddEntry(h4,"p_{T}[130,2500]GeV","lp");
    legend->SetTextSize(0.02);
    legend->Draw("same");
    c1->SaveAs("Plots/"+outName+".pdf");

    return;
}

void PlotRate2D(TH2F *hRates,TString outName) {

    TCanvas *c1 = new TCanvas;
    c1->cd();
    //c1->SetLogx(0);
    //c1->SetLogz(1);
    c1->SetRightMargin(0.15);
    gStyle->SetPaintTextFormat("4.12f");
    gStyle->SetTextFont(62);
    gStyle->SetOptLogy(1);
    //TPad *pad1 = new TPad();
    //pad1->cd();
    //pad1->SetRightMargin(0.15);
    //pad1->SetLogy(1);
    //pad1->SetLogz(1);
    //pad1->Draw();
    MiniTreeAnalyzer newanalyzer;
    hRates->GetXaxis()->SetTitle("#eta");
    //hRates->GetYaxis()->SetMoreLogLabels();
    //hRates->GetYaxis()->SetRange(0.1,2500);
    hRates->GetYaxis()->SetTitle("p_{T} [GeV]");
    hRates->GetYaxis()->SetTitleOffset(1.2);    
    hRates->GetZaxis()->SetTitleOffset(1);
    hRates->GetZaxis()->SetTitle("Chargeflip rates");
    hRates->SetMarkerColor(kRed+1);
    hRates->DrawCopy("COLZ1,TEXT45");
    //c1->SetLogy(1);
    //c1->Update(); 
    newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.04);
    newanalyzer.GetLabel(0.186,0.81,"Chargeflip Rates",0.03);
    c1->SaveAs("Plots/"+outName+".pdf");
}
