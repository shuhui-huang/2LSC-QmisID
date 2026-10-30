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

void PlotRate1D(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TString outName,TString label);
void PlotRate2D(TH2F *hRates,TString outName);
void GetHistPlot(TString FileName, TString TreeName, TString Version,TString label);

void PlotRates(){
    SetMyStyle();

    TString path = "../run/results/Likelihood_Data_v5/";
    //TString mcpath = "OutputOverlap212560_Weighted_wBDT_Full/";
    //std::vector<TString> data={"Data","Data2015","Data2016","Data2017","Data2018"};
    std::vector<TString> data={"Data"};
    //std::vector<TString> data={"Data2015"};
    //std::vector<TString> mcversion={"mc","mc16a","mc16d","mc16e"}; 
    //std::vector<TString> mcversion={"mc16a","mc16d","mc16e"}; 
    // std::vector<TString> mcversion={"mc"}; 
    // std::vector<TString> mcversion={"mc"}; 
    // std::vector<TString> mcversion={"Zjets","ttbar"}; 
    std::vector<TString> mcversion={"Zjets"}; 
  
    TString pathData="../run/results/Likelihood_Data_v5/Data";
    TString pathmc="../run/results/Likelihood_MC_v5/";
    TString pathmctruth = "../run/results/Truth_MC_Unweighted_Full/";

    for (unsigned int i=0; i<data.size(); i++){
          GetHistPlot(pathData+"/"+data[i]+"_Default_nJets_le1.root",data[i]+"_Default_nJets_le1_all_num",data[i],"Likelihood");  
    }


  for (unsigned int i=0; i<mcversion.size(); i++){
    // TString mcdir = "/mc/";
    // if (mcversion[i]!="mc") mcdir = mcversion[i]+"/";
    // GetHistPlot(pathmc+mcdir+"mc_Default.root","mc_Default_all_num", mcversion[i]);
    // for MC truth
    if (mcversion[i]=="ttbar") {
        GetHistPlot(pathmctruth+mcversion[i]+"_Default_nJets_le1.root",mcversion[i]+"_Default_nJets_le1_all_num", mcversion[i], "mcTruth");
    }
    else if (mcversion[i]=="Zjets") {
        GetHistPlot(pathmc+mcversion[i]+"_Default_nJets_le1.root",mcversion[i]+"_Default_nJets_le1_all_num", mcversion[i], "mcLikelihood");
        // GetHistPlot(pathmctruth+mcversion[i]+"_Default_nJets_le1.root",mcversion[i]+"_Default_nJets_le1_all_num", mcversion[i], "mcTruth");
    }
  }

return;
}

void GetHistPlot(TString FileName, TString TreeName, TString Version, TString label){
    TFile *f0 = new TFile(FileName);
    std::cout << "Openfile: " << FileName << std::endl;
    TH2F  *h0 = (TH2F*)f0->Get(TreeName);
    std::cout << "Opentree: " << TreeName << std::endl;
    
    TString CRflag="";
    if      (FileName.Contains("CR0")) CRflag="CR0";
    else if (FileName.Contains("CR1")) CRflag="CR1";
    else if (FileName.Contains("CR2")) CRflag="CR2";
    else if (FileName.Contains("CR3")) CRflag="CR3";
    else if (FileName.Contains("CR4")) CRflag="CR4";
    else if (FileName.Contains("CR5")) CRflag="CR5";
    Int_t   NEtaBin = 4;
    Int_t   NPtBin  = 4;
    // Float_t EtaBinning[]= {0,0.6,1.1,1.52,1.7,2.3,2.5};
    // Float_t PtBinning[] = {0,60,90,130,2500};
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
                Rates=1e-13;
                eRates=1e-13;
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


    PlotRate1D(hPt1,hPt2,hPt3,hPt4,Version+CRflag, label);
    PlotRate2D(hRates,Version+CRflag+label);

    return ;
}

void PlotRate1D(TH1F *h1,TH1F *h2,TH1F *h3,TH1F *h4,TString outName, TString label){
    TString CRlabel="";
    if      (outName.Contains("CR0")) CRlabel="el_mll_atPV (-inf,0)";
    else if (outName.Contains("CR1")) CRlabel="el_mll_atPV (0,0.05)";
    else if (outName.Contains("CR2")) CRlabel="el_mll_atPV (0.05,0.1)";
    else if (outName.Contains("CR3")) CRlabel="el_mll_atPV (0.1,0.2)";
    else if (outName.Contains("CR4")) CRlabel="el_mll_atPV (0.2,0.3)";
    else if (outName.Contains("CR5")) CRlabel="el_mll_atPV (0.3,+inf)";

    TCanvas *c1 = new TCanvas;
    c1->cd();
    c1->SetRightMargin(0.15);
    //gStyle->SetPaintTextFormat("4.12f");
    //gStyle->SetTextFont(62);i
    gStyle->SetErrorX(0.5);
    MiniTreeAnalyzer newanalyzer;
    h1->GetXaxis()->SetTitle("|#eta|");
    h1->GetYaxis()->SetTitle("#epsilon_{QmisID}");
    h1->GetYaxis()->SetTitleOffset(1.2);
    //if (outName=="Data2015") h1->SetAxisRange(1e-12, 1, "Y");
    //if (outName.Contains("mc16a")) h1->SetAxisRange(1e-13, 1, "Y");
    if (outName.Contains("mc")) h1->SetAxisRange(1e-5, 1, "Y");
    else h1->SetAxisRange(1e-5, 1, "Y");
    h1->SetMarkerColor(kGray+3);
    h1->SetMarkerStyle(kFullTriangleUp);
    h1->SetMarkerSize(1);
    h1->SetLineColor(kGray+3);
    h1->SetLineStyle(1);
    h1->SetLineWidth(1);
    h1->Draw("E1p");    

    h2->SetMarkerColor(kRed+1);
    h2->SetMarkerStyle(kFullTriangleUp);
    h2->SetMarkerSize(1);
    h2->SetLineColor(kRed+1);
    h2->SetLineStyle(1);
    h2->SetLineWidth(1);
    h2->Draw("SAMEE1p");

    h3->SetMarkerColor(kGreen+3);
    h3->SetMarkerStyle(kFullTriangleUp);
    h3->SetMarkerSize(1);
    h3->SetLineColor(kGreen+3);
    h3->SetLineStyle(1);
    h3->SetLineWidth(1);
    h3->Draw("SAMEE1p");

    h4->SetMarkerColor(kBlue);
    h4->SetMarkerStyle(kFullTriangleUp);
    h4->SetMarkerSize(1);
    h4->SetLineColor(kBlue);
    h4->SetLineStyle(1);
    h4->SetLineWidth(1);
    h4->Draw("SAMEE1p");
    // h4->SetMarkerColor(kOrange-3);
    // h4->SetMarkerStyle(kFullTriangleUp);
    // h4->SetMarkerSize(1);
    // h4->SetLineColor(kOrange-3);
    // h4->SetLineStyle(1);
    // h4->SetLineWidth(1);
    // h4->Draw("SAMEE1p");
    c1->SetLogy(1);
    newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.04);
    newanalyzer.GetLabel(0.186,0.84,"13 TeV, 140 fb^{-1}",0.03);
    if (outName.Contains("Data")) newanalyzer.GetLabel(0.188,0.80,"Data 2015-2018",0.03);
    if (outName.Contains("Zjets")) newanalyzer.GetLabel(0.188,0.80,"Zjets mc20",0.03);
    if (outName.Contains("ttbar")) newanalyzer.GetLabel(0.188,0.80,"ttbar mc20",0.03);
    // newanalyzer.GetLabel(0.188,0.76,"w/ECIDS",0.03);
    newanalyzer.GetLabel(0.188,0.72,CRlabel,0.03);
    
    //TLegend *legend = new TLegend(0.65,0.18,0.83,0.38, outName+" Truth-matched(Zjets)");
    //TLegend *legend = new TLegend(0.55,0.75,0.8,0.88, outName+" Truth-matched(Zjets)");
    //TLegend *legend = new TLegend(0.45,0.8,0.75,0.93, outName+" Likelihood(212750)");
    // TLegend *legend = new TLegend(0.45,0.8,0.75,0.93, outName+" Likelihood");
    TString methodlabel = " Likelihood";
    if (label.Contains("Truth")) methodlabel=" Truth-matched";
    TLegend *legend = new TLegend(0.45,0.75,0.75,0.93, outName+methodlabel);

    legend->AddEntry(h1,"p_{T} #in [20, 50] GeV","lp");
    legend->AddEntry(h2,"p_{T} #in [50, 100] GeV","lp");
    legend->AddEntry(h3,"p_{T} #in [100, 200]GeV","lp");
    legend->AddEntry(h4,"p_{T} > 200 GeV","lp");
    legend->SetTextSize(0.03);
    legend->Draw("same");
    c1->SaveAs("Plots/"+outName+"_"+label+".pdf");
    c1->SaveAs("Plots/"+outName+"_"+label+".root");
    c1->SaveAs("Plots/"+outName+"_"+label+".png");
    c1->SaveAs("Plots/"+outName+"_"+label+".eps");

    return;
}

void PlotRate2D(TH2F *hRates,TString outName) {

    TCanvas *c1 = new TCanvas;
    c1->cd();
    //c1->SetLogx(0);
    //c1->SetLogz(1);
    c1->SetRightMargin(0.15);
    gStyle->SetPaintTextFormat("4.8f");
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
    c1->SetLogy(1);
    //c1->Update(); 
    newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.04);
    newanalyzer.GetLabel(0.186,0.81,"Chargeflip Rates",0.03);
    c1->SaveAs("Plots/"+outName+"_2D.pdf");
    c1->SaveAs("Plots/"+outName+"_2D.png");
    c1->SaveAs("Plots/"+outName+"_2D.root");
}
