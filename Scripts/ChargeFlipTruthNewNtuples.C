#include "TString.h"
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <fstream>
#include <algorithm> // for using copy (library function)


class qMisID{

   public:
      qMisID(TString name,TString latex,TString file,TString tree,TString extra="");
      ~qMisID();

      inline void PrintLog(bool val){printLog=val;};

      void PrintSummary();

      void InitializeSelections();
      void InitializeHistos();
      void Execute();
      TH2F *GetRatio(TH2F *num,TH2F* den); 

      inline TString GetName(){return m_name;};
      inline TString GetLatex(){return m_latex;};

      inline void SetPath(TString path){m_path=path;};
      inline void SetBaselineSelection(TString base){m_baseline_sel="("+base+")";}; 
    
      void SetBaselineSelectionSplit(TString preselection_electron_1,TString preselection_electron_2,TString preselection_positron_1,TString preselection_positron_2);
      inline void SplitRatesElecPosi(bool doSplit){m_doSplit=doSplit;};

      inline void SetExtraSelection(TString extra){m_extra_sel=extra;};
      inline void SetSelectionLepton1(TString wrong_track,TString photon_conversion){m_sel_wrongTrack_1="("+wrong_track+")"; m_sel_photonConversion_1="("+photon_conversion+")";};
      inline void SetSelectionLepton2(TString wrong_track,TString photon_conversion){m_sel_wrongTrack_2="("+wrong_track+")"; m_sel_photonConversion_2="("+photon_conversion+")";};

      inline void SetSelectionRealLepton(TString sel_real_1,TString sel_real_2){m_sel_real_1="("+sel_real_1+")"; m_sel_real_2="("+sel_real_2+")";};
      inline void SetSelectionFakeLepton(TString sel_fake_1,TString sel_fake_2,TString ss){m_sel_fake_1="("+sel_fake_1+")"; m_sel_fake_2="("+sel_fake_2+")"; m_sel_SS=ss;};
      inline void SetSelectionNoReal(TString selection_ss_lep1,TString selection_ss_lep2){m_sel_noReal_lep1="("+selection_ss_lep1+")"; m_sel_noReal_lep2="("+selection_ss_lep2+")";};

      inline void SetWeight(TString weight){m_weight=weight;};

      inline TString GetFullSelectionWrongTrack1(){return m_full_sel_wrongTrack_1;};
      inline TString GetFullSelectionWrongTrack2(){return m_full_sel_wrongTrack_2;};     
      inline TString GetFullSelectionPhotonConversion1(){return m_full_sel_photonConversion_1;};
      inline TString GetFullSelectionPhotonConversion2(){return m_full_sel_photonConversion_2;};
      inline TString GetTotalSelectionLep1(){return m_full_sel_1;};
      inline TString GetTotalSelectionLep2(){return m_full_sel_2;};

      inline TString Get_qFlip_electron_1(){return m_qflip_ele_1;};
      inline TString Get_qFlip_electron_2(){return m_qflip_ele_2;};
      inline TString Get_qFlip_positron_1(){return m_qflip_pos_1;};
      inline TString Get_qFlip_positron_2(){return m_qflip_pos_2;};
  
      inline TString Get_selection_electron_1(){return m_sel_ele_1;};
      inline TString Get_selection_electron_2(){return m_sel_ele_2;};
      inline TString Get_selection_positron_1(){return m_sel_pos_1;};
      inline TString Get_selection_positron_2(){return m_sel_pos_2;};

      //get selection per event type
      inline TString GetSelection_RR(){return m_sel_RR;};

      inline TString GetSelection_QQ(){return m_sel_QQ;};

      inline TString GetSelection_FF(){return m_sel_FF;};

      inline TString GetSelection_RnoR_1(){return m_sel_RnoR_1;};
      inline TString GetSelection_RnoR_2(){return m_sel_RnoR_2;};

      inline TString GetSelection_RQ_1(){return m_sel_RQ_1;};
      inline TString GetSelection_RQ_2(){return m_sel_RQ_2;};

      inline TString GetSelection_RF_1(){return m_sel_RF_1;};
      inline TString GetSelection_RF_2(){return m_sel_RF_2;};

      inline TString GetSelection_QF_1(){return m_sel_QF_1;};
      inline TString GetSelection_QF_2(){return m_sel_QF_2;};

      inline TString GetSelection_FF_SS(){return m_sel_FF_SS;};
      inline TString GetSelection_RF_1_SS(){return m_sel_RF_1_SS;};
      inline TString GetSelection_RF_2_SS(){return m_sel_RF_2_SS;};

      inline TString GetSelection_QF_1_SS(){return m_sel_QF_1_SS;};
      inline TString GetSelection_QF_2_SS(){return m_sel_QF_2_SS;};
      //end sel per event type...

      inline TString GetFullSelectionWrongTrack(){return m_full_sel_wrongTrack_all;};
      inline TString GetFullSelectionPhotonConversion(){return m_full_sel_photonConversion_all;};
      inline TString GetTotalSelectionLep(){return m_full_sel_all;}
 
      inline TString GetTotalBaselineSelection(){return m_full_baseline_all;};


      inline void SetVarEta(TString var1,TString var2,Float_t *binning,int n){m_var_eta1=var1; m_var_eta2=var2; m_eta_binning=binning; m_nbins_eta=n;};
      inline void SetVarPt(TString var1,TString var2,Float_t *binning,int n){m_var_pt1=var1; m_var_pt2=var2;m_pt_binning=binning; m_nbins_pt=n;};
             
      std::vector<TString> FillVector(TString list);
   
      std::vector<TH2F*> pH_Ratios;
      std::vector<TH2F*> pH_Stat;
  
      void ComputeRatios();
      void GetPlots();

   private:

      TH2F* CreateEtaPtHisto(TString label);
      void ReadInputFile(TString file);
      bool printLog;
      bool m_doSplit;


      TString m_sel_RR;
      TString m_sel_QQ;
      TString m_sel_FF;
      TString m_sel_RQ_1;
      TString m_sel_RQ_2;
      TString m_sel_RF_1;
      TString m_sel_RF_2;
      TString m_sel_QF_1;
      TString m_sel_QF_2;
      TString m_sel_FF_SS;
      TString m_sel_RF_1_SS;
      TString m_sel_RF_2_SS;
      TString m_sel_QF_1_SS;
      TString m_sel_QF_2_SS;
      TString m_sel_SS;
      TString m_sel_noReal_lep1;
      TString m_sel_noReal_lep2;
      TString m_sel_RnoR_1;
      TString m_sel_RnoR_2;


      TString m_qflip_ele_1;
      TString m_qflip_ele_2;
      TString m_qflip_pos_1;
      TString m_qflip_pos_2;

      TString m_sel_ele_1;
      TString m_sel_ele_2;
      TString m_sel_pos_1;
      TString m_sel_pos_2;

      TString m_sel_real_1;
      TString m_sel_real_2;
      TString m_sel_fake_1;
      TString m_sel_fake_2;
  

      TString m_pre_electron_1;
      TString m_pre_electron_2;

      TString m_pre_positron_1;
      TString m_pre_positron_2;

      TString m_path;
      TString m_name;
      TString m_latex; 
      TString m_file;
      TString m_tree;
      TString m_baseline_sel;
      TString m_extra_sel;
      TString m_sel_wrongTrack_1;
      TString m_sel_photonConversion_1;
      TString m_sel_wrongTrack_2;
      TString m_sel_photonConversion_2;
      TString m_var_eta1;
      TString m_var_pt1;
      TString m_var_eta2;
      TString m_var_pt2;

      TString m_weight;
      TString m_full_sel_wrongTrack_1;
      TString m_full_sel_wrongTrack_2; 
      TString m_full_sel_photonConversion_1;
      TString m_full_sel_photonConversion_2;
      TString m_full_sel_1;
      TString m_full_sel_2;
      TString m_full_sel_wrongTrack_all;
      TString m_full_sel_photonConversion_all;
      TString m_full_sel_all;
      TString m_full_baseline_all;

      Float_t *m_eta_binning;
      Float_t *m_pt_binning;
      int m_nbins_eta;
      int m_nbins_pt;

      TH2F *m_histo_all_den;
      TH2F *m_histo_all_den1;
      TH2F *m_histo_all_den2;

      TH2F *m_histo_wrongTrack_num;
      TH2F *m_histo_photonConversion_num;
      TH2F *m_histo_all_num;

      TH2F *m_histo_wrongTrack_num1;
      TH2F *m_histo_wrongTrack_num2;
      TH2F *m_histo_photonConversion_num1;
      TH2F *m_histo_photonConversion_num2;     
      TH2F *m_histo_total_num1;
      TH2F *m_histo_total_num2; 

      TH2F *m_histo_ratio_wrongTrack_1;
      TH2F *m_histo_ratio_wrongTrack_2;
      TH2F *m_histo_ratio_photonConversion_1;
      TH2F *m_histo_ratio_photonConversion_2;
      TH2F *m_histo_ratio_total_1;
      TH2F *m_histo_ratio_total_2;

      TH2F *m_histo_ratio_wrongTrack;
      TH2F *m_histo_ratio_photonConversion;
      TH2F *m_histo_all_ratio;
      TH2F *m_histo_stat_wrongTrack;
      TH2F *m_histo_stat_photonConversion;
      TH2F *m_histo_all_stat;

      TH2F *m_histo_electron_num;
      TH2F *m_histo_positron_num;

      TH2F *m_histo_electron_den;
      TH2F *m_histo_positron_den;

      TH2F *m_histo_ratio_electron;
      TH2F *m_histo_ratio_positron;
   
      TH2F *m_histo_ratio_RnoR;
      TH2F *m_histo_ratio_RF;
      TH2F *m_histo_ratio_RQ;
      TH2F *m_histo_ratio_RF_SS;


      TH2F *m_histo_RnoR_num;
      TH2F *m_histo_RR_num;
      TH2F *m_histo_QQ_num;
      TH2F *m_histo_FF_num;
      TH2F *m_histo_RF_num;  
      TH2F *m_histo_RQ_num;
      TH2F *m_histo_QF_forQflip_num;
      TH2F *m_histo_QF_forFakes_num;
      TH2F *m_histo_FF_SS_num; 
      TH2F *m_histo_RF_SS_num;
      TH2F *m_histo_QF_forFakes_SS_num;

};

qMisID MeasureChargeFlip(TString name,TString file,TString extra_cut,bool useWeight,bool useBDT,bool doToy);


TString OutputDir="Truth_MC";
TString MCversion="";
TString CRtag="";

void ChargeFlipTruth(TString process,TString binning,bool useWeight,bool useBDT,TString label,TString cut,TString period,bool doToy){

   // if (!CRtag) std::cout << "--- Please indicate CRtag, e.g. CR0" << std::endl;

   if (!useWeight) OutputDir+="_Unweighted";
   else OutputDir+="_Weighted";

   if (useBDT) OutputDir+="_wBDT";

   if (doToy) OutputDir+="_Toy";
   else OutputDir+="_Full";

   TString InputPath = "../Inputlists/";
   MCversion=period; 
   TString sample = "";
   if (MCversion.Contains("all")) sample="MC";
 
   qMisID rates= MeasureChargeFlip(process+binning+label,InputPath+process+sample+".list",cut,useWeight,useBDT,doToy);


 return;
}

qMisID MeasureChargeFlip(TString name,TString file,TString extra_cut,bool useWeight,bool useBDT,bool doToy){


 TString ori="-1"; 

 if (name.Contains("Zjets")) ori="13";
 else if (name.Contains("vjets")) ori="13";
 else if (name.Contains("ttbar")) ori="10";
 else {
   std::cout << "-- Your name does not correspond to any of the programmed processes: Zjets and ttbar  --> Aborting!"  << std::endl;
   exit(1);
 }

//  TString preselection="fabs(lep0_eta)<2.47 && fabs(lep1_eta)<2.47 && (fabs(lep0_eta)<1.37  || fabs(lep0_eta)>1.52) && (fabs(lep1_eta)<1.37  || fabs(lep1_eta)>1.52)";
 TString preselection="1";

 TString preselection_electron_1="1";
 TString preselection_positron_1="1";

 TString preselection_electron_2="1";
 TString preselection_positron_2="1";




 if (useBDT) {
      preselection+=" && (el_ECIDS[0] && el_ECIDS[1])";

      preselection_electron_1+=" && el_ECIDS[0]";
      preselection_electron_2+=" && el_ECIDS[1]";
  
      preselection_positron_1+=" && el_ECIDS[0]";
      preselection_positron_2+=" && el_ECIDS[1]";
 } 

/*
Physics if all leptons are prompt_sameQ
QmisID if 2L and exact one promt_sameQ and one trident_flipQ or prompt_flipQ
*/

 TString selection_wrongTrack_lep1="(l1_IFF==3)"; // leading lep, FlipQ=3
 TString selection_wrongTrack_lep2="(l2_IFF==3)"; // sub-leading lep, FlipQ=3

 TString selection_photonConversion_lep1="l1_pt<0"; // shuhui: dummy selection
 TString selection_photonConversion_lep2="l2_pt<0"; // shuhui: dummy selection

 TString selection_real_lep1="(l1_IFF==2)"; //leading lep, PromptEle=2
 TString selection_real_lep2="(l2_IFF==2)"; //sub-leading lep, PromptEle=2

 TString selection_fake_lep1="!("+ selection_wrongTrack_lep1+") && !("+selection_photonConversion_lep1+") && !("+selection_real_lep1+")";
 TString selection_fake_lep2="!("+ selection_wrongTrack_lep2+") && !("+selection_photonConversion_lep2+") && !("+selection_real_lep2+")"; 

 TString selection_ss_lep1="("+selection_real_lep1+")"+" && !("+selection_real_lep2+") && same_charge==1";
 TString selection_ss_lep2="("+selection_real_lep2+")"+" && !("+selection_real_lep1+") && same_charge==1";


 TString TreeName="start";
 if (name.Contains("_nJets_le1")) TreeName="qmisid_cr";
 
 TString varEta0="fabs(l1_eta)";
 TString varEta1="fabs(l2_eta)";
 TString varPt0="l1_pt*0.001"; //converted to GeV
 TString varPt1="l2_pt*0.001"; //converted to GeV

 TString path="/eos/home-s/shudong/workspace/HHML/Run3/HHMLAnalysisCode/output/ntuple_QMisID/hhml_v2/";


 TString weight="weight";

 if (useBDT) weight+="*1";
 //if (useBDT) weight+="*weight_indiv_SF_EL_ChargeMisID"; // add egamma SFs




 int nEta=4;
 //Float_t EtaBinning[]={0,0.6,1.1,1.52,1.7,2.3,2.5};
//  Float_t EtaBinning[]={0,2.5};
 Float_t EtaBinning[]={0, 1.37, 1.52, 2.0, 2.6};

 //Float_t PtBinning_Temp[]={0,60,90,130,200,2500};
 Float_t PtBinning_Temp[]={0,60,130,2600};

 //int nPt=5;
 int nPt=4;
 
//  Float_t PtBinning_Default[]={0,60,90,130,2500};
 Float_t PtBinning_Default[]={20, 50, 100, 200, 2500};

 if (name.Contains("Default")) nPt=4;
 
//  std::cout << "name= "<< name << "nPt=" << nPt << std::endl;
//  std::cout << "--Shuhui debug PtBinning_Default[2]="<< PtBinning_Default[2] << std::endl;
//  std::cout << "--Shuhui debug sizeof(PtBinning_Default)="<< sizeof(PtBinning_Default) << std::endl;
//  std::cout << "--Shuhui debug sizeof(PtBinning_Default[0])="<< sizeof(PtBinning_Default[0]) << std::endl;

 Float_t PtBinning[]={20, 50, 100, 200, 2500};

//  if (name.Contains("Default")) {
//    std::cout << "--Shuhui debug start copying" <<  std::endl;
//    // std::copy(PtBinning_Default, PtBinning_Default+sizeof(PtBinning_Default)/sizeof(PtBinning_Default[0]), PtBinning);
//    std::copy(PtBinning_Default, PtBinning_Default+2, PtBinning); // shuhui: to be fixed

//  }
//  else std::copy(PtBinning_Temp,PtBinning_Temp+sizeof(PtBinning_Temp)/sizeof(PtBinning_Temp[0]), PtBinning);

 std::cout << "EtaBinning: ";
 for (unsigned int i=0; i<nEta+1; i++) std::cout << EtaBinning[i] << ", ";
 std::cout << " " << std::endl;
 std::cout << "PtBinning: ";
 for (unsigned int i=0; i<nPt+1; i++) std::cout << PtBinning[i] << ", ";
 std::cout << " " << std::endl;

 qMisID qRates(name,name,file,TreeName);
 
 qRates.SetPath(path);

 qRates.SetBaselineSelection(preselection);

 qRates.SetBaselineSelectionSplit(preselection_electron_1,preselection_electron_2,preselection_positron_1,preselection_positron_2);
 qRates.SplitRatesElecPosi(false),

 qRates.SetSelectionLepton1(selection_wrongTrack_lep1,selection_photonConversion_lep1);
 qRates.SetSelectionLepton2(selection_wrongTrack_lep2,selection_photonConversion_lep2);

 qRates.SetSelectionRealLepton(selection_real_lep1,selection_real_lep2);
 qRates.SetSelectionFakeLepton(selection_fake_lep1,selection_fake_lep2,"same_charge==1");
 qRates.SetSelectionNoReal(selection_ss_lep1,selection_ss_lep2);

 qRates.SetExtraSelection(extra_cut);
 
 qRates.SetVarEta(varEta0,varEta1,EtaBinning,nEta);
 qRates.SetVarPt(varPt0,varPt1,PtBinning,nPt);
 qRates.PrintLog(true);

 if (useWeight) qRates.SetWeight(weight);

 qRates.Execute();



 return qRates;  
}
qMisID::qMisID(TString name,TString latex,TString file,TString tree,TString extra=""):

 m_name(name),
 m_latex(latex),
 m_file(file),
 m_tree(tree),
 m_sel_RR(""),
 m_sel_QQ(""),
 m_sel_FF(""),
 m_sel_RQ_1(""),
 m_sel_RnoR_1(""),
 m_sel_RnoR_2(""),
 m_sel_RQ_2(""),
 m_sel_RF_1(""),
 m_sel_RF_2(""),
 m_sel_QF_1(""),
 m_sel_QF_2(""),
 m_sel_FF_SS(""),
 m_sel_RF_1_SS(""),
 m_sel_RF_2_SS(""),
 m_sel_QF_1_SS(""),
 m_sel_QF_2_SS(""),
 m_sel_noReal_lep1(""),
 m_sel_noReal_lep2(""),
 m_sel_SS(""),
 m_path(""),
 m_weight(""),
 m_baseline_sel(""),
 m_extra_sel(extra),
 m_qflip_ele_1(""),
 m_qflip_ele_2(""),
 m_qflip_pos_1(""),
 m_qflip_pos_2(""),
 m_pre_electron_1(""),
 m_pre_electron_2(""),
 m_pre_positron_1(""),
 m_pre_positron_2(""),
 m_sel_ele_1(""),
 m_sel_ele_2(""),
 m_sel_pos_1(""),
 m_sel_pos_2(""),
 m_sel_wrongTrack_1(""),
 m_sel_photonConversion_1(""),
 m_sel_wrongTrack_2(""),
 m_sel_photonConversion_2(""),
 m_full_sel_wrongTrack_1(""),
 m_full_sel_wrongTrack_2(""), 
 m_full_sel_photonConversion_1(""),
 m_full_sel_photonConversion_2(""),
 m_full_sel_1(""),
 m_full_sel_2(""),
 m_full_sel_wrongTrack_all(""),
 m_full_sel_photonConversion_all(""),
 m_full_sel_all(""),
 m_full_baseline_all(""),
 m_var_eta1(""),
 m_var_eta2(""),
 m_var_pt1(""),
 m_var_pt2(""),
 m_nbins_eta(0),
 m_nbins_pt(0),
 printLog(false),
 m_histo_RR_num(0),
 m_histo_RnoR_num(0),
 m_histo_QQ_num(0),
 m_histo_FF_num(0),
 m_histo_RF_num(0),
 m_histo_RQ_num(0),
 m_histo_QF_forQflip_num(0),
 m_histo_QF_forFakes_num(0),
 m_histo_FF_SS_num(0),
 m_histo_RF_SS_num(0),
 m_histo_QF_forFakes_SS_num(0),
 m_histo_electron_num(0),
 m_histo_positron_num(0),
 m_histo_electron_den(0),
 m_histo_positron_den(0),
 m_histo_ratio_electron(0),
 m_histo_ratio_positron(0),
 m_histo_ratio_RnoR(0),
 m_histo_ratio_RQ(0),
 m_histo_ratio_RF(0),
 m_histo_ratio_RF_SS(0),
 m_histo_all_den(0),
 m_histo_all_den1(0),
 m_histo_all_den2(0),
 m_histo_wrongTrack_num(0),
 m_histo_photonConversion_num(0),
 m_histo_all_num(0),
 m_histo_ratio_wrongTrack(0),
 m_histo_ratio_photonConversion(0),
 m_histo_all_ratio(0),
 m_histo_stat_wrongTrack(0),
 m_histo_stat_photonConversion(0),
 m_histo_all_stat(0),
 m_histo_wrongTrack_num1(0),
 m_histo_wrongTrack_num2(0),
 m_histo_photonConversion_num1(0),
 m_histo_photonConversion_num2(0),
 m_histo_total_num1(0),
 m_histo_total_num2(0),
 m_histo_ratio_wrongTrack_1(0),
 m_histo_ratio_wrongTrack_2(0),
 m_histo_ratio_photonConversion_1(0),
 m_histo_ratio_photonConversion_2(0),
 m_histo_ratio_total_1(0),
 m_histo_ratio_total_2(0),
 pH_Ratios(0),
 pH_Stat(0),
 m_doSplit(false)
{}

qMisID::~qMisID()
{

 //delete m_histo_num;
 //delete m_histo_den;
 //delete m_histo_ratio;

};

void qMisID::ReadInputFile(TString file){

  TString local_file=m_path+file;
 
  TFile *pInput = new TFile(local_file);
  if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
  }
  else {
      if (printLog) cout << "-- Opening file :" << m_path+file <<  endl;
  }


 TTree* local_Tree = (TTree*)pInput->Get(m_tree);
 if (printLog) std::cout<< "-- Tree entries:" << local_Tree->GetEntries() << std::endl;


 //Creating local histos
 TH2F *local_wrongTrack_num=CreateEtaPtHisto("local_wrongTrack_num");
 TH2F *local_photonConversion_num=CreateEtaPtHisto("local_photonConversion_num");
 TH2F *local_all_num=CreateEtaPtHisto("local_all_num");

 TH2F *local_wrongTrack_num1=CreateEtaPtHisto("local_wrongTrack_num1");
 TH2F *local_wrongTrack_num2=CreateEtaPtHisto("local_wrongTrack_num2");
 TH2F *local_photonConversion_num1=CreateEtaPtHisto("local_photonConversion_num1");
 TH2F *local_photonConversion_num2=CreateEtaPtHisto("local_photonConversion_num2");
 TH2F *local_total_num1=CreateEtaPtHisto("local_total_num1");
 TH2F *local_total_num2=CreateEtaPtHisto("local_total_num2");

 TH2F *local_den1=CreateEtaPtHisto("local_all_den1");
 TH2F *local_den2=CreateEtaPtHisto("local_all_den2");


 //separating in event types: RR, RQ, QQ, QF, FF, RF
 TH2F *local_RR_num1=CreateEtaPtHisto("local_RR_num1");
 TH2F *local_RR_num2=CreateEtaPtHisto("local_RR_num2");

 TH2F *local_RQ_num1=CreateEtaPtHisto("local_RQ_num1");
 TH2F *local_RQ_num2=CreateEtaPtHisto("local_RQ_num2");

 TH2F *local_QQ_num1=CreateEtaPtHisto("local_QQ_num1");
 TH2F *local_QQ_num2=CreateEtaPtHisto("local_QQ_num2");

 TH2F *local_QF_forQflip_num1=CreateEtaPtHisto("local_QF_forQflip_num1");
 TH2F *local_QF_forQflip_num2=CreateEtaPtHisto("local_QF_forQflip_num2");

 TH2F *local_RnoR_num1=CreateEtaPtHisto("local_RnoR_num1");
 TH2F *local_RnoR_num2=CreateEtaPtHisto("local_RnoR_num2");




 //fakes must be separated in different categories: All and SS
 TH2F *local_RF_num1=CreateEtaPtHisto("local_RF_num1");
 TH2F *local_RF_num2=CreateEtaPtHisto("local_RF_num2");

 TH2F *local_FF_num1=CreateEtaPtHisto("local_FF_num1");
 TH2F *local_FF_num2=CreateEtaPtHisto("local_FF_num2");

 TH2F *local_QF_forFakes_num1=CreateEtaPtHisto("local_QF_forFakes_num1");
 TH2F *local_QF_forFakes_num2=CreateEtaPtHisto("local_QF_forFakes_num2");

 //now for ss
 TH2F *local_RF_SS_num1=CreateEtaPtHisto("local_RF_SS_num1");
 TH2F *local_RF_SS_num2=CreateEtaPtHisto("local_RF_SS_num2");

 TH2F *local_FF_SS_num1=CreateEtaPtHisto("local_FF_SS_num1");
 TH2F *local_FF_SS_num2=CreateEtaPtHisto("local_FF_SS_num2");

 TH2F *local_QF_forFakes_SS_num1=CreateEtaPtHisto("local_QF_forFakes_SS_num1");
 TH2F *local_QF_forFakes_SS_num2=CreateEtaPtHisto("local_QF_forFakes_SS_num2");




 //separing in electrons and positrons
 TH2F *local_electron_num=CreateEtaPtHisto("local_electron_num");
 TH2F *local_positron_num=CreateEtaPtHisto("local_positron_num");

 TH2F *local_electron_den=CreateEtaPtHisto("local_electron_den");
 TH2F *local_positron_den=CreateEtaPtHisto("local_positron_den");

 TH2F *local_electron_num1=CreateEtaPtHisto("local_electron_num1");
 TH2F *local_electron_num2=CreateEtaPtHisto("local_electron_num2");
 TH2F *local_positron_num1=CreateEtaPtHisto("local_positron_num1");
 TH2F *local_positron_num2=CreateEtaPtHisto("local_positron_num2");

 TH2F *local_electron_den1=CreateEtaPtHisto("local_electron_den1");
 TH2F *local_electron_den2=CreateEtaPtHisto("local_electron_den2");
 TH2F *local_positron_den1=CreateEtaPtHisto("local_positron_den1");
 TH2F *local_positron_den2=CreateEtaPtHisto("local_positron_den2");



 //Gettting selection
 TString qflip_electron_1=Get_qFlip_electron_1();
 TString qflip_electron_2=Get_qFlip_electron_2();
 TString qflip_positron_1=Get_qFlip_positron_1();
 TString qflip_positron_2=Get_qFlip_positron_2();

 TString sel_electron_1=Get_selection_electron_1();
 TString sel_electron_2=Get_selection_electron_2();
 TString sel_positron_1=Get_selection_positron_1();
 TString sel_positron_2=Get_selection_positron_2();


 TString wrongTrack1=GetFullSelectionWrongTrack1();
 TString wrongTrack2=GetFullSelectionWrongTrack2();
 
 TString photon1=GetFullSelectionPhotonConversion1();
 TString photon2=GetFullSelectionPhotonConversion2();

 TString total1=GetTotalSelectionLep1();
 TString total2=GetTotalSelectionLep2();


 //selection par event type:
 TString sel_RR = GetSelection_RR();

 TString sel_QQ = GetSelection_QQ();

 TString sel_FF = GetSelection_FF();

 TString sel_RQ_1 = GetSelection_RQ_1();
 TString sel_RQ_2 = GetSelection_RQ_2();

 TString sel_RF_1 = GetSelection_RF_1();
 TString sel_RF_2 = GetSelection_RF_2();

 TString sel_QF_1 = GetSelection_QF_1();
 TString sel_QF_2 = GetSelection_QF_2();


 TString sel_RnoR_1=GetSelection_RnoR_1();
 TString sel_RnoR_2=GetSelection_RnoR_2();

 TString sel_FF_SS = GetSelection_FF_SS();
 TString sel_RF_SS_1 = GetSelection_RF_1_SS();
 TString sel_RF_SS_2 = GetSelection_RF_2_SS();

 TString sel_QF_SS_1 = GetSelection_QF_1_SS();
 TString sel_QF_SS_2 = GetSelection_QF_2_SS();
 //end of selection per event type


 TString base=GetTotalBaselineSelection();

 if (m_weight!=""){
    wrongTrack1=m_weight+"*("+wrongTrack1+")";
    wrongTrack2=m_weight+"*("+wrongTrack2+")";
    photon1=m_weight+"*("+photon1+")";
    photon2=m_weight+"*("+photon2+")";
    total1=m_weight+"*("+total1+")";
    total2=m_weight+"*("+total2+")";
    base=m_weight+"*("+base+")";

    qflip_electron_1=m_weight+"*("+qflip_electron_1+")";
    qflip_electron_2=m_weight+"*("+qflip_electron_2+")";
    qflip_positron_1=m_weight+"*("+qflip_positron_1+")";
    qflip_positron_2=m_weight+"*("+qflip_positron_2+")";

    sel_electron_1=m_weight+"*("+sel_electron_1+")";
    sel_electron_2=m_weight+"*("+sel_electron_2+")";
    sel_positron_1=m_weight+"*("+sel_positron_1+")";
    sel_positron_2=m_weight+"*("+sel_positron_2+")";

    sel_RR = m_weight+"*("+sel_RR+")";
    sel_QQ = m_weight+"*("+sel_QQ+")";
    sel_FF = m_weight+"*("+sel_FF+")";
    sel_RQ_1 = m_weight+"*("+sel_RQ_1+")";
    sel_RQ_2 = m_weight+"*("+sel_RQ_2+")";
   
    sel_RnoR_1 = m_weight+"*("+sel_RnoR_1+")";
    sel_RnoR_2 = m_weight+"*("+sel_RnoR_2+")";

    sel_RF_1 = m_weight+"*("+sel_RF_1+")";
    sel_RF_2 = m_weight+"*("+sel_RF_2+")";
    sel_QF_1 = m_weight+"*("+sel_QF_1+")";
    sel_QF_2 = m_weight+"*("+sel_QF_2+")";
    sel_FF_SS = m_weight+"*("+sel_FF_SS+")";
    sel_RF_SS_1 = m_weight+"*("+sel_RF_SS_1+")";
    sel_RF_SS_2 = m_weight+"*("+sel_RF_SS_2+")";
    sel_QF_SS_1 = m_weight+"*("+sel_QF_SS_1+")";
    sel_QF_SS_2 = m_weight+"*("+sel_QF_SS_2+")";
 }


 //Filling histos per event category:
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_RR_num1->GetName(),sel_RR,"goff"); //Fill RR
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_RR_num2->GetName(),sel_RR,"goff"); //Fill RR

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_QQ_num1->GetName(),sel_QQ,"goff"); //Fill QQ
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_QQ_num2->GetName(),sel_QQ,"goff"); //Fill QQ

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_FF_num1->GetName(),sel_FF,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_FF_num2->GetName(),sel_FF,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_RQ_num1->GetName(),sel_RQ_1,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_RQ_num2->GetName(),sel_RQ_2,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_RnoR_num1->GetName(),sel_RnoR_1,"goff"); //Fill RnoR
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_RnoR_num2->GetName(),sel_RnoR_2,"goff"); //Fill RnoR

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_QF_forQflip_num1->GetName(),sel_QF_1,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_QF_forQflip_num2->GetName(),sel_QF_2,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_QF_forFakes_num1->GetName(),sel_QF_2,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_QF_forFakes_num2->GetName(),sel_QF_1,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_RF_num1->GetName(),sel_RF_1,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_RF_num2->GetName(),sel_RF_2,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_FF_SS_num1->GetName(),sel_FF_SS,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_FF_SS_num2->GetName(),sel_FF_SS,"goff"); //Fill FF
 
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_QF_forFakes_SS_num1->GetName(),sel_QF_SS_2,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_QF_forFakes_SS_num2->GetName(),sel_QF_SS_1,"goff"); //Fill FF

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_RF_SS_num1->GetName(),sel_RF_SS_1,"goff"); //Fill FF
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_RF_SS_num2->GetName(),sel_RF_SS_2,"goff"); //Fill FF


 m_histo_RR_num->Add(local_RR_num1);
 m_histo_RR_num->Add(local_RR_num2);

 m_histo_RnoR_num->Add(local_RnoR_num1);
 m_histo_RnoR_num->Add(local_RnoR_num2);


 m_histo_QQ_num->Add(local_QQ_num1);
 m_histo_QQ_num->Add(local_QQ_num2);

 m_histo_FF_num->Add(local_FF_num1);
 m_histo_FF_num->Add(local_FF_num2);

 m_histo_RQ_num->Add(local_RQ_num1);
 m_histo_RQ_num->Add(local_RQ_num2);

 m_histo_RF_num->Add(local_RF_num1);
 m_histo_RF_num->Add(local_RF_num2);

 m_histo_QF_forQflip_num->Add(local_QF_forQflip_num1);
 m_histo_QF_forQflip_num->Add(local_QF_forQflip_num2);

 m_histo_QF_forFakes_num->Add(local_QF_forFakes_num1);
 m_histo_QF_forFakes_num->Add(local_QF_forFakes_num2);

 m_histo_FF_SS_num->Add(local_FF_SS_num1);
 m_histo_FF_SS_num->Add(local_FF_SS_num2); 

 m_histo_RF_SS_num->Add(local_RF_SS_num1);
 m_histo_RF_SS_num->Add(local_RF_SS_num2);

 m_histo_QF_forFakes_SS_num->Add(local_QF_forFakes_SS_num1);
 m_histo_QF_forFakes_SS_num->Add(local_QF_forFakes_SS_num2);


 
 


 //Filling histos
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_wrongTrack_num1->GetName(),wrongTrack1,"goff"); //Fill charge flipped electrons just for wrong track
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_wrongTrack_num2->GetName(),wrongTrack2,"goff"); //Fill charge flipped electrons just for wrong track

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_photonConversion_num1->GetName(),photon1,"goff"); //Fill charge flipped electrons just for photon conversion
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_photonConversion_num2->GetName(),photon2,"goff"); //Fill charge flipped electrons just for photon conversion

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_total_num1->GetName(),total1,"goff"); 
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_total_num2->GetName(),total2,"goff");     

 local_wrongTrack_num->Add(local_wrongTrack_num1); //total wrong track 1+2
 local_wrongTrack_num->Add(local_wrongTrack_num2); //total wrong track

 local_photonConversion_num->Add(local_photonConversion_num1); //total photon conversion 1+2
 local_photonConversion_num->Add(local_photonConversion_num2);

 local_all_num->Add(local_wrongTrack_num);  //adding wrong track + photon conversion
 local_all_num->Add(local_photonConversion_num);

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_den1->GetName(),base,"goff"); //baseline selection
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_den2->GetName(),base,"goff");
 

 m_histo_wrongTrack_num1->Add(local_wrongTrack_num1);
 m_histo_wrongTrack_num2->Add(local_wrongTrack_num2);
 m_histo_photonConversion_num1->Add(local_photonConversion_num1); 
 m_histo_photonConversion_num2->Add(local_photonConversion_num2);

 m_histo_total_num1->Add(local_total_num1); 
 m_histo_total_num2->Add(local_total_num2);

 m_histo_wrongTrack_num->Add(local_wrongTrack_num);
 m_histo_photonConversion_num->Add(local_photonConversion_num);

 m_histo_all_num->Add(local_all_num);

 m_histo_all_den1->Add(local_den1);
 m_histo_all_den2->Add(local_den2);
 m_histo_all_den->Add(local_den1);
 m_histo_all_den->Add(local_den2);


 //filling now for separated electron and positrons:

 if (m_doSplit){
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_electron_num1->GetName(),qflip_electron_1,"goff");   // Charge flipped e- for leading
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_electron_num2->GetName(),qflip_electron_2,"goff");   // Charge flipped e- for sub-leading

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_positron_num1->GetName(),qflip_positron_1,"goff");   // Charge flipped e+ for leading
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_positron_num2->GetName(),qflip_positron_2,"goff");   // Charge flipped e+ for sub-leading

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_electron_den1->GetName(),sel_electron_1,"goff");   // all e- for leading
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_electron_den2->GetName(),sel_electron_2,"goff");   // all e- for sub-leading
 
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_positron_den1->GetName(),sel_positron_1,"goff");   // all e+ for leading
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_positron_den2->GetName(),sel_positron_2,"goff");   // all e+ for sub-leading
 

 local_electron_num->Add(local_electron_num1);
 local_electron_num->Add(local_electron_num2);
 local_positron_num->Add(local_positron_num1);
 local_positron_num->Add(local_positron_num2);

 local_electron_den->Add(local_electron_den1);
 local_electron_den->Add(local_electron_den2);
 local_positron_den->Add(local_positron_den1);
 local_positron_den->Add(local_positron_den2);


 m_histo_electron_num->Add(local_electron_num);
 m_histo_positron_num->Add(local_positron_num);
 m_histo_electron_den->Add(local_electron_den);
 m_histo_positron_den->Add(local_positron_den); 
}

 
 std::cout << "###--------- DEBUG : For e- and e+: ------" << std::endl;
 std::cout << "### DEBUG : qFlip: leading e- : " << local_electron_num1->Integral() << std::endl;
 std::cout << "### DEBUG : qFlip: sub-leading e- : " << local_electron_num2->Integral() << std::endl;
 std::cout << "### DEBUG : qFlip: e- : " << local_electron_num->Integral() << std::endl;
 std::cout << "### DEBUG : qFlip: leading e+ : " << local_positron_num1->Integral() << std::endl;
 std::cout << "### DEBUG : qFlip: sub-leading e+ : " << local_positron_num2->Integral() << std::endl;
 std::cout << "### DEBUG : qFlip: e+ : " << local_positron_num->Integral() << std::endl; 
 std::cout << "............................................................"<< std::endl;
 std::cout << "### DEBUG : All: leading e- : " << local_electron_den1->Integral() << std::endl;
 std::cout << "### DEBUG : All: sub-leading e- : " << local_electron_den2->Integral() << std::endl;
 std::cout << "### DEBUG : All: e- : " << local_electron_den->Integral() << std::endl;
 std::cout << "### DEBUG : All: leading e+ : " << local_positron_den1->Integral() << std::endl;
 std::cout << "### DEBUG : All: sub-leading e+ : " << local_positron_den2->Integral() << std::endl;
 std::cout << "### DEBUG : All: e+ : " << local_positron_den->Integral() << std::endl;
 



 std::cout << "###--------- DEBUG : For event category: --------" << std::endl;
 std::cout << "### DEBUG : RR: " << m_histo_RR_num->Integral() << std::endl;
 std::cout << "### DEBUG : QQ: " << m_histo_QQ_num->Integral() << std::endl;
 std::cout << "### DEBUG : RQ: " << m_histo_RQ_num->Integral() << std::endl;
 std::cout << "### DEBUG : QF (for qFlip): " << m_histo_QF_forQflip_num->Integral() << std::endl;
 std::cout << "### DEBUG : QF (for Fakes): " << m_histo_QF_forFakes_num->Integral() << std::endl;
 std::cout << "### DEBUG : RF: " << m_histo_RF_num->Integral() << std::endl;
 std::cout << "### DEBUG : FF: " << m_histo_FF_num->Integral() << std::endl;
 std::cout << "### DEBUG : QF SS (for Fakes): " << m_histo_QF_forFakes_SS_num->Integral() << std::endl;
 std::cout << "### DEBUG : RF SS: " << m_histo_RF_SS_num->Integral() << std::endl;
 std::cout << "### DEBUG : FF SS: " << m_histo_FF_SS_num->Integral() << std::endl;
 std::cout << "### DEBUG : RnoR: " << m_histo_RnoR_num->Integral() << std::endl;

 std::cout << "###--------- DEBUG : For ALL: --------" << std::endl;
 std::cout << "### DEBUG : track 1: " << local_wrongTrack_num1->Integral() << std::endl;
 std::cout << "### DEBUG : track 2: " << local_wrongTrack_num2->Integral() << std::endl;
 std::cout << "### DEBUG : photon 1: " << local_photonConversion_num1->Integral() << std::endl;
 std::cout << "### DEBUG : photon 2: " << local_photonConversion_num2->Integral() << std::endl;
 std::cout << "### DEBUG : total 1: " << local_total_num1->Integral() << std::endl;
 std::cout << "### DEBUG : total 2: " << local_total_num2->Integral() << std::endl;
 std::cout << "### DEBUG : total track: " << local_wrongTrack_num->Integral() << std::endl;
 std::cout << "### DEBUG : total photon: " << local_photonConversion_num->Integral() << std::endl;
 std::cout << "### DEBUG : total all num: " << local_all_num->Integral() << std::endl;
 std::cout << "### DEBUG : total all den1: " << local_den1->Integral() << std::endl;
 std::cout << "### DEBUG : total all den2: " << local_den2->Integral() << std::endl;
 std::cout << "### DEBUG : total all den: " << local_den1->Integral()+local_den2->Integral() << std::endl;



 //TCanvas *pC= new TCanvas;
 //pC->cd(1);
 //m_histo_num->Draw("COLZ,TEXT");


  return;
}




TH2F* qMisID::CreateEtaPtHisto(TString label){
 
  TH2F *hist=new TH2F(m_name+"_"+label,m_name+"_"+label,m_nbins_eta,m_eta_binning,m_nbins_pt,m_pt_binning);
  hist->GetYaxis()->SetTitle("Lepton p_{T} [GeV]");
  hist->GetXaxis()->SetTitle("Lepton #eta");
  hist->Reset();

  hist->Sumw2(); 

 return hist;
}

void qMisID::Execute(){


  std::vector<TString> samples=FillVector(m_file);


  InitializeHistos();  

  InitializeSelections();

  if (printLog) PrintSummary();
  for (unsigned int i=0; i<samples.size(); i++) ReadInputFile(samples[i]);

  ComputeRatios();
  
  GetPlots();   
  



 return;
}

void qMisID::GetPlots(){


 //  m_histo_ratio->Draw();

 //for (unsigned int i=0; i<m_nbins_pt; i++){ 


 //   TH1D *pH=m_histo_ratio->ProjectionX("bin_pt",1,)


// }


 return;
}
TH2F *qMisID::GetRatio(TH2F *num,TH2F* den){

 TH2F *pH=(TH2F*)num->Clone();

 pH->Divide(den);

 return pH;
}
void qMisID::ComputeRatios(){

  m_histo_ratio_wrongTrack_1=GetRatio(m_histo_wrongTrack_num1,m_histo_all_den1);
  m_histo_ratio_wrongTrack_2=GetRatio(m_histo_wrongTrack_num2,m_histo_all_den2);
 
  m_histo_ratio_photonConversion_1=GetRatio(m_histo_photonConversion_num1,m_histo_all_den1);
  m_histo_ratio_photonConversion_2=GetRatio(m_histo_photonConversion_num2,m_histo_all_den2);

  m_histo_ratio_total_1=GetRatio(m_histo_total_num1,m_histo_all_den1);
  m_histo_ratio_total_2=GetRatio(m_histo_total_num2,m_histo_all_den2);

 //******
 /* TH2F *m_histo_all_den_noPhotonConv=(TH2F*)m_histo_all_den->Clone();
  TH2F *m_histo_all_den_noTrack=(TH2F*)m_histo_all_den->Clone();

  m_histo_all_den_noPhotonConv->Add(m_histo_photonConversion_num,-1);
  m_histo_all_den_noTrack->Add(m_histo_wrongTrack_num,-1);

  m_histo_ratio_wrongTrack=GetRatio(m_histo_wrongTrack_num,m_histo_all_den_noPhotonConv);
  m_histo_ratio_photonConversion=GetRatio(m_histo_photonConversion_num,m_histo_all_den_noTrack); 
 */
 //*****

   // try to correct with correct denominator: below included all SS: it has to include on SS from the respective source as just above!!
  m_histo_stat_wrongTrack=(TH2F*)m_histo_wrongTrack_num->Clone("m_histo_wrongTrack_num");
  m_histo_stat_photonConversion=(TH2F*)m_histo_photonConversion_num->Clone("m_histo_photonConversion_num");
  m_histo_ratio_wrongTrack=GetRatio(m_histo_wrongTrack_num,m_histo_all_den);
  m_histo_ratio_photonConversion=GetRatio(m_histo_photonConversion_num,m_histo_all_den); 

  m_histo_all_stat=(TH2F*)m_histo_all_num->Clone("m_histo_all_num");
  m_histo_all_ratio=GetRatio(m_histo_all_num,m_histo_all_den);

  
  //ratio for e- and e+
  m_histo_ratio_electron=GetRatio(m_histo_electron_num,m_histo_electron_den);
  m_histo_ratio_positron=GetRatio(m_histo_positron_num,m_histo_positron_den);


  //ratio per event category:
  m_histo_ratio_RnoR=GetRatio(m_histo_RnoR_num,m_histo_all_den);
  m_histo_ratio_RQ=GetRatio(m_histo_RQ_num,m_histo_all_den);
  m_histo_ratio_RF=GetRatio(m_histo_RF_num,m_histo_all_den);
  m_histo_ratio_RF_SS=GetRatio(m_histo_RF_SS_num,m_histo_all_den);
 
  pH_Stat.push_back(m_histo_stat_wrongTrack);
  pH_Stat.push_back(m_histo_stat_photonConversion);
  pH_Stat.push_back(m_histo_all_stat);

  pH_Ratios.push_back(m_histo_ratio_wrongTrack_1);
  pH_Ratios.push_back(m_histo_ratio_wrongTrack_2);
  pH_Ratios.push_back(m_histo_ratio_photonConversion_1);
  pH_Ratios.push_back(m_histo_ratio_photonConversion_2);
  pH_Ratios.push_back(m_histo_ratio_total_1);
  pH_Ratios.push_back(m_histo_ratio_total_2);
  pH_Ratios.push_back(m_histo_ratio_wrongTrack);
  pH_Ratios.push_back(m_histo_ratio_photonConversion);
  pH_Ratios.push_back(m_histo_all_ratio);
  pH_Ratios.push_back(m_histo_ratio_electron);
  pH_Ratios.push_back(m_histo_ratio_positron);
  pH_Ratios.push_back(m_histo_ratio_RnoR);
  pH_Ratios.push_back(m_histo_ratio_RQ);
  pH_Ratios.push_back(m_histo_ratio_RF);
  pH_Ratios.push_back(m_histo_ratio_RF_SS);


  TString fileOutput=OutputDir+"/"+GetName()+".root";

  system("mkdir -p "+OutputDir);
  if (MCversion.Contains("mc16")) {
      system("mkdir -p "+OutputDir+"/"+MCversion);
      fileOutput=OutputDir+"/"+MCversion+"/"+GetName()+".root";
  }




  TFile *pOutput = new TFile(fileOutput,"RECREATE");


  for (unsigned int i=0; i<pH_Ratios.size(); i++){

    TCanvas *pC = new TCanvas;
    pC->cd();
    pH_Ratios[i]->Draw("COLZ,TEXT");
    pH_Ratios[i]->Write();

  }
  for (unsigned int i=0; i<pH_Stat.size(); i++){

    TCanvas *pC = new TCanvas;
    pC->cd();
    pH_Stat[i]->Draw("COLZ,TEXT");
    pH_Stat[i]->Write();

  }

  pOutput->Close();


 /* 
  m_histo_all_ratio;
*/



 /*
 
 TCanvas *pC = new TCanvas;
 pC->Divide(2,2,0,0);
 pC->cd(1);
 m_histo_num->Draw("COLZ,TEXT");

 pC->cd(2);
 m_histo_den->Draw("COLZ,TEXT");

 pC->cd(3);

 m_histo_ratio->Draw("COLZ,TEXT");

 std::cout << "-- Integral rates: " << m_histo_ratio->Integral() << std::endl;

 */
 return;
}
void qMisID::InitializeHistos(){

  m_histo_all_den=CreateEtaPtHisto("all_den"); 
  m_histo_all_den1=CreateEtaPtHisto("all_den1");
  m_histo_all_den2=CreateEtaPtHisto("all_den2");

  m_histo_electron_num=CreateEtaPtHisto("electron_num");
  m_histo_positron_num=CreateEtaPtHisto("positron_num");
  m_histo_electron_den=CreateEtaPtHisto("electron_den");
  m_histo_positron_den=CreateEtaPtHisto("positron_den");


  m_histo_wrongTrack_num=CreateEtaPtHisto("wrongTrack_num");
  m_histo_photonConversion_num=CreateEtaPtHisto("photon_conversion_num");
  m_histo_all_num=CreateEtaPtHisto("all_num");

  m_histo_ratio_wrongTrack=CreateEtaPtHisto("ratio_wrongTrack");
  m_histo_ratio_photonConversion=CreateEtaPtHisto("ratio_photonConversion");
  m_histo_all_ratio=CreateEtaPtHisto("all_ratio");

  m_histo_wrongTrack_num1=CreateEtaPtHisto("wrongTrack_num1");
  m_histo_wrongTrack_num2=CreateEtaPtHisto("wrongTrack_num2");
  m_histo_photonConversion_num1=CreateEtaPtHisto("photonConversion_num1");
  m_histo_photonConversion_num2=CreateEtaPtHisto("photonConversion_num2");
  m_histo_total_num1=CreateEtaPtHisto("total_num1");
  m_histo_total_num2=CreateEtaPtHisto("total_num2");

  m_histo_ratio_wrongTrack_1=CreateEtaPtHisto("ratio_wrongTrack_1");
  m_histo_ratio_wrongTrack_2=CreateEtaPtHisto("ratio_wrongTrack_2");
  m_histo_ratio_photonConversion_1=CreateEtaPtHisto("ratio_photonConversion_1");
  m_histo_ratio_photonConversion_2=CreateEtaPtHisto("ratio_photonConversion_2");
  m_histo_ratio_total_1=CreateEtaPtHisto("ratio_total_1");
  m_histo_ratio_total_2=CreateEtaPtHisto("tatio_total_2");


  m_histo_ratio_electron=CreateEtaPtHisto("ratio_electron");
  m_histo_ratio_positron=CreateEtaPtHisto("ratio_positron");

  m_histo_ratio_RnoR=CreateEtaPtHisto("ratio_RnoR");
  m_histo_ratio_RQ=CreateEtaPtHisto("ratio_RQ");
  m_histo_ratio_RF=CreateEtaPtHisto("ratio_RF");
  m_histo_ratio_RF_SS=CreateEtaPtHisto("ratio_RF_SS");

  //histos per event category
  m_histo_RR_num=CreateEtaPtHisto("histo_RR_num");
  m_histo_QQ_num=CreateEtaPtHisto("histo_QQ_num");
  m_histo_RQ_num=CreateEtaPtHisto("histo_RQ_num");
  m_histo_QF_forQflip_num=CreateEtaPtHisto("histo_QF_forQflip_num");
  m_histo_RnoR_num=CreateEtaPtHisto("histo_RnoR_num"); 
 
  m_histo_FF_num=CreateEtaPtHisto("histo_FF_num");
  m_histo_RF_num=CreateEtaPtHisto("histo_RF_num");
  m_histo_QF_forFakes_num=CreateEtaPtHisto("histo_QF_forFakes_num");
  
  m_histo_FF_SS_num=CreateEtaPtHisto("histo_FF_SS_num");
  m_histo_RF_SS_num=CreateEtaPtHisto("histo_RF_SS_num");
  m_histo_QF_forFakes_SS_num=CreateEtaPtHisto("histo_QF_forFakes_SS_num");
  
  




 return;
}

void qMisID::SetBaselineSelectionSplit(TString preselection_electron_1,TString preselection_electron_2,TString preselection_positron_1,TString preselection_positron_2){

  m_pre_electron_1="("+preselection_electron_1+")"; 
  m_pre_electron_2="("+preselection_electron_2+")";

  m_pre_positron_1="("+preselection_positron_1+")";
  m_pre_positron_2="("+preselection_positron_2+")";


 return;
}
void qMisID::InitializeSelections(){


 m_full_baseline_all=m_baseline_sel;
 if (m_extra_sel != "") {
     m_full_baseline_all="("+m_baseline_sel+" && ("+m_extra_sel+"))";

     m_sel_ele_1="("+m_pre_electron_1+" && ("+m_extra_sel+"))";
     m_sel_ele_2="("+m_pre_electron_2+" && ("+m_extra_sel+"))";
     m_sel_pos_1="("+m_pre_positron_1+" && ("+m_extra_sel+"))";
     m_sel_pos_2="("+m_pre_positron_2+" && ("+m_extra_sel+"))";
 }


     
 m_full_sel_wrongTrack_1=m_full_baseline_all+" && "+m_sel_wrongTrack_1;
 m_full_sel_wrongTrack_2=m_full_baseline_all+" && "+m_sel_wrongTrack_2;

 m_full_sel_photonConversion_1=m_full_baseline_all+" && "+m_sel_photonConversion_1;
 m_full_sel_photonConversion_2=m_full_baseline_all+" && "+m_sel_photonConversion_2;

 m_full_sel_1=m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+")";
 m_full_sel_2=m_full_baseline_all+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";


 //for electrons and positrons:
  m_qflip_ele_1=m_sel_ele_1+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+")";
  m_qflip_ele_2=m_sel_ele_2+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";

  m_qflip_pos_1=m_sel_pos_1+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+")";
  m_qflip_pos_2=m_sel_pos_2+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";

  m_sel_RR = m_full_baseline_all+" && ("+m_sel_real_1+" && "+m_sel_real_2+")";
  m_sel_QQ = m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+") && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";
  m_sel_FF = m_full_baseline_all+" && ("+m_sel_fake_1+" && "+m_sel_fake_2+")";

  m_sel_RQ_1 = m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+") && ("+m_sel_real_2+")";
  m_sel_RQ_2 = m_full_baseline_all+" && "+m_sel_real_1+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";  

  m_sel_RF_1 = m_full_baseline_all+" && ("+m_sel_fake_1+" && "+m_sel_real_2+")";
  m_sel_RF_2 = m_full_baseline_all+" && ("+m_sel_real_1+" && "+m_sel_fake_2+")";
  // shuhui: Qflip events
  m_sel_QF_1 = m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+") && ("+m_sel_fake_2+")";
  m_sel_QF_2 = m_full_baseline_all+" && "+m_sel_fake_1+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";

  m_sel_FF_SS = m_sel_FF+" && "+m_sel_SS;
  m_sel_RF_1_SS = m_sel_RF_1+" && "+m_sel_SS;
  m_sel_RF_2_SS = m_sel_RF_2+" && "+m_sel_SS;

  m_sel_QF_1_SS = m_sel_QF_1+" && "+m_sel_SS;
  m_sel_QF_2_SS = m_sel_QF_2+" && "+m_sel_SS;


  m_sel_RnoR_1= m_full_baseline_all+" && ("+m_sel_noReal_lep1+")";
  m_sel_RnoR_2= m_full_baseline_all+" && ("+m_sel_noReal_lep2+")";


   

 //m_full_sel_wrongTrack_all=m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_wrongTrack_2+")";
 //m_full_sel_photonConversion_all=m_full_baseline_all+" && ("+m_sel_photonConversion_1+" || "+m_sel_photonConversion_2+")";

 //m_full_sel_all=m_full_baseline_all+" && ("+m_sel_photonConversion_1+" || "+m_sel_photonConversion_2+" ||"+m_sel_wrongTrack_1+" || "+m_sel_wrongTrack_2+")";


  return;
}
void qMisID::PrintSummary(){

 std::cout << "### Printing settings: " << std::endl;
 std::cout << "- Preselection: " << GetTotalBaselineSelection() << std::endl;
 std::cout << "- Extra cut: " << m_extra_sel << std::endl;
 std::cout << "- Sel. lep 1: " << GetTotalSelectionLep1() << std::endl;
 std::cout << "- Sel. lep 2: " << GetTotalSelectionLep2() << std::endl;

 std::cout << "- Selection RnoR (1) :" << m_sel_RnoR_1 << std::endl;
 std::cout << "- Selection RnoR (2) :" << m_sel_RnoR_2 << std::endl;


 return;
}
std::vector<TString> qMisID::FillVector(TString file){

  std::vector<TString> samples;

  std::ifstream infile(file.Data());

  if (!infile){
        std::cout <<"-- File with the list of samples: '"<< file <<"' does not exist --- ABORTING " << std::endl;
        exit(1);
  }

  std::string line;


  while (std::getline(infile,line) ){

        if (line.empty()) continue;

        TString ss = line;

        samples.push_back(ss);
  }

 return samples;

}
