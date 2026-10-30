#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"


#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/Yields.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/Channel.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/EventCut.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/VariableDistr2D.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/VariableDistr.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/PhysicsSample.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/PhysicsProcess.C"

#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/MiniTreeAnalyzer.C"
#include "/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/MyStyle.C"

std::vector<int> col={13,kRed,kAzure+1,kOrange+1,kViolet+2,kGreen+3};
//std::vector<int> col={2,8,kOrange+1,6,46,13};
std::vector<int> marker{20,22,23,24,25,26,27,28};
//std::vector<int> fillStyle={0,0,0,0,0};
std::vector<int> fillStyle={3004,3005,3004,3005,3004,3005,3004};


class qHisto{

     public:
         qHisto(TH2F* pH);
         ~qHisto();
         inline double GetMin(){return min;};
         inline double GetMax(){return max;};
         inline TH2F* GetHisto(){return m_pH;};

         void SetValues();
         void SetMin(float m_min){min=m_min;};
         void SetMax(float m_max){max=m_max;};

     private:

         double GetValue(TString str,TString out);
         double min;
         double max;
         TH2F *m_pH;

};

class qWeight{

    public:
      qWeight(double w,double up=0,double d=0);
      ~qWeight(); 

      inline void SetUp(double w){w_up=w;};
      inline void SetDown(double w){w_down=w;};
      inline double GetUp(){return w_up;};
      inline double GetDown(){return w_down;};
      inline double GetNominal(){return w_nom;};     


    private:
      double w_nom;
      double w_up;
      double w_down;


};
class qValidate{

   public:

        qValidate(TString filelist,TString tree);
        ~qValidate();

        std::vector<TString> FillVector(TString file);
        void Execute();
        inline void SetTreeName(TString tree){m_tree=tree;};
        inline void SetPathToFile(TString path){m_path=path;};
        void AddDependence(TString name,TString latex,std::map<TString,std::vector<TString> > rates);
        void AddDependence(TString name,TString latex,std::vector<TString> rates);
        inline void AddVariable(VariableDistr *var){v_var.push_back(var);};
        inline void AddBkgComposition(std::map<TString,float> bkgComp){m_bkg_comp=bkgComp;};
        inline void SetTag(TString tag){m_tag=tag;};  
        inline void SetInputRates(TString input){File_wRates=input;};
        inline void SetMCversion(TString mcversion){m_mcversion=mcversion;};
        inline void SetApplyBDT(bool doBDT){m_doBDT=doBDT;};
        inline void SetProcess(TString process){m_process=process;};
        inline void SetOutputDir(TString output){m_outputDir=output;}; 
        inline void SetData(bool isData){m_isData=true;};

   private:

       void InitializeHistos();
       void ReadInputFile(TString input);
       void FillTGraphAsymm();
       double GetWeight(double e1,double e2);
       qWeight* ComputeWeight(double el_pt1,double el_pt2,double el_eta1,double el_eta2,std::vector<qHisto*> local_rates,double var=-1);
       qWeight* GetRateReweighting(double el_pt1,double el_pt2,double el_eta1,double el_eta2,TString dependence,double var=-1);
       std::vector<qHisto*> GetHistos2D(std::vector<TString> rates); 
       void PlotDistributions();
       TH2F* ChooseHisto(std::vector<qHisto*> histos,double var);
       void SetHistoStyle(TH1D* &pH,int col,int width,int style,int marker=1);
       void SetGraphStyle(TGraphAsymmErrors* &pG,int col,int width,int style,int marker,int fill);

       std::vector<double> GetEpsilon(TH2F *pH,double el_pt,double el_eta);
       void GetYieldsEstimation();
       TGraphAsymmErrors* GetTGraphAsym(TH1D *nom,TH1D *up,TH1D *down);
       TGraphAsymmErrors* GetRatio(TH1D *pNum,TGraphAsymmErrors *pDen);
       double GetError(double num,double num_error,double den,double den_error);

       std::vector<VariableDistr*> v_var;
       std::map<TString,TString> latex_dependence;
       std::map<TString,std::vector<qHisto*>> m_dependence;    
       std::map<TString,std::vector<qHisto*>> m_method_A;
       std::map<TString, float> m_bkg_comp;

       TString m_tag;
       TString m_list;
       TString m_tree;
       TString m_mcversion;
       bool m_isData;
       bool m_doBDT;
       TString m_path;
       TString m_outputDir;
       TString m_process;
       std::map<TString,TH1D*> m_histo_SS;
       std::map<TString,TH1D*> m_histo_OS;     
       std::map<TString,TH1D*> m_histo_OSUp;
       std::map<TString,TH1D*> m_histo_OSDown;
       std::map<TString,TGraphAsymmErrors*> m_graph_OS;
       float Zjets_comp;
       float ttbar_comp;
       TString File_wRates;
};

bool doRealSS = false;
bool doplus = false;
TString doDep="";

void ClosureTest(TString process,TString inputRates,TString MCversion,TString var,TString outputDir, bool testRealSS, bool testplus, TString TestDep=""){
    doRealSS = testRealSS;
    std::cout << "doRealSS: " << doRealSS << std::endl;
    doplus = testplus;
    std::cout << "doplus Ptbinning: " << doplus << std::endl;
    if (TestDep!="") doDep = TestDep;
    std::cout << "doDep: " << doDep << std::endl;
    gROOT->LoadMacro("/afs/cern.ch/work/s/shuhui/Run3/GTA/Tools/MyStyle.C");
    SetMyStyle();

//    VariableDistr *v_pt1= new VariableDistr("el_pt[0]*0.001","PtLep","Leading lepton p_{T} [GeV]","Events /",10,0,250,false,false,true);
    double PtBinningnew[]={20,60,90,130,200,500};// shuhui: for new pt binning
   //  double PtBinning[]={28,60,90,130,500};// shuhui: for pt binning
    double PtBinning[]={20,50,100,200,2600};// shuhui: for pt binning
    VariableDistr *v_pt1= new VariableDistr("lep0_pT*0.001","PtLep","Leading lepton p_{T} [GeV]","Events /",4, PtBinning);
    VariableDistr *v_pt1new= new VariableDistr("lep0_pT*0.001","PtLepnew","Leading lepton p_{T} [GeV]","Events /",20,0,500, false,false,false);
   //  VariableDistr *v_pt1Std= new VariableDistr("el_pt[0]*0.001","PtLepStd","Leading lepton p_{T} [GeV]","Events /",20,0,500,false,false,false);

    VariableDistr *v_pt2= new VariableDistr("lep1_pT*0.001","PtSublep","Subleading lepton p_{T} [GeV]","Events /",20,0,500,false,false,false);

    VariableDistr *v_Mll= new VariableDistr("m_l0l1*0.001","Mll","M_{ee} [GeV]","Events /",40,71,111,false,false,false);
    //VariableDistr *v_eta1= new VariableDistr("el_eta[0]","EtaLep","Leading lepton #eta","Events /",10,-2.5,2.5,false,false,false);
    double EtaBinning[]={0,1.37,1.52,2.0,2.5}; // shuhui
    VariableDistr *v_eta1= new VariableDistr("fabs(el_eta[0])","EtaLep","Leading lepton |#eta|","Events /", 4, EtaBinning); 
    VariableDistr *v_eta1new = new VariableDistr("fabs(el_eta[0])","EtaLepnew","Leading lepton |#eta|","Events /", 10,0,2.5,false,false,false); 
    

    //double EtaBinning[]={-2.5,-2.3,-1.7,-1.52,-1.1,-0.6,0,0.6,1.1,1.52,1.7,2.3,2.5}; // shuhui

    TString path="/eos/home-s/shudong/workspace/HHML/Run3/HHARD/gen_ntuples/HHML_2lSC_v2_QmisID_ntuples/";
   
    TString bin_confi="Default";
    if(doplus)  bin_confi="plus";


  std::vector<TString> Ht_dep;
  std::vector<TString> njets_dep;
  std::vector<TString> nBjets_dep;
  std::vector<TString> Met_dep;
  std::vector<TString> Htjets_dep;
  std::vector<TString> Htlep_dep;
  std::vector<TString> PV_dep;
  std::vector<TString> mu_dep;
  std::vector<TString> elMll_dep;

  std::vector<TString> standard;
  std::vector<TString> standard_new;
  std::vector<TString> standard_new_data;


  if (process=="Zjets"){

     //standard={"Zjets_Default_nJets_ge0"};
     standard={MCversion+"_"+bin_confi+"_nJets_le1"}; // shuhui: rates to be read in
     standard_new={"Zjets_plus_nJets_ge0"};
    
     Ht_dep={"Zjets_"+bin_confi+"_Ht_ge0_le300","Zjets_"+bin_confi+"_Ht_ge300_le600","Zjets_"+bin_confi+"_Ht_ge600_le900","Zjets_"+bin_confi+"_Ht_ge900"};
     njets_dep={"Zjets_"+bin_confi+"_nJets_le4","Zjets_"+bin_confi+"_nJets_ge5"};

  }

  else if (process=="ttbar"){
     standard_new={"ttbar_"+bin_confi+"_nJets_le1"};
     //standard={"ttbar_Default_nJets_ge0"};
     standard={MCversion+"_"+bin_confi+"_nJets_ge0"}; // shuhui: rates to be read in
     Ht_dep={"ttbar_"+bin_confi+"_Ht_ge0_le300","ttbar_"+bin_confi+"_Ht_ge300_le600","ttbar_"+bin_confi+"_Ht_ge600"};
  }
  else if (process.Contains("mc16") || process.Contains("Data201")) standard={MCversion+"_"+bin_confi+"_nJets_ge0"};
  //else if (process.Contains("mc16a") || process.Contains("Data201")) standard={"mc_"+bin_confi+"_nJets_ge0"}; // shuhui for test
  else if (process=="Data" || process=="mc"){
      standard={MCversion+"_"+bin_confi+"_nJets_le1"}; // shuhui: rates to be read in 
      elMll_dep={"Data_"+bin_confi+"_nJets_ge0_CR0","Data_"+bin_confi+"_nJets_ge0_CR1","Data_"+bin_confi+"_nJets_ge0_CR2","Data_"+bin_confi+"_nJets_ge0_CR3","Data_"+bin_confi+"_nJets_ge0_CR4","Data_"+bin_confi+"_nJets_ge0_CR5"}; // shuhui: rates to be read in
  }

     TString label="";
     if (bin_confi=="plus") label="(+)";
     TString tag="Z+jets";
     if (process=="ttbar") tag="t#bar{t}";

     TString process_list=process+".list";

     //if (process.Contains("mc")) process_list="ZjetsMC_mc16a.list";//shuhuio: for mc16a test
     if (process.Contains("Zjets")) process_list="ZjetsMC.list";
     if (process.Contains("Data")) process_list="Data.list";
     if (process.Contains("2015")) process_list="Data2015.list";
     if (process.Contains("2016")) process_list="Data2016.list";
     if (process.Contains("2017")) process_list="Data2017.list";
     if (process.Contains("2018")) process_list="Data2018.list";
     process_list="../Inputlists_v0/"+process_list;

     // Initialization fot ValidationPlots 
     qValidate ValidationPlots(process_list,"tree");
     ValidationPlots.SetPathToFile(path);
     ValidationPlots.SetTag(tag);
     ValidationPlots.SetInputRates(inputRates);
     ValidationPlots.SetMCversion(MCversion);
     ValidationPlots.SetProcess(process);
     ValidationPlots.SetOutputDir(outputDir);
     if (process.Contains("Data")) ValidationPlots.SetData(true);
     if (inputRates.Contains("wBDT")) ValidationPlots.SetApplyBDT(true);

     //lines below are just for closure in both data and MC.
     //if (process=="Data" || process=="Zjets" || process=="mc" || process=="ttbar" || process=="vjets"){
     if (process=="vjets"){
        TString meth = "Truth.";
        if (MCversion.Contains("mc") || MCversion.Contains("Data")|| inputRates.Contains("Likelihood")) {
            meth = " Like."; // shuhui: MCversion is the rates to be used
         ValidationPlots.AddDependence("Standard","",standard);
         //ValidationPlots.AddDependence("Standard","ttbar Truth.",standard);
          //ValidationPlots.AddDependence("Standard",MCversion+meth,standard);
        //  ValidationPlots.AddDependence("Standard"," Likelihood",standard);
        //ValidationPlots.AddDependence("elMll",  "elMeePVCO" ,elMll_dep);// for CRttbarCO rate
          }
          if (bin_confi!="plus" && inputRates.Contains("Truth")){
          ValidationPlots.AddDependence("Standard", "Standard" ,standard); // shuhui
          //ValidationPlots.AddDependence("Ht",  "Ht dep." ,Ht_dep); // shuhui
          //ValidationPlots.AddDependence("nBjets", "nBjets dep." ,nBjets_dep); // shuhui
          }
          else if (bin_confi=="plus" && inputRates.Contains("Truth")) {
            ValidationPlots.AddDependence("Standard", "" ,standard);
          //  ValidationPlots.AddDependence("Ht",  "Ht dep." ,Ht_dep);
          }
         //ValidationPlots.AddDependence("Standard","Data Lik.",standard);
         //ValidationPlots.AddDependence("Def_New","Data Lik.",standard_new); // shuhui

     }
     else if (process=="Data"){
        ValidationPlots.AddDependence("Standard","",standard); // for CRttW closure or nomial closure
     //   ValidationPlots.AddDependence("elMll",  "elMeePVCO" ,elMll_dep);// for CRttbarCO rate

     }
     else ValidationPlots.AddDependence("Standard","",standard); //shuhui

     if (var=="all"){
      ValidationPlots.AddVariable(v_pt1);
      ValidationPlots.AddVariable(v_pt2);
      ValidationPlots.AddVariable(v_Mll);     
      ValidationPlots.AddVariable(v_eta1);
      ValidationPlots.AddVariable(v_pt1new);
      ValidationPlots.AddVariable(v_eta1new);


     }
 
     if (var=="pt1") ValidationPlots.AddVariable(v_pt1);

    

     ValidationPlots.Execute();

 return;
}
qValidate::qValidate(TString filelist,TString tree):
  m_list(filelist),
  m_tree(tree),
  m_path(""),
  v_var(),
  m_tag(""),
  m_bkg_comp(),
  latex_dependence(),
  File_wRates(""),
  m_method_A(),
  m_mcversion(""),
  m_process(""),
  m_outputDir(""),
  m_doBDT(false),
  m_isData(false),
  m_dependence(),
  m_histo_SS(),
  m_histo_OS(),
  m_histo_OSUp(),
  m_histo_OSDown(),
  m_graph_OS()
{
}
qValidate::~qValidate()
{}
void qValidate::ReadInputFile(TString input){


  std::cout << "-- Reading input file: " << m_path+input << std::endl;

  TFile *pFile = new TFile(m_path+input);

  TTree *pTree=(TTree*)pFile->Get(m_tree);

  int nEvents=pTree->GetEntries();
   nEvents=200000; // shuhui: for quick debug

  std::cout << "-- Total events: " << nEvents << std::endl;

  //setting vars

 float el0_pt=-1;
 float el1_pt=-1;
 float el0_eta=-1;
 float el1_eta=-1;
 bool isSS=0;
 bool isOS=0;
 char isel0_Tight=0;
 char isel1_Tight=0;

 double toGeV=0.001;

 double Mll01=-1;

 double totweight=1;


  //getting branches
 pTree->SetBranchAddress("lep0_pT",&el0_pt);
 pTree->SetBranchAddress("lep1_pT",&el1_pt);
 pTree->SetBranchAddress("lep0_eta",&el0_eta);
 pTree->SetBranchAddress("lep1_eta",&el1_eta);
 pTree->SetBranchAddress("lep0_Tight",&isel0_Tight);
 pTree->SetBranchAddress("lep1_Tight",&isel1_Tight);
 pTree->SetBranchAddress("same_charge",&isSS);
 pTree->SetBranchAddress("opposite_charge",&isOS);
 pTree->SetBranchAddress("m_l0l1",&Mll01);

  //getting weights
  if (!m_isData){
   pTree->SetBranchAddress("weight_NOSYS",&totweight);  
 }


//weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise

  double nZevent = 0;
  double nZevent2 = 0;
  double nSScount = 0;
  double nOScount = 0;
  for (unsigned int i=0; i<nEvents; i++){

        pTree->GetEntry(i);
       
        //std::cout << "-- print pt: " << (*el_charge)[0] << ", ... "<< (*el_charge)[1] << std::endl;    
        if (i%100000 ==0) std::cout << "-- event: " << i << std::endl;

        float el_pt1=el0_pt*toGeV;
        float el_pt2=el1_pt*toGeV;
        float el_eta1=fabs(el0_eta);
        float el_eta2=fabs(el1_eta);
        int el_charge1=0;
        int el_charge2=0;

        float Mll=Mll01*toGeV;
        float event_weight=1.0;

         if((bool)isSS){
            el_charge1=1;
            el_charge2=1;
            // std::cout << "--debug: isSS pair" << std::endl;
         } 
         else if ((bool)isOS){
            el_charge1=1;
            el_charge2=-1;
            // std::cout << "--debug: isOS pair" << std::endl;

         }


        
     //std::cout << "-- print pt: " << OSee << ", ... "<< loose_SSee << std::endl;
   //   if (!isOS) continue;
   //   if (nMuons>0) continue;     

     //if (nJets<1) continue; //request at least 1 jet for closure. 

    if (m_process!="ttbar"){
        if (Mll<81 || Mll>101) continue; //DO CLOSURE ONLY FOR Z-events (remember that other data events are fully blinded!).
        nZevent++;
    }
  
//             if (!(el_eta1<2.47) || !(el_eta2<2.47)) continue;
//         if ( (el_eta1>1.37 && el_eta1 <1.52) ) continue;
//         if ( (el_eta2>1.37 && el_eta2 <1.52) ) continue;

    //if ( !(((Tlepton_0_ID==11 && (lepton_0_mll_atConV>0.1 || lepton_0_mll_atConV<0))||(Tlepton_0_ID==13)) && ((Tlepton_1_ID==11 && (lepton_1_mll_atConV>0.1 || lepton_1_mll_atConV<0)) || (Tlepton_1_ID==13))) ) continue; // pass if falls in CRttW selection
      float Ht = 1.0; // shuhui: set dummy value
     if (Ht>0){  //Standard
        double event_weight=1;

        if (!m_isData) {
               //event_weight=weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise*(36184.86*(runNumber == 284500) + 43587.3*(runNumber == 300000) + 45691.0*(runNumber == 310000));
               event_weight=totweight; 
               //event_weight=1; 

         } 

         nZevent2++;  
         // std::cout << "-- print pt: " << el_pt1 << ",  "<< el_pt2 << std::endl;
         //std::cout << "-- print eta: " << el_eta1 << ", ... "<< el_eta2 << std::endl; 
         //if (el_pt2>200) std::cout << "-- print pt: " << el_pt1 << ", ... "<< el_pt2 << std::endl;
         //std::cout << "-- Print mu: " << mu << std::endl;
         //if (Ht>500) std::cout <<"######################################## DEBUG ############## " << std::endl;

        for (unsigned int v=0; v<v_var.size(); v++){

                double var_to_fill=0;
               //  std::cout << "-- in var :" << v_var[v]->GetTitle() << std::endl;
               //  std::cout << "v_var[v]->GetTitle().Contains=" << v_var[v]->GetTitle().Contains("PtLep") << endl;

                if (v_var[v]->GetTitle().Contains("PtLep")) {var_to_fill=el_pt1;}
                else if (v_var[v]->GetTitle().Contains("PtSublep")) var_to_fill=el_pt2;
                else if (v_var[v]->GetTitle().Contains("EtaLep")) var_to_fill=fabs(el_eta1); //shuhui: for EtaBinning [0,2.5]
                else if (v_var[v]->GetTitle().Contains("Mll"))    var_to_fill=Mll;

                else {
                    std::cout << "## Not sure which variable has to be filled... the code is not automatic: the variables to fill have to be setted by hand! -> Aborting." << std::endl;
                    exit(1);
                }
               //  std::cout << "-- var to fill: " << var_to_fill << std::endl;        

                if ((el_charge1*el_charge2 > 0)){
                    // std::cout << "SS event -> Mll :" << Mll << std::endl; 
                     
                     m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Fill(var_to_fill,event_weight);
                    //nSScount+=event_weight;
                    nSScount++;
                }
                else if ( el_charge1*el_charge2 < 0){
                    for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){      
                                        
                        double var_to_pass=0;

                        if (it->first.Contains("Standard")) var_to_pass=var_to_fill;
                        //else if (it->first.Contains("methodA")) var_to_pass=var_to_fill;   --> methodA_Ht
                        else {
                            std::cout << "## Not sure which variable has to be passed: Dependence is probably unknown or not added... the code is not automatic: the variables to fill have to be setted by hand! -> Aborting." << std::endl;
                            exit(1);
                         }
                        qWeight* chWeight=0;
                        
                        // if (it->first.Contains("methodA")) chWeight=GetRateReweighting(el_MeePV1,el_MeePV2,el_pt1,el_pt2,it->first,var_to_pass);
                        // else 
                        chWeight=ComputeWeight(el_pt1,el_pt2,el_eta1,el_eta2,m_dependence[it->first],var_to_pass);

                        double nom=chWeight->GetNominal();
                        //double nom=1;
                        double up=chWeight->GetUp();
                        double down=chWeight->GetDown();

                        m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Fill(var_to_fill,event_weight*nom);
                        m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first]->Fill(var_to_fill,event_weight*up);                         
                        m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]->Fill(var_to_fill,event_weight*down);
                        if (i<=10) {
                           // std::cout << "---DEBUG: runNumber " << runNumber << std::endl;
                            std::cout << "--- Shuhui DEBUG: chWeight_nom: " << nom << std::endl;
                            std::cout << "-- print pt: " << el_pt1 << ", ... "<< el_pt2 << std::endl;
                            std::cout << "-- print eta: " << el_eta1 << ", ... "<< el_eta2 << std::endl; 
                           // std::cout << "-- print mllatPV: " << el_MeePV1 << ", ... "<< el_MeePV2 << std::endl; 
                        }
                        //std::cout << "---DEBUG: chWeight_nom: " << nom << std::endl;
                        nOScount+=event_weight*nom;
                     }//filling OS weighted for different dependencies
                }
                else {
                //     std::cout << "-- Electron charge has not be filled: " << el_charge1 << " : " << el_charge2 << std::endl;

                } 

        }//filling histos for different vars  


     }//events cut!

  }//end loop over events
  std::cout<< "nZevent is: " << nZevent << std::endl;
  std::cout<< "nZevent2 is: " << nZevent2 << std::endl;
  std::cout<< "nSSevent count is: " << nSScount << std::endl;
  std::cout<< "nOSevent count is: " << nOScount << std::endl;
  for (unsigned int v=0; v<v_var.size(); v++){
      double nSSevent = m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Integral();
      for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){
          double nOSevent = m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Integral();
          std::cout << "---shuhui---:" << "In " << v_var[v]->GetTitle() << ": SSevent=" << nSSevent << " ; " << "OSevent=" << nOSevent << std::endl;
      }
  } 

  // to include overflow and underflow events
  for (unsigned int v=0; v<v_var.size(); v++){
      int last=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetNbinsX();
      double over=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(last+1);
      double under=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(0);
      m_histo_SS[v_var[v]->GetTitle()+"_SS"]->SetBinContent(last,m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(last)+over);
      m_histo_SS[v_var[v]->GetTitle()+"_SS"]->SetBinContent(last+1,0);
      //m_histo_SS[v_var[v]->GetTitle()+"_SS"]->SetBinContent(1,m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(1)+under);
      //m_histo_SS[v_var[v]->GetTitle()+"_SS"]->SetBinContent(0,0);
  }

  for (unsigned int v=0; v<v_var.size(); v++){
      for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){
           int last=m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->GetNbinsX();
           double over=m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->GetBinContent(last+1);
           double under=m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->GetBinContent(0);
           m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->SetBinContent(last,m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->GetBinContent(last)+over);
           m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->SetBinContent(last+1,0);
           //m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->SetBinContent(1,m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->GetBinContent(1)+under);
           //m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->SetBinContent(0,0);
      }
  }
  for (unsigned int v=0; v<v_var.size(); v++){
      double nSSevent = m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Integral();
      for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){
          double nOSevent = m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Integral();
          //std::cout << "---shuhui---: After inlcude overflow and underflow" << "In " << v_var[v]->GetTitle() << ": SSevent=" << nSSevent << " ; " << "OSevent=" << nOSevent << std::endl;
          std::cout << "---shuhui---: After inlcude overflow" << "In " << v_var[v]->GetTitle() << ": SSevent=" << nSSevent << " ; " << "OSevent=" << nOSevent << std::endl;
      }
  } 

  //std::cout<< "Add weight to Zevent: " << nZevent*event_weight <<  std::endl;
  return;
}
void qValidate::Execute(){

   std::vector<TString> samples=FillVector(m_list);

   if (m_dependence.count("methodA")==1 && (Zjets_comp <0 || ttbar_comp <0)){
     std::cout << " Method A habilitated but background composition not given!! --> Aborting !"<< std::endl;
     exit(1);
    } 

   InitializeHistos();  

   for (unsigned int i=0; i<samples.size(); i++) ReadInputFile(samples[i]);   

   FillTGraphAsymm(); 
  
   PlotDistributions();
    
   GetYieldsEstimation();

 return;
}
void qValidate::FillTGraphAsymm(){

 for (unsigned int v=0; v<v_var.size(); v++){ 

       for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){


           m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]=GetTGraphAsym(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first],
                                                                           m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first],
                                                                           m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]);

       } //method
 }//end var

 return;
}
void qValidate::SetGraphStyle(TGraphAsymmErrors* &pG,int col,int width,int style,int marker,int fill){

  pG->SetLineColor(col);
  pG->SetLineWidth(width);
  pG->SetMarkerStyle(marker);
  pG->SetMarkerColor(col);
   pG->SetFillStyle(fill);
  pG->SetLineStyle(style);
  pG->SetFillColor(col);
 
  return;
}
void qValidate::SetHistoStyle(TH1D* &pH,int col,int width,int style,int marker=1){

 pH->SetLineColor(col);
 pH->SetLineStyle(style);
 pH->SetLineWidth(width);
 pH->SetMarkerColor(col);
 pH->SetMarkerStyle(marker);

 return;
}
void qValidate::GetYieldsEstimation(){

  std::cout<<"-## Getting estimated yieds: " << std::endl;

  for (unsigned int v=0; v<v_var.size(); v++){

	//if  (v_var[v]->GetYield()){

	  std::cout << " SS events: " << m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Integral() << std::endl;

                   for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

                     double nom=m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Integral();
                     double down=nom-m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]->Integral();
                     double up=m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first]->Integral()-nom;

		     std::cout<< "OS weighted with '"<< it->first << "' : " << nom << " +/-" << up << "/" << down << std::endl;                            

                   }//loop over dependencies
 
	//}// get yield for this distribution
 }// loop over var



  return;
}
TGraphAsymmErrors* qValidate::GetTGraphAsym(TH1D *nom,TH1D *up,TH1D *down){

  TGraphAsymmErrors *pG= new TGraphAsymmErrors();

  int nbins=nom->GetNbinsX();

  for (unsigned int i=1; i<=nbins; i++){ 

     double x=nom->GetBinCenter(i);
     double y=nom->GetBinContent(i);
     double ex=0.5*nom->GetBinWidth(i);
     double y_up=up->GetBinContent(i)-y;
     double y_down=y-down->GetBinContent(i);

     pG->SetPoint(i-1,x,y);
     pG->SetPointError(i-1,ex,ex,y_down,y_up);

  } 


 return pG;
}
void qValidate::PlotDistributions(){


   for (unsigned int v=0; v<v_var.size(); v++){

       TCanvas *pC = new TCanvas("h_"+v_var[v]->GetTitle(),"h_"+v_var[v]->GetTitle(),700,600);

       TLegend *pLeg= new TLegend(0.6,0.65,0.9,0.9); //TLegend(0.5,0.5,0.9,0.9);

       TPad* c_plot = new TPad("plot","",0.,0.3,1.0,1.0);
       c_plot->SetRightMargin(0.08);
       c_plot->SetLeftMargin(0.15);
       c_plot->SetBottomMargin(0.02);
       c_plot->SetTopMargin(0.05);
       c_plot->Draw();

       TPad* c_ratio= new TPad("ratio","",0.,0,1,0.3);
       c_ratio->SetRightMargin(0.08);
       c_ratio->SetLeftMargin(0.15);
       c_ratio->SetBottomMargin(0.35);
       c_ratio->SetTopMargin(0.03);
       c_ratio->Draw();
       c_ratio->SetGridy(kTRUE);

       c_plot->cd();
       double plot_min=0;
       double plot_max=0;

       int maxBin=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetMaximumBin();
       plot_max=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(maxBin);


       if (v_var[v]->plotLogY()){
            c_plot->SetLogy();  //forNormal
            plot_min=0.1;
            plot_max=plot_max*10000;

       }
       else plot_max=2.0*plot_max;
      
       m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetYaxis()->SetRangeUser(plot_min,plot_max);
       m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetXaxis()->SetLabelOffset(1.5);


       SetHistoStyle(m_histo_SS[v_var[v]->GetTitle()+"_SS"],1,2,1);


       m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Draw("hist");      
       pLeg->AddEntry(m_histo_SS[v_var[v]->GetTitle()+"_SS"],"SS","l");
  

       if (m_dependence.count("Standard")==1){

          SetGraphStyle(m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"],4,2,1,1,3003);
          SetHistoStyle(m_histo_OS[v_var[v]->GetTitle()+"_OS_Standard"],4,2,1);
         
          m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"]->Draw("2,p");
          m_histo_OS[v_var[v]->GetTitle()+"_OS_Standard"]->Draw("hist,same");
          //pLeg->AddEntry(m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"],"OS weighted("+latex_dependence["Standard"]+")","fl");
          pLeg->AddEntry(m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"],"OS weighted","fl");
       }

       std::vector<TH1D*> ratios;
       std::map<TString,TGraphAsymmErrors*> gRatios;      

       int counter=0;
       for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

             if (!it->first.Contains("Standard")){
                SetGraphStyle(m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first],col[counter],2,7,marker[counter],fillStyle[counter]);
	        SetHistoStyle(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first],col[counter],2,7,marker[counter]);

                m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Draw("2,p");
                m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Draw("hist,same");                
                pLeg->AddEntry(m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first],"OS weighted ("+latex_dependence[it->first]+")","lp");
             }

             TH1D *pratio=(TH1D*) m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Clone();
             pratio->Divide(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]);
             pratio->SetName("Ratio_"+v_var[v]->GetTitle()+"_OS_"+it->first);

             TGraphAsymmErrors *g_pratio=GetRatio(m_histo_SS[v_var[v]->GetTitle()+"_SS"],m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]);
             g_pratio->SetName("gRatio_"+v_var[v]->GetTitle()+"_OS_"+it->first);


             if (!it->first.Contains("Standard")) {
                    SetHistoStyle(pratio,col[counter],2,7,marker[counter]);
                    SetGraphStyle(g_pratio,col[counter],2,7,marker[counter],fillStyle[counter]);       
             }
             else{
                    SetHistoStyle(pratio,4,2,1);
                    SetGraphStyle(g_pratio,4,2,1,1,3003);
             }
             ratios.push_back(pratio);
             gRatios[it->first]=g_pratio;

             if (!it->first.Contains("Standard")) counter+=1;
        }//end of depen
        pLeg->Draw("same");

        MiniTreeAnalyzer newanalyzer;
        newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.055);


        //36184.86*(runNumber == 284500) + 43587.3*(runNumber == 300000) + 45691.0*(runNumber == 310000)
        //36207.7*(runNumber==284500)+44307.4*(runNumber==300000)+(runNumber==310000)*58450.1 //shuhui
        TString lumi="139 fb^{-1}"; //shuhui: for data
        //TString lumi="103 fb^{-1}"; //shuhui
        if (m_mcversion=="mc16a") lumi="36.1 fb^{-1}";
        else if (m_mcversion=="mc16d") lumi="43.6 fb^{-1}"; 
        else if (m_mcversion=="mc16e") lumi="46.7 fb^{-1}";
        else if (m_mcversion.Contains("201")) lumi=m_mcversion;

        newanalyzer.GetLabel(0.186,0.82,"13 TeV, "+lumi,0.045);
        TString local_tag=m_tag+" ";
        if (m_mcversion.Contains("Data")) local_tag="Data2015-2018";
        else if (m_mcversion.Contains("201")) local_tag=m_mcversion;

        //if (m_doBDT) local_tag+="w/ ECIDS";

        if (m_mcversion.Contains("mc") || m_mcversion.Contains("vjets")){
   
             TString mc_tag=m_mcversion;
             if (m_mcversion=="mc" || m_mcversion.Contains("vjets")) mc_tag="mc16a/d/e"; 
             newanalyzer.GetLabel(0.186,0.77,mc_tag+":"+local_tag,0.045); 
             if (m_doBDT) newanalyzer.GetLabel(0.186,0.72,"w/ ECIDS",0.045);

        }
        else {
          if (m_doBDT) newanalyzer.GetLabel(0.186,0.77,"w/ ECIDS",0.045);
        }
        c_ratio->cd();       

            ratios[0]->GetYaxis()->SetNdivisions(505);
            ratios[0]->GetYaxis()->SetRangeUser(0,2.001);
            ratios[0]->GetXaxis()->SetLabelSize(0.1);
            ratios[0]->GetXaxis()->SetLabelOffset(0.02);
            ratios[0]->GetXaxis()->SetTitleSize(0.1);
            ratios[0]->GetYaxis()->SetLabelSize(0.1);
            ratios[0]->GetYaxis()->SetTitleSize(0.1);
            ratios[0]->GetYaxis()->SetTitleOffset(0.5);
            ratios[0]->GetYaxis()->SetTitle("SS/OS");
            ratios[0]->GetYaxis()->CenterTitle();
            ratios[0]->SetLineColor(kWhite);
            ratios[0]->Draw("hist");
        int kk=0;

        if (gRatios.count("Standard")==1){

           // gRatios["Standard"]->GetYaxis()->SetNdivisions(505);
           // gRatios["Standard"]->GetYaxis()->SetRangeUser(0,2.001);
           // gRatios["Standard"]->GetXaxis()->SetLabelSize(0.1);
           // gRatios["Standard"]->GetXaxis()->SetTitleSize(0.1);
           // gRatios["Standard"]->GetYaxis()->SetLabelSize(0.1);
           // gRatios["Standard"]->GetYaxis()->SetTitleSize(0.1);
           // gRatios["Standard"]->GetYaxis()->SetTitleOffset(0.5);
           // gRatios["Standard"]->GetYaxis()->SetTitle("SS/OS");
           // gRatios["Standard"]->GetYaxis()->CenterTitle();

            //gRatios["Standard"]->Draw("Ap same");
            gRatios["Standard"]->Draw("p same");

            kk+=1;
        }


      for (std::map<TString,TGraphAsymmErrors*>::iterator it=gRatios.begin(); it!=gRatios.end(); ++it){

            gRatios[it->first]->GetYaxis()->SetNdivisions(505); //404
            gRatios[it->first]->GetYaxis()->SetRangeUser(0,2.0);
            gRatios[it->first]->GetXaxis()->SetLabelSize(0.1);
            gRatios[it->first]->GetXaxis()->SetTitleSize(0.1);
            gRatios[it->first]->GetYaxis()->SetLabelSize(0.1);
            gRatios[it->first]->GetYaxis()->SetTitleSize(0.1);
            gRatios[it->first]->GetYaxis()->SetTitleOffset(0.5);
            gRatios[it->first]->GetYaxis()->SetTitle("SS/OS");       
            gRatios[it->first]->GetYaxis()->CenterTitle();

            if (!it->first.Contains("Standard")){

                  //if (kk==0) gRatios[it->first]->Draw("Ap");
                  //else  gRatios[it->first]->Draw("p");
                   gRatios[it->first]->Draw("p same");
                  kk+=1;
            }

       }

       /*
        for (unsigned int k=0; k<ratios.size(); k++){

          ratios[k]->GetYaxis()->SetRangeUser(0.5,3);
       
          ratios[k]->GetXaxis()->SetLabelSize(0.1);
          ratios[k]->GetXaxis()->SetTitleSize(0.1);
          ratios[k]->GetYaxis()->SetLabelSize(0.1);
          ratios[k]->GetYaxis()->SetTitleSize(0.1);
          ratios[k]->GetYaxis()->SetTitleOffset(0.5);
          ratios[k]->GetYaxis()->SetTitle("Ratio");


          if (k==0) ratios[k]->Draw("l");
          else ratios[k]->Draw("l,same");

          gRatios[k]->Draw("2,p");

        }
        */

        TLine* line1 = new TLine(ratios[0]->GetXaxis()->GetXmin(),1.,ratios[0]->GetXaxis()->GetXmax(),1.);
        line1->SetLineWidth(2);
        line1->Draw("same");

        system("mkdir -p "+m_outputDir);
        TString local_output=m_outputDir+"/h_"+m_process+"_"+v_var[v]->GetTitle();


        if (m_mcversion.Contains("mc16") && !m_process.Contains("mc") && !m_process.Contains("Data")) {
            system("mkdir -p "+m_outputDir+"/"+m_mcversion);
            local_output=m_outputDir+"/"+m_mcversion+"/h_"+m_process+"_"+v_var[v]->GetTitle();
        }

        if (m_doBDT) local_output+="wBDT";

        pC->SaveAs(local_output+".pdf");
        pC->SaveAs(local_output+".root");

   }//end of var

 return;
}
TGraphAsymmErrors* qValidate::GetRatio(TH1D *pNum,TGraphAsymmErrors *pDen){

 TGraphAsymmErrors *pG=new TGraphAsymmErrors();

 if (pDen->GetN()!=pNum->GetNbinsX()){
 
    std::cout << "--> Will try to divide different things for the ratio plot.. aborting!" << std::endl;
    exit(1);
 }

 Double_t* x_d=pDen->GetX();
 Double_t* y_d=pDen->GetY();
 Double_t* y_u=pDen->GetEYhigh();
 Double_t* y_l=pDen->GetEYlow();
 Double_t* x_e=pDen->GetEXlow();

 for (unsigned int i=1; i<=pNum->GetNbinsX(); i++){
  
    //std::cout << "DEBUG: " << pNum->GetBinCenter(i) << ", " << x_d[i-1] << std::endl;
    //std::cout << "---shuhui: pNum->GetBinCenter(i)==x_d[i-1]: " << (pNum->GetBinCenter(i)==x_d[i-1]) << std::endl;
    //if (pNum->GetBinCenter(i)==x_d[i-1]){
    if (1){//shuhui

         //std::cout << "DEBUG: " << pNum->GetBinCenter(i) << ", " << x_d[i-1] << std::endl;

         double num_error=pNum->GetBinError(i);
         double den_error_up=y_u[i-1];
         double den_error_down=y_l[i-1];
         double num=pNum->GetBinContent(i);
         double den=y_d[i-1];         

         double error_up=GetError(num,num_error,den,den_error_up);
         double error_down=GetError(num,num_error,den,den_error_down);
         pG->SetPoint(i-1,x_d[i-1],num/den);
         pG->SetPointError(i-1,x_e[i-1],x_e[i-1],error_down,error_up);

    }
    else {
    
       std::cout <<"--> Something wrong happended... aborting!" << std::endl;
       exit(1); 
    }
 }

 pG->GetXaxis()->SetTitle(pNum->GetXaxis()->GetTitle());
 pG->GetXaxis()->SetRangeUser(pNum->GetBinLowEdge(1),pNum->GetBinLowEdge(pNum->GetNbinsX()+1));


 return pG;
}
double qValidate::GetError(double num,double num_error,double den,double den_error){

  double term1=pow(num_error/num,2);
  double term2=pow(den_error/den,2);

  double val=(num/den)*TMath::Sqrt(term1+term2);
 
  return val;
}
TH2F* qValidate::ChooseHisto(std::vector<qHisto*> histos,double var){

  TH2F* pH=0;  

 // std::cout << "var: " << var << std::endl;

  for (unsigned int i=0; i<histos.size(); i++){

     double min=histos[i]->GetMin();
     double max=histos[i]->GetMax();

     if (min >=0 && max >=0){
        if (var >= min && var <= max) pH=histos[i]->GetHisto();
     }
     else if (min>=0 && max <0){
        if (var >= min) pH=histos[i]->GetHisto();
     }
     else if (min<0 && max >=0){
        if (var <= max) pH=histos[i]->GetHisto();
     }
     else std::cout << "-- Undefined values for min and max  : " << min << ", " << max << std::endl;

  }

  //std::cout << "-- min: " << min << ", max=" << max << std::endl; 

  //std::cout << "-- Will use histo : " << pH->GetName() << std::endl;

  return pH;
}
std::vector<double> qValidate::GetEpsilon(TH2F *pH,double el_pt,double el_eta){

  int bin=pH->FindBin(el_eta,el_pt);
  int binx=pH->GetXaxis()->FindBin(el_eta);
  int biny=pH->GetYaxis()->FindBin(el_pt);

  double epsilon=pH->GetBinContent(binx,biny);

  double up=epsilon+pH->GetBinError(binx,biny);
  double down=epsilon-pH->GetBinError(binx,biny);

  if (down<0) down=0;
  if (up>1) up=1;

  std::vector<double> eps={epsilon,up,down};

  return eps;
}

qWeight* qValidate::ComputeWeight(double el_pt1,double el_pt2,double el_eta1,double el_eta2,std::vector<qHisto*> local_rates,double var=-1){

   double weight=-1;

 
   TH2F *pH=0;

   if (local_rates.size()==1) pH=local_rates[0]->GetHisto();
   else pH=ChooseHisto(local_rates,var);
  
   std::vector<double> epsilon_1=GetEpsilon(pH,el_pt1,el_eta1);
   std::vector<double> epsilon_2=GetEpsilon(pH,el_pt2,el_eta2);

   weight=GetWeight(epsilon_1[0],epsilon_2[0]);
   //std::cout << "--- Shuhui DEBUG: epsilon1: "<< epsilon_1[0] << ", epsilon2: " << epsilon_2[0] << std::endl;  
   double up=GetWeight(epsilon_1[1],epsilon_2[1]);
   double down=GetWeight(epsilon_1[2],epsilon_2[2]);
      
   qWeight *total=new qWeight(weight,up,down);

   return total;
}

qWeight* qValidate::GetRateReweighting(double el_pt1,double el_pt2,double el_eta1,double el_eta2,TString dependence,double var=-1){
   
   double nom=0;
   double up2=0;
   double down2=0;

   for (std::map<TString,std::vector<qHisto*> >::iterator it=m_method_A.begin(); it!=m_method_A.end(); ++it){

         qWeight *qW=ComputeWeight(el_pt1,el_pt2,el_eta1,el_eta2,it->second,var);         

         double f=m_bkg_comp[it->first]; 
         double local_nom=qW->GetNominal();
         double sigma_up=qW->GetUp()-local_nom;
         double sigma_down=local_nom-qW->GetDown();

         nom=nom+f*qW->GetNominal();
         up2=up2+pow(f*sigma_up,2);
         down2=down2+pow(sigma_down,2);
   }
   
  qWeight *weight=new qWeight(nom,nom+TMath::Sqrt(up2),nom-TMath::Sqrt(down2));

  return weight;
}
double qValidate::GetWeight(double e1,double e2){

  double num=e1+e2-2*e1*e2;
  double den=1-e1-e2+2*e1*e2;

  double weight=num/den;

  return weight;
}

void qValidate::InitializeHistos(){

  std::cout << "-- Initliazing histos..." << std::endl;

  MiniTreeAnalyzer analyzer;  
  analyzer.AddWeight("1"); //dummy just to avoid that the y-label show "Unweighted events"
 
  for (unsigned int i=0; i< v_var.size(); i++){

       TH1D *pH_SS=analyzer.CreateHisto(v_var[i],"SS");

       m_histo_SS[v_var[i]->GetTitle()+"_SS"]=(TH1D*)pH_SS;

            for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

                   TH1D *pH_OS=analyzer.CreateHisto(v_var[i],"OS_"+it->first);
                   m_histo_OS[v_var[i]->GetTitle()+"_OS_"+it->first]=(TH1D*)pH_OS;

                   TH1D *pH_OSUp=analyzer.CreateHisto(v_var[i],"OSUp_"+it->first);
                   m_histo_OSUp[v_var[i]->GetTitle()+"_OSUp_"+it->first]=(TH1D*)pH_OSUp;

                   TH1D *pH_OSDown=analyzer.CreateHisto(v_var[i],"OSDown_"+it->first);
                   m_histo_OSDown[v_var[i]->GetTitle()+"_OSDown_"+it->first]=(TH1D*)pH_OSDown;


            }//end depen
  }//end var


 return;
}
std::vector<TString> qValidate::FillVector(TString file){

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
void qValidate::AddDependence(TString name,TString latex,std::vector<TString> rates){

  std::vector<qHisto*> histos2D=GetHistos2D(rates);

  m_dependence[name]=histos2D;

  latex_dependence[name]=latex;

 return;
}
void qValidate::AddDependence(TString name,TString latex,std::map<TString,std::vector<TString> > rates){

  std::vector<qHisto*> histos2D=GetHistos2D(rates.begin()->second);

  for (std::map<TString,std::vector<TString> >::iterator it=rates.begin(); it!=rates.end(); ++it){ 
    
     m_method_A[it->first]=GetHistos2D(it->second);
  }

  m_dependence[name]=histos2D; //dummy: it will not be used->it has to be filled by consistency

  latex_dependence[name]=latex;

  return;
}
std::vector<qHisto*> qValidate::GetHistos2D(std::vector<TString> rates){

    std::vector<qHisto*> local;

    for (unsigned int i=0; i<rates.size(); i++){

      TString file=File_wRates+"/"+rates[i]+".root";
      TString histoname=rates[i]+"_all_num";

      std::cout << "-- Reading file with rates: " << file << std::endl;

      TFile *pFile = new TFile(file);

      TH2F *pH=(TH2F*)pFile->Get(histoname);
      pH->SetDirectory(0);

      qHisto *local_histo=new qHisto(pH);

      local.push_back(local_histo);   

      pFile->Close();
   }


   if (local.size()==0){

        std::cout << "-- Could not read file with rates from : " << File_wRates  << std::endl;
        exit(1);
    }

   return local;
}
qHisto::qHisto(TH2F* pH):
 m_pH(pH),
 min(-1),
 max(-1)
{

  SetValues();   

}
qHisto::~qHisto()
{}
void qHisto::SetValues(){

    
 std::string name(m_pH->GetName());
 std::cout << "---- in internal class : " << name << std::endl;  
 TString local_name = (TString)name;
 MiniTreeAnalyzer analyzer;
 std::vector<std::string> tokens;

 analyzer.tokenizeString(name,'_',tokens);  


 TString str_min="";
 TString str_max="";

 for (unsigned int i=0; i<tokens.size(); i++){

     TString l_s=tokens[i];

      TString str( l_s(0,2) );
      if (str=="ge") str_min=tokens[i];
      else if (str=="le") str_max=tokens[i];
    /*  TString str=tokens[i];
      if (str.Contains("ge")) str_min=str;
      else if (str.Contains("le")) str_max=str;      
     */

 }

 if (str_min!="") min=GetValue(str_min,"ge");
 if (str_max!="") max=GetValue(str_max,"le");
 
 if (local_name.Contains("CR0")&&!local_name.Contains("CR01")) {SetMin(-1);SetMax(0);}
 if (local_name.Contains("CR01")) {SetMin(-1);SetMax(0.05);}
 if (local_name.Contains("CR1")&&!local_name.Contains("CR12")) {SetMin(0);SetMax(0.05);}
 if (local_name.Contains("CR12")) {SetMin(0);SetMax(0.1);}
 if (local_name.Contains("CR2")) {SetMin(0.05);SetMax(0.1);}
 if (local_name.Contains("CR3")) {SetMin(0.1);SetMax(0.2);}
 if (local_name.Contains("CR4")) {SetMin(0.2);SetMax(0.3);}
 if (local_name.Contains("CR5")) {SetMin(0.3);SetMax(-1);}


 std::cout << "- min: " << min << ", max: " << max << std::endl;

 return;
}
double qHisto::GetValue(TString str,TString out){

   double val=-1;

   str.ReplaceAll(out,"");

   val=str.Atof();

   return val; 
}
qWeight::qWeight(double w,double up=0,double d=0):
 w_nom(w),
 w_up(up),
 w_down(d)
{}
qWeight::~qWeight()
{}

