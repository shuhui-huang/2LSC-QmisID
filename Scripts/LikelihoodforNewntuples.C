#include <TMinuit.h>
#include <TH2F.h>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

extern void NLL(int &npar, double *gin,double &f, double *par, Int_t iflag);


int NBINS_ETA=-1;
int NBINS_PT=-1;


std::vector< std::vector< std::vector< std::vector<float> > > > Array_SS;
std::vector< std::vector< std::vector< std::vector<float> > > > Array_OS;

class Like{


   public:

       Like(TString name,TString input,bool isData,bool doBDT);
       ~Like();

       void Execute();         
       inline void SetPath(TString path){m_path=path;};   
       inline void SetTreeName(TString tree){m_tree=tree;};
       inline void SubtractBkg(bool subtract){m_subtract_bkg=subtract;};
       inline void SetZwindow(float low,float low_Z,float up_Z,float up){m_low=low; m_up=up; m_low_Z=low_Z; m_up_Z=up_Z;};
       inline void SetEtaBinning(Float_t *binning,int n){m_eta_binning=binning; m_bins_eta=n;};
       inline void SetPtBinning(Float_t *binning,int n){m_pt_binning=binning; m_bins_pt=n;};
              void SetHtCut(TString cut,float val1,float val2=-1);
       inline void SetNjetsCut(int val){m_njets_cut=val;};
       inline TH2F *GetRates(){return m_TH2_Rates;};
       inline void SetReferenceHisto(TString path,std::vector<TString> hist){m_path_ref=path; m_ref_hist=hist;};
       inline void SetOutputDir(TString outputDir,TString period){m_outputDir=outputDir; m_period=period;};

   private:

       void CreateHistos();
       void InitializeVectors();
       void ReadInputFile(TString input);
       int ComputeRates(std::vector<double> reference);
       std::vector<TString> FillVector(TString file);    
       void FillVectorsEtaPt(float eta1,float eta2,float pt1,float pt2,float charge1,float charge2,float Mll,float weight);
       void PrintSetup();
       bool isZpeak(float M);
       TH2F* CreateEtaPtHisto(TString label);
       void PlotRates();
       void DebugVector(); // check that it is using the number of events that it has to!
       std::vector< std::vector<double> > GetReferenceHistos();
       std::vector<double> GetReferenceValues(TH2F *pH);
       void WriteRates();


       Float_t *m_eta_binning;
       Float_t *m_pt_binning;

       TString m_outputDir;
       TString m_period;

       int m_bins_eta;
       int m_bins_pt;
       TString m_file;
       bool m_isData;
       TString m_name;
       TString m_path;
       TString m_tree;
       bool m_doBDT;
       bool m_subtract_bkg;

       float m_low;
       float m_up;
       float m_low_Z;
       float m_up_Z;
       float m_low_Ht;
       float m_up_Ht;
       int m_njets_cut;
       TH2F *m_TH2_Rates;
       TH1F *m_Mll_SS;
       TH1F *m_Mll_OS;
       TString m_path_ref;
       std::vector<TString> m_ref_hist;
        

};

Like *ComputeLikelihood(TString name,TString file,bool isData,bool doBDT,bool doSubtractBkg,TString outputDir,TString period="",TString pathReference="",std::vector<TString> ref_hist={});

void Likelihood(TString samples){

    std::vector<TString> ref_hist={"Zjets_Default_nJets_le1.root"};
    TString refDir="../run_31March/Truth_MC_Unweighted_Full/";
    TString InputPath = "../Inputlists/";


    //std::vector<TString> data={"2015","2016","2017","2018"};
    std::vector<TString> data={""};
    //std::vector<TString> mcversion={"mc16a","mc16d","mc16e"};
   //  std::vector<TString> mcversion={"mc20a","mc20d","mc20e"};
    std::vector<TString> mcversion={"mc"};

    TString OutputDirMC="Likelihood_MC_test";
    TString OutputDirData="Likelihood_Data_test";

    std::vector<TString> Regions={"_nJets_ge0", "_nJets_le1"};


   if (samples.Contains("Data")){
      for (unsigned int i=0; i<data.size(); i++){
         for (unsigned int j=0; j<Regions.size(); j++){
            std::cout << "-- Region: " << Regions[j] << std::endl;
            Like *pLikeData=ComputeLikelihood("Data"+data[i]+"_Default"+Regions[j],InputPath+"Data.list",true,false,true,OutputDirData,"Data"+data[i],refDir,ref_hist);
         
         }
      }
   }
   else{
      for (unsigned int i=0; i<mcversion.size(); i++){
         TString prefix=mcversion[i];
         if (prefix=="") prefix="mc";
         for (unsigned int j=0; j<Regions.size(); j++){
          Like *pLikeMC=ComputeLikelihood(prefix+"_Default"+Regions[j],InputPath+"ZjetsMC.list",false,false,true,OutputDirMC,mcversion[i],refDir,ref_hist);
          }
       }
   }




 return;
}
Like *ComputeLikelihood(TString name,TString file,bool isData,bool doBDT,bool doSubtractBkg,TString outputDir,TString period="",TString pathReference="",std::vector<TString> ref_hist={}){

  //Float_t EtaBinning[]={0,0.6,1.1,1.52,1.7,2.3,2.5};
  //Float_t PtBinning[]={0,60,90,130,200,2500};
  //Float_t PtBinning[]={0,60,90,130,2500};

//   double PtBinning[]   = {20, 50, 100, 200, 1000};
//   Float_t EtaBinning[]  = {0, 1.37, 1.52, 1.7, 2.5};
  Float_t PtBinning[]   = {20, 50, 100, 200, 2600};
  Float_t EtaBinning[]  = {0, 1.37, 1.52, 1.7, 2.6};

  NBINS_ETA=sizeof(EtaBinning)/sizeof(EtaBinning[0])-1;
  NBINS_PT=sizeof(PtBinning)/sizeof(PtBinning[0])-1;

  TString path="/eos/home-s/shudong/workspace/HHML/Run3/HHMLAnalysisCode/output/ntuple_QMisID/hhml_v2/";


  if (doBDT) outputDir+="_wBDT";

  TString tree="start";
  if (name.Contains("_nJets_le1")) tree="qmisid_cr";

  Like *pLike= new Like(name,file,isData,doBDT); //bool isData, bool doBDT
  pLike->SetPath(path);
  pLike->SetTreeName(tree);

  pLike->SetEtaBinning(EtaBinning,NBINS_ETA);
  pLike->SetPtBinning(PtBinning,NBINS_PT);
  pLike->SetZwindow(71,81,101,111);
  pLike->SubtractBkg(doSubtractBkg);
  pLike->SetReferenceHisto(pathReference,ref_hist);
  pLike->SetOutputDir(outputDir,period);

  pLike->SetNjetsCut(0);
  //pLike->SetHtCut("range",60,100); //low,up, range

  pLike->Execute();
    
  
 return pLike;
}

Like::Like(TString name,TString input,bool isData,bool doBDT):
 m_name(name),
 m_outputDir(""),
 m_period(""),
 m_doBDT(doBDT),
 m_subtract_bkg(0),
 m_file(input),
 m_isData(isData),
 m_path(""),
 m_tree(""),
 m_path_ref(""),
 m_ref_hist(0),
 m_bins_eta(-1),
 m_bins_pt(-1),
 m_low(-1),
 m_up(-1),
 m_low_Ht(-1),
 m_up_Ht(-1),
 m_low_Z(-1),
 m_up_Z(1),
 m_TH2_Rates(0),
 m_Mll_SS(0),
 m_Mll_OS(0),
 m_njets_cut(-1)
{}

Like::~Like()
{

Array_SS.clear();
Array_OS.clear();

Array_SS.shrink_to_fit();
Array_OS.shrink_to_fit();

}

void Like::Execute(){

    PrintSetup();

    std::vector<TString> samples=FillVector(m_file);

    CreateHistos();
    InitializeVectors();

    std::cout << "Check vector normalization (should be zero):" << std::endl;
    DebugVector();

    for (unsigned int i=0; i<samples.size(); i++) ReadInputFile(samples[i]);

    std::cout << "Check vector normalization (should be different from zero):" << std::endl;
    DebugVector();
    

    int out=-1;


    std::vector< std::vector<double> > reference=GetReferenceHistos();


    
    //for (unsigned int i=0; i<reference.size(); i++){
    for (unsigned int i=1; i<reference.size(); i++){ // shuhui: debug
       std::cout << "--shuhui-- After GetReferenceHistos(), reference.size = " << reference.size() << std::endl; // shuhui: debug
       std::cout << "-- In loop: trying reference histo :" << i << std::endl;
     //      for (unsigned int j=0; j<reference[i].size(); j++){
    //           std::cout << "--shuhui-- reference[" << i << "][" << i << "] :" << &(reference[i][j]) << std::endl; // shuhui: debug
       //    }
       out=ComputeRates(reference[i]);
   
       if (out==3){
        
          std::cout << "-- Rates have been succesfully computed..." << std::endl;
          //break;
       }
    }


    if (out!=3) {
         std::cout << "-- Sorry... rates could not be computed: try with a different input reference rates..." << std::endl;
         //exit(1);
    }
    else {
       PlotRates();
       WriteRates();
    }
   //PlotRates();
  


  return;
}
std::vector< std::vector<double> > Like::GetReferenceHistos(){

  std::vector< std::vector<double> >  vec;

  std::vector<double> defa;

  vec.push_back(defa); 

  //std::vector<TString> period={"/"};///,"/mc16a/","/mc16d/","/mc16e/"};
  //std::vector<TString> period={m_period+"/"};
  std::vector<TString> period={""};


  if (m_ref_hist.size()>0){

        for (unsigned int i=0; i<period.size(); i++){

              for (unsigned int ref=0; ref<m_ref_hist.size();  ref++){
      
                     TString inputfile=m_path_ref+period[i]+m_ref_hist[ref];
                    

                     std::cout << "-- Will read file reference :" << inputfile << std::endl;

                     TFile *pInput=new TFile(inputfile);

                     TString name=m_ref_hist[ref];
 
                     name.ReplaceAll(".root","_all_num");
                     std::cout << "--shuhui-- GetHisto with name :" << name << std::endl; //shuhui: debug
                    
                     TH2F *pH=(TH2F*)pInput->Get(name);

                     std::vector<double> local=GetReferenceValues(pH);

                     vec.push_back(local);

              } 
        }
  }
  
  return vec;
}

std::vector<double> Like::GetReferenceValues(TH2F *pH){

  std::vector<double> vec;

      for (int j=1; j<=m_bins_eta; j++){

  for (int i=1; i<=m_bins_pt; i++){

           double val=pH->GetBinContent(j,i);
            std::cout << "--shuhui-- GetHisto Content :" << val << std::endl;   //shuhui: debug
           vec.push_back(val);

      }
  }


  return vec;
}
void Like::SetHtCut(TString cut,float val1,float val2=-1){
  

  if (cut=="low") m_low_Ht=val1;
  else if (cut=="range"){
     m_low_Ht=val1;
     m_up_Ht=val2;
  }
  else if (cut=="up") m_up_Ht=val1;



 return;
}
void Like::ReadInputFile(TString input){

  TString file=m_path+input;

  TFile *pInput = new TFile(file);
  if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
  }
  else {
       cout << "-- Opening file :" << file <<  endl;
  }


 TTree* pTree = (TTree*)pInput->Get(m_tree);
 int nEvents=pTree->GetEntries();
//  int nEvents=200000; // shuhui: for debug usage
 std::cout<< "-- Tree entries:" << nEvents << std::endl;

 float el0_pt=-1;
 float el1_pt=-1;
 float el0_eta=-1;
 float el1_eta=-1;
 int isSS=0;
 int isOS=0;
 int isel0_Tight=0;
 int isel1_Tight=0;

//  std::vector<Float_t> *el_pt=new std::vector<Float_t>();
//  std::vector<Float_t> *el_eta=new std::vector<Float_t>();
//  std::vector<Int_t> *el_charge=new std::vector<Int_t>();
//  std::vector<Char_t> *el_ECIDS=new std::vector<Char_t>();
//  std::vector<Float_t> *el_phi=new std::vector<Float_t>();
//  std::vector<Float_t> *el_e=new std::vector<Float_t>();


 double toGeV=0.001;

//  unsigned int runNumber=-1;
//  int nJets=-1;
//  float HT_all=-1;
 double Mll01=-1;

 double totweight=1;



 //getting branches
 pTree->SetBranchAddress("l1_pt",&el0_pt);
 pTree->SetBranchAddress("l2_pt",&el1_pt);
 pTree->SetBranchAddress("l1_eta",&el0_eta);
 pTree->SetBranchAddress("l2_eta",&el1_eta);
 pTree->SetBranchAddress("l1_IDTight",&isel0_Tight);
 pTree->SetBranchAddress("l2_IDTight",&isel1_Tight);
 pTree->SetBranchAddress("same_charge",&isSS);
 pTree->SetBranchAddress("opposite_charge",&isOS);
 pTree->SetBranchAddress("m_l0l1",&Mll01);

  //getting weights

 if (!m_isData){
    pTree->SetBranchAddress("weight_NOSYS",&totweight);
 }


  double out=0.;

  double extra=0.;

  for (unsigned int i=0; i<nEvents; i++){

        pTree->GetEntry(i);

        if (i%100000 ==0) std::cout << "-- event: " << i << std::endl;

        float el_pt1=el0_pt*toGeV;
        float el_pt2=el1_pt*toGeV;
        float el_eta1=fabs(el0_eta);
        float el_eta2=fabs(el1_eta);
        int el_charge1=0;
        int el_charge2=0;

        float Mll=Mll01*toGeV;
        float event_weight=1.0;

         // shuhui: define charge manually as I don't want to change the structure of  FillVectorsEtaPt() function, where making use of charg1*charge2 to select SS or OS events.
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

         //shuhui: select two tight electrons
         if (!isel0_Tight || !isel1_Tight) continue; 
         // if (!(el0_eta<2.5) || !(el1_eta<2.5)) continue; 
        
        if (!m_isData){
               //event_weight=weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise*(36184.86*(runNumber == 284500) + 43587.3*(runNumber == 300000) + 45691.0*(runNumber == 310000));
               //event_weight=(36207.7*(runNumber==284500)+44307.4*(runNumber==300000)+(runNumber==310000)*58450.1)*(1/138965.16)*weight_normalise*weight_pileup*weight_jvt*weight_mc*weight_leptonSF*weight_bTagSF_MV2c10_Continuous;
               event_weight=totweight;
                      
               // if (m_doBDT) event_weight*=weight_indiv_SF_EL_ChargeMisID*weight_indiv_SF_EL_ChargeID;
        }
        //if (event_weight!=0){
        //std::cout << "**DEBUG**: 1.runNumber = " << runNumber << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_normalise = " << weight_normalise << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_pileup = " << weight_pileup << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_jvt = " << weight_jvt << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_mc = " << weight_mc << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_leptonSF = " << weight_leptonSF << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.weight_bTagSF_MV2c10_Continuous = " << weight_bTagSF_MV2c10_Continuous << std::endl; // shuhui: debug 
        //std::cout << "**DEBUG**: 1.event_weight = " << event_weight << std::endl; // shuhui: debug 
        //}

        if (m_doBDT){
            //  if (!el_BDT1 || !el_BDT2) continue;
        }

  
        if (Mll < m_low || Mll > m_up) {
               out+=event_weight;
               continue; //use only events between 71 and 111 GeV 
   
       }


        
 
      //   if (m_njets_cut !=-1 && jets<m_njets_cut) continue; //pass if jets>m_njets_cut   

        FillVectorsEtaPt(el_eta1,el_eta2,el_pt1,el_pt2,el_charge1,el_charge2,Mll,event_weight);

        extra+=event_weight;
        //if (runNumber!=-1) std::cout << "**DEBUG**: 2.event_weight = " << event_weight << std::endl; // shuhui: debug

        if (el_charge1*el_charge2>0) m_Mll_SS->Fill(Mll,event_weight);
        else if (el_charge1*el_charge2<0) m_Mll_OS->Fill(Mll,event_weight);




 }//end for events


 std::cout << "***** Total events used from this sample : " << extra << std::endl;
 std::cout << "-- outside the interesting range: " << out << std::endl;

 pInput->Close();

 return;
}
void Like::FillVectorsEtaPt(float eta1,float eta2,float pt1,float pt2,float charge1,float charge2,float Mll,float weight){

//std::cout <<"set::  eta1:" << eta1 << ", eta2:" << eta2 << ", pt1: " << pt1 << ", pt2: " << pt2 << std::endl;  


 int m_eta1=-1,m_eta2=-1,m_pt1=-1,m_pt2=-1;

 for (int i=0;i<m_bins_eta; i++)
 {
   if ((eta1>m_eta_binning[i]) && (eta1<=m_eta_binning[i+1]))  m_eta1=i;
   if ((eta2>m_eta_binning[i]) && (eta2<=m_eta_binning[i+1]))  m_eta2=i;
 }

for (int j=0;j<m_bins_pt; j++)
 {
   if ((pt1>m_pt_binning[j]) && (pt1<=m_pt_binning[j+1]))  m_pt1=j;
   if ((pt2>m_pt_binning[j]) && (pt2<=m_pt_binning[j+1]))  m_pt2=j;
 }




if ((m_eta1==-1) || (m_eta2==-1) || (m_pt1==-1) || (m_pt2==-1)){
    std::cout << "Something is wrong with this event: Check eta and pT" << std::endl;
    std::cout << "eta1, eta2, pt1, pt2 :" << eta1 << ", " << eta2  << ", " << pt1 << ", " << pt2 << std::endl;
    if (m_pt1==-1) std::cout << "assigning : pt1 :" << m_bins_pt-1  << std::endl;
    if (m_pt2==-1) std::cout << "assigning : pt2 :" << m_bins_pt-1  << std::endl;
}


if (m_pt1==-1 && pt1>m_pt_binning[m_bins_pt]) m_pt1=m_bins_pt-1;
if (m_pt2==-1 && pt2>m_pt_binning[m_bins_pt]) m_pt2=m_bins_pt-1;


int lowE = m_eta1;
int upE = m_eta2;
int lowPt = m_pt1;
int upPt = m_pt2;

if (lowE>upE){
  lowE=m_eta2;
  upE=m_eta1;
  lowPt=m_pt2;
  upPt=m_pt1;
}


bool isZ=isZpeak(Mll);
int qq=charge1*charge2;

if (m_subtract_bkg && !isZ) weight=-0.5*weight;


if (qq>0){

    Array_SS[lowE][upE][lowPt][upPt]=Array_SS[lowE][upE][lowPt][upPt]+weight;

}
else if (qq<0){

   Array_OS[lowE][upE][lowPt][upPt]=Array_OS[lowE][upE][lowPt][upPt]+weight;
}
else{
 std::cout << "-- the charge for one of these leptons has not been defined : " << charge1 << ", " << charge2 << std::endl;

 exit(0);

}

/*std::cout << " bin eta 1:" << m_eta1 << std::endl;
std::cout << " bin eta 2:" << m_eta2 << std::endl;
std::cout << " bin pt 1:" << m_pt1 << std::endl;
std::cout << " bin pt 2:" << m_pt2 << std::endl;
*/


 return;
}
bool Like::isZpeak(float M){

 if (M > m_low_Z && M < m_up_Z) return true;
 else return false;

}
int Like::ComputeRates(std::vector<double> reference){

   const int NPAR=m_bins_eta*m_bins_pt;


   //init minuit
   TMinuit *pMinuit = new TMinuit(NPAR);
 
   //give function NLL
   pMinuit->SetFCN(NLL); 
 

   //defining and setting parameters:
   Double_t vstart[NPAR];
   Double_t step[NPAR];
   Int_t ierflg = 1;
   double arglist[10];


   float init=0.00005;
  
   std::cout << "--shuhui-- In ComputeRates() reference.size = " << reference.size() << std::endl; // shuhui: debug

  if (reference.size()>0){
  
     if (reference.size()==NPAR){

        for (unsigned int i=0; i<reference.size(); i++){ vstart[i]=reference[i];
            std::cout << "--shuhui-- reference[" << i << "] = " << vstart[i] << std::endl; }// shuhui: debug

     }
     else {

        std::cout << "-- Something went wrong when reading the reference rates... ABORTING!" << std::endl;
        std::cout << "-- n reference rates: " << reference.size() << std::endl;
     }
    /*
     vstart[0]=vstart[0]
     
     for (int j=1; j<NPAR;++j){
        vstart[j]=vstart[j-1]*2;
     }

     */
  }
  else{
   for (unsigned int pt=0; pt<m_bins_pt; pt++){


      float local_start=(pt+1)*init;
      if (pt==2) {
         local_start=0.0005;
         
      }  
      if (pt>=3) {
         local_start=0.001;

      }

      int level=pt*m_bins_eta;

      vstart[level]=local_start;

      for (unsigned int eta=1; eta<m_bins_eta; eta++){

           int bin=level + eta;
           std::cout << "bin : " << bin << std::endl;
           vstart[bin]=vstart[bin-1]*2;


      }
   }
  }

   //float init=0.0005;
   //if (m_doBDT) init=0.1*init;

   for (unsigned int i=0; i<NPAR; i++){

       step[i]=1e-12;  

       if (m_doBDT) step[i]=0.5*step[i];
 
       const TString name=TString::Format("epsilon_%d",i+1);              

       float in=0.0000000000005;
       float end=0.99;
       if (m_name.Contains("Data") && NPAR==24){
            if (name.Contains("epsilon_9")){
                    in=0.0001; 
                    vstart[i]=0.0006;
                    end=0.001;

            }
            if (m_name.Contains("2015")){
                   step[i]=1e-12;
                   //in=0.0000000000005;

                   if (name.Contains("epsilon_2")) {
                        in=0.00001;
                        vstart[i]=0.000045;

                        end=0.1;
                   }
                   if (name.Contains("epsilon_7")) {
                        in=0.0003; //5*x10^12
                        vstart[i]=2*in;
                        end=0.0005;
                   }


            }
      }
   /*   else if (m_name.Contains("mc") && NPAR==24){
            if (name.Contains("epsilon_9")){
                    in=0.0001; 
                    vstart[i]=0.0006;
                    end=0.001;

            }
            if (m_name.Contains("16")){
                   step[i]=1e-12;
                   //in=0.0000000000005;

                   if (name.Contains("epsilon_2")) {
                        in=0.00001;
                        vstart[i]=0.000045;

                        end=0.1;
                   }
                   if (name.Contains("epsilon_7")) {
                        in=0.0003; //5*x10^12
                        vstart[i]=2*in;
                        end=0.0005;
                   }


            }
      }
    */
       pMinuit->mnparm(i,name,vstart[i],step[i],in,end,ierflg);

    }

   //TMinuit::mnexcm(const char *command, Double_t *plist, Int_t llist, Int_t &ierflg)

   //set output log
   pMinuit->SetPrintLevel(1); //-1 no output, 1 standard output

   //set error definition
   //arglist[0]=1;  //1 for Chi square, 0.5 for negative log likelihood
   //pMinuit->mnexcm("SET ERRDEF",arglist,1,ierrflg);
   pMinuit->SetErrorDef(1); //1 for Chi square, 0.5 for negative log likelihood


   //setting strategy: 1 standard, 2: try to improve the minimum (will take more time)
   arglist[0]=1;
   pMinuit->mnexcm("SET STR",arglist,1,ierflg);


   //Minimization itself:
   arglist[0]=10000000;
   pMinuit->mnexcm("MIGRAD",arglist,1,ierflg);


   //print statistic
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   pMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

   if(icstat < 3) cout<<"There is  a problem with Minuit! "<<endl;

   std::cout << std::endl;
   std::cout << " Minimum fonction = " << amin <<std:: endl;
   std::cout << " Estimated vert. distance to min. = " << edm <<std:: endl;
   std::cout << " Number of variable parameters = " << nvpar << std::endl;
   std::cout << " Highest number of parameters defined by user = " << nparx << std::endl;
   std::cout << " Status of covariance matrix = " << icstat << std::endl;
   std::cout << std::endl;


   double pvalue, perror, plow,phigh;
   TString pname;

   if(amin < 99999) {

   for (int i=0; i<m_bins_eta; i++){
       for (int j=0; j<m_bins_pt; j++){
         
           pMinuit->mnpout(m_bins_pt*i+j,pname,pvalue,perror,plow,phigh,ierflg);

           std::cout << " pname :" << pname << ", val: " << pvalue << std::endl;
           std::cout << "-------Eta : " << i+1 << ", pt : " << j+1 << std::endl;

           //ArrayLikelihoodEta[j]->SetBinContent(lu+1,pvalue);
	   // ArrayLikelihoodEta[j]->SetBinError(lu+1,perror);

 
           m_TH2_Rates->SetBinContent(i+1,j+1,pvalue);
           m_TH2_Rates->SetBinError(i+1,j+1,perror);
      }
    }
  }
  
 delete pMinuit;

 return icstat;
}
void Like::WriteRates(){


  TString path=m_outputDir+"/"+m_period+"/";
  if (m_period=="") path=m_outputDir+"/"; //shuhui: for all mc 

  system("mkdir -p "+path);

  TString filename=path+m_name+".root";


  TFile *pFile = new TFile(filename,"RECREATE");

  m_TH2_Rates->Write();
  m_Mll_OS->Write();
  m_Mll_SS->Write();


  pFile->Close();


  return;
}
void Like::PlotRates(){


  TCanvas *pR = new TCanvas;
  pR->cd();
  //pR->cd()->SetLogy();
  pR->cd()->SetLogz();

  m_TH2_Rates->Draw("COLZ,TEXT");

 
  TCanvas *pMll = new TCanvas();
  pMll->cd();
  pMll->SetLogy();
  m_Mll_OS->SetLineWidth(2);
  m_Mll_OS->GetXaxis()->SetTitle("M_{ee} [GeV]");
  m_Mll_OS->GetYaxis()->SetTitle("Events");

  m_Mll_SS->SetLineWidth(2);

  m_Mll_OS->SetLineColor(2);
  m_Mll_OS->SetMinimum(0.1); 

  m_Mll_OS->Draw("L");
  m_Mll_SS->Draw("L,same");


 return;
}
void Like::CreateHistos(){

 m_TH2_Rates=CreateEtaPtHisto(m_name);

 m_Mll_SS=new TH1F("SS","SS",m_up-m_low,m_low,m_up);
 m_Mll_OS=new TH1F("OS","OS",m_up-m_low,m_low,m_up);

 

 return;
}
void Like::PrintSetup(){

 std::cout << "## Initializing execution of Likelihood method with the following settings: " << std::endl;
 std::cout << "-- isData :" << m_isData << std::endl;
 std::cout << "-- doBDT :" << m_doBDT << std::endl;
 std::cout << "-- subtractBkg : " << m_subtract_bkg << std::endl;
 std::cout << "-- Eta binning : " << m_bins_eta << std::endl;
 
 //be carefull: number of elements is number of bins +1!
 for (unsigned int i=0; i<m_bins_eta+1; i++) std::cout << m_eta_binning[i] << ", ";
 std::cout << " " << std::endl;

 std::cout << "-- Pt binning : " << m_bins_pt << std::endl;

 for (unsigned int i=0; i<m_bins_pt+1; i++) std::cout << m_pt_binning[i] << ", ";
 std::cout << " " << std::endl;

 std::cout << "-- Apply nJets cut: " << m_njets_cut << std::endl;
 std::cout << "-- Apply Ht cut: " << m_low_Ht <<", " << m_up_Ht << std::endl;


 return;
}
void Like::InitializeVectors(){

 std::cout << "## Initializing 4D vectors..." << std::endl;

 Array_SS.clear();
 Array_OS.clear();

 Array_SS.shrink_to_fit();
 Array_OS.shrink_to_fit();

 std::vector<float> pt2(m_bins_pt,0); 

 std::vector<std::vector<float>> pt1_pt2(m_bins_pt,pt2);

 std::vector<std::vector<std::vector<float> > > eta_2(m_bins_eta,pt1_pt2);

 for (unsigned int i=0; i<m_bins_eta; i++) {

      Array_SS.push_back(eta_2);

      Array_OS.push_back(eta_2);
 }

 



 return;
}
void Like::DebugVector(){

  float local=0;

  for (unsigned int eta1=0; eta1<Array_SS.size(); eta1++){

      for (unsigned int eta2=0; eta2<Array_SS[eta1].size(); eta2++){   

          for (unsigned int pt1=0; pt1<Array_SS[eta1][eta2].size(); pt1++){

               for (unsigned int pt2=0; pt2<Array_SS[eta1][eta2][pt1].size(); pt2++){

                     float ss=Array_SS[eta1][eta2][pt1][pt2];       
                     float os=Array_OS[eta1][eta2][pt1][pt2];
          

                     local+=ss+os;
               }
          } // end pt1
      } // end eta2
  }//end eta1


 std::cout << "**** --> Compare this number with the one above (it no bkg subtraction then they should match perfectly) --> Total : " << local << std::endl;

 return;
}
void NLL(int &npar, double *gin, double &f, double *par, Int_t iflag){

  int maxEta = NBINS_ETA;
  int maxPt = NBINS_PT;
   
  float ss = 0;
  float all = 0;
  int par1=-1;
  int par2=-1;

  double function=0;


  //Will loop over eta bins then over pt bins
  for (unsigned int eta_1=0; eta_1<maxEta; eta_1++){
    for (unsigned int eta_2 = eta_1; eta_2<maxEta; eta_2++){

        for (unsigned int pt_1=0; pt_1<maxPt; pt_1++){
          for (unsigned int pt_2=0; pt_2<maxPt; pt_2++){


                  ss=Array_SS[eta_1][eta_2][pt_1][pt_2];
                  all=ss+Array_OS[eta_1][eta_2][pt_1][pt_2];

                  par1=maxPt*eta_1 + pt_1;
                  par2=maxPt*eta_2 + pt_2;

    //              std::cout << "eta_1:" << eta_1 << ", eta_2:" << eta_2 << ", pt_1: " << pt_1 << ", pt_2:" << pt_2 << " --> par 1: "<< par1 << ",  par2 : "<< par2 << std::endl;

                  if (all<=0){
                      all=.000000000000000001;
                      //std::cout <<"-- This configuration has a negative number of events... setting to e-10" << std::endl;
                  }


                  if (ss!=0) function += ss*log(all*(par[par1]+par[par2])) - all*(par[par1]+par[par2]);
                  else function += -all*(par[par1]+par[par2]);


    
           } // end pt_2
        }// end pt_1

   }// end eta_2
  }// end eta_1


  f = -1*function;


  return;
}
std::vector<TString> Like::FillVector(TString file){

  std::vector<TString> samples;

  std::ifstream infile(file.Data());

  if (!infile){
        std::cout <<"## File with the list of samples: '"<< file <<"' does not exist --- ABORTING " << std::endl;
        exit(1);
  }
  else std::cout << "## Reading file list from : " << file << std::endl;

  std::string line;


  while (std::getline(infile,line) ){

        if (line.empty()) continue;

        TString ss = line;

        samples.push_back(ss);
  }

 return samples;

}
TH2F* Like::CreateEtaPtHisto(TString label){

  TH2F *hist=new TH2F(m_name+"_all_num",m_name+"_all_num",m_bins_eta,m_eta_binning,m_bins_pt,m_pt_binning);
  hist->GetYaxis()->SetTitle("Lepton p_{T} [GeV]");
  hist->GetXaxis()->SetTitle("Lepton #eta");
  hist->Reset();

  hist->Sumw2();

 return hist;
}
