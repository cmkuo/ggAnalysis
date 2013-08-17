//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 20 12:23:43 2012 by ROOT version 5.27/06b
// from TTree st/a Tree
// found on file: ../DataForTraining/EmanueleDATANtuples/newtmva_Em_19Mar_MC.Real.All.root
//////////////////////////////////////////////////////////

#ifndef MakeTestMVAMacro_000_h
#define MakeTestMVAMacro_000_h
#include "TClass.h"
#include "TError.h"
#include "TLorentzVector.h"

#include <string>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class MakeTestMVAMacro_000 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         eleEoPout;
   Float_t         EoPout;
   Float_t         EoP;
   Float_t         IoEmIoP;
   Float_t         HoE;
   Float_t         eledeta;
   Float_t         deta;
   Float_t         dphi;
   Float_t         detacalo;
   Float_t         dphicalo;
   Float_t         s9s25;
   Float_t         phiwidth;
   Float_t         etawidth;
   Float_t         see;
   Float_t         sep;
   Float_t         spp;
   Float_t         fbrem;
   Int_t           nbrem;
   Int_t           missHits;
   Float_t         dist;
   Float_t         dcot;
   Float_t         d0;
   Float_t         dz;
   Float_t         ip3d;
   Float_t         ip3ds;
   Int_t           kfhits;
   Float_t         kfchi2;
   Float_t         gsfchi2;
   Float_t         e1x5e5x5;
   Float_t         SeedEMaxOverE;
   Float_t         SeedETopOverE;
   Float_t         SeedEBottomOverE;
   Float_t         SeedELeftOverE;
   Float_t         SeedERightOverE;
   Float_t         SeedE2ndOverE;
   Float_t         SeedE2x5RightOverE;
   Float_t         SeedE2x5LeftOverE;
   Float_t         SeedE2x5TopOverE;
   Float_t         SeedE2x5BottomOverE;
   Float_t         SeedE2x5MaxOverE;
   Float_t         SeedE1x5OverE;
   Float_t         SeedE2x2OverE;
   Float_t         SeedE3x3OverE;
   Float_t         SeedE5x5OverE;
   Float_t         R9;
   Int_t           matchConv;
   Int_t           ecaldriven;
   Float_t         scenergy;
   Float_t         scrawenergy;
   Float_t         scesenergy;
   Float_t         pt;
   Float_t         eta;
   Float_t         phi;
   Int_t           charge;
   Int_t           DenomFake;
   Int_t           DenomFakeSmurf;
   Float_t         trkIso;
   Float_t         ecalIso;
   Float_t         hcalIso;
   Float_t         combPFIsoHWW;
   Float_t         chaPFIso;
   Float_t         neuPFIso;
   Float_t         phoPFIso;
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Float_t         npu;
   Float_t         mcmatch;
   Float_t         bdthww;
   Float_t         bdthwwnoip;
   Float_t         bdthzz;
   Float_t         bdthzznoip;
   Float_t         bdthzzmc;
   Float_t         lh;
   Float_t         pfmva;
   Float_t         vertices;
   Float_t         rho;
   Float_t         PreShowerOverRaw;

   // List of branches
   TBranch        *b_eleEoPout;   //!
   TBranch        *b_EoPout;   //!
   TBranch        *b_EoP;   //!
   TBranch        *b_IoEmIoP;   //!
   TBranch        *b_HoE;   //!
   TBranch        *b_eledeta;   //!
   TBranch        *b_deta;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_detacalo;   //!
   TBranch        *b_dphicalo;   //!
   TBranch        *b_s9s25;   //!
   TBranch        *b_phiwidth;   //!
   TBranch        *b_etawidth;   //!
   TBranch        *b_see;   //!
   TBranch        *b_sep;   //!
   TBranch        *b_spp;   //!
   TBranch        *b_fbrem;   //!
   TBranch        *b_nbrems;   //!
   TBranch        *b_missHits;   //!
   TBranch        *b_dist;   //!
   TBranch        *b_dcot;   //!
   TBranch        *b_d0;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_ip3d;   //!
   TBranch        *b_ip3ds;   //!
   TBranch        *b_kfhits;   //!
   TBranch        *b_kfchi2;   //!
   TBranch        *b_gsfchi2;   //!
   TBranch        *b_e1x5e5x5;   //!
   TBranch        *b_SeedEMaxOverE;   //!
   TBranch        *b_SeedETopOverE;   //!
   TBranch        *b_SeedEBottomOverE;   //!
   TBranch        *b_SeedELeftOverE;   //!
   TBranch        *b_SeedERightOverE;   //!
   TBranch        *b_SeedE2ndOverE;   //!
   TBranch        *b_SeedE2x5RightOverE;   //!
   TBranch        *b_SeedE2x5LeftOverE;   //!
   TBranch        *b_SeedE2x5TopOverE;   //!
   TBranch        *b_SeedE2x5BottomOverE;   //!
   TBranch        *b_SeedE2x5MaxOverE;   //!
   TBranch        *b_SeedE1x5OverE;   //!
   TBranch        *b_SeedE2x2OverE;   //!
   TBranch        *b_SeedE3x3OverE;   //!
   TBranch        *b_SeedE5x5OverE;   //!
   TBranch        *b_R9;   //!
   TBranch        *b_matchConv;   //!
   TBranch        *b_ecaldriven;   //!
   TBranch        *b_scenergy;   //!
   TBranch        *b_scrawenergy;   //!
   TBranch        *b_scesenergy;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_DenomFake;   //!
   TBranch        *b_DenomFakeSmurf;   //!
   TBranch        *b_trkIso;   //!
   TBranch        *b_ecalIso;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_combPFIsoHWW;   //!
   TBranch        *b_chPFIso;   //!
   TBranch        *b_neuPFIso;   //!
   TBranch        *b_phoPFIso;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_mcmatch;   //!
   TBranch        *b_bdthww;   //!
   TBranch        *b_bdthwwnoip;   //!
   TBranch        *b_bdthzz;   //!
   TBranch        *b_bdthzznoip;   //!
   TBranch        *b_bdthzzmc;   //!
   TBranch        *b_lh;   //!
   TBranch        *b_pfmva;   //!
   TBranch        *b_vertices;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_PreShowerOverRaw;   //!

   string myfile_,myoutputfile_;
   bool sig_;


   MakeTestMVAMacro_000(TTree *tree=0, string filename = "", bool sig = false);
   virtual ~MakeTestMVAMacro_000();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MakeTestMVAMacro_000_cxx
MakeTestMVAMacro_000::MakeTestMVAMacro_000(TTree *tree, string filename, bool sig)
{
  
 
  myfile_ = filename;
  
  sig_ = sig;
  string sigback = "histo_back.root";
  if(sig)
    sigback = "histo_sig.root";
  
  myoutputfile_ = sigback;
  
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(myfile_.c_str());
    if (!f) {
      f = new TFile(myfile_.c_str());
    }
    tree = (TTree*)gDirectory->Get("ss");
    
  }
  Init(tree);
}
MakeTestMVAMacro_000::~MakeTestMVAMacro_000()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MakeTestMVAMacro_000::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MakeTestMVAMacro_000::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MakeTestMVAMacro_000::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eleEoPout", &eleEoPout, &b_eleEoPout);
   fChain->SetBranchAddress("EoPout", &EoPout, &b_EoPout);
   fChain->SetBranchAddress("EoP", &EoP, &b_EoP);
   fChain->SetBranchAddress("IoEmIoP", &IoEmIoP, &b_IoEmIoP);
   fChain->SetBranchAddress("HoE", &HoE, &b_HoE);
   fChain->SetBranchAddress("eledeta", &eledeta, &b_eledeta);
   fChain->SetBranchAddress("deta", &deta, &b_deta);
   fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
   fChain->SetBranchAddress("detacalo", &detacalo, &b_detacalo);
   fChain->SetBranchAddress("dphicalo", &dphicalo, &b_dphicalo);
   fChain->SetBranchAddress("s9s25", &s9s25, &b_s9s25);
   fChain->SetBranchAddress("phiwidth", &phiwidth, &b_phiwidth);
   fChain->SetBranchAddress("etawidth", &etawidth, &b_etawidth);
   fChain->SetBranchAddress("see", &see, &b_see);
   fChain->SetBranchAddress("sep", &sep, &b_sep);
   fChain->SetBranchAddress("spp", &spp, &b_spp);
   fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
   fChain->SetBranchAddress("nbrem", &nbrem, &b_nbrems);
   fChain->SetBranchAddress("missHits", &missHits, &b_missHits);
   fChain->SetBranchAddress("dist", &dist, &b_dist);
   fChain->SetBranchAddress("dcot", &dcot, &b_dcot);
   fChain->SetBranchAddress("d0", &d0, &b_d0);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("ip3d", &ip3d, &b_ip3d);
   fChain->SetBranchAddress("ip3ds", &ip3ds, &b_ip3ds);
   fChain->SetBranchAddress("kfhits", &kfhits, &b_kfhits);
   fChain->SetBranchAddress("kfchi2", &kfchi2, &b_kfchi2);
   fChain->SetBranchAddress("gsfchi2", &gsfchi2, &b_gsfchi2);
   fChain->SetBranchAddress("e1x5e5x5", &e1x5e5x5, &b_e1x5e5x5);
   fChain->SetBranchAddress("SeedEMaxOverE", &SeedEMaxOverE, &b_SeedEMaxOverE);
   fChain->SetBranchAddress("SeedETopOverE", &SeedETopOverE, &b_SeedETopOverE);
   fChain->SetBranchAddress("SeedEBottomOverE", &SeedEBottomOverE, &b_SeedEBottomOverE);
   fChain->SetBranchAddress("SeedELeftOverE", &SeedELeftOverE, &b_SeedELeftOverE);
   fChain->SetBranchAddress("SeedERightOverE", &SeedERightOverE, &b_SeedERightOverE);
   fChain->SetBranchAddress("SeedE2ndOverE", &SeedE2ndOverE, &b_SeedE2ndOverE);
   fChain->SetBranchAddress("SeedE2x5RightOverE", &SeedE2x5RightOverE, &b_SeedE2x5RightOverE);
   fChain->SetBranchAddress("SeedE2x5LeftOverE", &SeedE2x5LeftOverE, &b_SeedE2x5LeftOverE);
   fChain->SetBranchAddress("SeedE2x5TopOverE", &SeedE2x5TopOverE, &b_SeedE2x5TopOverE);
   fChain->SetBranchAddress("SeedE2x5BottomOverE", &SeedE2x5BottomOverE, &b_SeedE2x5BottomOverE);
   fChain->SetBranchAddress("SeedE2x5MaxOverE", &SeedE2x5MaxOverE, &b_SeedE2x5MaxOverE);
   fChain->SetBranchAddress("SeedE1x5OverE", &SeedE1x5OverE, &b_SeedE1x5OverE);
   fChain->SetBranchAddress("SeedE2x2OverE", &SeedE2x2OverE, &b_SeedE2x2OverE);
   fChain->SetBranchAddress("SeedE3x3OverE", &SeedE3x3OverE, &b_SeedE3x3OverE);
   fChain->SetBranchAddress("SeedE5x5OverE", &SeedE5x5OverE, &b_SeedE5x5OverE);
   fChain->SetBranchAddress("R9", &R9, &b_R9);
   fChain->SetBranchAddress("matchConv", &matchConv, &b_matchConv);
   fChain->SetBranchAddress("ecaldriven", &ecaldriven, &b_ecaldriven);
   fChain->SetBranchAddress("scenergy", &scenergy, &b_scenergy);
   fChain->SetBranchAddress("scrawenergy", &scrawenergy, &b_scrawenergy);
   fChain->SetBranchAddress("scesenergy", &scesenergy, &b_scesenergy);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("DenomFake", &DenomFake, &b_DenomFake);
   fChain->SetBranchAddress("DenomFakeSmurf", &DenomFakeSmurf, &b_DenomFakeSmurf);
   fChain->SetBranchAddress("trkIso", &trkIso, &b_trkIso);
   fChain->SetBranchAddress("ecalIso", &ecalIso, &b_ecalIso);
   fChain->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso);
   fChain->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW, &b_combPFIsoHWW);
   fChain->SetBranchAddress("chaPFIso", &chaPFIso, &b_chPFIso);
   fChain->SetBranchAddress("neuPFIso", &neuPFIso, &b_neuPFIso);
   fChain->SetBranchAddress("phoPFIso", &phoPFIso, &b_phoPFIso);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("mcmatch", &mcmatch, &b_mcmatch);
   fChain->SetBranchAddress("bdthww", &bdthww, &b_bdthww);
   fChain->SetBranchAddress("bdthwwnoip", &bdthwwnoip, &b_bdthwwnoip);
   fChain->SetBranchAddress("bdthzz", &bdthzz, &b_bdthzz);
   fChain->SetBranchAddress("bdthzznoip", &bdthzznoip, &b_bdthzznoip);
   fChain->SetBranchAddress("bdthzzmc", &bdthzzmc, &b_bdthzzmc);
   fChain->SetBranchAddress("lh", &lh, &b_lh);
   fChain->SetBranchAddress("pfmva", &pfmva, &b_pfmva);
   fChain->SetBranchAddress("vertices", &vertices, &b_vertices);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("PreShowerOverRaw", &PreShowerOverRaw, &b_PreShowerOverRaw);
   Notify();
}

Bool_t MakeTestMVAMacro_000::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MakeTestMVAMacro_000::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef MakeTestMVAMacro_000_cxx
