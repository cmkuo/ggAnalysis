#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"
#endif

void TrainElectronMVA_Trig_BDTG_V0() {
  
  TMVA::Tools::Instance();
  TFile* outputFile = TFile::Open("ElectronMVA_Pt5To35_Trig_BDTG.root", "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory("DanieleMVA", outputFile, "!V:!Silent");
  



  // tracking
  factory->AddVariable("fbrem", 'F');
  factory->AddVariable("kfchi2", 'F');
  factory->AddVariable("kfhits", 'I');
  factory->AddVariable("gsfchi2", 'F');

  // geom. matching
  factory->AddVariable("deta", 'F');
  factory->AddVariable("dphi", 'F');
  factory->AddVariable("detacalo", 'F');
  //  factory->AddVariable("dphicalo", 'F');  // pruned
  
  // pure ecal
  factory->AddVariable("see", 'F');
  factory->AddVariable("spp", 'F');
  factory->AddVariable("etawidth", 'F');
  factory->AddVariable("phiwidth", 'F');
  factory->AddVariable("e1x5e5x5", 'F');
  factory->AddVariable("R9", 'F');
  //  factory->AddVariable("nbrems", 'F');   //pruned 

  // energy matching
  factory->AddVariable("HoE", 'F');
  factory->AddVariable("EoP", 'F');
  factory->AddVariable("IoEmIoP", 'F');
  factory->AddVariable("eleEoPout", 'F');
  //  factory->AddVariable("EoPout", 'F');    //pruned
  
  factory->AddVariable("PreShowerOverRaw", 'F');

  // ip
  factory->AddVariable("d0", 'F');
  factory->AddVariable("ip3d", 'F');
  
  factory->AddSpectator( "eta");  
  factory->AddSpectator( "pt");
  //  factory->AddSpectator( "matchConv");


  //Split using parity of event number
  TFile* inputSignalAll = TFile::Open("./newtmva_Em_26Mar_Zmc.LooseTrig.Real.All.root");   
  TFile* inputBkgAll = TFile::Open("./newtmva_Em_26Mar_Trig.LooseTrig.Fakes.All.root");
 


  TTree *signalTraining     = (TTree*)inputSignalAll->Get("ss");
  TTree *backgroundTraining = (TTree*)inputBkgAll->Get("ss");
  TTree *signalTesting     = (TTree*)inputSignalAll->Get("st");
  TTree *backgroundTesting = (TTree*)inputBkgAll->Get("st");  
  factory->AddSignalTree    (signalTraining,1.0,TMVA::Types::kTraining);
  factory->AddBackgroundTree(backgroundTraining,1.0,TMVA::Types::kTraining);
  factory->AddSignalTree    (signalTesting,1.0,TMVA::Types::kTesting);
  factory->AddBackgroundTree(backgroundTesting,1.0,TMVA::Types::kTesting);
  factory->PrepareTrainingAndTestTree("pt > 10 && abs(dz) < 0.1",
                                      "pt > 10 && pt < 35 && abs(dz) < 0.1",
				      "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );





  //*****************************************************************
  //BDTG TrigV0
  //*****************************************************************

  TMVA::MethodBase* bdtCat_BDTG_TrigV0 = factory->BookMethod( TMVA::Types::kCategory, "BDTCat_BDTG_TrigV0","" );
  TMVA::MethodCategory* category_BDTG_TrigV0 = dynamic_cast<TMVA::MethodCategory*>(bdtCat_BDTG_TrigV0);
  category_BDTG_TrigV0->AddMethod("pt < 20 && abs(eta) <= 0.8",
		      "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:",
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_1",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  
  
  
  category_BDTG_TrigV0->AddMethod("pt < 20 && abs(eta) > 0.8 && abs(eta)<=1.485",
		       "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:",
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_2",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  
  
  
  category_BDTG_TrigV0->AddMethod("pt < 20 && abs(eta)> 1.485",
		       "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:PreShowerOverRaw:",
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_3",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  


  
  category_BDTG_TrigV0->AddMethod("pt >= 20 && abs(eta) <= 0.8",
		      "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:",	
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_4",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  
  
  
  category_BDTG_TrigV0->AddMethod("pt >= 20 && abs(eta) > 0.8 && abs(eta)<=1.485",
		      "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:",	
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_5",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  
  
  
  category_BDTG_TrigV0->AddMethod("pt >= 20 && abs(eta)> 1.485",
		       "fbrem:deta:dphi:see:spp:e1x5e5x5:R9:etawidth:phiwidth:detacalo:HoE:EoP:eleEoPout:IoEmIoP:gsfchi2:kfchi2:kfhits:d0:ip3d:PreShowerOverRaw:",   
				TMVA::Types::kBDT, 
				"Category_BDTG_TrigV0_6",
				"!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );
  




  

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  outputFile->Close();

  delete factory;
}
