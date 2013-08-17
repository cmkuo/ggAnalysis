{
  gROOT->Reset();
  gSystem->Load("../../../../../lib/slc5_amd64_gcc434/libEGammaEGammaAnalysisTools.so");  // Attention the lib name can change with CMSSW version
  gROOT->ProcessLine(".L MakeTestMVAMacro_000.C++");
  MakeTestMVAMacro_000 t(0,"/afs/cern.ch/user/b/benedet/public/MyTestNtupla/newtmva_Em_3Apr_MC.Real.root",true);
  t.Loop();
  gROOT->ProcessLine(".q");
}
