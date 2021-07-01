#ifndef MIXMACHINE_H
#define MIXMACHINE_H

//cpp dependencies
#include <string>
#include <vector>

//ROOT dependencies
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class mixMachine{
 public:
  enum mixMode{NONE=0,//Dummy mode - if on things should fail
	       INCLUSIVE= 1,
	       MULTI=2	  
  };

  mixMachine(){m_isInit = false;};
  ~mixMachine(){};

  mixMachine(std::string inMixMachineName, mixMachine::mixMode inMixMode, TEnv* inParams_p);
  bool Init(std::string inMixMachineName, mixMachine::mixMode inMixMode, TEnv* inParams_p);
  bool FillXY(Float_t fillX, Float_t fillY, Float_t fillWeight, std::string mixName);  
  bool FillXYRaw(Float_t fillX, Float_t fillY, Float_t fillWeight);
  bool FillXYMix(Float_t fillX, Float_t fillY, Float_t fillWeight);
  bool FillXYMixCorrection(Float_t fillX, Float_t fillY, Float_t fillWeight);
  bool FillXYTruth(Float_t fillX, Float_t fillY, Float_t fillWeight);
  bool FillXYTruthMatchedReco(Float_t fillX, Float_t fillY, Float_t fillWeight);
  bool FillX(Float_t fillX, Float_t fillWeight, std::string mixName);
  bool FillXRaw(Float_t fillX, Float_t fillWeight);
  bool FillXMix(Float_t fillX, Float_t fillWeight);
  bool FillXMixCorrection(Float_t fillX, Float_t fillWeight);
  bool FillXTruth(Float_t fillX, Float_t fillWeight);
  bool FillXTruthMatchedReco(Float_t fillX, Float_t fillWeight);

  void ComputeSub();
  void WriteToFile(TFile* inFile_p);
  void WriteToDirectory(TDirectoryFile* inDir_p);

  void Print();
  void Clean();

 private:
  bool m_isInit;
  std::string m_mixMachineName;
  mixMachine::mixMode m_mixMode;
  std::vector<std::string> m_noneMixNames = {"RAW", "SUB", "TRUTH", "TRUTHMATCHEDRECO"};
  std::vector<std::string> m_inclusiveMixNames = {"RAW", "MIX", "SUB", "TRUTH", "TRUTHMATCHEDRECO"};
  std::vector<std::string> m_multiMixNames = {"RAW", "MIX", "MIXCORRECTION", "MIXCORRECTED", "SUB", "TRUTH", "TRUTHMATCHEDRECO"};

  std::vector<TH1F*> m_hists1D;
  std::vector<TH2F*> m_hists2D;
  TEnv m_env;
  bool m_is2DUnfold;
  bool m_isMC;

  const static Int_t nMaxBins = 200;
  Int_t m_nBinsX;
  Double_t m_binsX[nMaxBins+1];

  Int_t m_nBinsY;
  Double_t m_binsY[nMaxBins+1];

  std::string m_titleX, m_titleY;
};

#endif
