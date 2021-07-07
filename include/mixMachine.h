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

  bool CheckMachinesMatchMode(mixMachine* machineToCheck);
  bool CheckMachinesMatchMC(mixMachine* machineToCheck);
  bool CheckMachinesMatch2D(mixMachine* machineToCheck);
  bool CheckMachinesMatchBins(mixMachine* machineToCheck, double precision = 0.0001);

  mixMode GetMixMode(){return m_mixMode;}
  bool GetIsMC(){return m_isMC;}
  bool GetIs2DUnfold(){return m_is2DUnfold;}
  std::string GetMixMachineName(){return m_mixMachineName;}

  TH1F* GetTH1FPtr(std::string histType);
  TH2F* GetTH2FPtr(std::string histType);

  Int_t GetNBinsX(){return m_nBinsX;}
  std::vector<double> GetBinsX(){return m_binsXVect;}
  Int_t GetNBinsY(){return m_nBinsY;}
  std::vector<double> GetBinsY(){return m_binsYVect;}

  std::vector<TH1F*> GetTH1F(){return m_hists1D;}
  std::vector<TH2F*> GetTH2F(){return m_hists2D;}

  bool Add(mixMachine* machineToAdd, double precision = 0.0001);
  bool Add(mixMachine* machineToAdd1, mixMachine* machineToAdd2, double precision = 0.0001);
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
  std::vector<double> m_binsXVect;

  Int_t m_nBinsY;
  Double_t m_binsY[nMaxBins+1];
  std::vector<double> m_binsYVect;

  std::string m_titleX, m_titleY;
};

#endif
