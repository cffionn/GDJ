#ifndef MIXMACHINE_H
#define MIXMACHINE_H

//cpp dependencies
#include <string>
#include <vector>

//ROOT dependencies
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

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
  bool FillXY(Double_t fillX, Double_t fillY, Double_t fillWeight, std::string mixName);  
  bool FillXYRaw(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYMix(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYMixCorrection(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYTruth(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYTruthWithRecoMatch(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYTruthNoRecoMatch(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYRawWithTruthMatch(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYRawNoTruthMatch(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYSingleTruthToMultiFake(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillXYSingleTruthToMultiFakeMix(Double_t fillX, Double_t fillY, Double_t fillWeight);
  bool FillX(Double_t fillX, Double_t fillWeight, std::string mixName);
  bool FillXRaw(Double_t fillX, Double_t fillWeight);
  bool FillXMix(Double_t fillX, Double_t fillWeight);
  bool FillXMixCorrection(Double_t fillX, Double_t fillWeight);
  bool FillXTruth(Double_t fillX, Double_t fillWeight);
  bool FillXTruthMatchedReco(Double_t fillX, Double_t fillWeight);
  bool FillXSingleTruthToMultiFake(Double_t fillX, Double_t fillWeight);
  bool FillXSingleTruthToMultiFakeMix(Double_t fillX, Double_t fillWeight);

  bool CheckMachinesMatchMode(mixMachine* machineToCheck);
  bool CheckMachinesMatchMC(mixMachine* machineToCheck);
  bool CheckMachinesMatch2D(mixMachine* machineToCheck);
  bool CheckMachinesMatchBins(mixMachine* machineToCheck, double precision = 0.0001);

  mixMode GetMixMode(){return m_mixMode;}
  bool GetIsMC(){return m_isMC;}
  bool GetIs2DUnfold(){return m_is2DUnfold;}
  std::string GetMixMachineName(){return m_mixMachineName;}

  TH1D* GetTH1DPtr(std::string histType);
  TH2D* GetTH2DPtr(std::string histType);

  Int_t GetNBinsX(){return m_nBinsX;}
  std::vector<double> GetBinsX(){return m_binsXVect;}
  Int_t GetNBinsY(){return m_nBinsY;}
  std::vector<double> GetBinsY(){return m_binsYVect;}

  std::vector<TH1D*> GetTH1D(){return m_hists1D;}
  std::vector<TH2D*> GetTH2D(){return m_hists2D;}

  bool Add(mixMachine* machineToAdd, double precision = 0.0001);
  bool Add(mixMachine* machineToAdd1, mixMachine* machineToAdd2, double precision = 0.0001);
  void ComputeSub();
  void WriteToFile(TFile* inFile_p);
  void WriteToDirectory(TDirectoryFile* inDir_p);

  void PushTrackingMap(int inVal);
  bool IsInTrackingMap(int inVal);
  void CleanTrackingMap();
  
  void Print();
  void Clean();

 private:
  bool m_isInit;
  std::string m_mixMachineName;
  mixMachine::mixMode m_mixMode;
  std::vector<std::string> m_noneMixNames = {"RAW", "SUB", "TRUTH", "TRUTHWITHRECOMATCH", "TRUTHNORECOMATCH", "RAWWITHTRUTHMATCH", "RAWNOTRUTHMATCH"};
  std::vector<std::string> m_inclusiveMixNames = {"RAW", "MIX", "SUB", "TRUTH", "TRUTHWITHRECOMATCH", "TRUTHNORECOMATCH", "RAWWITHTRUTHMATCH", "RAWNOTRUTHMATCH"};
  std::vector<std::string> m_multiMixNames = {"RAW", "MIX", "MIXCORRECTION", "MIXCORRECTED", "SUB", "TRUTH", "TRUTHWITHRECOMATCH", "TRUTHNORECOMATCH", "RAWWITHTRUTHMATCH", "RAWNOTRUTHMATCH", "SINGLETRUTHTOMULTIFAKE", "SINGLETRUTHTOMULTIFAKEMIX"};

  std::vector<TH1D*> m_hists1D;
  std::vector<TH2D*> m_hists2D;
  TEnv m_env;
  bool m_is2DUnfold;
  bool m_isMC;

  const static Int_t nMaxBins = 2000;
  Int_t m_nBinsX;
  Double_t m_binsX[nMaxBins+1];
  std::vector<double> m_binsXVect;

  Int_t m_nBinsY;
  Double_t m_binsY[nMaxBins+1];
  std::vector<double> m_binsYVect;

  std::string m_titleX, m_titleY;

  std::map<int,bool> m_trackingMap;
};

#endif
