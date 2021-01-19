//Author: Chris McGinn (2021.01.11)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TTree.h"

//RooUnfold
//https://gitlab.cern.ch/RooUnfold/RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
//#include "include/configParser.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"

int gdjHistToUnfold(std::string inConfigFileName)
{
  cppWatch globalTimer;
  globalTimer.start();

  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101             
  TRandom3* randGen_p = new TRandom3(randSeed);

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  //  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INRESPONSEFILENAME",
					      "INUNFOLDFILENAME",
                                              "OUTFILENAME"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
 
  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  
  TH2F* gammaPtJetXJReco_HalfResponse_p = nullptr;
  TH2F* gammaPtJetXJTruth_HalfResponse_p = nullptr;
  TH2F* gammaPtJetXJReco_HalfToUnfold_p = nullptr;
  TH2F* gammaPtJetXJTruth_HalfToUnfold_p = nullptr;

  TH2F* gammaPtJetXJReco_p = nullptr;
  TH2F* gammaPtJetXJTruth_p = nullptr;  
  
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");  
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  TFile* inUnfoldFile_p = nullptr;
  TEnv* inUnfoldFileConfig_p = nullptr;
  TEnv* inUnfoldFileLabels_p = nullptr;

  const bool isSameInFile = isStrSame(inResponseFileName, inUnfoldFileName);
  
  if(isSameInFile){
    inUnfoldFile_p = inResponseFile_p;
    inUnfoldFileConfig_p = inResponseFileConfig_p;
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }
  else{
    inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");  
    inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");    
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }

  std::vector<std::string> necessaryInFileParams = {"ISMC",
						    "NGAMMAPTBINS",
						    "GAMMAPTBINSLOW",
						    "GAMMAPTBINSHIGH",
						    "GAMMAPTBINSDOLOG",
						    "ISPP"};

 
  std::vector<std::string> pbpbParams = {"CENTBINS"};
  std::vector<std::string> mcParams = {"ASSOCGENMINPT"};

  std::vector<std::string> paramsToCompare = {"JETR",
					      "ISPP",
					      "CENTBINS",
					      "GAMMAEXCLUSIONDR",
					      "MIXJETEXCLUSIONDR",
					      "NGAMMAPTBINS",
					      "GAMMAPTBINSLOW",
					      "GAMMAPTBINSHIGH",
					      "GAMMAPTBINSDOLOG",
					      "NJTPTBINS",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "JTPTBINSDOLOG",
					      "GAMMAJTDPHI",
					      "NXJBINS",
					      "XJBINSLOW",
					      "XJBINSHIGH",
					      "XJBINSDOCUSTOM",
					      "XJBINSCUSTOM"};
    
  const std::string gammaJtDPhiStr = "DPhi0";
  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";
  
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  //Grab first set of parameters defined in necessaryInFileParams
  const Int_t nMaxPtBins = 200;
  //  const Int_t nMaxCentBins = 10;
  const bool isMC = inResponseFileConfig_p->GetValue("ISMC", 1);
  //  const bool isPP = inResponseFileConfig_p->GetValue("ISPP", 0);

  const Int_t nGammaPtBins = inResponseFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  Double_t gammaPtBins[nMaxPtBins+1];
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  const Int_t nXJBins = inResponseFileConfig_p->GetValue("NXJBINS", 20);
  const Float_t xjBinsLow = inResponseFileConfig_p->GetValue("XJBINSLOW", 0.0);
  const Float_t xjBinsHigh = inResponseFileConfig_p->GetValue("XJBINSHIGH", 2.0);
  const Bool_t xjBinsDoCustom = inResponseFileConfig_p->GetValue("XJBINSDOCUSTOM", 0);
  Double_t xjBins[nMaxPtBins+1];
  if(!xjBinsDoCustom) getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
  else{
    std::string xjBinsCustomStr = inResponseFileConfig_p->GetValue("XJBINSCUSTOM", "");
    std::vector<float> xjBinsCustomVect = strToVectF(xjBinsCustomStr);
    for(unsigned int i = 0; i < xjBinsCustomVect.size(); ++i){
      xjBins[i] = xjBinsCustomVect[i];
    }
  }
  
  //Enforce additional requirements on config depending on first set
  if(isMC){
    if(!checkEnvForParams(inResponseFileConfig_p, mcParams)) return 1;
  }

  /*
  int nCentBins = 1;
  std::vector<std::string> centBinsStr = {"PP"};
  if(!isPP){
    if(!checkEnvForParams(inResponseFileConfig_p, pbpbParams)) return 1;

    centBinsStr = strToVect(inResponseFileConfig_p->GetValue("CENTBINS", ""));
    nCentBins = centBinsStr.size()-1;
  }
  */
  //Enforce matching over specified parameters
  if(isSameInFile){
    if(!compEnvParams(inResponseFileConfig_p, inUnfoldFileConfig_p, paramsToCompare)) return 1;
  }

  //Construct matrices
  outFile_p->cd();
  gammaPtJetXJReco_HalfResponse_p = new TH2F("gammaPtJetXJReco_HalfResponse_h", ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
  gammaPtJetXJTruth_HalfResponse_p = new TH2F("gammaPtJetXJTruth_HalfResponse_h", ";Truth x_{J};Truth Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
  gammaPtJetXJReco_HalfToUnfold_p = new TH2F("gammaPtJetXJReco_HalfToUnfold_h", ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
  gammaPtJetXJTruth_HalfToUnfold_p = new TH2F("gammaPtJetXJTruth_HalfToUnfold_h", ";Truth x_{J};Truth Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

  gammaPtJetXJReco_p = new TH2F("gammaPtJetXJReco_h", ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
  gammaPtJetXJTruth_p = new TH2F("gammaPtJetXJTruth_h", ";Truth x_{J};Truth Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

  RooUnfoldResponse* rooRes_Half_p = nullptr;
  RooUnfoldResponse* rooRes_p = nullptr;
  
  setSumW2({gammaPtJetXJReco_HalfResponse_p, gammaPtJetXJTruth_HalfResponse_p, gammaPtJetXJReco_HalfToUnfold_p, gammaPtJetXJTruth_HalfToUnfold_p, gammaPtJetXJReco_p, gammaPtJetXJTruth_p});
  
  //Grab the unfolding ttree for matrix filling
  Float_t recoGammaPt_, truthGammaPt_;
  Float_t recoXJ_, truthXJ_;
  Float_t unfoldWeight_;

  TTree* unfoldXJTree_p = (TTree*)inResponseFile_p->Get("unfoldXJTree_p");
  unfoldXJTree_p->SetBranchAddress("recoGammaPt", &recoGammaPt_);
  unfoldXJTree_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldXJTree_p->SetBranchAddress("recoXJ", &recoXJ_);
  unfoldXJTree_p->SetBranchAddress("truthXJ", &truthXJ_);
  unfoldXJTree_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  
  //Construct our MC hists
  const ULong64_t nEntries = unfoldXJTree_p->GetEntries();
 
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    unfoldXJTree_p->GetEntry(entry);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;
    if(fillResponseHalf){
      gammaPtJetXJReco_HalfResponse_p->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
      gammaPtJetXJTruth_HalfResponse_p->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    else{
      gammaPtJetXJReco_HalfToUnfold_p->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
      gammaPtJetXJTruth_HalfToUnfold_p->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    
    gammaPtJetXJReco_p->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
    gammaPtJetXJTruth_p->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
  }

  outFile_p->cd();
  rooRes_Half_p = new RooUnfoldResponse(gammaPtJetXJReco_HalfResponse_p, gammaPtJetXJTruth_HalfResponse_p, "rooRes_Half", "");
  rooRes_p = new RooUnfoldResponse(gammaPtJetXJReco_p, gammaPtJetXJTruth_p, "rooRes", "");

  delete randGen_p;
  randGen_p = new TRandom3(randSeed); //Resetting to get the correct response matrix fills
  
  
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    unfoldXJTree_p->GetEntry(entry);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;
    if(fillResponseHalf) rooRes_Half_p->Fill(recoXJ_, recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);
    rooRes_p->Fill(recoXJ_, recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);
  }

  inResponseFile_p->Close();
  delete inResponseFile_p;

  if(!isSameInFile){
    inUnfoldFile_p->Close();
    delete inUnfoldFile_p;
  }

  outFile_p->cd();

  const int nIter = 21;

  for(int gI = 0; gI < nGammaPtBins; ++gI){
    TH1D* xjTruthProjection_p = new TH1D(("xjTruth_HalfToUnfold_GammaPt" + std::to_string(gI) + "_h").c_str(), ";x_{J};Arbitrary", nXJBins, xjBins);
    TH1D* xjRecoProjection_p = new TH1D(("xjReco_HalfToUnfold_GammaPt" + std::to_string(gI) + "_h").c_str(), ";x_{J};Arbitrary", nXJBins, xjBins);

    for(int xI = 0; xI < nXJBins; ++xI){
      xjTruthProjection_p->SetBinContent(xI+1, gammaPtJetXJTruth_HalfToUnfold_p->GetBinContent(xI+1, gI+1));
      xjTruthProjection_p->SetBinError(xI+1, gammaPtJetXJTruth_HalfToUnfold_p->GetBinError(xI+1, gI+1));      

      xjRecoProjection_p->SetBinContent(xI+1, gammaPtJetXJReco_HalfToUnfold_p->GetBinContent(xI+1, gI+1));
      xjRecoProjection_p->SetBinError(xI+1, gammaPtJetXJReco_HalfToUnfold_p->GetBinError(xI+1, gI+1));      
    }

    xjTruthProjection_p->Write("", TObject::kOverwrite);
    delete xjTruthProjection_p;

    xjRecoProjection_p->Write("", TObject::kOverwrite);
    delete xjRecoProjection_p;
  }


  for(unsigned int i = 1; i < nIter; ++ i){
    RooUnfoldBayes* rooBayes_Half_p = new RooUnfoldBayes(rooRes_Half_p, gammaPtJetXJReco_HalfToUnfold_p, i);  

    TH2D* unfolded_p = (TH2D*)rooBayes_Half_p->Hreco();//Very annoyingly Hreco is the name of the function that returns unfolded distribution

    unfolded_p->Write(("gammaPtJetXJReco_HalfToUnfold_Iter" + std::to_string(i) + "_h").c_str(), TObject::kOverwrite);

    for(int gI = 0; gI < nGammaPtBins; ++gI){
      TH1D* xjUnfoldedProjection_p = new TH1D(("xjReco_HalfToUnfold_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_h").c_str(), ";x_{J};Arbitrary", nXJBins, xjBins);

      for(int xI = 0; xI < nXJBins; ++xI){
	xjUnfoldedProjection_p->SetBinContent(xI+1, unfolded_p->GetBinContent(xI+1, gI+1));
	xjUnfoldedProjection_p->SetBinError(xI+1, unfolded_p->GetBinError(xI+1, gI+1));
      }

      outFile_p->cd();

      xjUnfoldedProjection_p->Write("", TObject::kOverwrite);
      delete xjUnfoldedProjection_p;
    }
    
    delete rooBayes_Half_p;
  }
  
  gammaPtJetXJReco_HalfResponse_p->Write("", TObject::kOverwrite);
  gammaPtJetXJTruth_HalfResponse_p->Write("", TObject::kOverwrite);
  gammaPtJetXJReco_HalfToUnfold_p->Write("", TObject::kOverwrite);
  gammaPtJetXJTruth_HalfToUnfold_p->Write("", TObject::kOverwrite);

  gammaPtJetXJReco_p->Write("", TObject::kOverwrite);
  gammaPtJetXJTruth_p->Write("", TObject::kOverwrite);

  delete gammaPtJetXJReco_HalfResponse_p;
  delete gammaPtJetXJTruth_HalfResponse_p;
  delete gammaPtJetXJReco_HalfToUnfold_p;
  delete gammaPtJetXJTruth_HalfToUnfold_p;

  delete gammaPtJetXJReco_p;
  delete gammaPtJetXJTruth_p;

  delete rooRes_Half_p;
  delete rooRes_p;

  inUnfoldFileConfig_p->SetValue("NITER", nIter);
  inUnfoldFileConfig_p->Write("config", TObject::kOverwrite);
  inUnfoldFileLabels_p->Write("label", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;
  
  std::cout << "GDJHISTTOUNFOLD COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHistToUnfold.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjHistToUnfold(argv[1]);
  return retVal;
}
