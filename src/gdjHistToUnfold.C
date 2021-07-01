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

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INRESPONSEFILENAME",
					      "INUNFOLDFILENAME",
					      "VARNAME",
                                              "OUTFILENAME",
					      "NITER",
					      "DOREWEIGHT",
					      "UNFOLDERRTYPE",
					      "NTOYS"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
 
  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string varName = config_p->GetValue("VARNAME", "");
  std::string varNameLower = returnAllLowercaseString(varName);

  const int nIter = config_p->GetValue("NITER", -1);
  const bool doReweight = config_p->GetValue("DOREWEIGHT", 0);

  const int unfoldErrType = config_p->GetValue("UNFOLDERRTYPE", 0);
  const int nToys = config_p->GetValue("NTOYS", 0);
  if(unfoldErrType < 0 || unfoldErrType > 1){//Valid error types are 0: default kErrors bin-by-bin errors (diagonal covariance matrix), 1: toys
    std::cout << "Given UNFOLDERRTYPE, " << unfoldErrType << ", is invalid. Please give valid type, return 1" << std::endl;    
    return 1;
  }
  if(nToys <= 0 && unfoldErrType == 1){
    std::cout << "NTOYS given \'" << nToys << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  
  const int nIterMax = 100;
  if(nIter > nIterMax){
    std::cout << "Requested number of iterations \'" << nIter << "\' exceeds maximum allowed \'" << nIterMax << "\'" << std::endl;
  }
  
  if(outFileName.find(".root") != std::string::npos){
    outFileName.replace(outFileName.rfind(".root"), 5, "");
    outFileName = outFileName + "_" + dateStr + ".root";
  }

  if(outFileName.find("/") == std::string::npos){
    outFileName = "output/" + dateStr + "/" + outFileName;
  }
  
  std::vector<std::string> validUnfoldVar = {"Pt", "XJ", "XJJ"};
  std::vector<std::string> validUnfoldVarTree = {"XJ", "XJ", "XJJ"};
  std::vector<std::string> validUnfoldVarStyle = {"Jet p_{T}", "x_{J}", "#vec{x}_{JJ#gamma}"};

  const int vectPos = vectContainsStrPos(varName, &validUnfoldVar);
  if(vectPos < 0){
    std::cout << "GIVEN VARNAME \'" << varName << "\' is invalid. Please pick from:" << std::endl;
    for(unsigned int i = 0; i < validUnfoldVar.size(); ++i){
      std::cout << " \'" << validUnfoldVar[i] << "\'" << std::endl;
    }
    std::cout << "return 1" << std::endl;
    return 1;
  }
  std::string varNameTree = validUnfoldVarTree[vectPos];
  std::string varNameStyle = validUnfoldVarStyle[vectPos];
  
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;
  
  const Int_t nMaxCentBins = 10;
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  

  TH1F* photonPtReco_HalfResponse_p[nMaxCentBins];
  TH1F* photonPtTruth_HalfResponse_p[nMaxCentBins];
  TH1F* photonPtReco_HalfToUnfold_p[nMaxCentBins];
  TH1F* photonPtTruth_HalfToUnfold_p[nMaxCentBins];

  TH2F* photonPtJetXJReco_HalfResponse_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_HalfResponse_p[nMaxCentBins];
  TH2F* photonPtJetXJReco_HalfToUnfold_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_HalfToUnfold_p[nMaxCentBins];

  TH1F* photonPtReco_p[nMaxCentBins];
  TH1F* photonPtTruth_p[nMaxCentBins];  

  TH2F* photonPtJetXJReco_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_p[nMaxCentBins];  
  
  const Int_t nBarrelAndEC = 2;
  const std::string barrelAndECStr[nBarrelAndEC] = {"Barrel", "EC"};

  TH1F* photonPtReco_PURCORR_p[nMaxCentBins][nBarrelAndEC];  
  TH1F* photonPtReco_PURCORR_COMBINED_p[nMaxCentBins];  
  TH1F* photonPtReco_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH1F* photonPtReco_TRUTH_COMBINED_p[nMaxCentBins];

  TH2F* photonPtJetXJReco_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJetXJReco_PURCORR_COMBINED_p[nMaxCentBins];
  TH2F* photonPtJetXJReco_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJetXJReco_TRUTH_COMBINED_p[nMaxCentBins];
  
  RooUnfoldResponse* rooResGamma_Half_p[nMaxCentBins];
  RooUnfoldResponse* rooResGamma_p[nMaxCentBins];
  
  RooUnfoldResponse* rooResGammaJet_Half_p[nMaxCentBins];
  RooUnfoldResponse* rooResGammaJet_p[nMaxCentBins];

  std::cout << inResponseFileName << std::endl;
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");  
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  TFile* inUnfoldFile_p = nullptr;
  TEnv* inUnfoldFileConfig_p = nullptr;
  TEnv* inUnfoldFileLabels_p = nullptr;

  std::cout << "INUNFOLDFILENAME: " << inUnfoldFileName << std::endl;
  
  const bool isSameInFile = isStrSame(inResponseFileName, inUnfoldFileName);
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<std::string> necessaryInFileParams = {"ISMC",
						    "NGAMMAPTBINS",
						    "GAMMAPTBINSLOW",
						    "GAMMAPTBINSHIGH",
						    "GAMMAPTBINSDOLOG",
						    "GAMMAPTBINSDOCUSTOM",
						    "NJTPTBINS",
						    "JTPTBINSLOW",
						    "JTPTBINSHIGH",
						    "JTPTBINSDOLOG",
						    "JTPTBINSDOCUSTOM",
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string gammaJtDPhiStr = "DPhi0";
  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  //  inResponseFileConfig_p->Print();
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  const bool isMCResponse = inResponseFileConfig_p->GetValue("ISMC", 1);
  //make sure response file given is MC

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!isMCResponse){
    std::cout << "Given responseFile \'" << inResponseFileName << "\' is not an MC file. return 1" << std::endl;
    return 1;
  }  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Grab first set of parameters defined in necessaryInFileParams
  const Int_t nMaxPtBins = 200;
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const bool isPP = inResponseFileConfig_p->GetValue("ISPP", 0);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nGammaPtBins = inResponseFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nGammaFinePtBins = inResponseFileConfig_p->GetValue("NGAMMAFINEPTBINS", 0);
  std::string gammaFinePtBinsStr = inResponseFileConfig_p->GetValue("GAMMAFINEPTBINS", "");
  std::vector<float> gammaFinePtBinsVect = strToVectF(gammaFinePtBinsStr);
  Double_t gammaFinePtBins[nMaxPtBins+1];

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int i = 0; i < gammaFinePtBinsVect.size(); ++i){
    gammaFinePtBins[i] = gammaFinePtBinsVect[i];
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    std::string gammaPtBinsStr = inResponseFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "gammaPtBinsDoCustom is true but gammaPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBinsTemp != nGammaPtBins){
      std::cout << "Number of gamma pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nJtPtBins = inResponseFileConfig_p->GetValue("NJTPTBINS", 0);
  const Float_t jtPtBinsLow = inResponseFileConfig_p->GetValue("JTPTBINSLOW", -1.0);
  const Float_t jtPtBinsHigh = inResponseFileConfig_p->GetValue("JTPTBINSHIGH", -1.0);
  const Bool_t jtPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOLOG", 1);
  const Bool_t jtPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOCUSTOM", 1);
  Double_t jtPtBins[nMaxPtBins+1];
  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else if(jtPtBinsDoCustom){
    std::string jtPtBinsStr = inResponseFileConfig_p->GetValue("JTPTBINSCUSTOM", "");
    if(jtPtBinsStr.size() == 0){
      std::cout << "jtPtBinsDoCustom is true but jtPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> jtPtBinsTemp = strToVectF(jtPtBinsStr);
    Int_t nJtPtBinsTemp = ((Int_t)jtPtBinsTemp.size()) - 1;
    if(nJtPtBinsTemp != nJtPtBins){
      std::cout << "Number of jt pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < jtPtBinsTemp.size(); ++gI){
      jtPtBins[gI] = jtPtBinsTemp[gI];
    }
  }
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
 
  //  const Float_t assocGenMinPt = inResponseFileConfig_p->GetValue("ASSOCGENMINPT", -1.0);
  
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
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  //Enforce additional requirements on config depending on first set
  if(isMC){
    if(!checkEnvForParams(inResponseFileConfig_p, mcParams)) return 1;
  }

  int nCentBins = 1;
  //  Float_t centBins[nMaxCentBins];
  std::vector<std::string> centBinsStr = {"PP"};
  std::vector<float> centBinsF = {};
  if(!isPP){
    if(!checkEnvForParams(inResponseFileConfig_p, pbpbParams)) return 1;

    centBinsStr.clear();
    centBinsF = strToVectF(inResponseFileConfig_p->GetValue("CENTBINS", ""));
    nCentBins = centBinsF.size()-1;
    std::vector<std::string> tempCentBinsStr = strToVect(inResponseFileConfig_p->GetValue("CENTBINS", ""));
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
    }

    nCentBins = centBinsStr.size();
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  //Enforce matching over specified parameters
  if(isSameInFile){
    if(!compEnvParams(inResponseFileConfig_p, inUnfoldFileConfig_p, paramsToCompare)) return 1;
  }

  //Gotta do this bit manually
  Int_t nVarBins = nXJBins;
  Double_t varBins[nMaxPtBins+1];
  if(isStrSame(varName, "Pt")){
    nVarBins = nJtPtBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = jtPtBins[vI];
    }
  }
  else{
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = xjBins[vI];
    }
  }
  
  //Construct matrices
  outFile_p->cd();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;					 
  

  for(Int_t cI = 0; cI < nCentBins; ++ cI){    
    photonPtReco_p[cI] = new TH1F(("photonPtReco_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);
    photonPtTruth_p[cI] = new TH1F(("photonPtTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);

    photonPtReco_HalfResponse_p[cI] = new TH1F(("photonPtReco_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);
    photonPtTruth_HalfResponse_p[cI] = new TH1F(("photonPtTruth_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), ";Truth Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);

    photonPtReco_HalfToUnfold_p[cI] = new TH1F(("photonPtReco_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);
    photonPtTruth_HalfToUnfold_p[cI] = new TH1F(("photonPtTruth_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), ";Truth Photon p_{T};Counts", nGammaFinePtBins, gammaFinePtBins);
    

    photonPtJetXJReco_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    photonPtJetXJReco_HalfResponse_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_HalfResponse_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJReco_HalfToUnfold_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_HalfToUnfold_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
   
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      photonPtReco_PURCORR_p[cI][eI] = (TH1F*)inUnfoldFile_p->Get((centBinsStr[cI] + "/photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_PURCORR_h").c_str());     

      photonPtJetXJReco_PURCORR_p[cI][eI] = (TH2F*)inUnfoldFile_p->Get((centBinsStr[cI] + "/photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str());     

      if(isMC){
	photonPtReco_TRUTH_p[cI][eI] = (TH1F*)inUnfoldFile_p->Get((centBinsStr[cI] + "/photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] +  "_TRUTH_h").c_str());
	photonPtJetXJReco_TRUTH_p[cI][eI] = (TH2F*)inUnfoldFile_p->Get((centBinsStr[cI] + "/photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str());
      }
      else{
	//	std::cout << "HISTOGRAMNAME: " << centBinsStr[cI] + "/photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] +  "_TRUTH_h" << std::endl; 
	photonPtReco_TRUTH_p[cI][eI] = (TH1F*)inResponseFile_p->Get((centBinsStr[cI] + "/photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] +  "_TRUTH_h").c_str());
	photonPtJetXJReco_TRUTH_p[cI][eI] = (TH2F*)inResponseFile_p->Get((centBinsStr[cI] + "/photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str());
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    photonPtReco_PURCORR_COMBINED_p[cI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), ";Reco. Photon p_{T};Counts (Purity Corrected)", nGammaFinePtBins, gammaFinePtBins);
    
    photonPtJetXJReco_PURCORR_COMBINED_p[cI] = new TH2F(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

    photonPtReco_TRUTH_COMBINED_p[cI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nGammaFinePtBins, gammaFinePtBins);   

    photonPtJetXJReco_TRUTH_COMBINED_p[cI] = new TH2F(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);   


    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t bIX = 0; bIX < photonPtReco_TRUTH_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      Float_t content = 0.0;
      Float_t error = 0.0;
    
      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	Float_t tempContent = photonPtReco_TRUTH_p[cI][eI]->GetBinContent(bIX+1);
	Float_t tempError = photonPtReco_TRUTH_p[cI][eI]->GetBinError(bIX+1);
	
	content += tempContent;
	error = TMath::Sqrt(error*error + tempError*tempError);
      }

      photonPtReco_TRUTH_COMBINED_p[cI]->SetBinContent(bIX+1, content);
      photonPtReco_TRUTH_COMBINED_p[cI]->SetBinError(bIX+1, error);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t bIX = 0; bIX < photonPtReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      Float_t content = 0.0;
      Float_t error = 0.0;

      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	Float_t tempContent = photonPtReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1);
	Float_t tempError = photonPtReco_PURCORR_p[cI][eI]->GetBinError(bIX+1);
	
	content += tempContent;
	error = TMath::Sqrt(error*error + tempError*tempError);
      }

      photonPtReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, content);
      photonPtReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, error);
    }
    
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t bIX = 0; bIX < photonPtJetXJReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      for(Int_t bIY = 0; bIY < photonPtJetXJReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins(); ++bIY){

	Float_t content = 0.0;
	Float_t error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetXJReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetXJReco_PURCORR_p[cI][eI]->GetBinError(bIX+1, bIY+1);

	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetXJReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetXJReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);	

	content = 0.0;
	error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetXJReco_TRUTH_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetXJReco_TRUTH_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	  
	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}
	  
	photonPtJetXJReco_TRUTH_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetXJReco_TRUTH_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);		
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    std::vector<TH1*> hists_p = {photonPtJetXJReco_p[cI], photonPtJetXJTruth_p[cI]};
    hists_p.push_back(photonPtJetXJReco_HalfResponse_p[cI]);
    hists_p.push_back(photonPtJetXJTruth_HalfResponse_p[cI]);
    hists_p.push_back(photonPtJetXJReco_HalfToUnfold_p[cI]);
    hists_p.push_back(photonPtJetXJTruth_HalfToUnfold_p[cI]);
    hists_p.push_back(photonPtJetXJReco_PURCORR_COMBINED_p[cI]);
    hists_p.push_back(photonPtJetXJReco_TRUTH_COMBINED_p[cI]);
    
    setSumW2(hists_p);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  //Grab the unfolding ttree for matrix filling
  const Int_t nJESSys = 18;
  const Int_t nJERSys = 9;
  const Int_t nJetSysAndNom = 1 + nJESSys + nJERSys;


  Float_t recoGammaPtPrev_ = -1.0;
  Float_t truthGammaPtPrev_ = -1.0;
  Float_t recoGammaPt_;
  Float_t truthGammaPt_;
  Float_t recoXJ_[nJetSysAndNom];
  Float_t recoJt1Pt_[nJetSysAndNom];
  Float_t recoJt2Pt_[nJetSysAndNom];
  Float_t truthXJ_;
  Float_t truthJt1Pt_;
  Float_t truthJt2Pt_;
  Float_t unfoldWeight_;
  Float_t unfoldCent_;
  
  TTree* unfoldTree_ResponseFile_p = (TTree*)inResponseFile_p->Get(("unfold" + varNameTree + "Tree_p").c_str());
  unfoldTree_ResponseFile_p->SetBranchAddress("recoGammaPt", &recoGammaPt_);
  unfoldTree_ResponseFile_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldTree_ResponseFile_p->SetBranchAddress(("reco" + varNameTree).c_str(), recoXJ_);
  unfoldTree_ResponseFile_p->SetBranchAddress(("truth" + varNameTree).c_str(), &truthXJ_);
  unfoldTree_ResponseFile_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  unfoldTree_ResponseFile_p->SetBranchAddress("unfoldCent", &unfoldCent_);

  //Some things have to be handled for a variable manually
  if(isStrSame("XJ", varName) || isStrSame("Pt", varName)){
    unfoldTree_ResponseFile_p->SetBranchAddress("recoJtPt", recoJt1Pt_);
    unfoldTree_ResponseFile_p->SetBranchAddress("truthJtPt", &truthJt1Pt_);
  }
  else if(isStrSame("XJJ", varName)){
    unfoldTree_ResponseFile_p->SetBranchAddress("recoJt1Pt", recoJt1Pt_);
    unfoldTree_ResponseFile_p->SetBranchAddress("truthJt1Pt", &truthJt1Pt_);
    unfoldTree_ResponseFile_p->SetBranchAddress("recoJt2Pt", recoJt2Pt_);
    unfoldTree_ResponseFile_p->SetBranchAddress("truthJt2Pt", &truthJt2Pt_);
  }
  
  //Construct our MC hists
  const ULong64_t nEntriesResponse = unfoldTree_ResponseFile_p->GetEntries();
 
  for(ULong64_t entry = 0; entry < nEntriesResponse; ++entry){
    unfoldTree_ResponseFile_p->GetEntry(entry);
    
    Int_t centPos = 0;
    if(!isPP) centPos = ghostPos(centBinsF, unfoldCent_, true);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;

    //Photon handling: Make sure its a unique photon
    if(TMath::Abs(recoGammaPt_ - recoGammaPtPrev_) > 0.01 && TMath::Abs(truthGammaPt_ - truthGammaPtPrev_) > 0.01){
      if(fillResponseHalf){
	photonPtReco_HalfResponse_p[centPos]->Fill(recoGammaPt_, unfoldWeight_);
	photonPtTruth_HalfResponse_p[centPos]->Fill(truthGammaPt_, unfoldWeight_);
      }
      else{
	photonPtReco_HalfToUnfold_p[centPos]->Fill(recoGammaPt_, unfoldWeight_);
	photonPtTruth_HalfToUnfold_p[centPos]->Fill(truthGammaPt_, unfoldWeight_);
      }      

      photonPtReco_p[centPos]->Fill(recoGammaPt_, unfoldWeight_);
      photonPtTruth_p[centPos]->Fill(truthGammaPt_, unfoldWeight_);

      recoGammaPtPrev_ = recoGammaPt_;
      truthGammaPtPrev_ = truthGammaPt_;
    }

    if(isStrSame(varName, "Pt")){
      recoXJ_[0] = recoJt1Pt_[0];
      truthXJ_ = truthJt1Pt_;   
    }    
    
    if(recoXJ_[0] < 0.0) continue; // Misses are only passed to the RooUnfoldResponse Object
    
    if(truthXJ_ < varBins[0]){//underflow
      truthXJ_ = (varBins[0] + varBins[1])/2.;
    }
    else if(truthXJ_ > varBins[nVarBins]){//overflow
      truthXJ_ = (varBins[nVarBins] + varBins[nVarBins-1])/2.;
    }
        
    if(fillResponseHalf){
      photonPtJetXJReco_HalfResponse_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, unfoldWeight_);
      photonPtJetXJTruth_HalfResponse_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    else{
      photonPtJetXJReco_HalfToUnfold_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, unfoldWeight_);
      photonPtJetXJTruth_HalfToUnfold_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    
    photonPtJetXJReco_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, unfoldWeight_);
    photonPtJetXJTruth_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
  }

  
  TH2F* reweightingHist_p[nMaxCentBins];
  if(!isSameInFile && doReweight){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      reweightingHist_p[cI] = new TH2F(("reweightingHist_" + centBinsStr[cI] + "_h").c_str(), (";Reco. Jet " + varNameStyle + ";Gamma p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

      for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
	for(Int_t bIY = 0; bIY < nGammaPtBins; ++bIY){
	
	  Float_t numeratorData = photonPtJetXJReco_PURCORR_COMBINED_p[cI]->GetBinContent(bIX+1, bIY+1);
	  Float_t denominatorMC = photonPtJetXJReco_p[cI]->GetBinContent(bIX+1, bIY+1);

	  Float_t ratioVal = 0.0;
	  if(denominatorMC > TMath::Power(10, -50)) ratioVal = numeratorData/denominatorMC;
	  
	  reweightingHist_p[cI]->SetBinContent(bIX+1, bIY+1, ratioVal);
	}
      }
    }

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      delete photonPtJetXJReco_p[cI];
      delete photonPtJetXJTruth_p[cI];


      photonPtJetXJReco_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
      photonPtJetXJTruth_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    }
    
    for(ULong64_t entry = 0; entry < nEntriesResponse; ++entry){
      unfoldTree_ResponseFile_p->GetEntry(entry);
    
      Int_t centPos = 0;
      if(!isPP) centPos = ghostPos(centBinsF, unfoldCent_, true);

      if(isStrSame(varName, "Pt")){
	recoXJ_[0] = recoJt1Pt_[0];
	truthXJ_ = truthJt1Pt_;   
      }

      if(recoXJ_[0] < 0.0) continue; // Misses are only passed to the RooUnfoldResponse Object
    
      if(truthXJ_ < varBins[0]){//underflow
	truthXJ_ = (varBins[0] + varBins[1])/2.;
      }
      else if(truthXJ_ > varBins[nVarBins]){//overflow
	truthXJ_ = (varBins[nVarBins] + varBins[nVarBins-1])/2.;
      }

      Int_t gammaPos = ghostPos(nGammaPtBins, gammaPtBins, recoGammaPt_);
      Int_t jtPos = ghostPos(nVarBins, varBins, recoXJ_[0]);

      Float_t reweightVal = reweightingHist_p[centPos]->GetBinContent(jtPos+1, gammaPos+1);
      
      photonPtJetXJReco_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, unfoldWeight_*reweightVal);
      photonPtJetXJTruth_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_*reweightVal);
    }
  }
  
  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    rooResGammaJet_Half_p[cI] = new RooUnfoldResponse(photonPtJetXJReco_HalfResponse_p[cI], photonPtJetXJTruth_HalfResponse_p[cI], ("rooResGammaJet_Half_" + centBinsStr[cI]).c_str(), "");
    rooResGammaJet_p[cI] = new RooUnfoldResponse(photonPtJetXJReco_p[cI], photonPtJetXJTruth_p[cI], ("rooResGammaJet_" + centBinsStr[cI]).c_str(), "");

    rooResGamma_Half_p[cI] = new RooUnfoldResponse(photonPtReco_HalfResponse_p[cI], photonPtTruth_HalfResponse_p[cI], ("rooResGamma_Half_" + centBinsStr[cI]).c_str(), "");
    rooResGamma_p[cI] = new RooUnfoldResponse(photonPtReco_p[cI], photonPtTruth_p[cI], ("rooResGamma_" + centBinsStr[cI]).c_str(), "");
  }
  
  delete randGen_p;
  randGen_p = new TRandom3(randSeed); //Resetting to get the correct response matrix fills

  recoGammaPtPrev_ = -1.0;
  truthGammaPtPrev_ = -1.0;
  
  for(ULong64_t entry = 0; entry < nEntriesResponse; ++entry){
    unfoldTree_ResponseFile_p->GetEntry(entry);
    
    Int_t centPos = 0;
    if(!isPP) centPos = ghostPos(centBinsF, unfoldCent_, true);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;
    //Photon handling: Make sure its a unique photon
    if(TMath::Abs(recoGammaPt_ - recoGammaPtPrev_) > 0.01 && TMath::Abs(truthGammaPt_ - truthGammaPtPrev_) > 0.01){
      if(fillResponseHalf){
	if(recoGammaPt_ >= gammaPtBins[0] && recoGammaPt_ <= gammaPtBins[nGammaPtBins]) rooResGamma_Half_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_);
	else rooResGamma_Half_p[centPos]->Miss(truthGammaPt_, unfoldWeight_);
      }

      if(recoGammaPt_ >= gammaPtBins[0] && recoGammaPt_ <= gammaPtBins[nGammaPtBins]) rooResGamma_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_);
      else rooResGamma_p[centPos]->Miss(truthGammaPt_, unfoldWeight_);

      recoGammaPtPrev_ = recoGammaPt_;
      truthGammaPtPrev_ = truthGammaPt_;
    }

    
    if(isStrSame(varName, "Pt")){
      recoXJ_[0] = recoJt1Pt_[0];
      truthXJ_ = truthJt1Pt_;
    }

    if(fillResponseHalf){
      if(recoXJ_[0] >= varBins[0] && recoXJ_[0] < varBins[nVarBins]) rooResGammaJet_Half_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);
      else rooResGammaJet_Half_p[centPos]->Miss(truthXJ_, truthGammaPt_, unfoldWeight_);
    }

    if(!isSameInFile && doReweight){
      Int_t gammaPos = ghostPos(nGammaPtBins, gammaPtBins, recoGammaPt_);
      Int_t jtPos = ghostPos(nVarBins, varBins, recoXJ_[0]);

      Float_t reweightVal = reweightingHist_p[centPos]->GetBinContent(jtPos+1, gammaPos+1);
      unfoldWeight_ *= reweightVal;
    }
    
    if(recoXJ_[0] >= varBins[0] && recoXJ_[0] < varBins[nVarBins]) rooResGammaJet_p[centPos]->Fill(recoXJ_[0], recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);
    else rooResGammaJet_p[centPos]->Miss(truthXJ_, truthGammaPt_, unfoldWeight_);
  }
  
  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    photonPtTruth_HalfToUnfold_p[cI]->Write(("photonPtTruth_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), TObject::kOverwrite);
    photonPtReco_HalfToUnfold_p[cI]->Write(("photonPtReco_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), TObject::kOverwrite);

    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);
    
    for(int gI = 0; gI < nGammaPtBins; ++gI){
      TH1D* xjTruthProjection_p = new TH1D((varNameLower + "Truth_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
      TH1D* xjRecoProjection_p = new TH1D((varNameLower + "Reco_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
      TH1D* xjRecoProjection_PURCORR_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
      TH1D* xjRecoProjection_TRUTH_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);

      for(int xI = 0; xI < nVarBins; ++xI){
	xjTruthProjection_p->SetBinContent(xI+1, photonPtJetXJTruth_HalfToUnfold_p[cI]->GetBinContent(xI+1, gI+1));
	xjTruthProjection_p->SetBinError(xI+1, photonPtJetXJTruth_HalfToUnfold_p[cI]->GetBinError(xI+1, gI+1));      

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

	xjRecoProjection_p->SetBinContent(xI+1, photonPtJetXJReco_HalfToUnfold_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_p->SetBinError(xI+1, photonPtJetXJReco_HalfToUnfold_p[cI]->GetBinError(xI+1, gI+1));

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

	xjRecoProjection_PURCORR_COMBINED_p->SetBinContent(xI+1, photonPtJetXJReco_PURCORR_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_PURCORR_COMBINED_p->SetBinError(xI+1, photonPtJetXJReco_PURCORR_COMBINED_p[cI]->GetBinError(xI+1, gI+1));

	xjRecoProjection_TRUTH_COMBINED_p->SetBinContent(xI+1, photonPtJetXJReco_TRUTH_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_TRUTH_COMBINED_p->SetBinError(xI+1, photonPtJetXJReco_TRUTH_COMBINED_p[cI]->GetBinError(xI+1, gI+1));
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
      
      xjTruthProjection_p->Write("", TObject::kOverwrite);
      delete xjTruthProjection_p;

      xjRecoProjection_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_p;

      xjRecoProjection_PURCORR_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_PURCORR_COMBINED_p;

      xjRecoProjection_TRUTH_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_TRUTH_COMBINED_p;
    }

    //Half unfold test gamma only
    for(int i = 1; i < nIter; ++i){
      RooUnfoldBayes* rooBayes_Half_p = new RooUnfoldBayes(rooResGamma_Half_p[cI], photonPtReco_HalfToUnfold_p[cI], i);
      TH1D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH1D*)rooBayes_Half_p->Hreco()->Clone(("photonPtReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str());
      else if(unfoldErrType == 1){
	rooBayes_Half_p->SetNToys(nToys);
	unfolded_p = (TH1D*)rooBayes_Half_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str());
      }
      TH1D* refolded_p = (TH1D*)rooResGamma_Half_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_Refolded_h").c_str());

      unfolded_p->Write(("photonPtReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_Refolded_h").c_str(), TObject::kOverwrite);

      delete refolded_p;
      delete unfolded_p;
      delete rooBayes_Half_p;
    }
    
    //Half unfold test gamma+jet
    for(int i = 1; i < nIter; ++i){
      RooUnfoldBayes* rooBayes_Half_p = new RooUnfoldBayes(rooResGammaJet_Half_p[cI], photonPtJetXJReco_HalfToUnfold_p[cI], i);  
      TH2D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH2D*)rooBayes_Half_p->Hreco()->Clone(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      else if(unfoldErrType == 1){
	rooBayes_Half_p->SetNToys(nToys);
	unfolded_p = (TH2D*)rooBayes_Half_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      }
      TH2D* refolded_p = (TH2D*)rooResGammaJet_Half_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_Refolded_h").c_str());
      
      unfolded_p->Write(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_Refolded_h").c_str(), TObject::kOverwrite);

      for(int gI = 0; gI < nGammaPtBins; ++gI){
	TH1D* xjUnfoldedProjection_p = new TH1D((varNameLower + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
	TH1D* xjRefoldedProjection_p = new TH1D((varNameLower + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_Refolded_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
	
	for(int xI = 0; xI < nVarBins; ++xI){
	  xjUnfoldedProjection_p->SetBinContent(xI+1, unfolded_p->GetBinContent(xI+1, gI+1));
	  xjUnfoldedProjection_p->SetBinError(xI+1, unfolded_p->GetBinError(xI+1, gI+1));

	  xjRefoldedProjection_p->SetBinContent(xI+1, refolded_p->GetBinContent(xI+1, gI+1));
	  xjRefoldedProjection_p->SetBinError(xI+1, refolded_p->GetBinError(xI+1, gI+1));
	}
	
	outFile_p->cd();
	
	xjUnfoldedProjection_p->Write("", TObject::kOverwrite);
	delete xjUnfoldedProjection_p;

	xjRefoldedProjection_p->Write("", TObject::kOverwrite);
	delete xjRefoldedProjection_p;
      }
      
      delete rooBayes_Half_p;
    }

    //Write the full histograms at reco level to file
    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    photonPtJetXJReco_PURCORR_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtJetXJReco_TRUTH_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    //Full unfold gamma only
    for(int i = 1; i < nIter; ++ i){
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGamma_p[cI], photonPtReco_PURCORR_COMBINED_p[cI], i);

      TH1D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH1D*)rooBayes_p->Hreco()->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      else if(unfoldErrType == 1){
        rooBayes_p->SetNToys(nToys);
        unfolded_p = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      }
      TH1D* refolded_p = (TH1D*)rooResGamma_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str());

      unfolded_p->Write(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str(), TObject::kOverwrite);

      delete unfolded_p;
      delete refolded_p;
      delete rooBayes_p;
    }

    //Full unfold gamma+jet, should be instant when the unfold file is also mc    
    for(int i = 1; i < nIter; ++ i){
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGammaJet_p[cI], photonPtJetXJReco_PURCORR_COMBINED_p[cI], i);  
      
      TH2D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH2D*)rooBayes_p->Hreco()->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      else if(unfoldErrType == 1){
	rooBayes_p->SetNToys(nToys);
	unfolded_p = (TH2D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      }
      TH2D* refolded_p = (TH2D*)rooResGammaJet_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str());

      unfolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str(), TObject::kOverwrite);

      for(int gI = 0; gI < nGammaPtBins; ++gI){
	TH1D* xjUnfoldedProjection_p = new TH1D((varNameLower + "Reco_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
	TH1D* xjRefoldedProjection_p = new TH1D((varNameLower + "Reco_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
	
	for(int xI = 0; xI < nVarBins; ++xI){
	  xjUnfoldedProjection_p->SetBinContent(xI+1, unfolded_p->GetBinContent(xI+1, gI+1));
	  xjUnfoldedProjection_p->SetBinError(xI+1, unfolded_p->GetBinError(xI+1, gI+1));

	  xjRefoldedProjection_p->SetBinContent(xI+1, refolded_p->GetBinContent(xI+1, gI+1));
	  xjRefoldedProjection_p->SetBinError(xI+1, refolded_p->GetBinError(xI+1, gI+1));
	}
	
	outFile_p->cd();
	
	xjUnfoldedProjection_p->Write("", TObject::kOverwrite);
	delete xjUnfoldedProjection_p;
	
	xjRefoldedProjection_p->Write("", TObject::kOverwrite);
	delete xjRefoldedProjection_p;
      }
      
      delete rooBayes_p;
    }

    
    photonPtJetXJReco_HalfResponse_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_HalfResponse_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJReco_HalfToUnfold_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_HalfToUnfold_p[cI]->Write("", TObject::kOverwrite);
    
    photonPtJetXJReco_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_p[cI]->Write("", TObject::kOverwrite);

    photonPtJetXJReco_PURCORR_COMBINED_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJReco_TRUTH_COMBINED_p[cI]->Write("", TObject::kOverwrite);
    
    if(!isSameInFile && doReweight){
      reweightingHist_p[cI]->Write("", TObject::kOverwrite);
    }
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;  

    delete photonPtReco_HalfResponse_p[cI];
    delete photonPtTruth_HalfResponse_p[cI];
    delete photonPtReco_HalfToUnfold_p[cI];
    delete photonPtTruth_HalfToUnfold_p[cI];

    delete photonPtReco_p[cI];
    delete photonPtTruth_p[cI];
    
    
    delete photonPtJetXJReco_HalfResponse_p[cI];
    delete photonPtJetXJTruth_HalfResponse_p[cI];
    delete photonPtJetXJReco_HalfToUnfold_p[cI];
    delete photonPtJetXJTruth_HalfToUnfold_p[cI];
    
    delete photonPtJetXJReco_p[cI];
    delete photonPtJetXJTruth_p[cI];

    delete photonPtJetXJReco_PURCORR_COMBINED_p[cI];
    delete photonPtJetXJReco_TRUTH_COMBINED_p[cI];
    
    delete rooResGammaJet_Half_p[cI];
    delete rooResGammaJet_p[cI];

    if(!isSameInFile && doReweight){
      delete reweightingHist_p[cI];
    }
  }
  
  inUnfoldFileConfig_p->SetValue("VARNAME", varName.c_str());
  inUnfoldFileConfig_p->SetValue("NITER", nIter);

  inUnfoldFileConfig_p->Write("config", TObject::kOverwrite);
  inUnfoldFileLabels_p->Write("label", TObject::kOverwrite);

  inResponseFile_p->Close();
  delete inResponseFile_p;  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  if(!isSameInFile){
    inUnfoldFile_p->Close();
    delete inUnfoldFile_p;
  }

  outFile_p->Close();
  delete outFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
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
