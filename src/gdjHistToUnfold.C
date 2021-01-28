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
					      "NITER"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
 
  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string varName = config_p->GetValue("VARNAME", "");
  std::string varNameLower = returnAllLowercaseString(varName);

  const int nIter = config_p->GetValue("NITER", -1);
  const int nIterMax = 100;
  if(nIter > nIterMax){
    std::cout << "Requested number of iterations \'" << nIter << "\' exceeds maximum allowed \'" << nIterMax << "\'" << std::endl;
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(outFileName.find(".root") != std::string::npos){
    outFileName.replace(outFileName.rfind(".root"), 5, "");
    outFileName = outFileName + "_" + dateStr + ".root";
  }

  if(outFileName.find("/") == std::string::npos){
    outFileName = "output/" + dateStr + "/" + outFileName;
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  std::vector<std::string> validUnfoldVar = {"XJ", "XJJ"};
  std::vector<std::string> validUnfoldVarStyle = {"x_{J}", "#vec{x}_{JJ#gamma}"};
  const int vectPos = vectContainsStrPos(varName, &validUnfoldVar);
  if(vectPos < 0){
    std::cout << "GIVEN VARNAME \'" << varName << "\' is invalid. Please pick from:" << std::endl;
    for(unsigned int i = 0; i < validUnfoldVar.size(); ++i){
      std::cout << " \'" << validUnfoldVar[i] << "\'" << std::endl;
    }
    std::cout << "return 1" << std::endl;
    return 1;
  }
  std::string varNameStyle = validUnfoldVarStyle[vectPos];
  
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;
  
  const Int_t nMaxCentBins = 10;
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  
  TH2F* photonPtJetXJReco_HalfResponse_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_HalfResponse_p[nMaxCentBins];
  TH2F* photonPtJetXJReco_HalfToUnfold_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_HalfToUnfold_p[nMaxCentBins];

  TH2F* photonPtJetXJReco_p[nMaxCentBins];
  TH2F* photonPtJetXJTruth_p[nMaxCentBins];  

  TH2F* photonPtJetXJReco_SUBData_p[nMaxCentBins];

  RooUnfoldResponse* rooRes_Half_p[nMaxCentBins];
  RooUnfoldResponse* rooRes_p[nMaxCentBins];
  
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");  
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  TFile* inUnfoldFile_p = nullptr;
  TEnv* inUnfoldFileConfig_p = nullptr;
  TEnv* inUnfoldFileLabels_p = nullptr;

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
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  //Grab first set of parameters defined in necessaryInFileParams
  const Int_t nMaxPtBins = 200;
  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  const bool isPP = inResponseFileConfig_p->GetValue("ISPP", 0);

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

  //Construct matrices
  outFile_p->cd();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++ cI){    
    photonPtJetXJReco_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);

    photonPtJetXJReco_HalfResponse_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_HalfResponse_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_HalfResponse_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJReco_HalfToUnfold_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);
    photonPtJetXJTruth_HalfToUnfold_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_HalfToUnfold_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nXJBins, xjBins, nGammaPtBins, gammaPtBins);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    photonPtJetXJReco_SUBData_p[cI] = (TH2F*)inUnfoldFile_p->Get((centBinsStr[cI] + "/photonPtJtXJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str());

    //    std::cout << centBinsStr[cI] + "/photonPtJtXJReco_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_SUB_h" << std::endl;
    
    std::vector<TH1*> hists_p = {photonPtJetXJReco_p[cI], photonPtJetXJTruth_p[cI]};
    hists_p.push_back(photonPtJetXJReco_HalfResponse_p[cI]);
    hists_p.push_back(photonPtJetXJTruth_HalfResponse_p[cI]);
    hists_p.push_back(photonPtJetXJReco_HalfToUnfold_p[cI]);
    hists_p.push_back(photonPtJetXJTruth_HalfToUnfold_p[cI]);
    hists_p.push_back(photonPtJetXJReco_SUBData_p[cI]);
    
    setSumW2(hists_p);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  //Grab the unfolding ttree for matrix filling
  Float_t recoGammaPt_, truthGammaPt_;
  Float_t recoXJ_, truthXJ_;
  Float_t unfoldWeight_;
  Float_t unfoldCent_;
  
  TTree* unfoldTree_ResponseFile_p = (TTree*)inResponseFile_p->Get(("unfold" + varName + "Tree_p").c_str());
  unfoldTree_ResponseFile_p->SetBranchAddress("recoGammaPt", &recoGammaPt_);
  unfoldTree_ResponseFile_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldTree_ResponseFile_p->SetBranchAddress(("reco" + varName).c_str(), &recoXJ_);
  unfoldTree_ResponseFile_p->SetBranchAddress(("truth" + varName).c_str(), &truthXJ_);
  unfoldTree_ResponseFile_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  unfoldTree_ResponseFile_p->SetBranchAddress("unfoldCent", &unfoldCent_);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  //Construct our MC hists
  const ULong64_t nEntriesResponse = unfoldTree_ResponseFile_p->GetEntries();
 
  for(ULong64_t entry = 0; entry < nEntriesResponse; ++entry){
    unfoldTree_ResponseFile_p->GetEntry(entry);

    Int_t centPos = 0;
    if(!isPP) centPos = ghostPos(centBinsF, unfoldCent_, true);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;
    if(fillResponseHalf){
      photonPtJetXJReco_HalfResponse_p[centPos]->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
      photonPtJetXJTruth_HalfResponse_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    else{
      photonPtJetXJReco_HalfToUnfold_p[centPos]->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
      photonPtJetXJTruth_HalfToUnfold_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
    }
    
    photonPtJetXJReco_p[centPos]->Fill(recoXJ_, recoGammaPt_, unfoldWeight_);
    photonPtJetXJTruth_p[centPos]->Fill(truthXJ_, truthGammaPt_, unfoldWeight_);
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    rooRes_Half_p[cI] = new RooUnfoldResponse(photonPtJetXJReco_HalfResponse_p[cI], photonPtJetXJTruth_HalfResponse_p[cI], ("rooRes_Half_" + centBinsStr[cI]).c_str(), "");
    rooRes_p[cI] = new RooUnfoldResponse(photonPtJetXJReco_p[cI], photonPtJetXJTruth_p[cI], ("rooRes_" + centBinsStr[cI]).c_str(), "");
  }
  
  delete randGen_p;
  randGen_p = new TRandom3(randSeed); //Resetting to get the correct response matrix fills

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

  for(ULong64_t entry = 0; entry < nEntriesResponse; ++entry){
    unfoldTree_ResponseFile_p->GetEntry(entry);
    
    Int_t centPos = 0;
    if(!isPP) centPos = ghostPos(centBinsF, unfoldCent_, true);

    bool fillResponseHalf = randGen_p->Uniform() < 0.5;
    if(fillResponseHalf) rooRes_Half_p[centPos]->Fill(recoXJ_, recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);

    rooRes_p[centPos]->Fill(recoXJ_, recoGammaPt_, truthXJ_, truthGammaPt_, unfoldWeight_);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  outFile_p->cd();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(int gI = 0; gI < nGammaPtBins; ++gI){
      TH1D* xjTruthProjection_p = new TH1D((varNameLower + "Truth_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nXJBins, xjBins);
      TH1D* xjRecoProjection_p = new TH1D((varNameLower + "Reco_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nXJBins, xjBins);
      TH1D* xjRecoProjection_SUBData_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_SUBData_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nXJBins, xjBins);

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

      for(int xI = 0; xI < nXJBins; ++xI){
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

	xjTruthProjection_p->SetBinContent(xI+1, photonPtJetXJTruth_HalfToUnfold_p[cI]->GetBinContent(xI+1, gI+1));
	xjTruthProjection_p->SetBinError(xI+1, photonPtJetXJTruth_HalfToUnfold_p[cI]->GetBinError(xI+1, gI+1));      

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

	xjRecoProjection_p->SetBinContent(xI+1, photonPtJetXJReco_HalfToUnfold_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_p->SetBinError(xI+1, photonPtJetXJReco_HalfToUnfold_p[cI]->GetBinError(xI+1, gI+1));

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

	xjRecoProjection_SUBData_p->SetBinContent(xI+1, photonPtJetXJReco_SUBData_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_SUBData_p->SetBinError(xI+1, photonPtJetXJReco_SUBData_p[cI]->GetBinError(xI+1, gI+1));

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
      
      xjTruthProjection_p->Write("", TObject::kOverwrite);
      delete xjTruthProjection_p;

      xjRecoProjection_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_p;

      xjRecoProjection_SUBData_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_SUBData_p;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
    //Half unfold test
    for(int i = 1; i < nIter; ++ i){
      RooUnfoldBayes* rooBayes_Half_p = new RooUnfoldBayes(rooRes_Half_p[cI], photonPtJetXJReco_HalfToUnfold_p[cI], i);  

      TH2D* unfolded_p = (TH2D*)rooBayes_Half_p->Hreco();//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      
      unfolded_p->Write(("photonPtJetXJReco_HalfToUnfold_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_h").c_str(), TObject::kOverwrite);

      for(int gI = 0; gI < nGammaPtBins; ++gI){
	TH1D* xjUnfoldedProjection_p = new TH1D((varNameLower + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nXJBins, xjBins);
	
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

    //Full unfold, should be instant when the unfold file is also mc    
    for(int i = 1; i < nIter; ++ i){
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooRes_p[cI], photonPtJetXJReco_SUBData_p[cI], i);  

      TH2D* unfolded_p = (TH2D*)rooBayes_p->Hreco();//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      
      unfolded_p->Write(("photonPtJetXJReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_SUBData_h").c_str(), TObject::kOverwrite);

      for(int gI = 0; gI < nGammaPtBins; ++gI){
	TH1D* xjUnfoldedProjection_p = new TH1D((varNameLower + "Reco_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_SUBData_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nXJBins, xjBins);
	
	for(int xI = 0; xI < nXJBins; ++xI){
	  xjUnfoldedProjection_p->SetBinContent(xI+1, unfolded_p->GetBinContent(xI+1, gI+1));
	  xjUnfoldedProjection_p->SetBinError(xI+1, unfolded_p->GetBinError(xI+1, gI+1));
	}
	
	outFile_p->cd();
	
	xjUnfoldedProjection_p->Write("", TObject::kOverwrite);
	delete xjUnfoldedProjection_p;
      }
      
      delete rooBayes_p;
    }
    
    photonPtJetXJReco_HalfResponse_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_HalfResponse_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJReco_HalfToUnfold_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_HalfToUnfold_p[cI]->Write("", TObject::kOverwrite);
    
    photonPtJetXJReco_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetXJTruth_p[cI]->Write("", TObject::kOverwrite);

    photonPtJetXJReco_SUBData_p[cI]->Write("", TObject::kOverwrite);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;  

    delete photonPtJetXJReco_HalfResponse_p[cI];
    delete photonPtJetXJTruth_HalfResponse_p[cI];
    delete photonPtJetXJReco_HalfToUnfold_p[cI];
    delete photonPtJetXJTruth_HalfToUnfold_p[cI];
    
    delete photonPtJetXJReco_p[cI];
    delete photonPtJetXJTruth_p[cI];

    delete photonPtJetXJReco_SUBData_p[cI];
    
    delete rooRes_Half_p[cI];
    delete rooRes_p[cI];
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
