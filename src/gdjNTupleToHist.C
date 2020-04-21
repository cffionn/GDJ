//Author: Chris McGinn (2020.02.19)
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

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"

void fillTH1(TH1F* inHist_p, Float_t fillVal, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal, weight);
  }
  return;
}

void fillTH2(TH2F* inHist_p, Float_t fillVal1, Float_t fillVal2, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal1, fillVal2);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal1, fillVal2, weight);
  }
  return;
}

int gdjNTupleToHist(std::string inConfigFileName)
{
  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101
  TRandom3* randGen_p = new TRandom3(randSeed);

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  configParser config(inConfigFileName);

  std::vector<std::string> necessaryParams = {"INFILENAME",
                                              "OUTFILENAME",
                                              "CENTFILENAME",
					      "MIXFILENAME",
					      "DOMIX",
                                              "ISPP",
					      "ISMC",
					      "CENTBINS",
					      "NGAMMAPTBINS",
					      "GAMMAPTBINSLOW",
					      "GAMMAPTBINSHIGH",
					      "GAMMAPTBINSDOLOG",
					      "NETABINS",
					      "ETABINSLOW",
					      "ETABINSHIGH",
					      "ETABINSDOABS",
					      "NPHIBINS",
					      "NGAMMAPTBINSSUB",
					      "GAMMAPTBINSSUBLOW",
					      "GAMMAPTBINSSUBHIGH",
					      "GAMMAPTBINSDOLOG",
					      "NETABINSSUB",
					      "ETABINSSUBLOW",
					      "ETABINSSUBHIGH",
					      "ETABINSSUBDOABS",
					      "NJTPTBINS",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "JTPTBINSDOLOG",
					      "NDPHIBINS",
					      "GAMMAJTDPHI",  
					      "GAMMAJTDPHI",  
					      "NXJBINS",
					      "XJBINSLOW",
					      "XJBINSHIGH"};

  std::vector<std::string> mixParams = {"DOMIXCENT",
					"NMIXCENTBINS",
					"MIXCENTBINSLOW",
					"MIXCENTBINSHIGH",
					"DOMIXPSI2",
					"NMIXPSI2BINS",
					"MIXPSI2BINSLOW",
					"MIXPSI2BINSHIGH",
					"DOMIXVZ",
					"NMIXVZBINS",
					"MIXVZBINSLOW",
					"MIXVZBINSHIGH"};
				       
  
  if(!config.ContainsParamSet(necessaryParams)) return 1;  
 
  std::string inROOTFileName = config.GetConfigVal("INFILENAME");
  std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
  std::string outFileName = config.GetConfigVal("OUTFILENAME");
  std::string inGRLFileName = config.GetConfigVal("GRLFILENAME");
  std::string inMixFileName = config.GetConfigVal("MIXFILENAME");
  const bool doMix = std::stoi(config.GetConfigVal("DOMIX"));

  if(doMix && !config.ContainsParamSet(mixParams)) return 1;
  
  //Mixing categories, temp hardcoding
  //Centrality, percent level

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxMixBins = 200;
  bool doMixCent = false;
  Int_t nMixCentBins = -1;
  Float_t mixCentBinsLow = 0;
  Float_t mixCentBinsHigh = 100;
  Double_t mixCentBins[nMaxMixBins+1];

  bool doMixPsi2 = false;
  Int_t nMixPsi2Bins = -1;
  Float_t mixPsi2BinsLow = -TMath::Pi()/2.;
  Float_t mixPsi2BinsHigh = TMath::Pi()/2.;
  Double_t mixPsi2Bins[nMaxMixBins+1];

  bool doMixVz = false;
  Int_t nMixVzBins = -1;
  Float_t mixVzBinsLow = -15.;
  Float_t mixVzBinsHigh = 15.;
  Double_t mixVzBins[nMaxMixBins+1];

  std::vector<std::vector<unsigned long long> > mixVect;
  std::vector<std::vector<unsigned long long> > keyVect;
  
  if(doMix){
    if(!check.checkFileExt(inMixFileName, ".root")) return 1;

    doMixCent = (bool)std::stoi(config.GetConfigVal("DOMIXCENT"));

    if(doMixCent){      
      nMixCentBins = std::stoi(config.GetConfigVal("NMIXCENTBINS"));
      if(nMixCentBins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixCentBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixCentBinsLow = std::stof(config.GetConfigVal("MIXCENTBINSLOW"));
      mixCentBinsHigh = std::stof(config.GetConfigVal("MIXCENTBINSHIGH"));    
      getLinBins(mixCentBinsLow, mixCentBinsHigh, nMixCentBins, mixCentBins);

      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixCentBins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);
	keyVect.push_back({(unsigned long long)mI});
      }
    }

    doMixPsi2 = (bool)std::stoi(config.GetConfigVal("DOMIXPSI2"));

    if(doMixPsi2){
      nMixPsi2Bins = std::stoi(config.GetConfigVal("NMIXPSI2BINS"));
      if(nMixPsi2Bins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixPsi2Bins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixPsi2BinsLow = std::stof(config.GetConfigVal("MIXPSI2BINSLOW"));
      mixPsi2BinsHigh = std::stof(config.GetConfigVal("MIXPSI2BINSHIGH"));    
      getLinBins(mixPsi2BinsLow, mixPsi2BinsHigh, nMixPsi2Bins, mixPsi2Bins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixPsi2Bins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);

	if(keyVect.size() != 0){
	  for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
	    std::vector<unsigned long long> tempVect = keyVect[vI];
	    tempVect.push_back((unsigned long long)mI);
	    tempKeyVect.push_back(tempVect);
	  }
	}
	else tempKeyVect.push_back({(unsigned long long)mI});
      }

      keyVect = tempKeyVect;
    }

    doMixVz = (bool)std::stoi(config.GetConfigVal("DOMIXVZ"));

    if(doMixVz){
      nMixVzBins = std::stoi(config.GetConfigVal("NMIXVZBINS"));
      if(nMixVzBins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixVzBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixVzBinsLow = std::stof(config.GetConfigVal("MIXVZBINSLOW"));
      mixVzBinsHigh = std::stof(config.GetConfigVal("MIXVZBINSHIGH"));    
      getLinBins(mixVzBinsLow, mixVzBinsHigh, nMixVzBins, mixVzBins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixVzBins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);

	if(keyVect.size() != 0){
	  for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
	    std::vector<unsigned long long> tempVect = keyVect[vI];
	    tempVect.push_back((unsigned long long)mI);
	    tempKeyVect.push_back(tempVect);
	  }
	}
	else tempKeyVect.push_back({(unsigned long long)mI});
      }

      keyVect = tempKeyVect;
    }

  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  keyHandler runLumiKey("runLumiHandler");//for lumi estimation
  runLumiKey.Init({1000000, 10000});//runnumbers, then lumi
  
  keyHandler keyBoy("mixingHandler");//For Mixing
  std::map<unsigned long long, std::vector<std::vector<TLorentzVector> > > mixingMap;
  if(doMix){
    std::vector<unsigned long long> sizes;
    for(unsigned int vI = 0; vI < mixVect.size(); ++vI){
      sizes.push_back(mixVect[vI].size()+1);
    }

    //    return 0;
    
    keyBoy.Init(sizes);    

    for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
      unsigned long long key = keyBoy.GetKey(keyVect[vI]);
      mixingMap[key] = {};
    }
  }  
  
  const bool isMC = std::stoi(config.GetConfigVal("ISMC"));

  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
  if(!isMC && !check.checkFileExt(inGRLFileName, "xml")) return 1; // GRL File xml
  if(doMix && !check.checkFileExt(inMixFileName, "root")) return 1; // Check mixed file exists if mixing is requested
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();
  
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const bool isPP = std::stoi(config.GetConfigVal("ISPP"));
  const Int_t nMaxSubBins = 10;
  const Int_t nMaxCentBins = 10;
  Int_t nCentBins = 1;

 std::string systStr = "PP";
  if(!isPP) systStr = "PbPb";
  
  std::vector<int> centBins;
  std::vector<std::string> centBinsStr = {systStr};
  std::map<std::string, std::string> binsToLabelStr;
  if(!isPP){
    centBins = strToVectI(config.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;

    if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

    centBinsStr.clear();
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));

      binsToLabelStr[centBinsStr[cI]] = std::to_string(centBins[cI]) + "-" + std::to_string(centBins[cI+1]) + "%";
    }
  }
  else binsToLabelStr[centBinsStr[0]] = "pp";

  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = std::stoi(config.GetConfigVal("NGAMMAPTBINS"));
  if(!goodBinning(inConfigFileName, nMaxPtBins, nGammaPtBins, "NGAMMAPTBINS")) return 1;
  const Float_t gammaPtBinsLow = std::stof(config.GetConfigVal("GAMMAPTBINSLOW"));
  const Float_t gammaPtBinsHigh = std::stof(config.GetConfigVal("GAMMAPTBINSHIGH"));
  const Bool_t gammaPtBinsDoLog = std::stof(config.GetConfigVal("GAMMAPTBINSDOLOG"));
  Double_t gammaPtBins[nMaxPtBins+1];
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  std::vector<std::string> genGammaPtBinsStr, recoGammaPtBinsStr;
  for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
    genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(pI));
    recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(pI));

    binsToLabelStr[genGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[recoGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
  }
  genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(nGammaPtBins));
  recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(nGammaPtBins));

  binsToLabelStr[genGammaPtBinsStr[genGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[recoGammaPtBinsStr[recoGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);

  const Int_t nMaxEtaPhiBins = 100;
  const Int_t nEtaBins = std::stoi(config.GetConfigVal("NETABINS"));
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nEtaBins, "NETABINS")) return 1;
  const Float_t etaBinsLow = std::stof(config.GetConfigVal("ETABINSLOW"));
  const Float_t etaBinsHigh = std::stof(config.GetConfigVal("ETABINSHIGH"));
  const Bool_t etaBinsDoAbs = std::stoi(config.GetConfigVal("ETABINSDOABS"));
  if(etaBinsDoAbs && etaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsLow \'" << etaBinsLow << "\' less than 0 despite requested etaBinsDoAbs \'" << etaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  Double_t etaBins[nMaxEtaPhiBins+1];
  getLinBins(etaBinsLow, etaBinsHigh, nEtaBins, etaBins);

  const Int_t nPhiBins = std::stoi(config.GetConfigVal("NPHIBINS"));
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nPhiBins, "NPHIBINS")) return 1;
  Double_t phiBins[nMaxEtaPhiBins+1];
  getLinBins(-TMath::Pi()+0.01, TMath::Pi()+0.01, nPhiBins, phiBins);
 
 //Pt sub bins handling
  const Int_t nGammaPtBinsSub = std::stoi(config.GetConfigVal("NGAMMAPTBINSSUB"));
  if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
  const Float_t gammaPtBinsSubLow = std::stof(config.GetConfigVal("GAMMAPTBINSSUBLOW"));
  const Float_t gammaPtBinsSubHigh = std::stof(config.GetConfigVal("GAMMAPTBINSSUBHIGH"));
  if(gammaPtBinsSubLow < gammaPtBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubLow \'" << gammaPtBinsSubLow << "\' less than gammaPtBinsLow \'" << gammaPtBinsLow << "\'. return 1" << std::endl;
  if(gammaPtBinsSubHigh > gammaPtBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubHigh \'" << gammaPtBinsSubHigh << "\' greater than gammaPtBinsHigh \'" << gammaPtBinsHigh << "\'. return 1" << std::endl;
  if(gammaPtBinsSubLow < gammaPtBinsLow || gammaPtBinsSubHigh > gammaPtBinsHigh) return 1;
  const Bool_t gammaPtBinsSubDoLog = std::stof(config.GetConfigVal("GAMMAPTBINSDOLOG"));
  Double_t gammaPtBinsSub[nMaxSubBins+1];
  if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
  else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
  std::vector<std::string> gammaPtBinsSubStr;
  for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
    gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(pI));

    binsToLabelStr[gammaPtBinsSubStr[pI]] = prettyString(gammaPtBinsSub[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[pI+1], 1, false);
  }
  gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(nGammaPtBinsSub));

  binsToLabelStr[gammaPtBinsSubStr[gammaPtBinsSubStr.size()-1]] = prettyString(gammaPtBinsSub[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false);

  //Eta sub bins handling
  const Int_t nEtaBinsSub = std::stoi(config.GetConfigVal("NETABINSSUB"));
  if(!goodBinning(inConfigFileName, nMaxSubBins, nEtaBinsSub, "NETABINSSUB")) return 1;
  const Float_t etaBinsSubLow = std::stof(config.GetConfigVal("ETABINSSUBLOW"));
  const Float_t etaBinsSubHigh = std::stof(config.GetConfigVal("ETABINSSUBHIGH"));
  if(etaBinsSubLow < etaBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than etaBinsLow \'" << etaBinsLow << "\'. return 1" << std::endl;
  if(etaBinsSubHigh > etaBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubHigh \'" << etaBinsSubHigh << "\' greater than etaBinsHigh \'" << etaBinsHigh << "\'. return 1" << std::endl;
  if(etaBinsSubLow < etaBinsLow || etaBinsSubHigh > etaBinsHigh) return 1;
  const	Bool_t etaBinsSubDoAbs = std::stoi(config.GetConfigVal("ETABINSSUBDOABS"));
  if(etaBinsSubDoAbs && etaBinsSubLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than 0 despite requested etaBinsSubDoAbs \'" << etaBinsSubDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t etaBinsSub[nMaxSubBins+1];
  getLinBins(etaBinsSubLow, etaBinsSubHigh, nEtaBinsSub, etaBinsSub);
  std::vector<std::string> etaBinsSubStr;
  std::string preStr = "";
  if(etaBinsSubDoAbs) preStr = "Abs";
  for(Int_t eI = 0; eI < nEtaBinsSub; ++eI){
    etaBinsSubStr.push_back(preStr + "Eta" + std::to_string(eI));

    if(etaBinsSubDoAbs) binsToLabelStr[etaBinsSubStr[eI]] = prettyString(etaBinsSub[eI], 2, false) + "<|#eta_{#gamma}|<" + prettyString(etaBinsSub[eI+1], 2, false);
    else binsToLabelStr[etaBinsSubStr[eI]] = prettyString(etaBinsSub[eI], 2, false) + " < #eta_{#gamma} < " + prettyString(etaBinsSub[eI+1], 2, false);
  }
  etaBinsSubStr.push_back(preStr + "Eta" + std::to_string(nEtaBinsSub));
  if(etaBinsSubDoAbs) binsToLabelStr[etaBinsSubStr[etaBinsSubStr.size()-1]] = prettyString(etaBinsSub[0], 2, false) + "<|#eta_{#gamma}|<" + prettyString(etaBinsSub[nEtaBinsSub], 2, false);
  else binsToLabelStr[etaBinsSubStr[etaBinsSubStr.size()-1]] = prettyString(etaBinsSub[0], 2, false) + " < #eta_{#gamma} < " + prettyString(etaBinsSub[nEtaBinsSub], 2, false);
  
  const Int_t nJtPtBins = std::stoi(config.GetConfigVal("NJTPTBINS"));
  if(!goodBinning(inConfigFileName, nMaxPtBins, nJtPtBins, "NJTPTBINS")) return 1;
  const Float_t jtPtBinsLow = std::stof(config.GetConfigVal("JTPTBINSLOW"));
  const Float_t jtPtBinsHigh = std::stof(config.GetConfigVal("JTPTBINSHIGH"));
  const Bool_t jtPtBinsDoLog = std::stof(config.GetConfigVal("JTPTBINSDOLOG"));
  Double_t jtPtBins[nMaxPtBins+1];
  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);

  std::string jtPtBinsStr = "JtPt0";
  std::string jtPtBinsLabel = prettyString(jtPtBinsLow,1,false) + " < p_{T,jet} < " + prettyString(jtPtBinsHigh,1,false);
  binsToLabelStr[jtPtBinsStr] = jtPtBinsLabel;   
    
  const Int_t nDPhiBins = std::stoi(config.GetConfigVal("NDPHIBINS"));
  const Double_t gammaJtDPhiCut = mathStringToNum(config.GetConfigVal("GAMMAJTDPHI"));  
  const std::string gammaJtDPhiStr = "DPhi0";
  std::string gammaJtDPhiLabel = returnAllCapsString(config.GetConfigVal("GAMMAJTDPHI"));  
  if(gammaJtDPhiLabel.find("PI") != std::string::npos) gammaJtDPhiLabel.replace(gammaJtDPhiLabel.find("PI"), 2, "#pi");
  gammaJtDPhiLabel = "|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiLabel;
  binsToLabelStr[gammaJtDPhiStr] = gammaJtDPhiLabel;
  
  const Int_t nXJBins = std::stoi(config.GetConfigVal("NXJBINS"));
  const Float_t xjBinsLow = std::stof(config.GetConfigVal("XJBINSLOW"));
  const Float_t xjBinsHigh = std::stof(config.GetConfigVal("XJBINSHIGH"));
  Double_t xjBins[nMaxPtBins+1];
  getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* runNumber_p = nullptr;
  TH1F* lumiFractionPerRun_p = nullptr;
  TH1F* pthat_p = nullptr;
  TH1F* pthat_Unweighted_p = nullptr;
  TH1F* centrality_p = nullptr;
  TH1F* centrality_Unweighted_p = nullptr;
  TH1F* photonPtVCentEta_p[nMaxCentBins][nMaxSubBins+1];
  TH2F* photonGenResVCentEta_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH2F* photonEtaPt_p[nMaxCentBins];
  TH2F* photonEtaPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  
  TH1F* photonJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH2F* photonJtGenResVCentGenPtRecoPt_p[nMaxCentBins][nMaxPtBins][nMaxPtBins];
  
  TH1F* photonMixJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1F* photonSubJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1F* photonGenJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1F* photonGenMatchedJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  
  TH1F* photonJtFakeVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(isMC){
    pthat_p = new TH1F(("pthat_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    pthat_Unweighted_p = new TH1F(("pthat_Unweighted_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    centerTitles({pthat_p, pthat_Unweighted_p});
  }
  
  if(!isPP){
    centrality_p = new TH1F(("centrality_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
    centerTitles(centrality_p);

    if(isMC){
      centrality_Unweighted_p = new TH1F(("centrality_Unweighted_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
      centerTitles(centrality_Unweighted_p);
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      photonPtVCentEta_p[cI][eI] = new TH1F(("photonPtVCentEta_" + centBinsStr[cI] + "_" + etaBinsSubStr[eI]+ "_h").c_str(), ";#gamma p_{T} [GeV];Counts", nGammaPtBins, gammaPtBins);
      centerTitles(photonPtVCentEta_p[cI][eI]);

      if(isMC){
	photonGenResVCentEta_p[cI][eI] = new TH2F(("photonGenResVCentEta_" + centBinsStr[cI] + "_" + etaBinsSubStr[eI] + "_h").c_str(), ";Reco. Photon p_{T};Gen. Photon p_{T}", nGammaPtBins, gammaPtBins, nGammaPtBins, gammaPtBins);
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonEtaVCentPt_p[cI][pI] = new TH1F(("photonEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";#gamma #eta;Counts", nEtaBins, etaBins);
      photonPhiVCentPt_p[cI][pI] = new TH1F(("photonPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";#gamma #phi;Counts", nPhiBins, phiBins);

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";#Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonJtPtVCentPt_p[cI][pI] = new TH1F(("photonJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
      photonJtEtaVCentPt_p[cI][pI] = new TH1F(("photonJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
      photonJtXJVCentPt_p[cI][pI] = new TH1F(("photonJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);

      photonJtMultVCentPt_p[cI][pI] = new TH1F(("photonJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 10, -0.5, 9.5);

      centerTitles({photonEtaVCentPt_p[cI][pI], photonPhiVCentPt_p[cI][pI], photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI], photonJtMultVCentPt_p[cI][pI]});
      setSumW2({photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI]});
  
      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonMixJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";Mixed event #Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);	
	photonMixJtPtVCentPt_p[cI][pI] = new TH1F(("photonMixJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
	photonMixJtEtaVCentPt_p[cI][pI] = new TH1F(("photonMixJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
	photonMixJtXJVCentPt_p[cI][pI] = new TH1F(("photonMixJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);

	photonMixJtMultVCentPt_p[cI][pI] = new TH1F(("photonMixJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Mixed Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 10, -0.5, 9.5);

	photonSubJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonSubJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";Subtracted event #Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);	
	photonSubJtPtVCentPt_p[cI][pI] = new TH1F(("photonSubJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
	photonSubJtEtaVCentPt_p[cI][pI] = new TH1F(("photonSubJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
	photonSubJtXJVCentPt_p[cI][pI] = new TH1F(("photonSubJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);

	photonSubJtMultVCentPt_p[cI][pI] = new TH1F(("photonSubJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Subtracted Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 10, -0.5, 9.5);

	centerTitles({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI], photonSubJtDPhiVCentPt_p[cI][pI], photonSubJtPtVCentPt_p[cI][pI], photonSubJtEtaVCentPt_p[cI][pI], photonSubJtXJVCentPt_p[cI][pI], photonSubJtMultVCentPt_p[cI][pI]});
	setSumW2({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI], photonSubJtDPhiVCentPt_p[cI][pI], photonSubJtPtVCentPt_p[cI][pI], photonSubJtEtaVCentPt_p[cI][pI], photonSubJtXJVCentPt_p[cI][pI], photonSubJtMultVCentPt_p[cI][pI]});
      }
  
      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonGenJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";Gen. #Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);	
	photonGenJtPtVCentPt_p[cI][pI] = new TH1F(("photonGenJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
	photonGenJtEtaVCentPt_p[cI][pI] = new TH1F(("photonGenJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
	photonGenJtXJVCentPt_p[cI][pI] = new TH1F(("photonGenJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);
	photonGenJtMultVCentPt_p[cI][pI] = new TH1F(("photonGenJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen. Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 10, -0.5, 9.5);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";Gen.-matched #Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);	
	photonGenMatchedJtPtVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
	photonGenMatchedJtEtaVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
	photonGenMatchedJtXJVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);
	photonGenMatchedJtMultVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen.-matched Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 10, -0.5, 9.5);

	photonJtFakeVCentPt_p[cI][pI] = new TH1F(("photonJtFakeVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} ;#frac{N_{Fake jets}}{N_{All jets}}", nJtPtBins, jtPtBins);
	
	centerTitles({photonGenJtDPhiVCentPt_p[cI][pI], photonGenJtPtVCentPt_p[cI][pI], photonGenJtEtaVCentPt_p[cI][pI], photonGenJtXJVCentPt_p[cI][pI], photonGenJtMultVCentPt_p[cI][pI], photonGenMatchedJtDPhiVCentPt_p[cI][pI], photonGenMatchedJtPtVCentPt_p[cI][pI], photonGenMatchedJtEtaVCentPt_p[cI][pI], photonGenMatchedJtXJVCentPt_p[cI][pI], photonGenMatchedJtMultVCentPt_p[cI][pI], photonJtFakeVCentPt_p[cI][pI]});

      	setSumW2({photonGenJtDPhiVCentPt_p[cI][pI], photonGenJtPtVCentPt_p[cI][pI], photonGenJtEtaVCentPt_p[cI][pI], photonGenJtXJVCentPt_p[cI][pI], photonGenJtMultVCentPt_p[cI][pI], photonGenMatchedJtDPhiVCentPt_p[cI][pI], photonGenMatchedJtPtVCentPt_p[cI][pI], photonGenMatchedJtEtaVCentPt_p[cI][pI], photonGenMatchedJtXJVCentPt_p[cI][pI], photonGenMatchedJtMultVCentPt_p[cI][pI], photonJtFakeVCentPt_p[cI][pI]});
      }
    }

    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2] = new TH2F(("photonJtGenResVCentGenPtRecoPt_" + centBinsStr[cI] + "_" + genGammaPtBinsStr[ptI] + "_" + recoGammaPtBinsStr[ptI2] + "_h").c_str(), "", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	}
      }
    }
    
    photonEtaPt_p[cI] = new TH2F(("photonEtaPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;#gamma p_{T} [GeV]", nEtaBins, etaBins, nGammaPtBins, gammaPtBins);

    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI] = new TH2F(("photonEtaPhi_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";#gamma #eta;#gamma #phi", nEtaBins, etaBins, nPhiBins, phiBins);
      centerTitles(photonEtaPhiVCentPt_p[cI][pI]);
    }
    
    centerTitles(photonEtaPt_p[cI]);
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TFile* inFile_p = new TFile(inROOTFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("runNumber", 1);
  Int_t runMin = inTree_p->GetMinimum("runNumber");
  Int_t runMax = inTree_p->GetMaximum("runNumber");
  Int_t nRunBins = runMax - runMin;
  Float_t runMinF = ((Float_t)runMin) - 0.5;
  Float_t runMaxF = ((Float_t)runMax) + 0.5;

  outFile_p->cd();
  runNumber_p = new TH1F(("runNumber_" + systStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);
  if(!isMC) lumiFractionPerRun_p = new TH1F(("lumiFractionPerRun_" + systStr + "_h").c_str(), ";Run;Fraction of Lumiblocks", nRunBins+1, runMinF, runMaxF);
  centerTitles(runNumber_p);
  if(!isMC) centerTitles(lumiFractionPerRun_p);

  std::map<unsigned long long, bool> runLumiIsFired;
  std::map<int, int> runLumiCounter;
  std::map<int, int> runLumiTotal;
  if(!isMC){
    std::ifstream inFile(inGRLFileName.c_str());
    std::string tempStr;

    std::string currRunStr = "";
    
    while(std::getline(inFile, tempStr)){
      if(tempStr.find("<Run") != std::string::npos){
	tempStr.replace(0, tempStr.find(">")+1, "");
	tempStr.replace(tempStr.rfind("<"), tempStr.size(), "");
	
	currRunStr = tempStr;

	if(runLumiCounter.count(std::stoi(currRunStr)) != 0) std::cout << "Warning - counts found already for run \'" << currRunStr << "\'" << std::endl;
	else{
	  runLumiCounter[std::stoi(currRunStr)] = 0;
	  runLumiTotal[std::stoi(currRunStr)] = 0;
	}
      }
      else if(currRunStr.size() != 0 && tempStr.find("<LB") != std::string::npos){
	tempStr.replace(0, tempStr.find("\"")+1, "");
	tempStr.replace(tempStr.rfind("\""), tempStr.size(), "");
	std::string firstNumStr = tempStr.substr(0, tempStr.find("\""));
	std::string secondNumStr = tempStr;
	while(secondNumStr.find("\"") != std::string::npos){
	  secondNumStr.replace(0, secondNumStr.find("\"")+1, "");
	}
	while(firstNumStr.size() < 3){firstNumStr = "0" + firstNumStr;}
	while(secondNumStr.size() < 3){secondNumStr = "0" + secondNumStr;}

	Int_t firstNum = std::stoi(firstNumStr);
	Int_t secondNum = std::stoi(secondNumStr);
	
	for(Int_t nI = firstNum; nI <= secondNum; ++nI){
	  ++(runLumiTotal[std::stoi(currRunStr)]);
	  runLumiIsFired[runLumiKey.GetKey({(unsigned long long)std::stoi(currRunStr), (unsigned long long)nI})] = false;
	}
      }
    }
    
    inFile.close();
  }
  
  inFile_p->cd();
  
  //Grab the hltbranches for some basic prescale checks
  std::vector<std::string> listOfBranches = getVectBranchList(inTree_p);
  std::vector<std::string> hltList;
  std::vector<std::string> hltListPres;
  std::string hltStr = "HLT_";
  std::string prescaleStr = "_prescale";

  //HLTLists are built
  for(auto const & branchStr : listOfBranches){
    if(branchStr.size() < hltStr.size()) continue;
    if(!isStrSame(branchStr.substr(0, hltStr.size()), hltStr)) continue;

    if(branchStr.find(prescaleStr) != std::string::npos) hltListPres.push_back(branchStr);
    else hltList.push_back(branchStr);
  }

  bool allHLTPrescalesFound = true;
  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    std::string modStr = hltList[hI] + prescaleStr;
    if(!vectContainsStr(modStr, &hltListPres)){
      allHLTPrescalesFound = false;
      std::cout << "HLT " << hltList[hI] << " has no prescale. return 1" << std::endl;
    }
  }
  if(!allHLTPrescalesFound) return 1;

  float hltPrescaleDelta = 0.01;
  std::vector<bool*> hltVect;
  std::vector<float*> hltPrescaleVect;
  Int_t runNumber;
  UInt_t lumiBlock;
  Float_t pthat;
  Float_t sampleWeight;
  Float_t ncollWeight;
  Float_t fullWeight;
  Float_t fcalA_et, fcalC_et;
  Float_t evtPlane2Phi;
  std::vector<float>* vert_z_p=nullptr;
  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
  
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  
  
  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<int>* akt4hi_truthpos_p=nullptr;

  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;
  
  TFile* mixFile_p = nullptr;
  TTree* mixTree_p = nullptr;
  if(doMix){
    mixFile_p = new TFile(inMixFileName.c_str(), "READ");
    mixTree_p = (TTree*)mixFile_p->Get("gammaJetTree_p");

    mixTree_p->SetBranchStatus("*", 0);
    mixTree_p->SetBranchStatus("vert_z", 1);

    
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_pt", 1);
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_eta", 1);
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_phi", 1);

    mixTree_p->SetBranchAddress("vert_z", &vert_z_p);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);

    if(!isPP){
      mixTree_p->SetBranchStatus("fcalA_et", 1);
      mixTree_p->SetBranchStatus("fcalC_et", 1);
      
      mixTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
      mixTree_p->SetBranchAddress("fcalC_et", &fcalC_et);

      if(doMixPsi2){
	mixTree_p->SetBranchStatus("evtPlane2Phi", 1);	
	mixTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
      }
    }

    
    const ULong64_t nMixEntries = mixTree_p->GetEntries();

    for(ULong64_t entry = 0; entry < nMixEntries; ++entry){
      mixTree_p->GetEntry(entry);

      double vert_z = vert_z_p->at(0);
      vert_z /= 1000.;
      if(vert_z <= -15. || vert_z >= 15.) continue;      
      //      if(vert_z <= vzMixBinsLow || vert_z >= vzMixBinsHigh) continue;
      
      Double_t cent = -1;
      unsigned long long centPos = 0;
      unsigned long long psi2Pos = 0;
      if(!isPP){
	cent = centTable.GetCent(fcalA_et + fcalC_et);
	if(cent < mixCentBinsLow || cent >= mixCentBinsHigh) continue;
	if(doMixCent) centPos = ghostPos(nMixCentBins, mixCentBins, cent);

	if(doMixPsi2){
	  if(evtPlane2Phi > TMath::Pi()/2) evtPlane2Phi -= TMath::Pi();
	  else if(evtPlane2Phi < -TMath::Pi()/2) evtPlane2Phi += TMath::Pi();

	  psi2Pos = ghostPos(nMixPsi2Bins, mixPsi2Bins, evtPlane2Phi);
	}	
      }      

      unsigned long long vzPos = 0;
      if(doMixVz) vzPos = ghostPos(nMixVzBins, mixVzBins, vert_z);
      
      std::vector<unsigned long long> eventKeyVect;
      if(doMixCent) eventKeyVect.push_back(centPos);
      if(doMixPsi2) eventKeyVect.push_back(psi2Pos);
      if(doMixVz) eventKeyVect.push_back(vzPos);
      
      unsigned long long key = keyBoy.GetKey(eventKeyVect);//, vzPos, evtPlanePos});
      
      std::vector<TLorentzVector> jets;
      for(unsigned int jI = 0; jI < akt4hi_em_xcalib_jet_pt_p->size(); ++jI){
	if(akt4hi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	if(akt4hi_em_xcalib_jet_eta_p->at(jI) <= etaBinsLow) continue;
	if(akt4hi_em_xcalib_jet_eta_p->at(jI) >= etaBinsHigh) continue;

	TLorentzVector temp;
	temp.SetPtEtaPhiM(akt4hi_em_xcalib_jet_pt_p->at(jI), akt4hi_em_xcalib_jet_eta_p->at(jI), akt4hi_em_xcalib_jet_phi_p->at(jI), 0.0);

	jets.push_back(temp);
      }

      mixingMap[key].push_back(jets);
    }
    
    mixFile_p->Close();
    delete mixFile_p;
    inFile_p->cd();
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inTree_p->SetBranchStatus("*", 0);

  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    hltVect.push_back(new bool(false));
    hltPrescaleVect.push_back(new float(0.0));

    inTree_p->SetBranchStatus(hltList[hI].c_str(), 1);
    inTree_p->SetBranchStatus(hltListPres[hI].c_str(), 1);
  }  

  inTree_p->SetBranchStatus("runNumber", 1);
  inTree_p->SetBranchStatus("lumiBlock", 1);

  if(isMC){
    inTree_p->SetBranchStatus("pthat", 1);
    inTree_p->SetBranchStatus("sampleWeight", 1);
    inTree_p->SetBranchStatus("ncollWeight", 1);
    inTree_p->SetBranchStatus("fullWeight", 1);

    inTree_p->SetBranchStatus("truth_pt", 1);
    inTree_p->SetBranchStatus("truth_eta", 1);
    inTree_p->SetBranchStatus("truth_phi", 1);
    inTree_p->SetBranchStatus("truth_pdg", 1);

    inTree_p->SetBranchStatus("truthPhotonPt", 1);
    inTree_p->SetBranchStatus("truthPhotonEta", 1);
    inTree_p->SetBranchStatus("truthPhotonPhi", 1);
  }
  
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
    if(doMixPsi2) inTree_p->SetBranchStatus("evtPlane2Phi", 1);
  }

  inTree_p->SetBranchStatus("vert_z", 1);
  
  inTree_p->SetBranchStatus("photon_pt", 1);
  inTree_p->SetBranchStatus("photon_eta", 1);
  inTree_p->SetBranchStatus("photon_phi", 1);
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_pt", 1);
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_eta", 1);
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_phi", 1);
  
  if(isMC){
    inTree_p->SetBranchStatus("akt4hi_truthpos", 1);
    
    inTree_p->SetBranchStatus("akt4_truth_jet_pt", 1);
    inTree_p->SetBranchStatus("akt4_truth_jet_eta", 1);
    inTree_p->SetBranchStatus("akt4_truth_jet_phi", 1);
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    inTree_p->SetBranchAddress(hltList[hI].c_str(), hltVect[hI]);
    inTree_p->SetBranchAddress(hltListPres[hI].c_str(), hltPrescaleVect[hI]);
  }

  inTree_p->SetBranchAddress("runNumber", &runNumber);
  inTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  if(isMC){
    inTree_p->SetBranchAddress("pthat", &pthat);
    inTree_p->SetBranchAddress("sampleWeight", &sampleWeight);
    inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
    inTree_p->SetBranchAddress("fullWeight", &fullWeight);

    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
    inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);

    inTree_p->SetBranchAddress("truthPhotonPt", &truthPhotonPt);
    inTree_p->SetBranchAddress("truthPhotonEta", &truthPhotonEta);
    inTree_p->SetBranchAddress("truthPhotonPhi", &truthPhotonPhi);
  }

  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
    if(doMixPsi2) inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
  }

  inTree_p->SetBranchAddress("vert_z", &vert_z_p);
  
  inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
  inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
  inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);

  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);

  if(isMC){
    inTree_p->SetBranchAddress("akt4hi_truthpos", &akt4hi_truthpos_p);    

    inTree_p->SetBranchAddress("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
  }

  const ULong64_t nEntries = inTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::vector<std::vector<Double_t> > gammaCountsPerPtCent;
  for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
    gammaCountsPerPtCent.push_back({});

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      gammaCountsPerPtCent[pI].push_back(0.0);
    }
  }
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
  std::cout << "Processing " << nEntries << " events..." << std::endl;

  bool didOneFireMiss = false;
  std::vector<int> skippedCent;
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
  
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    double vert_z = vert_z_p->at(0);
    vert_z /= 1000.;
    if(vert_z <= -15. || vert_z >= 15.) continue;      
    //    if(vert_z <= vzMixBinsLow || vert_z >= vzMixBinsHigh) continue;

    if(!didOneFireMiss && !isMC){//only check this once per input
      //check at least one of the purported selection triggers fired
      bool oneFire = false;
      for(unsigned int hI = 0; hI < hltVect.size(); ++hI){
	if(*(hltVect[hI])){
	  oneFire = true;
	  break;
	}
      }

      if(!oneFire){
	std::cout << "WARNING - YOU HAVE EVENTS w/ NO TRIGGERS!!!" << std::endl;
	didOneFireMiss = true;
      }
    }
    
    //Check prescale is 1 in case i made a mistake on first unprescaled

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    for(unsigned int hI = 0; hI < hltPrescaleVect.size(); ++hI){
      if(TMath::Abs((*(hltPrescaleVect[hI])) - 1.0) > hltPrescaleDelta){
	std::cout << "WARNING - prescale for \'" << hltList[hI] << "\' has non-unity value, \'" << (*(hltPrescaleVect[hI])) << "\'." << std::endl;
      }
    }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 


    Int_t centPos = -1;
    Double_t cent = -1;
    if(!isPP){
      cent = centTable.GetCent(fcalA_et + fcalC_et);
      centPos = ghostPos(centBins, cent, true, doGlobalDebug);
    }
    else centPos = 0;

    if(centPos < 0){
      bool vectContainsCent = vectContainsInt((Int_t)cent, &skippedCent);

      if(!vectContainsCent){
	std::cout << "gdjNTupleToHist Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
	skippedCent.push_back((Int_t)cent);
      }
      
      continue;
    }

    if(!isMC) fullWeight = -1.0;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    unsigned long long tempKey = runLumiKey.GetKey({(unsigned long long)runNumber, (unsigned long long)lumiBlock});    
    runLumiIsFired[tempKey] = true;
      
    fillTH1(runNumber_p, runNumber, fullWeight);	
    if(!isPP){
      fillTH1(centrality_p, cent, fullWeight);
      if(isMC) centrality_Unweighted_p->Fill(cent);
    }

    if(isMC){
      fillTH1(pthat_p, pthat, fullWeight);
      pthat_Unweighted_p->Fill(pthat);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      if(photon_pt_p->at(pI) < gammaPtBins[0]) continue;
      if(photon_pt_p->at(pI) >= gammaPtBins[nGammaPtBins]) continue;

      Float_t etaValMain = photon_eta_p->at(pI);
      Float_t etaValSub = etaValMain;
      if(etaBinsDoAbs) etaValMain = TMath::Abs(etaValMain);
      if(etaBinsSubDoAbs) etaValSub = TMath::Abs(etaValSub);

      if(etaValMain <= etaBins[0]) continue;
      if(etaValMain >= etaBins[nEtaBins]) continue;
      
      Int_t ptPos = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(pI), true, doGlobalDebug);
      Int_t etaPos = ghostPos(nEtaBinsSub, etaBinsSub, etaValSub, true, doGlobalDebug);
    
      if(etaPos >= 0){
	fillTH1(photonPtVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), fullWeight);
	fillTH1(photonPtVCentEta_p[centPos][nEtaBinsSub], photon_pt_p->at(pI), fullWeight);


	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	if(isMC && truthPhotonPt > 0){
	  if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) < 0.2){
	    fillTH2(photonGenResVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), truthPhotonPt, fullWeight);
	    fillTH2(photonGenResVCentEta_p[centPos][nEtaBinsSub], photon_pt_p->at(pI), truthPhotonPt, fullWeight);
	  }
	}
      }
      
      if(ptPos >= 0){
	if(!isMC){
	  ++(gammaCountsPerPtCent[ptPos][centPos]);
	  ++(gammaCountsPerPtCent[nGammaPtBinsSub][centPos]);
	}
	else{
	  gammaCountsPerPtCent[ptPos][centPos] += fullWeight;
	  gammaCountsPerPtCent[nGammaPtBinsSub][centPos] += fullWeight;
	}
	
	fillTH1(photonEtaVCentPt_p[centPos][ptPos], etaValMain, fullWeight);
	fillTH1(photonPhiVCentPt_p[centPos][ptPos], photon_phi_p->at(pI), fullWeight);
	fillTH2(photonEtaPhiVCentPt_p[centPos][ptPos], etaValMain, photon_phi_p->at(pI), fullWeight);
	
      	fillTH1(photonEtaVCentPt_p[centPos][nGammaPtBinsSub], etaValMain, fullWeight);
	fillTH1(photonPhiVCentPt_p[centPos][nGammaPtBinsSub], photon_phi_p->at(pI), fullWeight);
	fillTH2(photonEtaPhiVCentPt_p[centPos][nGammaPtBinsSub], etaValMain, photon_phi_p->at(pI), fullWeight);

	int multCounter = 0;
	int multCounterGen = 0;
	int multCounterGenMatched = 0;
	
	for(unsigned int jI = 0; jI < akt4hi_em_xcalib_jet_pt_p->size(); ++jI){
	  if(akt4hi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	  if(akt4hi_em_xcalib_jet_eta_p->at(jI) <= etaBinsLow) continue;
	  if(akt4hi_em_xcalib_jet_eta_p->at(jI) >= etaBinsHigh) continue;
	  
	  Float_t dR = getDR(akt4hi_em_xcalib_jet_eta_p->at(jI), akt4hi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(pI), photon_phi_p->at(pI));
	  if(dR < 0.3) continue;
	  
	  Float_t dPhi = TMath::Abs(getDPHI(akt4hi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(pI)));
	  
	  fillTH1(photonJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	  fillTH1(photonJtDPhiVCentPt_p[centPos][nGammaPtBinsSub], dPhi, fullWeight);

	  if(isMC){
	    if(akt4hi_truthpos_p->at(jI) >= 0){
	      if(dPhi < TMath::Pi()/3.){
		std::cout << "Event to check: " << entry << std::endl;
	      }

	      fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	      fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][nGammaPtBinsSub], dPhi, fullWeight);
	    }
	  }
	  
	  if(dPhi >= gammaJtDPhiCut){
	    fillTH1(photonJtPtVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtPtVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);

	    ++multCounter;
	    
	    if(isMC){
	      if(akt4hi_truthpos_p->at(jI) < 0){
		fillTH1(photonJtFakeVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	      }
	      else{
 		fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
		fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
		fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
		fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
		fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
		fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
		
		++multCounterGenMatched;
		
		if(truthPhotonPt > 0){
		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		  if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) < 0.2){
		    if(truthPhotonPt >= gammaPtBins[0] && truthPhotonPt < gammaPtBins[nGammaPtBins]){
		      Int_t genPtPos = ghostPos(nGammaPtBins, gammaPtBins, truthPhotonPt, true, doGlobalDebug);
		      Int_t recoPtPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(pI), true, doGlobalDebug);		
		      
		      fillTH2(photonJtGenResVCentGenPtRecoPt_p[centPos][genPtPos][recoPtPos], akt4hi_em_xcalib_jet_pt_p->at(jI), akt4_truth_jet_pt_p->at(akt4hi_truthpos_p->at(jI)), fullWeight);
		    }
		  }
		}
	      }
	    }
	  }
	}
      
	if(isMC){
	  for(unsigned int tI = 0; tI < akt4_truth_jet_pt_p->size(); ++tI){
	    if(akt4_truth_jet_pt_p->at(tI) < jtPtBinsLow) continue;
	    if(akt4_truth_jet_eta_p->at(tI) <= etaBinsLow) continue;
	    if(akt4_truth_jet_eta_p->at(tI) >= etaBinsHigh) continue;
	  
	    Float_t dR = getDR(akt4_truth_jet_eta_p->at(tI), akt4_truth_jet_phi_p->at(tI), photon_eta_p->at(pI), photon_phi_p->at(pI));
	    if(dR < 0.3) continue;
	  
	    Float_t dPhi = TMath::Abs(getDPHI(akt4_truth_jet_phi_p->at(tI), photon_phi_p->at(pI)));
	  
	    fillTH1(photonGenJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	    fillTH1(photonGenJtDPhiVCentPt_p[centPos][nGammaPtBinsSub], dPhi, fullWeight);

	    if(dPhi >= gammaJtDPhiCut){
	      fillTH1(photonGenJtPtVCentPt_p[centPos][ptPos], akt4_truth_jet_pt_p->at(tI), fullWeight);
	      fillTH1(photonGenJtPtVCentPt_p[centPos][nGammaPtBinsSub], akt4_truth_jet_pt_p->at(tI), fullWeight);
	      fillTH1(photonGenJtEtaVCentPt_p[centPos][ptPos], akt4_truth_jet_eta_p->at(tI), fullWeight);
	      fillTH1(photonGenJtEtaVCentPt_p[centPos][nGammaPtBinsSub], akt4_truth_jet_eta_p->at(tI), fullWeight);
	      fillTH1(photonGenJtXJVCentPt_p[centPos][ptPos], akt4_truth_jet_pt_p->at(tI)/photon_pt_p->at(pI), fullWeight);
	      fillTH1(photonGenJtXJVCentPt_p[centPos][nGammaPtBinsSub], akt4_truth_jet_pt_p->at(tI)/photon_pt_p->at(pI), fullWeight);
	      ++multCounterGen;
	    }
	  }
	}

	fillTH1(photonJtMultVCentPt_p[centPos][ptPos], multCounter, fullWeight);
	fillTH1(photonJtMultVCentPt_p[centPos][nGammaPtBinsSub], multCounter, fullWeight);

	if(isMC){
	  fillTH1(photonGenJtMultVCentPt_p[centPos][ptPos], multCounterGen, fullWeight);
	  fillTH1(photonGenJtMultVCentPt_p[centPos][nGammaPtBinsSub], multCounterGen, fullWeight);

	  fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][ptPos], multCounterGenMatched, fullWeight);
	  fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][nGammaPtBinsSub], multCounterGenMatched, fullWeight);
	}
	
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	if(doMix){

	  unsigned long long mixCentPos = 0;
	  unsigned long long mixPsi2Pos = 0;
	  if(!isPP){
	    if(doMixCent) mixCentPos = ghostPos(nMixCentBins, mixCentBins, cent);
	    if(doMixPsi2){
	      if(evtPlane2Phi > TMath::Pi()/2) evtPlane2Phi -= TMath::Pi();
	      else if(evtPlane2Phi < -TMath::Pi()/2) evtPlane2Phi += TMath::Pi();
	      mixPsi2Pos = ghostPos(nMixPsi2Bins, mixPsi2Bins, evtPlane2Phi);
	    }
	  }

	  unsigned long long mixVzPos = 0;
	  if(doMixVz) mixVzPos = ghostPos(nMixVzBins, mixVzBins, vert_z);
	  
	  std::vector<unsigned long long> eventKeyVect;
	  if(doMixCent) eventKeyVect.push_back(mixCentPos);
	  if(doMixPsi2) eventKeyVect.push_back(mixPsi2Pos);
	  if(doMixVz) eventKeyVect.push_back(mixVzPos);
	  
	  unsigned long long key = keyBoy.GetKey(eventKeyVect);
	  unsigned long long maxPos = mixingMap[key].size();
	  if(maxPos == 0){
	    std::cout << "WHOOPS NO AVAILABLE MIXED EVENT. bailing" << std::endl;
	    std::cout << key << ", " << mixCentPos << ", " << cent << std::endl;
	    return 1;
	  }

	  unsigned long long jetPos = maxPos;
	  while(jetPos == maxPos){jetPos = randGen_p->Uniform(0, maxPos-1);}
	  
	  std::vector<TLorentzVector> jets = mixingMap[key][jetPos];
	       
	  multCounter = 0;
	  for(unsigned int jI = 0; jI < jets.size(); ++jI){
	    if(jets[jI].Pt()< jtPtBinsLow) continue;
	    if(jets[jI].Eta() <= etaBinsLow) continue;
	    if(jets[jI].Eta() >= etaBinsHigh) continue;
	  
	    Float_t dR = getDR(jets[jI].Eta(), jets[jI].Phi(), photon_eta_p->at(pI), photon_phi_p->at(pI));
	    if(dR < 0.3) continue;
	    
	    Float_t dPhi = TMath::Abs(getDPHI(jets[jI].Phi(), photon_phi_p->at(pI)));
	  
	    fillTH1(photonMixJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	    fillTH1(photonMixJtDPhiVCentPt_p[centPos][nGammaPtBinsSub], dPhi, fullWeight);
	    
	    if(dPhi >= gammaJtDPhiCut){
	      fillTH1(photonMixJtPtVCentPt_p[centPos][ptPos], jets[jI].Pt(), fullWeight);
	      fillTH1(photonMixJtPtVCentPt_p[centPos][nGammaPtBinsSub], jets[jI].Pt(), fullWeight);
	      fillTH1(photonMixJtEtaVCentPt_p[centPos][ptPos], jets[jI].Eta(), fullWeight);
	      fillTH1(photonMixJtEtaVCentPt_p[centPos][nGammaPtBinsSub], jets[jI].Eta(), fullWeight);
	      fillTH1(photonMixJtXJVCentPt_p[centPos][ptPos], jets[jI].Pt()/photon_pt_p->at(pI), fullWeight);
	      fillTH1(photonMixJtXJVCentPt_p[centPos][nGammaPtBinsSub], jets[jI].Pt()/photon_pt_p->at(pI), fullWeight);

	      ++multCounter;
	    }	    
	  }

	  fillTH1(photonMixJtMultVCentPt_p[centPos][ptPos], multCounter, fullWeight);
	  fillTH1(photonMixJtMultVCentPt_p[centPos][nGammaPtBinsSub], multCounter, fullWeight);	  
	}
      }
      
      fillTH2(photonEtaPt_p[centPos], etaValMain, photon_pt_p->at(pI), fullWeight);
    }
  }  
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  //Pre-write and delete some of these require some mods
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      if(isMC){
	photonJtFakeVCentPt_p[cI][pI]->Divide(photonJtPtVCentPt_p[cI][pI]);
	const double maxFakeVal = 0.45;
	double tempMax = getMax(photonJtFakeVCentPt_p[cI][pI]);

	if(tempMax > maxFakeVal) photonJtFakeVCentPt_p[cI][pI]->SetMaximum(maxFakeVal);
      }

      if(doMix){
	photonSubJtDPhiVCentPt_p[cI][pI]->Add(photonJtDPhiVCentPt_p[cI][pI], photonMixJtDPhiVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtPtVCentPt_p[cI][pI]->Add(photonJtPtVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtEtaVCentPt_p[cI][pI]->Add(photonJtEtaVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtXJVCentPt_p[cI][pI]->Add(photonJtXJVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtMultVCentPt_p[cI][pI]->Add(photonJtMultVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI], 1.0, -1.0);
      }
      
      photonJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

	photonSubJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      }
    }
  }

  runNumber_p->Write("", TObject::kOverwrite);

  if(!isPP){
    centrality_p->Write("", TObject::kOverwrite);
    if(isMC) centrality_Unweighted_p->Write("", TObject::kOverwrite);
  }

  if(isMC){
    pthat_p->Write("", TObject::kOverwrite);
    pthat_Unweighted_p->Write("", TObject::kOverwrite);
  }
  else{
    for(auto const & iter : runLumiIsFired){
      if(iter.second){
	unsigned long long run = (runLumiKey.InvertKey(iter.first))[0];
	
	++(runLumiCounter[run]);
      }
    }

    double num = 0.0;
    double denom = 0.0;
    for(auto const & iter : runLumiTotal){
      num += ((double)runLumiCounter[iter.first]);
      denom += ((double)iter.second);
      double val = ((double)runLumiCounter[iter.first])/((double)iter.second);
      int binVal = lumiFractionPerRun_p->FindBin(iter.first);
      lumiFractionPerRun_p->SetBinContent(binVal, val);
      lumiFractionPerRun_p->SetBinError(binVal, 0.0);
    }

    std::cout << "Total lumiblocks: " << num << "/" << denom << "=" << num/denom << std::endl;
    
    lumiFractionPerRun_p->Write("", TObject::kOverwrite);
  }
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TDirectoryFile* centDir_p = (TDirectoryFile*)outFile_p->mkdir(centBinsStr[cI].c_str());
    centDir_p->cd();
    
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      photonPtVCentEta_p[cI][eI]->Write("", TObject::kOverwrite);

      if(isMC){
	photonGenResVCentEta_p[cI][eI]->Write("", TObject::kOverwrite);
      }
    }    

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
 
    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      photonEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

      	photonSubJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      }
      
      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonJtFakeVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);	
      }
    }

    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2]->Write("", TObject::kOverwrite);
	}
      }
    }
    
    photonEtaPt_p[cI]->Write("", TObject::kOverwrite);

    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
    }
    centDir_p->Close();
    delete centDir_p;
  }

  //Before we do a bunch of deletions lets auto-generate some stats tables
  if(!isMC){
    std::string lStr = "l | ";
    for(int gI = 0; gI < nGammaPtBinsSub; ++gI){
      lStr = lStr + "l | ";
    }
    lStr = lStr + "l";
    std::string multiColumn = std::to_string(nGammaPtBinsSub+2);
    
    std::string jtPtStr = prettyString(jtPtBins[0], 1, false) + " < \\ptj < " + prettyString(jtPtBins[nJtPtBins], 1, false);
    
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
    std::cout << "\\vspace{-0.6cm}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
    std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With Exactly 2 Jets (" << jtPtStr << ")} \\\\" << std::endl;
    std::string columnStr = " Centrality &";
    for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
      columnStr = columnStr + " " + prettyString(gammaPtBinsSub[pI], 1, false) + " < \\ptg < " + prettyString(gammaPtBinsSub[pI+1], 1, false) + " &";
    }
    columnStr = columnStr + " " + prettyString(gammaPtBinsSub[0], 1, false) + " < \\ptg < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false) + " \\\\ \\hline" ;
    
    std::cout << columnStr << std::endl;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){ 
      std::string outStr = binsToLabelStr[centBinsStr[cI]];
      if(outStr.find("%") != std::string::npos){
	outStr.replace(outStr.find("%"), 1, "\\%");
      }
      
      for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
	outStr = outStr + " & " + std::to_string(((int)photonJtMultVCentPt_p[cI][pI]->GetBinContent(3)));
      }
      
      std::cout << outStr << " \\\\" << std::endl;
    }  
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    std::cout << std::endl;
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
    std::cout << "\\vspace{-0.6cm}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
    std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With 2 or More Jets (" << jtPtStr << ")} \\\\" << std::endl;
    std::cout << columnStr << std::endl;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){ 
      std::string outStr = binsToLabelStr[centBinsStr[cI]];
      if(outStr.find("%") != std::string::npos){
	outStr.replace(outStr.find("%"), 1, "\\%");
      }
      
      for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
	int countNJet = 0;
	
	for(Int_t bIX = 2; bIX < photonJtMultVCentPt_p[cI][pI]->GetXaxis()->GetNbins(); ++bIX){
	  countNJet += photonJtMultVCentPt_p[cI][pI]->GetBinContent(bIX+1);
	}
	
	outStr = outStr + " & " + std::to_string(countNJet);
      }

      std::cout << outStr << " \\\\" << std::endl;
    }  
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
  }
  
  delete runNumber_p;

  if(isMC){
    delete pthat_p;
    delete pthat_Unweighted_p;
  }
  else delete lumiFractionPerRun_p;

  if(!isPP){
    delete centrality_p;
    if(isMC) delete centrality_Unweighted_p;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      delete photonPtVCentEta_p[cI][eI];

      if(isMC){
	delete photonGenResVCentEta_p[cI][eI];
      }
    }

    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      delete photonEtaVCentPt_p[cI][pI];
      delete photonPhiVCentPt_p[cI][pI];

      delete photonJtDPhiVCentPt_p[cI][pI];
      delete photonJtPtVCentPt_p[cI][pI];
      delete photonJtEtaVCentPt_p[cI][pI];
      delete photonJtXJVCentPt_p[cI][pI];
      delete photonJtMultVCentPt_p[cI][pI];

      if(doMix){
	delete photonMixJtDPhiVCentPt_p[cI][pI];
	delete photonMixJtPtVCentPt_p[cI][pI];
	delete photonMixJtEtaVCentPt_p[cI][pI];
	delete photonMixJtXJVCentPt_p[cI][pI];
	delete photonMixJtMultVCentPt_p[cI][pI];

      	delete photonSubJtDPhiVCentPt_p[cI][pI];
	delete photonSubJtPtVCentPt_p[cI][pI];
	delete photonSubJtEtaVCentPt_p[cI][pI];
	delete photonSubJtXJVCentPt_p[cI][pI];
	delete photonSubJtMultVCentPt_p[cI][pI];
      }
      
      if(isMC){
	delete photonGenJtDPhiVCentPt_p[cI][pI];
	delete photonGenJtPtVCentPt_p[cI][pI];
	delete photonGenJtEtaVCentPt_p[cI][pI];
	delete photonGenJtXJVCentPt_p[cI][pI];
	delete photonGenJtMultVCentPt_p[cI][pI];

	delete photonGenMatchedJtDPhiVCentPt_p[cI][pI];
	delete photonGenMatchedJtPtVCentPt_p[cI][pI];
	delete photonGenMatchedJtEtaVCentPt_p[cI][pI];
	delete photonGenMatchedJtXJVCentPt_p[cI][pI];
	delete photonGenMatchedJtMultVCentPt_p[cI][pI];

	delete photonJtFakeVCentPt_p[cI][pI];	
      }
    }

    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  delete photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2];
	}
      }
    }
    
    delete photonEtaPt_p[cI];

    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
      delete photonEtaPhiVCentPt_p[cI][pI];
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form          
  TEnv configEnv; // We will convert to tenv and store in file
  for(auto const & val : configMap){
    configEnv.SetValue(val.first.c_str(), val.second.c_str()); //Fill out the map
  }
  configEnv.Write("config", TObject::kOverwrite);  

  TEnv labelEnv;
  for(auto const & lab : binsToLabelStr){
    labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
  }
  labelEnv.Write("label", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  std::cout << "GDJNTUPLETOHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNTupleToHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjNTupleToHist(argv[1]);
  return retVal;
}
