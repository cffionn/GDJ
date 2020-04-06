//Author: Chris McGinn (2020.02.19)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
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

  if(!config.ContainsParamSet(necessaryParams)) return 1;  

  std::string inROOTFileName = config.GetConfigVal("INFILENAME");
  std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
  std::string outFileName = config.GetConfigVal("OUTFILENAME");
  std::string inMixFileName = config.GetConfigVal("MIXFILENAME");
  const bool doMix = std::stoi(config.GetConfigVal("DOMIX"));

  //Mixing categories, temp hardcoding
  //Centrality, percent level
  const Int_t nCentMixBins = 100;
  Float_t centMixBinsLow = 0;
  Float_t centMixBinsHigh = 100;
  Double_t centMixBins[nCentMixBins+1];
  getLinBins(centMixBinsLow, centMixBinsHigh, nCentMixBins, centMixBins);

  const Int_t nVzMixBins = 30;
  Float_t vzMixBinsLow = -15.0;
  Float_t vzMixBinsHigh = 15.0;
  Double_t vzMixBins[nVzMixBins+1];
  getLinBins(vzMixBinsLow, vzMixBinsHigh, nVzMixBins, vzMixBins);

  //Currently not used but will be in the future
  const Int_t nEvtPlaneMixBins = 16;
  Float_t evtPlaneMixBinsLow = -TMath::Pi()/2.;
  Float_t evtPlaneMixBinsHigh = TMath::Pi()/2.;
  Double_t evtPlaneMixBins[nEvtPlaneMixBins+1];
  getLinBins(evtPlaneMixBinsLow, evtPlaneMixBinsHigh, nEvtPlaneMixBins, evtPlaneMixBins);

  keyHandler keyBoy;
  keyBoy.Init({nCentMixBins-1, nVzMixBins-1, nEvtPlaneMixBins-1});
  std::map<unsigned long long, std::vector<std::vector<TLorentzVector> > > mixingMap;
  if(doMix){
    for(Int_t cI = 0; cI < nCentMixBins; ++cI){
      for(Int_t vI = 0; vI < nVzMixBins; ++vI){
	for(Int_t eI = 0; eI < nEvtPlaneMixBins; ++eI){
	  unsigned long long key = keyBoy.GetKey({(unsigned long long)cI, (unsigned long long)vI, (unsigned long long)eI});
	  mixingMap[key] = {};
	}
      }
    }
  }  
  
  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
  if(doMix && !check.checkFileExt(inMixFileName, ".txt")) return 1; // Check mixed file exists if mixing is requested
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();
  
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const bool isPP = std::stoi(config.GetConfigVal("ISPP"));
  const bool isMC = std::stoi(config.GetConfigVal("ISMC"));
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
  TH2F* photonJtGenResVCentGenPtRecoPt_p[nMaxCentBins][nMaxPtBins][nMaxPtBins];
  
  TH1F* photonMixJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];

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
    
      centerTitles({photonEtaVCentPt_p[cI][pI], photonPhiVCentPt_p[cI][pI], photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI]});
      setSumW2({photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI]});
  
      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonMixJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_h").c_str(), ";#Delta#phi_{#gamma,jet};#frac{N_{#gamma,jet}}{N_{#gamma}}", nDPhiBins, 0, TMath::Pi() + 0.01);	
	photonMixJtPtVCentPt_p[cI][pI] = new TH1F(("photonMixJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
	photonMixJtEtaVCentPt_p[cI][pI] = new TH1F(("photonMixJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet #eta;#frac{N_{#gamma,jet}}{N_{#gamma}}", nEtaBins, etaBins);
	photonMixJtXJVCentPt_p[cI][pI] = new TH1F(("photonMixJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";x_{J,#gamma};#frac{N_{#gamma,jet}}{N_{#gamma}}", nXJBins, xjBins);

	centerTitles({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI]});
	setSumW2({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI]});
      }
  
      if(isMC){
	photonJtFakeVCentPt_p[cI][pI] = new TH1F(("photonJtFakeVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + jtPtBinsStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} ;#frac{N_{Fake jets}}{N_{All jets}}", nJtPtBins, jtPtBins);
	
	centerTitles(photonJtFakeVCentPt_p[cI][pI]);
	setSumW2(photonJtFakeVCentPt_p[cI][pI]);
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
  centerTitles(runNumber_p);
  
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
  Float_t pthat;
  Float_t sampleWeight;
  Float_t ncollWeight;
  Float_t fullWeight;
  Float_t fcalA_et, fcalC_et;
  std::vector<float>* vert_z_p=nullptr;
  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;
  
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
    mixTree_p->SetBranchStatus("fcalA_et", 1);
    mixTree_p->SetBranchStatus("fcalC_et", 1);
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_pt", 1);
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_eta", 1);
    mixTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_phi", 1);

    mixTree_p->SetBranchAddress("vert_z", &vert_z_p);
    mixTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    mixTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
    mixTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);

    const ULong64_t nMixEntries = mixTree_p->GetEntries();

    for(ULong64_t entry = 0; entry < nMixEntries; ++entry){
      mixTree_p->GetEntry(entry);

      double vert_z = vert_z_p->at(0);
      vert_z /= 1000.;
      if(vert_z <= vzMixBinsLow || vert_z >= vzMixBinsHigh) continue;
      
      Double_t cent = -1;
      unsigned long long centPos = 0;
      if(!isPP){
	cent = centTable.GetCent(fcalA_et + fcalC_et);
	centPos = ghostPos(nCentMixBins, centMixBins, cent);
      }

      if(!isPP){
	if(cent < centMixBinsLow || cent >= centMixBinsHigh) continue;
      }

      unsigned long long vzPos = ghostPos(nVzMixBins, vzMixBins, vert_z);
      unsigned long long evtPlanePos = 0; // to be added

      unsigned long long key = keyBoy.GetKey({centPos, vzPos, evtPlanePos});
      
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

  if(isMC){
    inTree_p->SetBranchStatus("pthat", 1);
    inTree_p->SetBranchStatus("sampleWeight", 1);
    inTree_p->SetBranchStatus("ncollWeight", 1);
    inTree_p->SetBranchStatus("fullWeight", 1);

    inTree_p->SetBranchStatus("truth_pt", 1);
    inTree_p->SetBranchStatus("truth_eta", 1);
    inTree_p->SetBranchStatus("truth_phi", 1);
    inTree_p->SetBranchStatus("truth_pdg", 1);
  }
  
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
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
  if(isMC){
    inTree_p->SetBranchAddress("pthat", &pthat);
    inTree_p->SetBranchAddress("sampleWeight", &sampleWeight);
    inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
    inTree_p->SetBranchAddress("fullWeight", &fullWeight);

    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
    inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);
  }

  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
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
    if(vert_z <= vzMixBinsLow || vert_z >= vzMixBinsHigh) continue;

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
    
      Int_t truthPos = -1;
      if(isMC){
	if(doGlobalDebug){
	  std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  std::cout << truth_eta_p->size() << ", " << truth_pt_p->size() << ", " << truth_phi_p->size() << ", " << truth_pdg_p->size() << std::endl;
	}
	
	for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	  if(truth_pdg_p->at(tI) != 22) continue;
	  if(getDR(truth_eta_p->at(tI), truth_phi_p->at(tI), photon_eta_p->at(pI), photon_phi_p->at(pI)) < 0.2){
	    truthPos = tI;
	    break;
	  }
	}
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
      }
      
      if(etaPos >= 0){
	fillTH1(photonPtVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), fullWeight);
	fillTH1(photonPtVCentEta_p[centPos][nEtaBinsSub], photon_pt_p->at(pI), fullWeight);


	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	if(isMC){
	  if(truthPos >= 0){
	    fillTH2(photonGenResVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), truth_pt_p->at(truthPos), fullWeight);
	    fillTH2(photonGenResVCentEta_p[centPos][nEtaBinsSub], photon_pt_p->at(pI), truth_pt_p->at(truthPos), fullWeight);
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
	
	for(unsigned int jI = 0; jI < akt4hi_em_xcalib_jet_pt_p->size(); ++jI){
	  if(akt4hi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	  if(akt4hi_em_xcalib_jet_eta_p->at(jI) <= etaBinsLow) continue;
	  if(akt4hi_em_xcalib_jet_eta_p->at(jI) >= etaBinsHigh) continue;
	  
	  Float_t dR = getDR(akt4hi_em_xcalib_jet_eta_p->at(jI), akt4hi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(pI), photon_phi_p->at(pI));
	  if(dR < 0.3) continue;
	  
	  Float_t dPhi = TMath::Abs(getDPHI(akt4hi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(pI)));
	  
	  fillTH1(photonJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	  fillTH1(photonJtDPhiVCentPt_p[centPos][nGammaPtBinsSub], dPhi, fullWeight);
	  
	  if(dPhi >= gammaJtDPhiCut){
	    fillTH1(photonJtPtVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtPtVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][nGammaPtBinsSub], akt4hi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);

	    if(isMC){
	      if(akt4hi_truthpos_p->at(jI) < 0){
		fillTH1(photonJtFakeVCentPt_p[centPos][ptPos], akt4hi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	      }
	      else{
		if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		if(truthPos >= 0){
		  if(truth_pt_p->at(truthPos) >= gammaPtBins[0] && truth_pt_p->at(truthPos) < gammaPtBins[nGammaPtBins]){
		    Int_t genPtPos = ghostPos(nGammaPtBins, gammaPtBins, truth_pt_p->at(truthPos), true, doGlobalDebug);
		    Int_t recoPtPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(pI), true, doGlobalDebug);		
		
		    fillTH2(photonJtGenResVCentGenPtRecoPt_p[centPos][genPtPos][recoPtPos], akt4hi_em_xcalib_jet_pt_p->at(jI), akt4_truth_jet_pt_p->at(akt4hi_truthpos_p->at(jI)), fullWeight);
		  }
		}
	      }
	    }
	  }
	}
      
	if(doMix){
	  unsigned long long centMixPos = 0;
	  if(!isPP) centMixPos = ghostPos(nCentMixBins, centMixBins, cent);
	  
	  unsigned long long vzPos = ghostPos(nVzMixBins, vzMixBins, vert_z);
	  unsigned long long evtPlanePos = 0; // to be added
	  
	  unsigned long long key = keyBoy.GetKey({centMixPos, vzPos, evtPlanePos});
	  unsigned long long maxPos = mixingMap[key].size();
	  unsigned long long jetPos = maxPos;
	  while(jetPos == maxPos){jetPos = randGen_p->Uniform(0, maxPos);}
	  
	  std::vector<TLorentzVector> jets = mixingMap[key][jetPos];
	  
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
	    }	    
	  }
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
      
      photonJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
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

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      }
      
      if(isMC){
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

  delete runNumber_p;

  if(isMC){
    delete pthat_p;
    delete pthat_Unweighted_p;
  }

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

      if(doMix){
	delete photonJtDPhiVCentPt_p[cI][pI];
	delete photonJtPtVCentPt_p[cI][pI];
	delete photonJtEtaVCentPt_p[cI][pI];
	delete photonJtXJVCentPt_p[cI][pI];
      }
      
      if(isMC){
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
