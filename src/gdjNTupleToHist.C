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
  cppWatch globalTimer;
  globalTimer.start();

  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101
  TRandom3* randGen_p = new TRandom3(randSeed);


  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INFILENAME",
                                              "OUTFILENAME",
                                              "CENTFILENAME",
					      "MIXFILENAME",
					      "JETR",
					      "ASSOCGENMINPT",
					      "GAMMAEXCLUSIONDR",
					      "MIXJETEXCLUSIONDR",
					      "DOMIX",
					      "MIXCAP",
					      "NMIXEVENTS",
					      "DOSTRICTMIX",
                                              "ISPP",
					      "ISMC",
					      "CENTBINS",
					      "NGAMMAPTBINS",
					      "GAMMAPTBINSLOW",
					      "GAMMAPTBINSHIGH",
					      "GAMMAPTBINSDOLOG",
					      "GAMMAPTBINSDOCUSTOM",
					      "GAMMAPTBINSCUSTOM",
					      "NGAMMAETABINS",
					      "GAMMAETABINSLOW",
					      "GAMMAETABINSHIGH",
					      "GAMMAETABINSDOABS",
					      "NJTETABINS",
					      "JTETABINSLOW",
					      "JTETABINSHIGH",
					      "JTETABINSDOABS",
					      "NJTETABINSSUB",
					      "JTETABINSSUBLOW",
					      "JTETABINSSUBHIGH",
					      "JTETABINSSUBDOABS",
					      "NPHIBINS",
					      "NGAMMAPTBINSSUB",
					      "GAMMAPTBINSSUBLOW",
					      "GAMMAPTBINSSUBHIGH",
					      "GAMMAPTBINSSUBDOLOG",
					      "NGAMMAETABINSSUB",
					      "GAMMAETABINSSUBLOW",
					      "GAMMAETABINSSUBHIGH",
					      "GAMMAETABINSSUBDOABS",
					      "NJTPTBINS",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "JTPTBINSDOLOG",
					      "JTPTBINSDOCUSTOM",
					      "JTPTBINSCUSTOM",
					      "NDPHIBINS",
					      "DPHIBINSLOW",
					      "DPHIBINSHIGH",
					      "GAMMAJTDPHI",  
					      "NXJBINS",
					      "XJBINSLOW",
					      "XJBINSHIGH",
					      "XJBINSDOCUSTOM",
					      "XJBINSCUSTOM"};

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

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  //  if(!config.ContainsParamSet(necessaryParams)) return 1;  
 
  std::string inROOTFileName = config_p->GetValue("INFILENAME", "");
  std::string inCentFileName = config_p->GetValue("CENTFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string inGRLFileName = config_p->GetValue("GRLFILENAME", "");
  std::string inMixFileName = config_p->GetValue("MIXFILENAME", "");

  const std::string nMaxEvtStr = config_p->GetValue("NEVT", "");
  ULong64_t nMaxEvt = 0;
  if(nMaxEvtStr.size() != 0){
    nMaxEvt = std::stol(nMaxEvtStr);
  }

  const int jetR = config_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);
  const double assocGenMinPt = config_p->GetValue("ASSOCGENMINPT", 15.0);
  const double gammaExclusionDR = config_p->GetValue("GAMMAEXCLUSIONDR", 0.5);
  const double mixJetExclusionDR = config_p->GetValue("MIXJETEXCLUSIONDR", 0.5);

  const bool doMix = config_p->GetValue("DOMIX", 0);

  if(doMix){
    if(!checkEnvForParams(config_p, mixParams)) return 1;
  } 
  
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
  std::vector<std::string> nameVect;
  std::vector<std::vector<unsigned long long> > keyVect;
  
  if(doMix){
    if(!check.checkFileExt(inMixFileName, ".root")) return 1;

    doMixCent = (bool)config_p->GetValue("DOMIXCENT", 0);

    if(doMixCent){      
      nMixCentBins = config_p->GetValue("NMIXCENTBINS", 10);
      if(nMixCentBins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixCentBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixCentBinsLow = config_p->GetValue("MIXCENTBINSLOW", 0.0);
      mixCentBinsHigh = config_p->GetValue("MIXCENTBINSHIGH", 80.0);
      getLinBins(mixCentBinsLow, mixCentBinsHigh, nMixCentBins, mixCentBins);

      mixVect.push_back({});
      nameVect.push_back("Centrality");
      for(Int_t mI = 0; mI < nMixCentBins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);
	keyVect.push_back({(unsigned long long)mI});
      }
    }

    doMixPsi2 = (bool)config_p->GetValue("DOMIXPSI2", 0);

    if(doMixPsi2){
      nMixPsi2Bins = config_p->GetValue("NMIXPSI2BINS", 16);
      if(nMixPsi2Bins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixPsi2Bins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixPsi2BinsLow = config_p->GetValue("MIXPSI2BINSLOW", 0.0);
      mixPsi2BinsHigh = config_p->GetValue("MIXPSI2BINSHIGH", TMath::Pi()/2.0);
      getLinBins(mixPsi2BinsLow, mixPsi2BinsHigh, nMixPsi2Bins, mixPsi2Bins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      nameVect.push_back("Psi2");
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

    doMixVz = (bool)config_p->GetValue("DOMIXVZ", 0);

    if(doMixVz){
      nMixVzBins = config_p->GetValue("NMIXVZBINS", 10);
      if(nMixVzBins > nMaxMixBins){
	std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixVzBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixVzBinsLow = config_p->GetValue("MIXVZBINSLOW", -15.0);
      mixVzBinsHigh = config_p->GetValue("MIXVZBINSHIGH", 15.0);    
      getLinBins(mixVzBinsLow, mixVzBinsHigh, nMixVzBins, mixVzBins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      nameVect.push_back("vZ");
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

    if(doMixCent){
      std::cout << "MIXING IN CENTRALItY: " << std::endl;
      for(Int_t cI = 0; cI < nMixCentBins; ++cI){
	std::cout << " BIN " << cI << ": " << mixCentBins[cI] << "-" << mixCentBins[cI+1] << std::endl;
      }      
    }

    if(doMixPsi2){
      std::cout << "MIXING IN PSI2: " << std::endl;
      for(Int_t cI = 0; cI < nMixPsi2Bins; ++cI){
	std::cout << " BIN " << cI << ": " << mixPsi2Bins[cI] << "-" << mixPsi2Bins[cI+1] << std::endl;
      }      
    }

    if(doMixVz){
      std::cout << "MIXING IN VZ: " << std::endl;
      for(Int_t cI = 0; cI < nMixVzBins; ++cI){
	std::cout << " BIN " << cI << ": " << mixVzBins[cI] << "-" << mixVzBins[cI+1] << std::endl;
      }      
    }

  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  keyHandler runLumiKey("runLumiHandler");//for lumi estimation
  runLumiKey.Init({1000000, 10000});//runnumbers, then lumi
  
  keyHandler keyBoy("mixingHandler");//For Mixing
  std::map<unsigned long long, std::vector<std::vector<TLorentzVector> > > mixingMap;
  std::map<unsigned long long, std::vector<unsigned long long> > mixingMapEvtUseCounter;
  std::map<unsigned long long, unsigned long long> mixingMapCounter, signalMapCounterPre, signalMapCounterPost;

  if(doMix){
    std::vector<unsigned long long> sizes;
    for(unsigned int vI = 0; vI < mixVect.size(); ++vI){
      sizes.push_back(mixVect[vI].size()+1);
    }

    //    return 0;
    
    keyBoy.Init(nameVect, sizes);    

    for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
      unsigned long long key = keyBoy.GetKey(keyVect[vI]);
      mixingMap[key] = {};     
      mixingMap[key].reserve(40);

      mixingMapEvtUseCounter[key] = {};

      mixingMapCounter[key] = 0;
      signalMapCounterPre[key] = 0;
      signalMapCounterPost[key] = 0;
    }
  }  
  
  const bool isMC = config_p->GetValue("ISMC", 0);

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

  const bool isPP = config_p->GetValue("ISPP", 1);
  const Int_t nMaxSubBins = 60;
  const Int_t nMaxCentBins = 10;
  Int_t nCentBins = 1;

 std::string systStr = "PP";
  if(!isPP) systStr = "PbPb";
  
  std::vector<int> centBins;
  std::vector<std::string> centBinsStr = {systStr};
  std::map<std::string, std::string> binsToLabelStr;
  if(!isPP){
    centBins = strToVectI(config_p->GetValue("CENTBINS", "0,10,30,80"));
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
  const Int_t nGammaPtBins = config_p->GetValue("NGAMMAPTBINS", 10);
  if(!goodBinning(inConfigFileName, nMaxPtBins, nGammaPtBins, "NGAMMAPTBINS")) return 1;
  const Float_t gammaPtBinsLow = config_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = config_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)config_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)config_p->GetValue("GAMMAPTBINSDOCUSTOM", 0);
  Double_t gammaPtBins[nMaxPtBins+1];
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    std::string gammaPtBinsStr = config_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "CUSTOM GAMMA PT BINS REQUESTED, BUT \'GAMMAPTBINSCUSTOM\' not defined in config \'" << inConfigFileName << "\'. return 1" << std::endl;
    }
    std::vector<float> gammaPtBinsCustom = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsCustom.size()) - 1;
    if(nGammaPtBinsTemp != nGammaPtBins){
      std::cout << "REQUESTED CUSTON BINS \'" << gammaPtBinsStr << "\' HAS SIZE \'" << nGammaPtBinsTemp << "\' does not match given nGammaPtBins \'" << nGammaPtBins << "\'. return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsCustom.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsCustom[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  std::vector<std::string> genGammaPtBinsStr, recoGammaPtBinsStr, gammaPtBinsStr;
  for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
    genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(pI));
    recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(pI));
    gammaPtBinsStr.push_back("GammaPt" + std::to_string(pI));

    binsToLabelStr[genGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[recoGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[gammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
  }
  genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(nGammaPtBins));
  recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(nGammaPtBins));
  gammaPtBinsStr.push_back("GammaPt" + std::to_string(nGammaPtBins));

  binsToLabelStr[genGammaPtBinsStr[genGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[recoGammaPtBinsStr[recoGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[gammaPtBinsStr[gammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);

  const Int_t nMaxEtaPhiBins = 100;

  const Int_t nGammaEtaBins = config_p->GetValue("NGAMMAETABINS", 12);
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nGammaEtaBins, "NGAMMAETABINS")) return 1;
  const Float_t gammaEtaBinsLow = config_p->GetValue("GAMMAETABINSLOW", 0.0);
  const Float_t gammaEtaBinsHigh = config_p->GetValue("GAMMAETABINSHIGH", 2.37);
  const Bool_t gammaEtaBinsDoAbs = (bool)config_p->GetValue("GAMMAETABINSDOABS", 1);
  if(gammaEtaBinsDoAbs && gammaEtaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaEtaBinsLow \'" << gammaEtaBinsLow << "\' less than 0 despite requested gammaEtaBinsDoAbs \'" << gammaEtaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  Double_t gammaEtaBins[nMaxEtaPhiBins+1];
  getLinBins(gammaEtaBinsLow, gammaEtaBinsHigh, nGammaEtaBins, gammaEtaBins);

  const Int_t nJtEtaBins = config_p->GetValue("NJTETABINS", 14);
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nJtEtaBins, "NJTETABINS")) return 1;
  const Float_t jtEtaBinsLow = config_p->GetValue("JTETABINSLOW", 0.0);
  const Float_t jtEtaBinsHigh = config_p->GetValue("JTETABINSHIGH", 2.8);
  const Bool_t jtEtaBinsDoAbs = config_p->GetValue("JTETABINSDOABS", 1);
  if(jtEtaBinsDoAbs && jtEtaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains jtEtaBinsLow \'" << jtEtaBinsLow << "\' less than 0 despite requested jtEtaBinsDoAbs \'" << jtEtaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t jtEtaBins[nMaxEtaPhiBins+1];
  getLinBins(jtEtaBinsLow, jtEtaBinsHigh, nJtEtaBins, jtEtaBins);

  const Int_t nJtEtaBinsSub = config_p->GetValue("NJTETABINSSUB", 14);
  if(!goodBinning(inConfigFileName, nMaxSubBins, nJtEtaBinsSub, "NJTETABINSSUB")) return 1;
  const Float_t jtEtaBinsSubLow = config_p->GetValue("JTETABINSSUBLOW", 0.0);
  const Float_t jtEtaBinsSubHigh = config_p->GetValue("JTETABINSSUBHIGH", 2.8);
  const Bool_t jtEtaBinsSubDoAbs = config_p->GetValue("JTETABINSSUBDOABS", 1);
  if(jtEtaBinsSubDoAbs && jtEtaBinsSubLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains jtEtaBinsSubLow \'" << jtEtaBinsSubLow << "\' less than 0 despite requested jtEtaBinsSubDoAbs \'" << jtEtaBinsSubDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t jtEtaBinsSub[nMaxSubBins+1];
  const Bool_t jtEtaBinsSubDoLin = config_p->GetValue("JTETABINSSUBDOLIN", 0);
  const Bool_t jtEtaBinsSubDoLog = config_p->GetValue("JTETABINSSUBDOLOG", 0);
  const Bool_t jtEtaBinsSubDoCustom = config_p->GetValue("JTETABINSSUBDOCUSTOM", 0);
  if(jtEtaBinsSubDoLin) getLinBins(jtEtaBinsSubLow, jtEtaBinsSubHigh, nJtEtaBinsSub, jtEtaBinsSub);
  else if(jtEtaBinsSubDoLog){
    if(jtEtaBinsSubLow <= 0){
      std::cout << "JTETABINSSUBDOLOG called despite JTETABINSSUBLOW \'" << jtEtaBinsSubLow << "\' less than or equal  to 0. return 1" << std::endl;
      return 1;
    }
    getLogBins(jtEtaBinsSubLow, jtEtaBinsSubHigh, nJtEtaBinsSub, jtEtaBinsSub);
  }
  else if(jtEtaBinsSubDoCustom){
    if(!checkEnvForParams(config_p, {"JTETABINSSUBCUSTOM"})) return 1;
    std::vector<float> jtEtaBinsSubCustom = strToVectF(config_p->GetValue("JTETABINSSUBCUSTOM", ""));
    if(nJtEtaBinsSub != ((int)jtEtaBinsSubCustom.size()) - 1){
      return 1;
    }

    for(unsigned int jI = 0; jI < jtEtaBinsSubCustom.size(); ++jI){
      jtEtaBinsSub[jI] = jtEtaBinsSubCustom[jI];
    }
  }
  else{
    std::cout << "None of 'JTETABINSSUBDOLIN', 'JTETABINSSUBDOLOG', or 'JTETABINSSUBDOCUSTOM' are specified as true. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> jtEtaBinsSubStr;
  for(Int_t pI = 0; pI < nJtEtaBinsSub; ++pI){
    if(pI < 10) jtEtaBinsSubStr.push_back("JtEta0" + std::to_string(pI));
    else jtEtaBinsSubStr.push_back("JtEta" + std::to_string(pI));

    binsToLabelStr[jtEtaBinsSubStr[pI]] = prettyString(jtEtaBinsSub[pI], 1, false) + " < #eta_{Jet} < " + prettyString(jtEtaBinsSub[pI+1], 1, false);
  }
  if(nJtEtaBinsSub < 10) jtEtaBinsSubStr.push_back("JtEta0" + std::to_string(nJtEtaBinsSub));
  else jtEtaBinsSubStr.push_back("JtEta" + std::to_string(nJtEtaBinsSub));

  binsToLabelStr[jtEtaBinsSubStr[jtEtaBinsSubStr.size()-1]] = prettyString(jtEtaBinsSub[0], 1, false) + " < #eta_{Jet} < " + prettyString(jtEtaBinsSub[nJtEtaBinsSub], 1, false);

  const Int_t nPhiBins = config_p->GetValue("NPHIBINS", 32);
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nPhiBins, "NPHIBINS")) return 1;
  Double_t phiBins[nMaxEtaPhiBins+1];
  getLinBins(-TMath::Pi()+0.01, TMath::Pi()+0.01, nPhiBins, phiBins);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
 
 //Pt sub bins handling
  //  const Int_t nGammaPtBinsSub = config_p->GetValue("NGAMMAPTBINSSUB", 4);
  //  if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
  //  const Float_t gammaPtBinsSubLow = config_p->GetValue("GAMMAPTBINSSUBLOW", 80.0);
  //  const Float_t gammaPtBinsSubHigh = config_p->GetValue("GAMMAPTBINSSUBHIGH", 280.0);
  //  if(gammaPtBinsSubLow < gammaPtBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubLow \'" << gammaPtBinsSubLow << "\' less than gammaPtBinsLow \'" << gammaPtBinsLow << "\'. return 1" << std::endl;
  //  if(gammaPtBinsSubHigh > gammaPtBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubHigh \'" << gammaPtBinsSubHigh << "\' greater than gammaPtBinsHigh \'" << gammaPtBinsHigh << "\'. return 1" << std::endl;
  //  if(gammaPtBinsSubLow < gammaPtBinsLow || gammaPtBinsSubHigh > gammaPtBinsHigh) return 1;
  //  const Bool_t gammaPtBinsSubDoLog = config_p->GetValue("GAMMAPTBINSSUBDOLOG", 1);
  //  Double_t gammaPtBinsSub[nMaxSubBins+1];
  //  if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
  //  else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
  //  std::vector<std::string> gammaPtBinsSubStr;
  //  for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
  //    gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(pI));

  //  binsToLabelStr[gammaPtBinsSubStr[pI]] = prettyString(gammaPtBinsSub[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[pI+1], 1, false);
  //  }
  //  gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(nGammaPtBinsSub));

  //  binsToLabelStr[gammaPtBinsSubStr[gammaPtBinsSubStr.size()-1]] = prettyString(gammaPtBinsSub[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Eta sub bins handling
  const Int_t nGammaEtaBinsSub = config_p->GetValue("NGAMMAETABINSSUB", 4);
  if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaEtaBinsSub, "NGAMMAETABINSSUB")) return 1;
  const Float_t gammaEtaBinsSubLow = config_p->GetValue("GAMMAETABINSSUBLOW", 0.0);
  const Float_t gammaEtaBinsSubHigh = config_p->GetValue("GAMMAETABINSSUBHIGH", 2.37);
  //  if(etaBinsSubLow < etaBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than etaBinsLow \'" << etaBinsLow << "\'. return 1" << std::endl;
  //  if(etaBinsSubHigh > etaBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubHigh \'" << etaBinsSubHigh << "\' greater than etaBinsHigh \'" << etaBinsHigh << "\'. return 1" << std::endl;
  //  if(etaBinsSubLow < etaBinsLow || etaBinsSubHigh > etaBinsHigh) return 1;

  const	Bool_t gammaEtaBinsSubDoAbs = config_p->GetValue("GAMMAETABINSSUBDOABS", 1);
  if(gammaEtaBinsSubDoAbs && gammaEtaBinsSubLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaEtaBinsSubLow \'" << gammaEtaBinsSubLow << "\' less than 0 despite requested gammaEtaBinsSubDoAbs \'" << gammaEtaBinsSubDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t gammaEtaBinsSub[nMaxSubBins+1];
  getLinBins(gammaEtaBinsSubLow, gammaEtaBinsSubHigh, nGammaEtaBinsSub, gammaEtaBinsSub);
  std::vector<std::string> gammaEtaBinsSubStr;
  std::string preStr = "";
  if(gammaEtaBinsSubDoAbs) preStr = "Abs";
  for(Int_t eI = 0; eI < nGammaEtaBinsSub; ++eI){
    gammaEtaBinsSubStr.push_back(preStr + "Eta" + std::to_string(eI));

    if(gammaEtaBinsSubDoAbs) binsToLabelStr[gammaEtaBinsSubStr[eI]] = prettyString(gammaEtaBinsSub[eI], 2, false) + "<|#eta_{#gamma}|<" + prettyString(gammaEtaBinsSub[eI+1], 2, false);
    else binsToLabelStr[gammaEtaBinsSubStr[eI]] = prettyString(gammaEtaBinsSub[eI], 2, false) + " < #eta_{#gamma} < " + prettyString(gammaEtaBinsSub[eI+1], 2, false);
  }
  gammaEtaBinsSubStr.push_back(preStr + "Eta" + std::to_string(nGammaEtaBinsSub));
  if(gammaEtaBinsSubDoAbs) binsToLabelStr[gammaEtaBinsSubStr[gammaEtaBinsSubStr.size()-1]] = prettyString(gammaEtaBinsSub[0], 2, false) + "<|#eta_{#gamma}|<" + prettyString(gammaEtaBinsSub[nGammaEtaBinsSub], 2, false);
  else binsToLabelStr[gammaEtaBinsSubStr[gammaEtaBinsSubStr.size()-1]] = prettyString(gammaEtaBinsSub[0], 2, false) + " < #eta_{#gamma} < " + prettyString(gammaEtaBinsSub[nGammaEtaBinsSub], 2, false);



  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nJtPtBins = config_p->GetValue("NJTPTBINS", 20);
  if(!goodBinning(inConfigFileName, nMaxPtBins, nJtPtBins, "NJTPTBINS")) return 1;
  const Float_t jtPtBinsLow = config_p->GetValue("JTPTBINSLOW", 30.0);
  const Float_t jtPtBinsHigh = config_p->GetValue("JTPTBINSHIGH", 230.0);
  const Bool_t jtPtBinsDoLog = config_p->GetValue("JTPTBINSDOLOG", 1);
  const Bool_t jtPtBinsDoCustom = config_p->GetValue("JTPTBINSDOCUSTOM", 1);
  Double_t jtPtBins[nMaxPtBins+1];
  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else if(jtPtBinsDoCustom){
    std::string jtPtBinsStr = config_p->GetValue("JTPTBINSCUSTOM", "");
    if(jtPtBinsStr.size() == 0){
      std::cout << "CUSTOM JT PT BINS REQUESTED, BUT \'JTPTBINSCUSTOM\' not defined in config \'" << inConfigFileName << "\'. return 1" << std::endl;
    }
    std::vector<float> jtPtBinsCustom = strToVectF(jtPtBinsStr);
    Int_t nJtPtBinsTemp = ((Int_t)jtPtBinsCustom.size()) - 1;
    if(nJtPtBinsTemp != nJtPtBins){
      std::cout << "REQUESTED CUSTON BINS \'" << jtPtBinsStr << "\' HAS SIZE \'" << nJtPtBinsTemp << "\' does not match given nJtPtBins \'" << nJtPtBins << "\'. return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < jtPtBinsCustom.size(); ++gI){
      jtPtBins[gI] = jtPtBinsCustom[gI];
    }
  }
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);

  std::string jtPtBinsGlobalStr = "GlobalJtPt0";
  std::string jtPtBinsGlobalLabel = prettyString(jtPtBinsLow,1,false) + " < p_{T,jet} < " + prettyString(jtPtBinsHigh,1,false);
  binsToLabelStr[jtPtBinsGlobalStr] = jtPtBinsGlobalLabel;   

  std::string multiJtCutGlobalStr = "MultiJt0";
  std::string multiJtCutGlobalLabel = "N_{Jet,Reco.} >= 2";
  binsToLabelStr[multiJtCutGlobalStr] = multiJtCutGlobalLabel;   

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<std::string> jtPtBinsStr;
  for(Int_t pI = 0; pI < nJtPtBins; ++pI){
    if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
    else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));

    binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
  }
  jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
  binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nDPhiBins = config_p->GetValue("NDPHIBINS", 16);
  const Float_t dPhiBinsLow = config_p->GetValue("DPHIBINSLOW", 0.0);
  const Float_t dPhiBinsHigh = config_p->GetValue("DPHIBINSHIGH", TMath::Pi());
  const Double_t gammaJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string gammaJtDPhiStr = "DPhi0";
  std::string gammaJtDPhiLabel = returnAllCapsString(config_p->GetValue("GAMMAJTDPHI", ""));
  if(gammaJtDPhiLabel.find("PI") != std::string::npos) gammaJtDPhiLabel.replace(gammaJtDPhiLabel.find("PI"), 2, "#pi");
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  gammaJtDPhiLabel = "|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiLabel;
  binsToLabelStr[gammaJtDPhiStr] = gammaJtDPhiLabel;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nXJBins = config_p->GetValue("NXJBINS", 20);
  const Float_t xjBinsLow = config_p->GetValue("XJBINSLOW", 0.0);
  const Float_t xjBinsHigh = config_p->GetValue("XJBINSHIGH", 2.0);
  const Bool_t xjBinsDoCustom = config_p->GetValue("XJBINSDOCUSTOM", 0);
  std::string xjBinsCustomStr = config_p->GetValue("XJBINSCUSTOM", "");
  Double_t xjBins[nMaxPtBins+1];
  if(!xjBinsDoCustom) getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
  else{
    if(xjBinsCustomStr.size() == 0){
      std::cout << "No custom bins supplied despite XJBINSDOCUSOM \'" << xjBinsDoCustom <<  "'\". Reverting to uniform binning." << std::endl;
      getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
    }
    else{
      std::vector<float> xjBinsCustomVect = strToVectF(xjBinsCustomStr);
      //Check that the bins are ordered and the number matches nXJBins
      Int_t nXJBinsCustomVect = xjBinsCustomVect.size();
      if(nXJBinsCustomVect-1 != nXJBins){
	std::cout << "Number of supplied NXJBINS, \'" << nXJBins << "\', does not match number from custom bin string \'" << nXJBinsCustomVect-1 << "\'. Reverting to uniform bins" << std::endl;
        getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
      }
      else{
        bool areBinsOrdered = true;
        for(unsigned int i = 1; i < xjBinsCustomVect.size(); ++i){
          if(xjBinsCustomVect[i] < xjBinsCustomVect[i-1]){
            areBinsOrdered = false;
            break;
          }
        }

        if(areBinsOrdered){
          for(unsigned int i = 0; i < xjBinsCustomVect.size(); ++i){
            xjBins[i] = xjBinsCustomVect[i];
          }
        }
        else{
	  std::cout << "XJBINSCUSTOM provided \'" << xjBinsCustomStr << "\'is not ordered. Reverting to uniform bins"  << std::endl;
          getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
        }
      }
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const unsigned long long mixCap = config_p->GetValue("MIXCAP", 100000000);
  const unsigned long long nMixEvents = config_p->GetValue("NMIXEVENTS", 1);
  const bool doStrictMix = config_p->GetValue("DOSTRICTMIX", true);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  Float_t recoGammaPt_, truthGammaPt_;
  Float_t recoXJ_, truthXJ_;
  Float_t recoJtPt_, truthJtPt_;
  Float_t recoXJJ_, truthXJJ_;
  Float_t recoJt1Pt_, truthJt1Pt_;
  Float_t recoJt2Pt_, truthJt2Pt_;
  Float_t unfoldWeight_;
  Float_t unfoldCent_;
  
  TTree* unfoldXJTree_p = nullptr;
  TTree* unfoldXJJTree_p = nullptr;
  if(isMC){
    unfoldXJTree_p = new TTree("unfoldXJTree_p", "");
    unfoldXJTree_p->Branch("recoGammaPt", &recoGammaPt_, "recoGammaPt/F"); 
    unfoldXJTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F"); 
    unfoldXJTree_p->Branch("recoXJ", &recoXJ_, "recoXJ/F"); 
    unfoldXJTree_p->Branch("truthXJ", &truthXJ_, "truthXJ/F"); 
    unfoldXJTree_p->Branch("recoJtPt", &recoJtPt_, "recoJtPt/F"); 
    unfoldXJTree_p->Branch("truthJtPt", &truthJtPt_, "truthJtPt/F"); 
    unfoldXJTree_p->Branch("unfoldWeight", &unfoldWeight_, "unfoldWeight/F"); 
    unfoldXJTree_p->Branch("unfoldCent", &unfoldCent_, "unfoldCent/F"); 

    unfoldXJJTree_p = new TTree("unfoldXJJTree_p", "");
    unfoldXJJTree_p->Branch("recoGammaPt", &recoGammaPt_, "recoGammaPt/F"); 
    unfoldXJJTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F"); 
    unfoldXJJTree_p->Branch("recoXJJ", &recoXJJ_, "recoXJJ/F"); 
    unfoldXJJTree_p->Branch("truthXJJ", &truthXJJ_, "truthXJJ/F"); 
    unfoldXJJTree_p->Branch("recoJt1Pt", &recoJt1Pt_, "recoJt1Pt/F"); 
    unfoldXJJTree_p->Branch("truthJt1Pt", &truthJt1Pt_, "truthJt1Pt/F"); 
    unfoldXJJTree_p->Branch("recoJt2Pt", &recoJt2Pt_, "recoJt2Pt/F"); 
    unfoldXJJTree_p->Branch("truthJt2Pt", &truthJt2Pt_, "truthJt2Pt/F"); 
    unfoldXJJTree_p->Branch("unfoldWeight", &unfoldWeight_, "unfoldWeight/F"); 
    unfoldXJJTree_p->Branch("unfoldCent", &unfoldCent_, "unfoldCent/F"); 
  }

  TH1F* mixingCentrality_p = nullptr;
  TH1F* mixingPsi2_p = nullptr;
  TH1F* mixingVz_p = nullptr;
  TH2F* mixingCentralityPsi2_p = nullptr;
  TH2F* mixingCentralityVz_p = nullptr;
  TH2F* mixingPsi2Vz_p = nullptr;
  TH2F* mixingCentralityPsi2_VzBinned_p[nMaxMixBins];

  if(doMix){
    if(doMixCent){
      mixingCentrality_p = new TH1F("mixingCentrality_h", ";Centrality;Counts", nMixCentBins, mixCentBins);
      if(doMixPsi2){
	mixingCentralityPsi2_p = new TH2F("mixingCentralityPsi2_h", ";Centrality;#Psi_{2}", nMixCentBins, mixCentBins, nMixPsi2Bins, mixPsi2Bins);

	if(doMixVz){
	  for(Int_t vI = 0; vI < nMixVzBins; ++vI){
	    mixingCentralityPsi2_VzBinned_p[vI] = new TH2F(("mixingCentralityPsi2_Vz" + std::to_string(vI) + "_h").c_str(), ";;", nMixCentBins, mixCentBins, nMixPsi2Bins, mixPsi2Bins);
	  }
	}
      }

      if(doMixVz) mixingCentralityVz_p = new TH2F("mixingCentralityVz_h", ";Centrality;v_{z}", nMixCentBins, mixCentBins, nMixVzBins, mixVzBins);
    }
    if(doMixPsi2){
      mixingPsi2_p = new TH1F("mixingPsi2_h", ";#Psi_{2};Counts", nMixPsi2Bins, mixPsi2Bins);
      if(doMixVz) mixingPsi2Vz_p = new TH2F("mixingPsi2Vz_h", ";#Psi_{2};v_{z}", nMixPsi2Bins, mixPsi2Bins, nMixVzBins, mixVzBins);
    }
    if(doMixVz) mixingVz_p = new TH1F("mixingVz_h", ";v_{z};Counts", nMixVzBins, mixVzBins);
  }

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
  TH1F* photonMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded
  TH2F* photonPtJtPtVCent_RAW_p[nMaxCentBins];
  TH2F* photonPtJtXJVCent_RAW_p[nMaxCentBins];
  TH2F* photonPtJtXJJVCent_RAW_p[nMaxCentBins];


  TH1F* photonJtRecoOverGenVCentJtPt_p[nMaxCentBins][nMaxPtBins];
  TH2F* photonJtCorrOverUncorrVCentJtEta_p[nMaxCentBins][nMaxSubBins];
  TH2F* photonJtGenResVCentGenPtRecoPt_p[nMaxCentBins][nMaxPtBins][nMaxPtBins];
  
  TH1F* photonMixJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixCorrectionMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixCorrectedMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixCorrectionMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixCorrectedMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonMixJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded
  TH2F* photonPtJtPtVCent_MIX_p[nMaxCentBins];
  TH2F* photonPtJtXJVCent_MIX_p[nMaxCentBins];
  TH2F* photonPtJtXJJVCent_MIX_p[nMaxCentBins];
  TH2F* photonPtJtXJJVCent_MIXCorrected_p[nMaxCentBins];
  TH2F* photonPtJtXJJVCent_MIXCorrection_p[nMaxCentBins];

  TH1F* photonSubJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonSubJtMultModVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded
  TH2F* photonPtJtPtVCent_SUB_p[nMaxCentBins];
  TH2F* photonPtJtXJVCent_SUB_p[nMaxCentBins];
  TH2F* photonPtJtXJJVCent_SUB_p[nMaxCentBins];

  TH1F* photonGenJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1F* photonGenMatchedJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonGenMatchedMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
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

    for(Int_t eI = 0; eI < nGammaEtaBinsSub+1; ++eI){
      photonPtVCentEta_p[cI][eI] = new TH1F(("photonPtVCentEta_" + centBinsStr[cI] + "_" + gammaEtaBinsSubStr[eI]+ "_h").c_str(), ";#gamma p_{T} [GeV];Counts", nGammaPtBins, gammaPtBins);
      centerTitles(photonPtVCentEta_p[cI][eI]);

      if(isMC){
	photonGenResVCentEta_p[cI][eI] = new TH2F(("photonGenResVCentEta_" + centBinsStr[cI] + "_" + gammaEtaBinsSubStr[eI] + "_h").c_str(), ";Reco. Photon p_{T};Gen. Photon p_{T}", nGammaPtBins, gammaPtBins, nGammaPtBins, gammaPtBins);
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonEtaVCentPt_p[cI][pI] = new TH1F(("photonEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #eta;Counts", nGammaEtaBins, gammaEtaBins);
      photonPhiVCentPt_p[cI][pI] = new TH1F(("photonPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #phi;Counts", nPhiBins, phiBins);

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";#Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      photonJtPtVCentPt_p[cI][pI] = new TH1F(("photonJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
      photonMultiJtPtVCentPt_p[cI][pI] = new TH1F(("photonMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
      photonJtEtaVCentPt_p[cI][pI] = new TH1F(("photonJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
      photonJtXJVCentPt_p[cI][pI] = new TH1F(("photonJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
      photonMultiJtXJVCentPt_p[cI][pI] = new TH1F(("photonMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

      photonMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

      photonMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

      photonJtMultVCentPt_p[cI][pI] = new TH1F(("photonJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

      setSumW2({photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonMultiJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI], photonMultiJtXJVCentPt_p[cI][pI], photonMultiJtXJJVCentPt_p[cI][pI], photonMultiJtDPhiJJVCentPt_p[cI][pI]});
  
      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonMixJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Mixed event #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonMixJtPtVCentPt_p[cI][pI] = new TH1F(("photonMixJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonMixMultiJtPtVCentPt_p[cI][pI] = new TH1F(("photonMixMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonMixJtEtaVCentPt_p[cI][pI] = new TH1F(("photonMixJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonMixJtXJVCentPt_p[cI][pI] = new TH1F(("photonMixJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonMixMultiJtXJVCentPt_p[cI][pI] = new TH1F(("photonMixMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

	photonMixMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonMixMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

	photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonMixCorrectionMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event correction #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

	photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonMixCorrectedMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event corrected #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);


	photonMixMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonMixMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonMixCorrectionMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event correction #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonMixCorrectedMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event corrected #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixJtMultVCentPt_p[cI][pI] = new TH1F(("photonMixJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Mixed Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	setSumW2({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixMultiJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], photonMixMultiJtXJVCentPt_p[cI][pI], photonMixMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI], photonMixMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI]});
      }

      if(doMix || isPP){
	photonSubJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonSubJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Subtracted event #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonSubJtPtVCentPt_p[cI][pI] = new TH1F(("photonSubJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonSubMultiJtPtVCentPt_p[cI][pI] = new TH1F(("photonSubMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);

	photonSubJtEtaVCentPt_p[cI][pI] = new TH1F(("photonSubJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonSubJtXJVCentPt_p[cI][pI] = new TH1F(("photonSubJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonSubMultiJtXJVCentPt_p[cI][pI] = new TH1F(("photonSubMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonSubMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonSubMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonSubMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonSubJtMultVCentPt_p[cI][pI] = new TH1F(("photonSubJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Subtracted Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	photonSubJtMultModVCentPt_p[cI][pI] = new TH1F(("photonSubJtMultModVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Subtracted Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);
      
	setSumW2({photonSubJtDPhiVCentPt_p[cI][pI], photonSubJtPtVCentPt_p[cI][pI], photonSubMultiJtPtVCentPt_p[cI][pI], photonSubJtEtaVCentPt_p[cI][pI], photonSubJtXJVCentPt_p[cI][pI], photonSubMultiJtXJVCentPt_p[cI][pI], photonSubMultiJtXJJVCentPt_p[cI][pI], photonSubMultiJtDPhiJJVCentPt_p[cI][pI], photonSubJtMultVCentPt_p[cI][pI], photonSubJtMultModVCentPt_p[cI][pI]});
      }
  
      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonGenJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Gen. #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonGenJtPtVCentPt_p[cI][pI] = new TH1F(("photonGenJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonGenJtEtaVCentPt_p[cI][pI] = new TH1F(("photonGenJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonGenJtXJVCentPt_p[cI][pI] = new TH1F(("photonGenJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonGenMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen. #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonGenMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen. #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonGenJtMultVCentPt_p[cI][pI] = new TH1F(("photonGenJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen. Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Gen.-matched #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonGenMatchedJtPtVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonGenMatchedJtEtaVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonGenMatchedJtXJVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMatchedMultiJtXJJVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen.-matched #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen.-matched #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonGenMatchedJtMultVCentPt_p[cI][pI] = new TH1F(("photonGenMatchedJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen.-matched Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);
      
	


	photonJtFakeVCentPt_p[cI][pI] = new TH1F(("photonJtFakeVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} ;#frac{N_{Fake jets}}{N_{All jets}}", nJtPtBins, jtPtBins);
	
      	setSumW2({photonGenJtDPhiVCentPt_p[cI][pI], photonGenJtPtVCentPt_p[cI][pI], photonGenJtEtaVCentPt_p[cI][pI], photonGenJtXJVCentPt_p[cI][pI], photonGenMultiJtXJJVCentPt_p[cI][pI], photonGenMultiJtDPhiJJVCentPt_p[cI][pI], photonGenJtMultVCentPt_p[cI][pI], photonGenMatchedJtDPhiVCentPt_p[cI][pI], photonGenMatchedJtPtVCentPt_p[cI][pI], photonGenMatchedJtEtaVCentPt_p[cI][pI], photonGenMatchedJtXJVCentPt_p[cI][pI], photonGenMatchedMultiJtXJJVCentPt_p[cI][pI], photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI], photonGenMatchedJtMultVCentPt_p[cI][pI], photonJtFakeVCentPt_p[cI][pI]});
      }
    }
	  
    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2] = new TH2F(("photonJtGenResVCentGenPtRecoPt_" + centBinsStr[cI] + "_" + genGammaPtBinsStr[ptI] + "_" + recoGammaPtBinsStr[ptI2] + "_" + gammaJtDPhiStr + "_h").c_str(), "", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	}
      }

      for(Int_t gI = 0; gI < nJtPtBins; ++gI){
	photonJtRecoOverGenVCentJtPt_p[cI][gI] = new TH1F(("photonJtRecoOverGenVCentJtPt_" + centBinsStr[cI] + "_" + jtPtBinsStr[gI] + "_" + gammaPtBinsStr[nGammaPtBins] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Reco./Gen.;Counts (Weighted)", 51, 0, 2.0);

	setSumW2(photonJtRecoOverGenVCentJtPt_p[cI][gI]);
      }
    }
  
    for(Int_t gI = 0; gI < nJtEtaBinsSub+1; ++gI){
      photonJtCorrOverUncorrVCentJtEta_p[cI][gI] = new TH2F(("photonJtCorrOverUncorrVCentJtEta_" + centBinsStr[cI] + "_" + jtEtaBinsSubStr[gI] + "_" + gammaPtBinsStr[nGammaPtBins] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Uncorrected Jet p_{T} [GeV];Corrected/Uncorrected", nJtPtBins, jtPtBins, 150, 0.5, 2.0);
      setSumW2(photonJtCorrOverUncorrVCentJtEta_p[cI][gI]);
    }      
  
    photonPtJtPtVCent_RAW_p[cI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
    photonPtJtXJVCent_RAW_p[cI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
    photonPtJtXJJVCent_RAW_p[cI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

    setSumW2({photonPtJtPtVCent_RAW_p[cI], photonPtJtXJVCent_RAW_p[cI], photonPtJtXJJVCent_RAW_p[cI]});

    if(doMix){
      photonPtJtPtVCent_MIX_p[cI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJVCent_MIX_p[cI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_MIX_p[cI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_MIXCorrection_p[cI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_MIXCorrection_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_MIXCorrected_p[cI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_MIXCorrected_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

      setSumW2({photonPtJtPtVCent_MIX_p[cI], photonPtJtXJVCent_MIX_p[cI], photonPtJtXJJVCent_MIX_p[cI], photonPtJtXJJVCent_MIXCorrection_p[cI], photonPtJtXJJVCent_MIXCorrected_p[cI]});
    }

    if(doMix || isPP){
      photonPtJtPtVCent_SUB_p[cI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJVCent_SUB_p[cI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_SUB_p[cI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

      setSumW2({photonPtJtPtVCent_SUB_p[cI], photonPtJtXJVCent_SUB_p[cI], photonPtJtXJJVCent_SUB_p[cI]});
    }

    photonEtaPt_p[cI] = new TH2F(("photonEtaPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;#gamma p_{T} [GeV]", nGammaEtaBins, gammaEtaBins, nGammaPtBins, gammaPtBins);

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI] = new TH2F(("photonEtaPhi_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #eta;#gamma #phi", nGammaEtaBins, gammaEtaBins, nPhiBins, phiBins);
    }
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
  const Double_t mmToCMDivFactor = 10.;

  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
  
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;  
  std::vector<float>* photon_etcone30_p=nullptr;
  
  std::vector<float>* aktRhi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_uncorrpt_p=nullptr;
  std::vector<float>* aktRhi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_uncorreta_p=nullptr;
  std::vector<float>* aktRhi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_phi_p=nullptr;
  std::vector<int>* aktRhi_truthpos_p=nullptr;

  std::vector<float>* aktR_truth_jet_pt_p=nullptr;
  std::vector<float>* aktR_truth_jet_eta_p=nullptr;
  std::vector<float>* aktR_truth_jet_phi_p=nullptr;
  std::vector<int>* aktR_truth_jet_recopos_p=nullptr;

  TFile* mixFile_p = nullptr;
  TTree* mixTree_p = nullptr;
  if(doMix){
    inFile_p->cd();
    inTree_p->SetBranchStatus("*", 0);

    if(doMixVz){
      inTree_p->SetBranchStatus("vert_z", 1);
      inTree_p->SetBranchAddress("vert_z", &vert_z_p);
    }

    if(doMixCent){
      inTree_p->SetBranchStatus("fcalA_et", 1);
      inTree_p->SetBranchStatus("fcalC_et", 1);
      
      inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
      inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
    }
    
    if(doMixPsi2){
      inTree_p->SetBranchStatus("evtPlane2Phi", 1);
      inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
    }

    ULong64_t nEntriesTemp = inTree_p->GetEntries();
    if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, (ULong64_t)nMaxEvt);
    const ULong64_t nSigEntries = nEntriesTemp;

    for(ULong64_t entry = 0; entry < nSigEntries; ++entry){
      inTree_p->GetEntry(entry);
    
      double vert_z = vert_z_p->at(0);
      vert_z /= mmToCMDivFactor;
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
         
      ++(signalMapCounterPre[key]);
    }    

    mixFile_p = new TFile(inMixFileName.c_str(), "READ");
    mixTree_p = (TTree*)mixFile_p->Get("gammaJetTree_p");

    mixTree_p->SetBranchStatus("*", 0);
    mixTree_p->SetBranchStatus("vert_z", 1);

    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);

    mixTree_p->SetBranchAddress("vert_z", &vert_z_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);

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

    nEntriesTemp = mixTree_p->GetEntries();
    if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, (ULong64_t)nMaxEvt*10);
    const ULong64_t nMixEntries = nEntriesTemp;

    for(ULong64_t entry = 0; entry < nMixEntries; ++entry){
      mixTree_p->GetEntry(entry);

      double vert_z = vert_z_p->at(0);
      vert_z /= mmToCMDivFactor;
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
      
      if(mixingMapCounter[key] >= mixCap) continue;

      
      std::vector<TLorentzVector> jets;
      for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
	if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;
	if(aktRhi_em_xcalib_jet_eta_p->at(jI) <= jtEtaBinsLow) continue;
	if(aktRhi_em_xcalib_jet_eta_p->at(jI) >= jtEtaBinsHigh) continue;

	TLorentzVector temp;
	temp.SetPtEtaPhiM(aktRhi_em_xcalib_jet_pt_p->at(jI), aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), 0.0);

	jets.push_back(temp);
      }
    
      if(doMixCent){
	mixingCentrality_p->Fill(cent);
	if(doMixPsi2){
	  mixingCentralityPsi2_p->Fill(cent, evtPlane2Phi);

	  if(doMixVz) mixingCentralityPsi2_VzBinned_p[vzPos]->Fill(cent, evtPlane2Phi);
	}
	if(doMixVz) mixingCentralityVz_p->Fill(cent, vert_z);
      }
      if(doMixPsi2){
	mixingPsi2_p->Fill(evtPlane2Phi);
	if(doMixVz) mixingPsi2Vz_p->Fill(evtPlane2Phi, vert_z);
      }
      if(doMixVz) mixingVz_p->Fill(vert_z);
    
      mixingMap[key].push_back(jets);
      ++(mixingMapCounter[key]);
      //      ++(signalMapCounter[key]);   
    }
    
    mixFile_p->Close();
    delete mixFile_p;
    inFile_p->cd();

    //Gonna dump out some statements on the statistics of the mixed events 
    unsigned long long maxKey = 0;
    unsigned long long maximumVal = 0;

    unsigned long long minKey = 0;
    unsigned long long minimumVal = 2147483647; //hardcoded max value of a signed integer

    double aveVal = 0;
    double tempCounter = 0;

    std::vector<unsigned long long> tempCounts;

    for(auto const & mixes : mixingMapCounter){
      
      if(doStrictMix){
	if(mixes.second < signalMapCounterPre[mixes.first]*nMixEvents){	
	  std::cout << "Mixing has less events than signal (" << signalMapCounterPre[mixes.first] << ") times nMixEvens (" << nMixEvents << "): " << mixes.second << "<" << signalMapCounterPre[mixes.first] << " * " << nMixEvents << "." << std::endl;
	  std::cout << "Key is: " << mixes.first << std::endl;
	  std::cout << "Key bin is: " << keyBoy.GetKeyStr(mixes.first) << std::endl;
	  std::cout << "return 1" << std::endl;
	  
	  inFile_p->Close();
	  delete inFile_p;
	  return 1;
	}
      }

      if(mixes.second < minimumVal){
	minKey = mixes.first;
	minimumVal = mixes.second;
      }
      else if(mixes.second > maximumVal){
	maxKey = mixes.first;
	maximumVal = mixes.second;
      }
      
      aveVal += mixes.second;
      tempCounts.push_back(mixes.second);
      ++tempCounter;
    }

    std::sort(std::begin(tempCounts), std::end(tempCounts));
  
    std::cout << "MINIMUM NUMBER TO MIX, CORRESPONDING KEY: " << minimumVal << " (" << keyBoy.GetKeyStr(minKey) << "), key=" << minKey << std::endl;
    std::cout << "MAXIMUM NUMBER TO MIX, CORRESPONDING KEY: " << maximumVal << " (" << keyBoy.GetKeyStr(maxKey) << "), key=" << maxKey << std::endl;
    std::cout << "AVERAGE NUMBER TO MIX: " << aveVal/tempCounter << std::endl;
    std::cout << "MEDIAN NUMBER TO MIX: " << tempCounts[tempCounts.size()/2] << std::endl;
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
    if(!isPP) inTree_p->SetBranchStatus("ncollWeight", 1);
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
  inTree_p->SetBranchStatus("photon_tight", 1);
  inTree_p->SetBranchStatus("photon_etcone30", 1);
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorrpt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorreta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);
  
  if(isMC){
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_truthpos").c_str(), 1);
    
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_pt").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_eta").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_phi").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_recopos").c_str(), 1);
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
    if(!isPP) inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
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
  inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
  inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);

  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorrpt").c_str(), &aktRhi_em_xcalib_jet_uncorrpt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_pt").c_str(), &aktRhi_constit_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorreta").c_str(), &aktRhi_em_xcalib_jet_uncorreta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_eta").c_str(), &aktRhi_constit_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);

  if(isMC){
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_truthpos").c_str(), &aktRhi_truthpos_p);    
    
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_pt").c_str(), &aktR_truth_jet_pt_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_eta").c_str(), &aktR_truth_jet_eta_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_phi").c_str(), &aktR_truth_jet_phi_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_recopos").c_str(), &aktR_truth_jet_recopos_p);
  }


  std::map<int, int> perEventPhotonCounts;
  for(unsigned int i = 0; i < 10; ++i){
    perEventPhotonCounts[i] = 0;
  }


  Double_t recoJtPtMin = 100000.;
  
  ULong64_t nEntriesTemp = inTree_p->GetEntries();
  if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, nMaxEvt);
  const ULong64_t nEntries = nEntriesTemp;
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::vector<std::vector<Double_t> > gammaCountsPerPtCent;
  for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
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
    vert_z /= mmToCMDivFactor;
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

    if(!isMC) fullWeight = 1.0;
  
    Float_t mixWeight = fullWeight/(double)nMixEvents;

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
  
    //First loop for the corrections plots
    for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_uncorrpt_p->size(); ++jI){
      double jtEtaForBin = aktRhi_em_xcalib_jet_uncorreta_p->at(jI);
      if(jtEtaBinsSubDoAbs) jtEtaForBin = TMath::Abs(jtEtaForBin);
      
      if(jtEtaForBin < jtEtaBinsSubLow) continue;
      if(jtEtaForBin >= jtEtaBinsSubHigh) continue;
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
      
      if(aktRhi_em_xcalib_jet_uncorrpt_p->at(jI) <= jtPtBinsLow) continue;
      if(aktRhi_em_xcalib_jet_uncorrpt_p->at(jI) > jtPtBinsHigh) continue;
      
      int etaPos = ghostPos(nJtEtaBinsSub, jtEtaBinsSub, jtEtaForBin);
      fillTH2(photonJtCorrOverUncorrVCentJtEta_p[centPos][etaPos], aktRhi_em_xcalib_jet_uncorrpt_p->at(jI), aktRhi_em_xcalib_jet_pt_p->at(jI)/aktRhi_em_xcalib_jet_uncorrpt_p->at(jI), fullWeight);	  
      fillTH2(photonJtCorrOverUncorrVCentJtEta_p[centPos][nJtEtaBinsSub], aktRhi_em_xcalib_jet_uncorrpt_p->at(jI), aktRhi_em_xcalib_jet_pt_p->at(jI)/aktRhi_em_xcalib_jet_uncorrpt_p->at(jI), fullWeight);	  
    }


    Int_t nGoodPhotons = 0;
    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      if(!isGoodPhoton(isPP, photon_tight_p->at(pI), photon_etcone30_p->at(pI), photon_eta_p->at(pI))) continue;
      //      if(!photon_tight_p->at(pI)) continue;
      // above now handled with photonutil.h
      if(photon_pt_p->at(pI) < gammaPtBins[0]) continue;
      if(photon_pt_p->at(pI) >= gammaPtBins[nGammaPtBins]) continue;


      //Isolation as taken from internal note of 2015 data analysis      
      //now handled by photonutil.h
      /*
      if(!isPP){
	if(photon_etcone30_p->at(pI) > 8.0) continue;
      }
      else{
	if(photon_etcone30_p->at(pI) > 3.0) continue;
      }
      */
      //Gap now handled by photonUtil.h
      //Float_t gammaAbsEta = TMath::Abs(photon_eta_p->at(pI));
      //if(gammaAbsEta >= 1.37 && gammaAbsEta < 1.52) continue;
      
      Float_t etaValMain = photon_eta_p->at(pI);
      Float_t etaValSub = etaValMain;
      if(gammaEtaBinsDoAbs) etaValMain = TMath::Abs(etaValMain);
      if(gammaEtaBinsSubDoAbs) etaValSub = TMath::Abs(etaValSub);

      if(etaValMain <= gammaEtaBins[0]) continue;
      if(etaValMain >= gammaEtaBins[nGammaEtaBins]) continue;

      ++nGoodPhotons;
      
      Int_t ptPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(pI), true, doGlobalDebug);

      Int_t etaPos = ghostPos(nGammaEtaBinsSub, gammaEtaBinsSub, etaValSub, true, doGlobalDebug);

      if(etaPos >= 0){
	fillTH1(photonPtVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), fullWeight);
	fillTH1(photonPtVCentEta_p[centPos][nGammaEtaBinsSub], photon_pt_p->at(pI), fullWeight);


	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	if(isMC && truthPhotonPt > 0){
	  if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) < 0.2){
	    fillTH2(photonGenResVCentEta_p[centPos][etaPos], photon_pt_p->at(pI), truthPhotonPt, fullWeight);
	    fillTH2(photonGenResVCentEta_p[centPos][nGammaEtaBinsSub], photon_pt_p->at(pI), truthPhotonPt, fullWeight);
	  }
	}

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
      }
      
      if(ptPos >= 0){

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	if(!isMC){
	  ++(gammaCountsPerPtCent[ptPos][centPos]);
	  ++(gammaCountsPerPtCent[nGammaPtBins][centPos]);
	}
	else{
	  gammaCountsPerPtCent[ptPos][centPos] += fullWeight;
	  gammaCountsPerPtCent[nGammaPtBins][centPos] += fullWeight;
	}
	
	fillTH1(photonEtaVCentPt_p[centPos][ptPos], etaValMain, fullWeight);
	fillTH1(photonPhiVCentPt_p[centPos][ptPos], photon_phi_p->at(pI), fullWeight);
	fillTH2(photonEtaPhiVCentPt_p[centPos][ptPos], etaValMain, photon_phi_p->at(pI), fullWeight);
	
      	fillTH1(photonEtaVCentPt_p[centPos][nGammaPtBins], etaValMain, fullWeight);
	fillTH1(photonPhiVCentPt_p[centPos][nGammaPtBins], photon_phi_p->at(pI), fullWeight);
	fillTH2(photonEtaPhiVCentPt_p[centPos][nGammaPtBins], etaValMain, photon_phi_p->at(pI), fullWeight);
	

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	int multCounter = 0;
	int multCounterGen = 0;
	int multCounterGenMatched = 0;	
      
	std::vector<TLorentzVector> goodJets, goodJetsDPhi;
	std::vector<int> goodJetsTPos, goodJetsDPhiTPos;
	std::vector<std::vector<TLorentzVector> > goodJetsMix, goodJetsDPhiMix;
	goodJetsMix.push_back({});
	goodJetsMix.push_back({});
	goodJetsDPhiMix.push_back({});
	goodJetsDPhiMix.push_back({});

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	//Start w/ unpartnered truth
	if(isMC){
	  for(unsigned int jI = 0; jI < aktR_truth_jet_pt_p->size(); ++jI){
	    if(aktR_truth_jet_recopos_p->at(jI) >= 0) continue;
	    
	    recoGammaPt_ = photon_pt_p->at(pI);
	    truthGammaPt_ = truthPhotonPt;
	    
	    recoXJ_ = -1;
	    recoJtPt_ = -1;
	    
	    truthXJ_ = aktR_truth_jet_pt_p->at(jI)/truthPhotonPt;
	    truthJtPt_ = aktR_truth_jet_pt_p->at(jI);
	    
	    unfoldWeight_ = fullWeight;
	    unfoldCent_ = cent;
	    
	    unfoldXJTree_p->Fill();
	  }
	}

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;           

	for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
	  if(aktRhi_em_xcalib_jet_eta_p->at(jI) <= jtEtaBinsLow) continue;
	  if(aktRhi_em_xcalib_jet_eta_p->at(jI) >= jtEtaBinsHigh) continue;

	  Float_t dR = getDR(aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(pI), photon_phi_p->at(pI));
	  if(dR < gammaExclusionDR) continue;

	  if(recoJtPtMin > aktRhi_em_xcalib_jet_pt_p->at(jI)) recoJtPtMin = aktRhi_em_xcalib_jet_pt_p->at(jI);
	  
	  if(isMC){
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	    int pos = aktRhi_truthpos_p->at(jI);
	    if(pos >= 0){
	      if(aktR_truth_jet_pt_p->at(pos) >= jtPtBinsLow && aktR_truth_jet_pt_p->at(pos) < jtPtBinsHigh){
		Int_t genJtPtPos = ghostPos(nJtPtBins, jtPtBins, aktR_truth_jet_pt_p->at(pos), true, doGlobalDebug);		
		fillTH1(photonJtRecoOverGenVCentJtPt_p[centPos][genJtPtPos], aktRhi_em_xcalib_jet_pt_p->at(jI)/aktR_truth_jet_pt_p->at(pos), fullWeight);
	      }
	    }
	  }
	  
	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	  
	  if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	  if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	  TLorentzVector tempJet;
	  tempJet.SetPtEtaPhiM(aktRhi_em_xcalib_jet_pt_p->at(jI), aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), 0.0);
	  goodJets.push_back(tempJet);
	  goodJetsTPos.push_back(-1);
	  if(isMC){
	    if(aktRhi_truthpos_p->at(jI) >= 0){
	    
	      if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
		goodJetsTPos[goodJetsTPos.size()-1] = aktRhi_truthpos_p->at(jI);
	      }
	      else goodJetsTPos[goodJetsTPos.size()-1] = -1;
	    }
	    else goodJetsTPos[goodJetsTPos.size()-1] = -1;
	  }

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	  Float_t dPhi = TMath::Abs(getDPHI(aktRhi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(pI)));
	  
	  fillTH1(photonJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	  fillTH1(photonJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);

	  if(isMC){
	    if(aktRhi_truthpos_p->at(jI) >= 0){
	      if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
		fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
		fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);
	      }
	    }
	  }

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	  
	  if(dPhi >= gammaJtDPhiCut){
	    goodJetsDPhi.push_back(tempJet);
	    goodJetsDPhiTPos.push_back(-1);
	    if(isMC){
	      if(aktRhi_truthpos_p->at(jI) >= 0){
		if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
		  goodJetsDPhiTPos[goodJetsDPhiTPos.size()-1] = aktRhi_truthpos_p->at(jI);
		}
		else goodJetsDPhiTPos[goodJetsDPhiTPos.size()-1] = -1;
	      }
	      else goodJetsDPhiTPos[goodJetsDPhiTPos.size()-1] = -1;
	    }

	    fillTH1(photonJtPtVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtPtVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtEtaVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_eta_p->at(jI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
	    fillTH1(photonJtXJVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
    
	    if(ptPos == 0 && aktRhi_em_xcalib_jet_pt_p->at(jI) > 300){
	      std::cout << std::endl;
	      std::cout << "WEIRD ASS EVENT ALERT (Entry=" << entry << ":" << std::endl;
	      std::cout << " gamma pt eta phi: " << photon_pt_p->at(pI) << ", " << photon_eta_p->at(pI) << ", " << photon_phi_p->at(pI) << std::endl;
	      std::cout << " jet pt eta phi: " << aktRhi_em_xcalib_jet_pt_p->at(jI) << ", " << aktRhi_em_xcalib_jet_eta_p->at(jI) << ", " << aktRhi_em_xcalib_jet_phi_p->at(jI) << std::endl;
	    }

	    fillTH2(photonPtJtPtVCent_RAW_p[centPos], aktRhi_em_xcalib_jet_pt_p->at(jI), photon_pt_p->at(pI), fullWeight);
	    fillTH2(photonPtJtXJVCent_RAW_p[centPos], aktRhi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), photon_pt_p->at(pI), fullWeight);
	     
	    if(isMC && truthPhotonPt > 0){
	      recoGammaPt_ = photon_pt_p->at(pI);
	      recoXJ_ = aktRhi_em_xcalib_jet_pt_p->at(jI)/recoGammaPt_;
	      recoJtPt_ = aktRhi_em_xcalib_jet_pt_p->at(jI);
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
	      if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) < 0.2){
		if(aktRhi_truthpos_p->at(jI) >= 0){

		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		  Float_t tempTruthJtPt = aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI));
		  Float_t tempTruthJtPhi= aktR_truth_jet_phi_p->at(aktRhi_truthpos_p->at(jI));
		  Float_t tempTruthJtEta = aktR_truth_jet_eta_p->at(aktRhi_truthpos_p->at(jI));
		  
		  TLorentzVector tempRecoJet1, tempRecoJet2;
		  TLorentzVector tempTruthJet1, tempTruthJet2;
		  tempTruthJet1.SetPtEtaPhiM(tempTruthJtPt, tempTruthJtEta, tempTruthJtPhi, 0.0);
		  tempRecoJet1.SetPtEtaPhiM(aktRhi_em_xcalib_jet_pt_p->at(jI), aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), 0.0);
	     
		  //		    if(tempTruthJtPt >= jtPtBins[0] && tempTruthJtPt < jtPtBins[nJtPtBins]){
		  truthGammaPt_ = truthPhotonPt;
		  truthXJ_ = tempTruthJtPt/truthGammaPt_;
		  truthJtPt_ = tempTruthJtPt;
		  unfoldWeight_ = fullWeight;
		  if(isPP) unfoldCent_ = 0;
		  else unfoldCent_ = cent;
						      
		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		  unfoldXJTree_p->Fill();
		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		  //Now do double XJJ
		  for(unsigned int jI2 = jI+1; jI2 < aktRhi_em_xcalib_jet_pt_p->size(); ++jI2){
		    if(aktRhi_truthpos_p->at(jI2) < 0) continue;
		    
		    Float_t dPhi = TMath::Abs(getDPHI(aktRhi_em_xcalib_jet_phi_p->at(jI2), photon_phi_p->at(pI)));
		    if(dPhi < gammaJtDPhiCut) continue;			
		    //			if(aktRhi_em_xcalib_jet_pt_p->at(jI2) < jtPtBins[0] && aktRhi_em_xcalib_jet_pt_p->at(jI2) >= jtPtBins[nJtPtBins]) continue;
						
		    Float_t tempTruthJtPt2 = aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI2));
		    Float_t tempTruthJtPhi2 = aktR_truth_jet_phi_p->at(aktRhi_truthpos_p->at(jI2));
		    Float_t tempTruthJtEta2 = aktR_truth_jet_eta_p->at(aktRhi_truthpos_p->at(jI2));
		    
		    //		    if(tempTruthJtPt < jtPtBins[0] && tempTruthJtPt >= jtPtBins[nJtPtBins]) continue;		    
		    tempTruthJet2.SetPtEtaPhiM(tempTruthJtPt2, tempTruthJtEta2, tempTruthJtPhi2, 0.0);
		    tempRecoJet2.SetPtEtaPhiM(aktRhi_em_xcalib_jet_pt_p->at(jI2), aktRhi_em_xcalib_jet_eta_p->at(jI2), aktRhi_em_xcalib_jet_phi_p->at(jI2), 0.0);


		    recoJt1Pt_ = tempRecoJet1.Pt();
		    recoJt2Pt_ = tempRecoJet2.Pt();

		    truthJt1Pt_ = tempTruthJet1.Pt();
		    truthJt2Pt_ = tempTruthJet2.Pt();


		    tempTruthJet2 += tempTruthJet1;
		    tempRecoJet2 += tempRecoJet1;
		    truthXJJ_ = tempTruthJet2.Pt()/truthGammaPt_;
		    recoXJJ_ = tempRecoJet2.Pt()/recoGammaPt_;
		    
		    unfoldXJJTree_p->Fill();		   				     	 
		  }
		}
	      }
	    }
	  
	    ++multCounter;
	  
	    if(isMC){
	      if(aktRhi_truthpos_p->at(jI) < 0){
		fillTH1(photonJtFakeVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
	      }
	      else{
		if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
		  fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
		  fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);

		  fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_eta_p->at(jI), fullWeight);
		  fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_eta_p->at(jI), fullWeight);
		  
		  fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
		  fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][nGammaPtBins], aktRhi_em_xcalib_jet_pt_p->at(jI)/photon_pt_p->at(pI), fullWeight);
		
		  ++multCounterGenMatched;
		
		  if(truthPhotonPt > 0){
		    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		    if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) < 0.2){
		      if(truthPhotonPt >= gammaPtBins[0] && truthPhotonPt < gammaPtBins[nGammaPtBins]){
			Int_t genPtPos = ghostPos(nGammaPtBins, gammaPtBins, truthPhotonPt, true, doGlobalDebug);
			Int_t recoPtPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(pI), true, doGlobalDebug);		
			
			fillTH2(photonJtGenResVCentGenPtRecoPt_p[centPos][genPtPos][recoPtPos], aktRhi_em_xcalib_jet_pt_p->at(jI), aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)), fullWeight); 
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
            
	if(isMC){
	  std::vector<TLorentzVector> goodTruthJets, goodTruthJetsDPhi;

	  for(unsigned int tI = 0; tI < aktR_truth_jet_pt_p->size(); ++tI){
	    if(aktR_truth_jet_pt_p->at(tI) < jtPtBinsLow) continue;
	    if(aktR_truth_jet_pt_p->at(tI) >= jtPtBinsHigh) continue;
	    if(aktR_truth_jet_eta_p->at(tI) <= jtEtaBinsLow) continue;
	    if(aktR_truth_jet_eta_p->at(tI) >= jtEtaBinsHigh) continue;
	  
	    Float_t dR = getDR(aktR_truth_jet_eta_p->at(tI), aktR_truth_jet_phi_p->at(tI), photon_eta_p->at(pI), photon_phi_p->at(pI));
	    if(dR < gammaExclusionDR) continue;
	  
	    Float_t dPhi = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(tI), photon_phi_p->at(pI)));
	  
	    fillTH1(photonGenJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	    fillTH1(photonGenJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);

	    TLorentzVector tempJet;
	    tempJet.SetPtEtaPhiM(aktR_truth_jet_pt_p->at(tI), aktR_truth_jet_eta_p->at(tI), aktR_truth_jet_phi_p->at(tI), 0.0);
	    goodTruthJets.push_back(tempJet);
	 	    
	    if(dPhi >= gammaJtDPhiCut){
	      fillTH1(photonGenJtPtVCentPt_p[centPos][ptPos], aktR_truth_jet_pt_p->at(tI), fullWeight);
	      fillTH1(photonGenJtPtVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_pt_p->at(tI), fullWeight);
	      fillTH1(photonGenJtEtaVCentPt_p[centPos][ptPos], aktR_truth_jet_eta_p->at(tI), fullWeight);
	      fillTH1(photonGenJtEtaVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_eta_p->at(tI), fullWeight);
	      fillTH1(photonGenJtXJVCentPt_p[centPos][ptPos], aktR_truth_jet_pt_p->at(tI)/photon_pt_p->at(pI), fullWeight);
	      fillTH1(photonGenJtXJVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_pt_p->at(tI)/photon_pt_p->at(pI), fullWeight);
	      ++multCounterGen;

	      goodTruthJetsDPhi.push_back(tempJet);
	    }
	  }

	  for(unsigned int gI = 0; gI < goodTruthJetsDPhi.size(); ++gI){
	    TLorentzVector jet1 = goodTruthJetsDPhi[gI];

	    for(unsigned int gI2 = gI+1; gI2 < goodTruthJetsDPhi.size(); ++gI2){
	      TLorentzVector jet2 = goodTruthJetsDPhi[gI2];
	      double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	      jet2 += jet1;

	      fillTH1(photonGenMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), fullWeight);
	      fillTH1(photonGenMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(pI), fullWeight);

	      fillTH1(photonGenMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
	      fillTH1(photonGenMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);


	    }
	  }
	}
      
	fillTH1(photonJtMultVCentPt_p[centPos][ptPos], multCounter, fullWeight);
	fillTH1(photonJtMultVCentPt_p[centPos][nGammaPtBins], multCounter, fullWeight);	
      
	if(goodJetsDPhi.size() >= 2){
	  for(auto const & jet : goodJetsDPhi){
	    fillTH1(photonMultiJtPtVCentPt_p[centPos][ptPos], jet.Pt(), fullWeight);
	    fillTH1(photonMultiJtPtVCentPt_p[centPos][nGammaPtBins], jet.Pt(), fullWeight);
	    fillTH1(photonMultiJtXJVCentPt_p[centPos][ptPos], jet.Pt()/photon_pt_p->at(pI), fullWeight);
	    fillTH1(photonMultiJtXJVCentPt_p[centPos][nGammaPtBins], jet.Pt()/photon_pt_p->at(pI), fullWeight);
	  }

	  for(unsigned int jI = 0; jI < goodJetsDPhi.size(); ++jI){
	    TLorentzVector jet1 = goodJetsDPhi[jI];

	    for(unsigned int jI2 = jI+1; jI2 < goodJetsDPhi.size(); ++jI2){
	      TLorentzVector jet2 = goodJetsDPhi[jI2];
	      
	      double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	      jet2 += jet1;

	      fillTH1(photonMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), fullWeight);
	      fillTH1(photonMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(pI), fullWeight);

	      fillTH1(photonMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
	      fillTH1(photonMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);

	      fillTH2(photonPtJtXJJVCent_RAW_p[centPos], jet2.Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), fullWeight);


	      if(isMC){		
		if(goodJetsDPhiTPos[jI] >= 0){
		  if(goodJetsDPhiTPos[jI2] >= 0){
		    fillTH1(photonGenMatchedMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), fullWeight);
		    fillTH1(photonGenMatchedMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(pI), fullWeight);

		    fillTH1(photonGenMatchedMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
		    fillTH1(photonGenMatchedMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);
		  }
		}
	      }

	    }
	  }
	}
      	
	if(isMC){
	  fillTH1(photonGenJtMultVCentPt_p[centPos][ptPos], multCounterGen, fullWeight);
	  fillTH1(photonGenJtMultVCentPt_p[centPos][nGammaPtBins], multCounterGen, fullWeight);

	  fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][ptPos], multCounterGenMatched, fullWeight);
	  fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][nGammaPtBins], multCounterGenMatched, fullWeight);
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
	  
	  if(doGlobalDebug) std::cout << "CENT: " << mixCentPos << ", " << cent << std::endl;
	  if(doGlobalDebug) std::cout << "PSI2: " << mixPsi2Pos << ", " << evtPlane2Phi << std::endl;
	  if(doGlobalDebug) std::cout << "VZ: " << mixVzPos << ", " << vert_z << std::endl;
      
	  unsigned long long key = keyBoy.GetKey(eventKeyVect);
	  unsigned long long maxPos = mixingMap[key].size();
	  if(maxPos == 0){
	    std::cout << "WHOOPS NO AVAILABLE MIXED EVENT. bailing" << std::endl;
	    std::cout << key << ", " << mixCentPos << ", " << cent << std::endl;
	    return 1;
	  }
	  unsigned long long nCurrentMixEvents = 0;
	
	  while(nCurrentMixEvents < nMixEvents){
	    unsigned long long jetPos = maxPos;
	    while(jetPos == maxPos){jetPos = randGen_p->Uniform(0, maxPos-1);}
	    ++(signalMapCounterPost[key]);
	    std::vector<TLorentzVector> jets = mixingMap[key][jetPos];
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, JETS SIZE, MAX, CHOSEN: " << __FILE__ << ", " << __LINE__ << ", " << jets.size() << ", " << maxPos << ", " << jetPos << std::endl; 
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << key << ", " << jetPos << ", " << std::endl; 

	    unsigned long long jetPos2 = maxPos;
	    while(jetPos2 == jetPos || jetPos2 == maxPos){jetPos2 = randGen_p->Uniform(0, maxPos-1);}
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	    std::vector<TLorentzVector> jets2 = mixingMap[key][jetPos2];
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	    int multCounterMix = 0;
	    goodJetsMix[0].clear();
	    goodJetsDPhiMix[0].clear();
	    goodJetsMix[1].clear();
	    goodJetsDPhiMix[1].clear();
	  
	    for(unsigned int jI = 0; jI < jets.size(); ++jI){
	      if(jets[jI].Pt()< jtPtBinsLow) continue;
	      
	      if(jets[jI].Eta() <= jtEtaBinsLow) continue;
	      if(jets[jI].Eta() >= jtEtaBinsHigh) continue;
	  
	      Float_t dR = getDR(jets[jI].Eta(), jets[jI].Phi(), photon_eta_p->at(pI), photon_phi_p->at(pI));
	      if(dR < gammaExclusionDR) continue;
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
	      goodJetsMix[0].push_back(jets[jI]);
	    
	      Float_t dPhi = TMath::Abs(getDPHI(jets[jI].Phi(), photon_phi_p->at(pI)));
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
	      fillTH1(photonMixJtDPhiVCentPt_p[centPos][ptPos], dPhi, mixWeight);
	      fillTH1(photonMixJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, mixWeight);
	      
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	      if(dPhi >= gammaJtDPhiCut){
		goodJetsDPhiMix[0].push_back(jets[jI]);
		fillTH1(photonMixJtPtVCentPt_p[centPos][ptPos], jets[jI].Pt(), mixWeight);
		fillTH1(photonMixJtPtVCentPt_p[centPos][nGammaPtBins], jets[jI].Pt(), mixWeight);
		fillTH1(photonMixJtEtaVCentPt_p[centPos][ptPos], jets[jI].Eta(), mixWeight);
		fillTH1(photonMixJtEtaVCentPt_p[centPos][nGammaPtBins], jets[jI].Eta(), mixWeight);
		fillTH1(photonMixJtXJVCentPt_p[centPos][ptPos], jets[jI].Pt()/photon_pt_p->at(pI), mixWeight);
		fillTH1(photonMixJtXJVCentPt_p[centPos][nGammaPtBins], jets[jI].Pt()/photon_pt_p->at(pI), mixWeight);
	     
		fillTH2(photonPtJtPtVCent_MIX_p[centPos], jets[jI].Pt(), photon_pt_p->at(pI), mixWeight);
		fillTH2(photonPtJtXJVCent_MIX_p[centPos], jets[jI].Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);

		++multCounterMix;
	      }	    
	    }
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	    for(unsigned int jI = 0; jI < jets2.size(); ++jI){
	      if(jets2[jI].Pt()< jtPtBinsLow) continue;
	      if(jets2[jI].Eta() <= jtEtaBinsLow) continue;
	      if(jets2[jI].Eta() >= jtEtaBinsHigh) continue;
	      
	      Float_t dR = getDR(jets2[jI].Eta(), jets2[jI].Phi(), photon_eta_p->at(pI), photon_phi_p->at(pI));
	      if(dR < gammaExclusionDR) continue;
	    	      
	      goodJetsMix[1].push_back(jets2[jI]);
	      
	      Float_t dPhi = TMath::Abs(getDPHI(jets2[jI].Phi(), photon_phi_p->at(pI)));
	      if(dPhi >= gammaJtDPhiCut){
		goodJetsDPhiMix[1].push_back(jets2[jI]);
	      }	    
	    }
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	
	    if(goodJetsDPhiMix[0].size() >= 2){
	      for(auto const & jet : goodJetsDPhiMix[0]){
		fillTH1(photonMixMultiJtPtVCentPt_p[centPos][ptPos], jet.Pt(), mixWeight);
		fillTH1(photonMixMultiJtPtVCentPt_p[centPos][nGammaPtBins], jet.Pt(), mixWeight);
		fillTH1(photonMixMultiJtXJVCentPt_p[centPos][ptPos], jet.Pt()/photon_pt_p->at(pI), mixWeight);
		fillTH1(photonMixMultiJtXJVCentPt_p[centPos][nGammaPtBins], jet.Pt()/photon_pt_p->at(pI), mixWeight);
	      }
	      
	      //First do pure background
	      for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
		TLorentzVector jet1 = goodJetsDPhiMix[0][bI];
		
		for(unsigned int bI2 = bI+1; bI2 < goodJetsDPhiMix[0].size(); ++bI2){
		  TLorentzVector jet2 = goodJetsDPhiMix[0][bI2];
		  
		  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
		  jet2 += jet1;
		  
		  fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), mixWeight);
		  fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(pI), mixWeight);
		  
		  fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
		  fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);

		  fillTH2(photonPtJtXJJVCent_MIX_p[centPos], jet2.Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);		  
		}
	      }
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	  
	      //Now do mixed background (signal jet + associated w/ fake jets)
	      for(unsigned int sI = 0; sI < goodJetsDPhi.size(); ++sI){
		TLorentzVector jet1 = goodJetsDPhi[sI];
		
		for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
		  TLorentzVector jet2 = goodJetsDPhiMix[0][bI];
	
		  Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		  if(dR < mixJetExclusionDR) continue;

		  
		  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
		  jet2 += jet1;
		  
		  fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), mixWeight);
		  fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][nGammaPtBins],jet2.Pt()/photon_pt_p->at(pI), mixWeight);
		  
		  fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
		  fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);

		  fillTH2(photonPtJtXJJVCent_MIX_p[centPos], jet2.Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);		  

		}
	      }
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
	      
	      //Now calculate the mixed event correction, for cases where your gamma + single jet embed accidentally picked a fake jet
	      for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
		TLorentzVector jet1 = goodJetsDPhiMix[0][bI];
		
		for(unsigned int bI2 = 0; bI2 < goodJetsDPhiMix[1].size(); ++bI2){
		  TLorentzVector jet2 = goodJetsDPhiMix[1][bI2];

		  Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		  if(dR < mixJetExclusionDR) continue;
		  
		  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
		  jet2 += jet1;
		  
		  fillTH1(photonMixCorrectionMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(pI), mixWeight);
		  fillTH1(photonMixCorrectionMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(pI), mixWeight);

		  fillTH1(photonMixCorrectionMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
		  fillTH1(photonMixCorrectionMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);

		  fillTH2(photonPtJtXJJVCent_MIXCorrection_p[centPos], jet2.Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);		  
		}
	      }	      	      
	    }
	    
	    fillTH1(photonMixJtMultVCentPt_p[centPos][ptPos], multCounterMix, mixWeight);
	    fillTH1(photonMixJtMultVCentPt_p[centPos][nGammaPtBins], multCounterMix, mixWeight);	  
	    
	    fillTH1(photonSubJtMultModVCentPt_p[centPos][ptPos], TMath::Max(0, multCounter - multCounterMix), mixWeight);
	    fillTH1(photonSubJtMultModVCentPt_p[centPos][nGammaPtBins], TMath::Max(0, multCounter - multCounterMix), mixWeight);

	    ++nCurrentMixEvents;
	  }
	}
	else if(isPP){
	  fillTH1(photonSubJtMultModVCentPt_p[centPos][ptPos], TMath::Max(0, multCounter), fullWeight);
	  fillTH1(photonSubJtMultModVCentPt_p[centPos][nGammaPtBins], TMath::Max(0, multCounter), fullWeight);
	}
      }
      
      fillTH2(photonEtaPt_p[centPos], etaValMain, photon_pt_p->at(pI), fullWeight);
    }

    if(nGoodPhotons >= 10) nGoodPhotons = 9;
    ++(perEventPhotonCounts[nGoodPhotons]);
  }

  if(doMix){
    //Dump info on the mixing
    unsigned long long maxRatioKey = 0;
    double maximumRatioVal = 0;
    unsigned long long maximumNumVal = 0;
    unsigned long long maximumDenomVal = 0;
    
    unsigned long long minRatioKey = 0;
    double minimumRatioVal = 2147483647; //hardcoded max value of a signed integer
    unsigned long long minimumNumVal = 0;
    unsigned long long minimumDenomVal = 0;
    
    std::vector<double> tempRatioVals;
    
    for(auto const & mixes : mixingMapCounter){
      if(signalMapCounterPost[mixes.first] == 0) continue;
      
      double ratio = ((double)mixes.second)/(double)signalMapCounterPost[mixes.first];
      
      if(ratio < minimumRatioVal){
	minRatioKey = mixes.first;
	minimumRatioVal = ratio;
	minimumDenomVal = signalMapCounterPost[mixes.first];
	minimumNumVal = mixes.second;
      }
      else if(ratio > maximumRatioVal){
	maxRatioKey = mixes.first;
	maximumRatioVal = ratio;
	maximumDenomVal = signalMapCounterPost[mixes.first];
	maximumNumVal = mixes.second;
      }
      
      tempRatioVals.push_back(ratio);
    }
    
    std::sort(std::begin(tempRatioVals), std::end(tempRatioVals));
    
    std::cout << "LOWEST RATIO MIX-PER-SIGNAL, CORRESPONDING KEY: " << minimumNumVal << "/" << minimumDenomVal << "=" << minimumRatioVal << " (" << keyBoy.GetKeyStr(minRatioKey) << "), key=" << minRatioKey << std::endl;
    std::cout << "HIGHEST RATIO MIX-PER-SIGNAL, CORRESPONDING KEY: " << maximumNumVal << "/" << maximumDenomVal << "=" << maximumRatioVal << " (" << keyBoy.GetKeyStr(maxRatioKey) << "), key=" << maxRatioKey << std::endl;
    std::cout << "MEDIAN RATIO MIX-PER-SIGNAL: " << tempRatioVals[tempRatioVals.size()/2] << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Number of events w/ N photons" << std::endl;
  for(auto const & val: perEventPhotonCounts){
    std::cout << " N=" << val.first << ": " << val.second << std::endl;
  }
  std::cout << std::endl;
  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  if(isMC){
    unfoldXJTree_p->Write("", TObject::kOverwrite);
    delete unfoldXJTree_p;

    unfoldXJJTree_p->Write("", TObject::kOverwrite);
    delete unfoldXJJTree_p;
  }


  //Pre-write and delete some of these require some mods
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      if(isMC){
	photonJtFakeVCentPt_p[cI][pI]->Divide(photonJtPtVCentPt_p[cI][pI]);
	const double maxFakeVal = 0.45;
	double tempMax = getMax(photonJtFakeVCentPt_p[cI][pI]);

	if(tempMax > maxFakeVal) photonJtFakeVCentPt_p[cI][pI]->SetMaximum(maxFakeVal);
      }

      if(doMix){
	photonSubJtDPhiVCentPt_p[cI][pI]->Add(photonJtDPhiVCentPt_p[cI][pI], photonMixJtDPhiVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtPtVCentPt_p[cI][pI]->Add(photonJtPtVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Add(photonMultiJtPtVCentPt_p[cI][pI], photonMixMultiJtPtVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtEtaVCentPt_p[cI][pI]->Add(photonJtEtaVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtXJVCentPt_p[cI][pI]->Add(photonJtXJVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Add(photonMultiJtXJVCentPt_p[cI][pI], photonMixMultiJtXJVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubJtMultVCentPt_p[cI][pI]->Add(photonJtMultVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI], 1.0, -1.0);

	//Correct the mixed event
	photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI]->Add(photonMixMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI], 1.0, -1.0);
	photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI]->Add(photonMixMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI], 1.0, -1.0);

	//Now subtract off the corrected mixed event
	photonSubMultiJtXJJVCentPt_p[cI][pI]->Add(photonMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI], 1.0, -1.0);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Add(photonMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI], 1.0, -1.0);
      }
      else if(isPP){
	photonSubJtDPhiVCentPt_p[cI][pI]->Add(photonJtDPhiVCentPt_p[cI][pI]);
	photonSubJtPtVCentPt_p[cI][pI]->Add(photonJtPtVCentPt_p[cI][pI]);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Add(photonMultiJtPtVCentPt_p[cI][pI]);
	photonSubJtEtaVCentPt_p[cI][pI]->Add(photonJtEtaVCentPt_p[cI][pI]);
	photonSubJtXJVCentPt_p[cI][pI]->Add(photonJtXJVCentPt_p[cI][pI]);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Add(photonMultiJtXJVCentPt_p[cI][pI]);
	photonSubJtMultVCentPt_p[cI][pI]->Add(photonJtMultVCentPt_p[cI][pI]);	


	photonSubMultiJtXJJVCentPt_p[cI][pI]->Add(photonMultiJtXJJVCentPt_p[cI][pI]);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Add(photonMultiJtDPhiJJVCentPt_p[cI][pI]);
      }
      
      photonJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonMultiJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonMultiJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
    
      photonMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      photonMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMultiJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMultiJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	
	photonGenMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      }

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixMultiJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixMultiJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

	photonMixMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

	photonSubJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

	photonSubMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
 	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      }
      else if(isPP){
	photonSubJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
      }
    }


    if(doMix){
      //Correct XJJ Mix Event
      photonPtJtXJJVCent_MIXCorrected_p[cI]->Add(photonPtJtXJJVCent_MIX_p[cI], photonPtJtXJJVCent_MIXCorrection_p[cI], 1.0, -1.0);

      //Subtract off the mixed event
      photonPtJtPtVCent_SUB_p[cI]->Add(photonPtJtPtVCent_RAW_p[cI], photonPtJtPtVCent_MIX_p[cI], 1.0, -1.0);
      photonPtJtXJVCent_SUB_p[cI]->Add(photonPtJtXJVCent_RAW_p[cI], photonPtJtXJVCent_MIX_p[cI], 1.0, -1.0);
      photonPtJtXJJVCent_SUB_p[cI]->Add(photonPtJtXJJVCent_RAW_p[cI], photonPtJtXJJVCent_MIXCorrected_p[cI], 1.0, -1.0);
    }
    else if(isPP){
      photonPtJtPtVCent_SUB_p[cI]->Add(photonPtJtPtVCent_RAW_p[cI]);
      photonPtJtXJVCent_SUB_p[cI]->Add(photonPtJtXJVCent_RAW_p[cI]);
      photonPtJtXJJVCent_SUB_p[cI]->Add(photonPtJtXJJVCent_RAW_p[cI]);
    }
  }


  if(doMix){
    TDirectoryFile* mixDir_p = (TDirectoryFile*)outFile_p->mkdir("mixDir");
    mixDir_p->cd();

    if(doMixCent){
      mixingCentrality_p->Write("", TObject::kOverwrite);
      delete mixingCentrality_p;

      if(doMixPsi2){
	mixingCentralityPsi2_p->Write("", TObject::kOverwrite);
	delete mixingCentralityPsi2_p;

	if(doMixVz){
	  for(Int_t vI = 0; vI < nMixVzBins; ++vI){
	    mixingCentralityPsi2_VzBinned_p[vI]->Write("", TObject::kOverwrite);
	    delete mixingCentralityPsi2_VzBinned_p[vI];
	  }
	}
      }
      if(doMixVz){
	mixingCentralityVz_p->Write("", TObject::kOverwrite);
	delete mixingCentralityVz_p;
      }
    }
    if(doMixPsi2){
      mixingPsi2_p->Write("", TObject::kOverwrite);
      delete mixingPsi2_p;

      if(doMixVz){
	mixingPsi2Vz_p->Write("", TObject::kOverwrite);
	delete mixingPsi2Vz_p;
      }      
    }
    if(doMixVz){
      mixingVz_p->Write("", TObject::kOverwrite);
      delete mixingVz_p;
    }

    mixDir_p->Close();
    delete mixDir_p;
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

    for(Int_t eI = 0; eI < nGammaEtaBinsSub+1; ++eI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      photonPtVCentEta_p[cI][eI]->Write("", TObject::kOverwrite);

      if(isMC){
	photonGenResVCentEta_p[cI][eI]->Write("", TObject::kOverwrite);
      }
    }    

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 
    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonMultiJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonMultiJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixMultiJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixMultiJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonMixJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

      	photonSubJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtMultModVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      }
      else if(isPP){
      	photonSubJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonSubJtMultModVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      }


      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedJtPtVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedJtEtaVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedJtXJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedMultiJtXJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
	photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtPtVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedJtXJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedMultiJtXJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
	photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonGenMatchedJtMultVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);

	photonJtFakeVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);	
      }
    }

    photonPtJtPtVCent_RAW_p[cI]->Write("", TObject::kOverwrite);
    photonPtJtXJVCent_RAW_p[cI]->Write("", TObject::kOverwrite);

    photonPtJtXJJVCent_RAW_p[cI]->Write("", TObject::kOverwrite);

    if(doMix){
      photonPtJtPtVCent_MIX_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJVCent_MIX_p[cI]->Write("", TObject::kOverwrite);

      photonPtJtXJJVCent_MIX_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_MIXCorrection_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_MIXCorrected_p[cI]->Write("", TObject::kOverwrite);

      photonPtJtPtVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
    }
    else if(isPP){
      photonPtJtPtVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_SUB_p[cI]->Write("", TObject::kOverwrite);
    }

    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2]->Write("", TObject::kOverwrite);
	}
      }

      for(Int_t gI = 0; gI < nJtPtBins; ++gI){
        photonJtRecoOverGenVCentJtPt_p[cI][gI]->Write("", TObject::kOverwrite);
      }
    }
    

    for(Int_t gI = 0; gI < nJtEtaBinsSub+1; ++gI){
      photonJtCorrOverUncorrVCentJtEta_p[cI][gI]->Write("", TObject::kOverwrite);
    }
    photonEtaPt_p[cI]->Write("", TObject::kOverwrite);

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
    }
    centDir_p->Close();
    delete centDir_p;
  }

  //Before we do a bunch of deletions lets auto-generate some stats tables
  if(!isMC){
    std::string lStr = "l | ";
    for(int gI = 0; gI < nGammaPtBins; ++gI){
      lStr = lStr + "l | ";
    }
    // lStr = lStr + "l";
    std::string multiColumn = std::to_string(nGammaPtBins+2);

    std::string tempDPhiStr = gammaJtDPhiLabel.substr(gammaJtDPhiLabel.find(">")+2, gammaJtDPhiLabel.size());
    tempDPhiStr.replace(tempDPhiStr.find("#"),1,"\\");
    std::string jtPtStr = "R=" +jetRStr + ", " + prettyString(jtPtBins[0], 1, false) + " < \\ptj < " + prettyString(jtPtBins[nJtPtBins], 1, false) + "; |\\dphigj| > " + tempDPhiStr;
    
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
    std::cout << "\\vspace{-0.6cm}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
    std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With Exactly 2 Jets (" << jtPtStr << ")} \\\\" << std::endl;
    std::string columnStr = " Centrality &";
    for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
      columnStr = columnStr + " " + prettyString(gammaPtBins[pI], 1, false) + " < \\ptg < " + prettyString(gammaPtBins[pI+1], 1, false) + " &";
    }
    columnStr = columnStr + " " + prettyString(gammaPtBins[0], 1, false) + " < \\ptg < " + prettyString(gammaPtBins[nGammaPtBins], 1, false) + " \\\\ \\hline" ;
    
    std::cout << columnStr << std::endl;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){ 
      std::string outStr = binsToLabelStr[centBinsStr[cI]];
      if(outStr.find("%") != std::string::npos){
	outStr.replace(outStr.find("%"), 1, "\\%");
      }
      
      for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
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
      
      for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
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

    std::cout << std::endl;
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
    std::cout << "\\vspace{-0.6cm}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
    std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With 1 or More Jets (" << jtPtStr << ")} \\\\" << std::endl;
    std::cout << columnStr << std::endl;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){ 
      std::string outStr = binsToLabelStr[centBinsStr[cI]];
      if(outStr.find("%") != std::string::npos){
	outStr.replace(outStr.find("%"), 1, "\\%");
      }
      
      for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
	int countNJet = 0;
	
	for(Int_t bIX = 1; bIX < photonJtMultVCentPt_p[cI][pI]->GetXaxis()->GetNbins(); ++bIX){
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
    for(Int_t eI = 0; eI < nGammaEtaBinsSub+1; ++eI){
      delete photonPtVCentEta_p[cI][eI];

      if(isMC){
	delete photonGenResVCentEta_p[cI][eI];
      }
    }

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      delete photonEtaVCentPt_p[cI][pI];
      delete photonPhiVCentPt_p[cI][pI];

      delete photonJtDPhiVCentPt_p[cI][pI];
      delete photonJtPtVCentPt_p[cI][pI];
      delete photonMultiJtPtVCentPt_p[cI][pI];
      delete photonJtEtaVCentPt_p[cI][pI];
      delete photonJtXJVCentPt_p[cI][pI];
      delete photonMultiJtXJVCentPt_p[cI][pI];
      delete photonMultiJtXJJVCentPt_p[cI][pI];
      delete photonMultiJtDPhiJJVCentPt_p[cI][pI];
      delete photonJtMultVCentPt_p[cI][pI];

      if(doMix){
	delete photonMixJtDPhiVCentPt_p[cI][pI];
	delete photonMixJtPtVCentPt_p[cI][pI];
	delete photonMixMultiJtPtVCentPt_p[cI][pI];
	delete photonMixJtEtaVCentPt_p[cI][pI];
	delete photonMixJtXJVCentPt_p[cI][pI];
	delete photonMixMultiJtXJVCentPt_p[cI][pI];
	delete photonMixMultiJtXJJVCentPt_p[cI][pI];
	delete photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI];
	delete photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI];
	delete photonMixMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonMixJtMultVCentPt_p[cI][pI];

      	delete photonSubJtDPhiVCentPt_p[cI][pI];
	delete photonSubJtPtVCentPt_p[cI][pI];
	delete photonSubMultiJtPtVCentPt_p[cI][pI];
	delete photonSubJtEtaVCentPt_p[cI][pI];
	delete photonSubJtXJVCentPt_p[cI][pI];
	delete photonSubMultiJtXJVCentPt_p[cI][pI];
	delete photonSubMultiJtXJJVCentPt_p[cI][pI];
	delete photonSubMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonSubJtMultVCentPt_p[cI][pI];
	delete photonSubJtMultModVCentPt_p[cI][pI];
      }
      else if(isPP){
      	delete photonSubJtDPhiVCentPt_p[cI][pI];
	delete photonSubJtPtVCentPt_p[cI][pI];
	delete photonSubMultiJtPtVCentPt_p[cI][pI];
	delete photonSubJtEtaVCentPt_p[cI][pI];
	delete photonSubJtXJVCentPt_p[cI][pI];
	delete photonSubMultiJtXJVCentPt_p[cI][pI];
	delete photonSubMultiJtXJJVCentPt_p[cI][pI];
	delete photonSubMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonSubJtMultVCentPt_p[cI][pI];
	delete photonSubJtMultModVCentPt_p[cI][pI];
      }

      if(isMC){
	delete photonGenJtDPhiVCentPt_p[cI][pI];
	delete photonGenJtPtVCentPt_p[cI][pI];
	delete photonGenJtEtaVCentPt_p[cI][pI];
	delete photonGenJtXJVCentPt_p[cI][pI];
	delete photonGenMultiJtXJJVCentPt_p[cI][pI];
	delete photonGenMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonGenJtMultVCentPt_p[cI][pI];

	delete photonGenMatchedJtDPhiVCentPt_p[cI][pI];
	delete photonGenMatchedJtPtVCentPt_p[cI][pI];
	delete photonGenMatchedJtEtaVCentPt_p[cI][pI];
	delete photonGenMatchedJtXJVCentPt_p[cI][pI];
	delete photonGenMatchedMultiJtXJJVCentPt_p[cI][pI];
	delete photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI];
	delete photonGenMatchedJtMultVCentPt_p[cI][pI];

	delete photonJtFakeVCentPt_p[cI][pI];	
      }
    }

    delete photonPtJtPtVCent_RAW_p[cI];
    delete photonPtJtXJVCent_RAW_p[cI];
    delete photonPtJtXJJVCent_RAW_p[cI];

    if(doMix){
      delete photonPtJtPtVCent_MIX_p[cI];
      delete photonPtJtXJVCent_MIX_p[cI];
      delete photonPtJtXJJVCent_MIX_p[cI];
      delete photonPtJtXJJVCent_MIXCorrection_p[cI];
      delete photonPtJtXJJVCent_MIXCorrected_p[cI];

      delete photonPtJtPtVCent_SUB_p[cI];
      delete photonPtJtXJVCent_SUB_p[cI];
      delete photonPtJtXJJVCent_SUB_p[cI];
    }
    else if(isPP){
      delete photonPtJtPtVCent_SUB_p[cI];
      delete photonPtJtXJVCent_SUB_p[cI];
      delete photonPtJtXJJVCent_SUB_p[cI];
    }

    if(isMC){
      for(Int_t ptI = 0; ptI < nGammaPtBins+1; ++ptI){
	for(Int_t ptI2 = 0; ptI2 < nGammaPtBins+1; ++ptI2){
	  delete photonJtGenResVCentGenPtRecoPt_p[cI][ptI][ptI2];
	}
      }

      for(Int_t gI = 0; gI < nJtPtBins; ++gI){
        delete photonJtRecoOverGenVCentJtPt_p[cI][gI];
      }
    }
    
    for(Int_t gI = 0; gI < nJtEtaBinsSub+1; ++gI){
      delete photonJtCorrOverUncorrVCentJtEta_p[cI][gI];
    }

    delete photonEtaPt_p[cI];

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      delete photonEtaPhiVCentPt_p[cI][pI];
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  /*
  std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form          
  configMap["RECOJTPTMIN"] = prettyString(recoJtPtMin, 1, false);
  
  TEnv configEnv; // We will convert to tenv and store in file
  for(auto const & val : configMap){
    configEnv.SetValue(val.first.c_str(), val.second.c_str()); //Fill out the map
  }
  configEnv.Write("config", TObject::kOverwrite);  
  */

  config_p->SetValue("RECOJTPTMIN", prettyString(recoJtPtMin, 1, false).c_str());
  config_p->Write("config", TObject::kOverwrite);

  TEnv labelEnv;
  for(auto const & lab : binsToLabelStr){
    labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
  }
  labelEnv.Write("label", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  globalTimer.stop();
  double globalTimeCPU = globalTimer.totalCPU();
  double globalTimeWall = globalTimer.totalWall();

  std::cout << "FINAL TIME: " << std::endl;
  std::cout << " WALL: " << globalTimeWall << std::endl;
  std::cout << " CPU:  " << globalTimeCPU << std::endl;
  
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
