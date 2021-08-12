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
#include "TF1.h"
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
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/mixMachine.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/sampleHandler.h"
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
  bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INFILENAME",
                                              "OUTFILENAME",
                                              "CENTFILENAME",
					      "MIXFILENAME",
					      "DOUNIFIEDPURITY",
					      "PURITYFILENAME",
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
					      "DOPTISOCORRECTION",
					      "DOCENTISOCORRECTION",
					      "ISOLATIONR",
					      "SIDEBANDTYPE",
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
					      "GAMMAMULTIJTDPHI",  
					      "NXJBINS",
					      "XJBINSLOW",
					      "XJBINSHIGH",
					      "XJBINSDOCUSTOM",
					      "XJBINSCUSTOM",
					      "NAJBINS",
					      "AJBINSLOW",
					      "AJBINSHIGH",
					      "AJBINSDOCUSTOM",
					      "AJBINSCUSTOM"};

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
  const bool doUnifiedPurity = config_p->GetValue("DOUNIFIEDPURITY", 0);
  std::string inPurityFileName = config_p->GetValue("PURITYFILENAME", "");

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
  if(!check.checkFileExt(inPurityFileName, "root")) return 1;
  
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
  const Int_t nBarrelAndEC = 3;
  const std::string barrelAndECStr[nBarrelAndEC] = {"Barrel", "EC", "BarrelAndEC"};
  const std::string barrelAndECFitStr[nBarrelAndEC] = {"Eta0p00to1p37", "Eta1p52to2p37", "Eta0p00to2p37"};
  
  Int_t barrelAndECComboPos = -1;
  for(Int_t bI = 0; bI < nBarrelAndEC; ++bI){
    if(isStrSame(barrelAndECStr[bI], "BarrelAndEC")){
      barrelAndECComboPos = bI;
      break;
    }
  }


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

  const bool doPtIsoCorrection = config_p->GetValue("DOPTISOCORRECTION", 1);
  const bool doCentIsoCorrection = config_p->GetValue("DOCENTISOCORRECTION", 1);
  const int isolationR = config_p->GetValue("ISOLATIONR", 3);
  const photonType sidebandType = (photonType)config_p->GetValue("SIDEBANDTYPE", 2);

  std::vector<int> validIsoR = {2,3,4};
  std::vector<int> validSidebandType = {1,2,3,4,5};
  //SIDEBAND TYPE 1: tight and non-isolated
  //SIDEBAND TYPE 2: non-tight and isolated
  //SIDEBAND TYPE 3: non-tight and non-isolated
  //SIDEBAND TYPE 4: non-isolated
  //SIDEBAND TYPE 5: non-tight
  if(!vectContainsInt(isolationR, &validIsoR)){
    std::cout << "Given ISOLATIONR \'" << isolationR << "\' is not valid. Please pick from:" << std::endl;
    for(auto const & isoRs : validIsoR){
      std::cout << " " << isoRs << std::endl;
    }
    std::cout << "return 1." << std::endl;
    return 1;
  }
  if(!vectContainsInt(sidebandType, &validSidebandType)){
    std::cout << "Given SIDEBANDTYPE \'" << sidebandType << "\' is not valid. Please pick from:" << std::endl;
    for(auto const & isoRs : validSidebandType){
      std::cout << " " << isoRs << std::endl;
    }
    std::cout << "return 1." << std::endl;
    return 1;
  }

  //Set the gamma pt bins
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
  std::string gammaPtBinsStrForConfig = "";
  std::vector<std::string> genGammaPtBinsStr, recoGammaPtBinsStr, gammaPtBinsStr;
  for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
    genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(pI));
    recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(pI));
    gammaPtBinsStr.push_back("GammaPt" + std::to_string(pI));

    gammaPtBinsStrForConfig = gammaPtBinsStrForConfig + std::to_string(gammaPtBins[pI]) + ", ";

    binsToLabelStr[genGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[recoGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[gammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
  }
  gammaPtBinsStrForConfig = gammaPtBinsStrForConfig + std::to_string(gammaPtBins[nGammaPtBins]);

  genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(nGammaPtBins));
  recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(nGammaPtBins));
  gammaPtBinsStr.push_back("GammaPt" + std::to_string(nGammaPtBins));

  binsToLabelStr[genGammaPtBinsStr[genGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[recoGammaPtBinsStr[recoGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[gammaPtBinsStr[gammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  //Also set the fine grain photon bins right now

  
  const Int_t nGammaFinePtBins = (gammaPtBins[nGammaPtBins-1] - gammaPtBins[0])/5 + 1;
  if(!goodBinning(inConfigFileName, nMaxPtBins, nGammaFinePtBins, "NGAMMAFINEPTBINS")) return 1;
  Double_t gammaFinePtBins[nMaxPtBins+1];
  std::string gammaFinePtBinsStr = "";

  for(Int_t gI = 0; gI < nGammaFinePtBins; ++gI){
    gammaFinePtBins[gI] = gammaPtBins[0]  + 5*gI;

    gammaFinePtBinsStr = gammaFinePtBinsStr + std::to_string((int)gammaFinePtBins[gI]) + ",";
  }
  gammaFinePtBins[nGammaFinePtBins] = gammaPtBins[nGammaPtBins];
  gammaFinePtBinsStr = gammaFinePtBinsStr + std::to_string((int)gammaPtBins[nGammaPtBins]);
  
 
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

  std::string jtPtBinsStrForConfig = "";
  std::vector<std::string> jtPtBinsStr;
  for(Int_t pI = 0; pI < nJtPtBins; ++pI){
    if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
    else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));

    jtPtBinsStrForConfig = jtPtBinsStrForConfig + std::to_string(jtPtBins[pI]) + ", ";

    binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
  }
  jtPtBinsStrForConfig = jtPtBinsStrForConfig + std::to_string(jtPtBins[nJtPtBins]);

  jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
  binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nDPhiBins = config_p->GetValue("NDPHIBINS", 16);
  const Float_t dPhiBinsLow = config_p->GetValue("DPHIBINSLOW", 0.0);
  const Float_t dPhiBinsHigh = config_p->GetValue("DPHIBINSHIGH", TMath::Pi());
  Double_t dphiBins[nMaxPtBins+1];
  getLinBins(dPhiBinsLow, dPhiBinsHigh, nDPhiBins, dphiBins);
  std::string dphiBinsStrForConfig;
  for(Int_t jI = 0; jI < nDPhiBins; ++jI){
    dphiBinsStrForConfig = dphiBinsStrForConfig + prettyString(dphiBins[jI],8,false) + ",";      
  }  
  dphiBinsStrForConfig = dphiBinsStrForConfig + prettyString(dphiBins[nDPhiBins],8,false);

  //Gotta do this because otherwise the bin limits are slightly off
  std::vector<float> dphiBinsStrForConfigVect = strToVectF(dphiBinsStrForConfig);
  for(unsigned int i = 0; i < dphiBinsStrForConfigVect.size(); ++i){
    dphiBins[i] = dphiBinsStrForConfigVect[i];
  }
  
  
  const Double_t gammaJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);  
  const Double_t gammaMultiJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAMULTIJTDPHI", ""), doGlobalDebug);  

  /*
  std::cout << "CHECK GAMMA JTDPHI, MULTIJTDPHI:" << std::endl;
  std::cout << "GammaJtDPhiCut: " << gammaJtDPhiCut << std::endl;
  std::cout << "GammaMultiJtDPhiCut: " << gammaMultiJtDPhiCut << std::endl;
  std::cout << std::endl; 
  return 1;
  */
  
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
  std::string xjBinsStrForConfig = "";
  for(Int_t jI = 0; jI < nXJBins; ++jI){
    xjBinsStrForConfig = xjBinsStrForConfig + std::to_string(xjBins[jI]) + ",";
  }
  xjBinsStrForConfig = xjBinsStrForConfig + std::to_string(xjBins[nXJBins]);

  const Int_t nAJBins = config_p->GetValue("NAJBINS", 20);
  const Float_t ajBinsLow = config_p->GetValue("AJBINSLOW", 0.0);
  const Float_t ajBinsHigh = config_p->GetValue("AJBINSHIGH", 2.0);
  const Bool_t ajBinsDoCustom = config_p->GetValue("AJBINSDOCUSTOM", 0);
  std::string ajBinsCustomStr = config_p->GetValue("AJBINSCUSTOM", "");
  Double_t ajBins[nMaxPtBins+1];
  if(!ajBinsDoCustom) getLinBins(ajBinsLow, ajBinsHigh, nAJBins, ajBins);
  else{
    if(ajBinsCustomStr.size() == 0){
      std::cout << "No custom bins supplied despite AJBINSDOCUSOM \'" << ajBinsDoCustom <<  "'\". Reverting to uniform binning." << std::endl;
      getLinBins(ajBinsLow, ajBinsHigh, nAJBins, ajBins);
    }
    else{
      std::vector<float> ajBinsCustomVect = strToVectF(ajBinsCustomStr);
      //Check that the bins are ordered and the number matches nAJBins
      Int_t nAJBinsCustomVect = ajBinsCustomVect.size();
      if(nAJBinsCustomVect-1 != nAJBins){
	std::cout << "Number of supplied NAJBINS, \'" << nAJBins << "\', does not match number from custom bin string \'" << nAJBinsCustomVect-1 << "\'. Reverting to uniform bins" << std::endl;
        getLinBins(ajBinsLow, ajBinsHigh, nAJBins, ajBins);
      }
      else{
        bool areBinsOrdered = true;
        for(unsigned int i = 1; i < ajBinsCustomVect.size(); ++i){
          if(ajBinsCustomVect[i] < ajBinsCustomVect[i-1]){
            areBinsOrdered = false;
            break;
          }
        }

        if(areBinsOrdered){
          for(unsigned int i = 0; i < ajBinsCustomVect.size(); ++i){
            ajBins[i] = ajBinsCustomVect[i];
          }
        }
        else{
	  std::cout << "AJBINSCUSTOM provided \'" << ajBinsCustomStr << "\'is not ordered. Reverting to uniform bins"  << std::endl;
          getLinBins(ajBinsLow, ajBinsHigh, nAJBins, ajBins);
        }
      }
    }
  }
  std::string ajBinsStrForConfig = "";
  for(Int_t jI = 0; jI < nAJBins; ++jI){
    ajBinsStrForConfig = ajBinsStrForConfig + std::to_string(ajBins[jI]) + ",";
  }
  ajBinsStrForConfig = ajBinsStrForConfig + std::to_string(ajBins[nAJBins]);

  const unsigned long long mixCap = config_p->GetValue("MIXCAP", 100000000);
  const unsigned long long nMixEvents = config_p->GetValue("NMIXEVENTS", 1);
  const bool doStrictMix = config_p->GetValue("DOSTRICTMIX", true);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nJESSys = 18;
  const Int_t nJERSys = 9;

  const Int_t nJetSysAndNom = 1 + nJESSys + nJERSys;

  std::string jetSysAndNomStr[nJetSysAndNom] = {"Nominal"}; 
  for(Int_t sI = 0; sI < nJESSys; ++sI){
    jetSysAndNomStr[1 + sI] = "JES" + std::to_string(sI);
  }
  for(Int_t sI = 0; sI < nJERSys; ++sI){
    jetSysAndNomStr[1 + nJESSys + sI] = "JER" + std::to_string(sI);
  }


  Float_t recoGammaPt_;
  Float_t recoGammaPhi_;
  Float_t truthGammaPt_;
  Float_t truthGammaPhi_;

  Float_t recoXJ_[nJetSysAndNom];
  Float_t recoJtPt_[nJetSysAndNom];
  Float_t recoJtPhi_;
  Float_t recoXJJ_[nJetSysAndNom];
  Float_t recoJt1Pt_[nJetSysAndNom];
  Float_t recoJt2Pt_[nJetSysAndNom];
  Float_t recoJt1Phi_;
  Float_t recoJt2Phi_;
  Float_t truthXJ_;
  Float_t truthJtPt_;
  Float_t truthJtPhi_;
  Float_t truthXJJ_;
  Float_t truthJt1Pt_;
  Float_t truthJt2Pt_;
  Float_t truthJt1Phi_;
  Float_t truthJt2Phi_;

  Float_t unfoldWeight_;
  Float_t unfoldCent_;
  
  TTree* unfoldPhotonPtTree_p = nullptr;
  TTree* unfoldXJTree_p = nullptr;
  TTree* unfoldXJJTree_p = nullptr;

  if(isMC){
    unfoldPhotonPtTree_p = new TTree("unfoldPhotonPtTree_p", "");
    unfoldPhotonPtTree_p->Branch("recoGammaPt", &recoGammaPt_, "recoGammaPt/F");
    unfoldPhotonPtTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F");
    unfoldPhotonPtTree_p->Branch("unfoldWeight", &unfoldWeight_, "unfoldWeight/F"); 
    unfoldPhotonPtTree_p->Branch("unfoldCent", &unfoldCent_, "unfoldCent/F"); 

    unfoldXJTree_p = new TTree("unfoldXJTree_p", "");
    unfoldXJTree_p->Branch("recoGammaPt", &recoGammaPt_, "recoGammaPt/F"); 
    unfoldXJTree_p->Branch("recoGammaPhi", &recoGammaPhi_, "recoGammaPhi/F"); 
    unfoldXJTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F"); 
    unfoldXJTree_p->Branch("truthGammaPhi", &truthGammaPhi_, "truthGammaPhi/F"); 
    unfoldXJTree_p->Branch("recoXJ", recoXJ_, ("recoXJ[" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldXJTree_p->Branch("recoJtPt", recoJtPt_, ("recoJtPt[" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldXJTree_p->Branch("recoJtPhi", &recoJtPhi_, "recoJtPhi/F"); 
    unfoldXJTree_p->Branch("truthXJ", &truthXJ_, "truthXJ/F"); 
    unfoldXJTree_p->Branch("truthJtPt", &truthJtPt_, "truthJtPt/F"); 
    unfoldXJTree_p->Branch("truthJtPhi", &truthJtPhi_, "truthJtPhi/F"); 
    unfoldXJTree_p->Branch("unfoldWeight", &unfoldWeight_, "unfoldWeight/F"); 
    unfoldXJTree_p->Branch("unfoldCent", &unfoldCent_, "unfoldCent/F"); 

    unfoldXJJTree_p = new TTree("unfoldXJJTree_p", "");
    unfoldXJJTree_p->Branch("recoGammaPt", &recoGammaPt_, "recoGammaPt/F"); 
    unfoldXJJTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F"); 

    unfoldXJJTree_p->Branch("recoXJJ", recoXJJ_, ("recoXJJ[" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldXJJTree_p->Branch("recoJt1Pt", recoJt1Pt_, ("recoJt1Pt[" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldXJJTree_p->Branch("recoJt2Pt", recoJt2Pt_, ("recoJt2Pt[" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldXJJTree_p->Branch("recoJt1Phi", &recoJt1Phi_, "recoJt1Phi/F"); 
    unfoldXJJTree_p->Branch("recoJt2Phi", &recoJt2Phi_, "recoJt2Phi/F"); 
  
    unfoldXJJTree_p->Branch("truthXJJ", &truthXJJ_, "truthXJJ/F"); 
    unfoldXJJTree_p->Branch("truthJt1Pt", &truthJt1Pt_, "truthJt1Pt/F"); 
    unfoldXJJTree_p->Branch("truthJt2Pt", &truthJt2Pt_, "truthJt2Pt/F"); 
    unfoldXJJTree_p->Branch("truthJt1Phi", &truthJt1Phi_, "truthJt1Phi/F"); 
    unfoldXJJTree_p->Branch("truthJt2Phi", &truthJt2Phi_, "truthJt2Phi/F"); 
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

  TH1F* truthPtOfMatches_Unweighted_p[nMaxPtBins];
  TH1F* truthPtOfMatches_Weighted_p[nMaxPtBins];
  TH2F* truthPtDeltaROfMatches_Unweighted_p[nMaxPtBins];
  TH2F* truthPtDeltaROfMatches_Weighted_p[nMaxPtBins];

  for(Int_t i = 0; i < nJtPtBins; ++i){
    truthPtOfMatches_Unweighted_p[i]= new TH1F(("truthPtOfMatches_Weighted_RecoPt" + std::to_string(i) + "_h").c_str(), ";Truth Jet p_{T};Counts (Weighted)", 50, -0.5, 99.5);
    truthPtOfMatches_Weighted_p[i] = new TH1F(("truthPtOfMatches_Unweighted_RecoPt" + std::to_string(i) + "_h").c_str(), ";Truth Jet p_{T};Counts (Unweighted)", 50, -0.5, 99.5);

    truthPtDeltaROfMatches_Unweighted_p[i]= new TH2F(("truthPtDeltaROfMatches_Weighted_RecoPt" + std::to_string(i) + "_h").c_str(), ";Truth Jet p_{T};#DeltaR_{Reco.,Gen.}", 50, -0.5, 99.5, 20, 0.0, 0.4);
    truthPtDeltaROfMatches_Weighted_p[i] = new TH2F(("truthPtDeltaROfMatches_Unweighted_RecoPt" + std::to_string(i) + "_h").c_str(), ";Truth Jet p_{T};#DeltaR_{Reco.,Gen.}", 50, -0.5, 99.5, 20, 0.0, 0.4);

    truthPtOfMatches_Weighted_p[i]->Sumw2();
    truthPtDeltaROfMatches_Weighted_p[i]->Sumw2();
  }

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
  //Make this a class, one for multijet one for single jet, so you dont have to carry stuff around
  
  TEnv photonPtJtPtVCent_Config;
  photonPtJtPtVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtPtVCent_Config.SetValue("ISMC", isMC);
  photonPtJtPtVCent_Config.SetValue("NBINSX", nJtPtBins);
  photonPtJtPtVCent_Config.SetValue("BINSX", jtPtBinsStrForConfig.c_str());
  photonPtJtPtVCent_Config.SetValue("TITLEX", "Jet p_{T}");  
  photonPtJtPtVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtPtVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtPtVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtXJVCent_Config;
  photonPtJtXJVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtXJVCent_Config.SetValue("ISMC", isMC);
  photonPtJtXJVCent_Config.SetValue("NBINSX", nXJBins);
  photonPtJtXJVCent_Config.SetValue("BINSX", xjBinsStrForConfig.c_str());
  photonPtJtXJVCent_Config.SetValue("TITLEX", "Jet x_{J#gamma}");  
  photonPtJtXJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtXJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtXJVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtDPhiVCent_Config;
  photonPtJtDPhiVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtDPhiVCent_Config.SetValue("ISMC", isMC);
  photonPtJtDPhiVCent_Config.SetValue("NBINSX", nDPhiBins);
  photonPtJtDPhiVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  photonPtJtDPhiVCent_Config.SetValue("TITLEX", "Jet x_{J#gamma}");  
  photonPtJtDPhiVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtDPhiVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtDPhiVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtXJJVCent_Config;
  photonPtJtXJJVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtXJJVCent_Config.SetValue("ISMC", isMC);
  photonPtJtXJJVCent_Config.SetValue("NBINSX", nXJBins);
  photonPtJtXJJVCent_Config.SetValue("BINSX", xjBinsStrForConfig.c_str());
  photonPtJtXJJVCent_Config.SetValue("TITLEX", "Jet #vec{x}_{JJ#gamma}");  
  photonPtJtXJJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtXJJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtXJJVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtAJJVCent_Config;
  photonPtJtAJJVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtAJJVCent_Config.SetValue("ISMC", isMC);
  photonPtJtAJJVCent_Config.SetValue("NBINSX", nAJBins);
  photonPtJtAJJVCent_Config.SetValue("BINSX", ajBinsStrForConfig.c_str());
  photonPtJtAJJVCent_Config.SetValue("TITLEX", "Jet A_{JJ#gamma}");  
  photonPtJtAJJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtAJJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtAJJVCent_Config.SetValue("TITLEY", "Photon p_{T}");
 
  TEnv photonPtJtDPhiJJGVCent_Config;
  photonPtJtDPhiJJGVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtDPhiJJGVCent_Config.SetValue("ISMC", isMC);
  photonPtJtDPhiJJGVCent_Config.SetValue("NBINSX", nDPhiBins);
  photonPtJtDPhiJJGVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  photonPtJtDPhiJJGVCent_Config.SetValue("TITLEX", "Jet #Delta#phi_{JJ#gamma}");  
  photonPtJtDPhiJJGVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtDPhiJJGVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtDPhiJJGVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtDPhiJJVCent_Config;
  photonPtJtDPhiJJVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtDPhiJJVCent_Config.SetValue("ISMC", isMC);
  photonPtJtDPhiJJVCent_Config.SetValue("NBINSX", nDPhiBins);
  photonPtJtDPhiJJVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  photonPtJtDPhiJJVCent_Config.SetValue("TITLEX", "Jet #Delta#phi_{JJ}");  
  photonPtJtDPhiJJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtDPhiJJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtDPhiJJVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv photonPtJtDRJJVCent_Config;
  photonPtJtDRJJVCent_Config.SetValue("IS2DUNFOLD", true);
  photonPtJtDRJJVCent_Config.SetValue("ISMC", isMC);
  photonPtJtDRJJVCent_Config.SetValue("NBINSX", nDPhiBins);
  photonPtJtDRJJVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  photonPtJtDRJJVCent_Config.SetValue("TITLEX", "Jet #DeltaR_{JJ}");  
  photonPtJtDRJJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  photonPtJtDRJJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  photonPtJtDRJJVCent_Config.SetValue("TITLEY", "Photon p_{T}");

  TEnv leadingJtPtJtDPhiJJGVCent_Config;
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("IS2DUNFOLD", true);
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("ISMC", isMC);
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("NBINSX", nDPhiBins);
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("TITLEX", "Jet #Delta#phi_{JJ#gamma}");  
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("NBINSY", nGammaPtBins);
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  leadingJtPtJtDPhiJJGVCent_Config.SetValue("TITLEY", "Leading Jet p_{T}");

  TEnv leadingJtPtJtDPhiJJVCent_Config;
  leadingJtPtJtDPhiJJVCent_Config.SetValue("IS2DUNFOLD", true);
  leadingJtPtJtDPhiJJVCent_Config.SetValue("ISMC", isMC);
  leadingJtPtJtDPhiJJVCent_Config.SetValue("NBINSX", nDPhiBins);
  leadingJtPtJtDPhiJJVCent_Config.SetValue("BINSX", dphiBinsStrForConfig.c_str());
  leadingJtPtJtDPhiJJVCent_Config.SetValue("TITLEX", "Jet #Delta#phi_{JJ}");  
  leadingJtPtJtDPhiJJVCent_Config.SetValue("NBINSY", nGammaPtBins);
  leadingJtPtJtDPhiJJVCent_Config.SetValue("BINSY", gammaPtBinsStrForConfig.c_str());
  leadingJtPtJtDPhiJJVCent_Config.SetValue("TITLEY", "Leading Jet p_{T}");

 
  //";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins

  mixMachine* photonPtJtPtVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtXJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtXJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtAJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiJJGVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDRJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];

  //Separating this declaration as its fundamentally different
  //we will unfold in leading jet pt and azimuthal angle (between vector jets and photon and between jets themselves)
  mixMachine* leadingJtPtJtDPhiJJGVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* leadingJtPtJtDPhiJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC];

  mixMachine* photonPtJtPtVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtXJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtXJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtAJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDPhiJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* photonPtJtDRJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];

  mixMachine* leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];
  mixMachine* leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC];

  TH1F* photonPtVCent_RAW_p[nMaxCentBins][nBarrelAndEC];//This is needed for the purity correction - it is a pure partner of JtPt and JtXJ to be filled on every one of their fills. unclear yet whether it also works for XJJ


  
  TH2F* photonPtJtPtVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_RAW_p[nMaxCentBins][nBarrelAndEC];
  
  TH1F* photonPtVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  
  
  
  TH2F* photonPtJtPtVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];
  

  TH1F* photonJtRecoOverGenVCentJtPt_p[nMaxCentBins][nMaxPtBins];
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
  
  TH2F* photonPtJtPtVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXPURE_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXMIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXMIXCorrected_p[nMaxCentBins][nBarrelAndEC];
  
  
  TH2F* photonPtJtAJJVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXPURE_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXMIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXMIXCorrected_p[nMaxCentBins][nBarrelAndEC];

  TH2F* photonPtJtDPhiJJGVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXPURE_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXMIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[nMaxCentBins][nBarrelAndEC];

  TH2F* photonPtJtDPhiJJVCent_MIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXPURE_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXMIX_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXMIXCorrected_p[nMaxCentBins][nBarrelAndEC];

  TH2F* photonPtJtPtVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXSidebandCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXSidebandCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXSidebandCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXSidebandCorrection_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[nMaxCentBins][nBarrelAndEC];

  TH2F* photonPtJtDPhiJJVCent_MIXSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[nMaxCentBins][nBarrelAndEC];
  

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
  TH2F* photonPtJtPtVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_SUB_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_SUB_p[nMaxCentBins][nBarrelAndEC];

  TH2F* photonPtJtPtVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_SUBSideband_p[nMaxCentBins][nBarrelAndEC];  

  TH1F* photonPtVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtPtVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDRJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];

  TH2F* leadingJtPtJtDPhiJJGVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2F* leadingJtPtJtDPhiJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  
  TH1F* photonPtVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtPtVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];

  TH1F* photonPtVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtPtVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  
  TH2F* photonPtJtPtVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtXJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtAJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJGVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2F* photonPtJtDPhiJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  

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

  TH1F* vzPassing_p = new TH1F("vzPassing_h", ";vz;Counts", 80, -20.0, 20.0);
  TH1F* centPassing_p = new TH1F("centPassing_h", ";Centrality;Counts", 80, -0.5, 79.5);
  TH1F* truPhoPtPassing_p = new TH1F("truPhoPt_h", ";True Photon p_{T};Counts", 255, -5.0, 505.0);
  TH1F* leadingPhoPtPassing_p = new TH1F("leadingPhoPt_h", ";Leading Photon p_{T};Counts", 255, -5.0, 505.0);
  TH1F* leadingPhoEtaPassing_p = new TH1F("leadingPhoEta_h", ";Leading Photon #eta;Counts", 84, -2.1, 2.1);
  TH1F* leadingPhoIsoPassing_p = new TH1F("leadingPhoIso_h", ";Leading Photon Iso;Counts", 80, -10, 30);
  TH1F* leadingPhoTightPassing_p = new TH1F("leadingPhoTight_h", ";Leading Photon Tight;Counts", 2, -0.5, 1.5);

  TH1F* jtEtaPassing_p = new TH1F("jtEtaPassing_h", ";Jet #eta;Counts", 124, -3.1, 3.1);
  TH1F* jtGammaDRPassing_p = new TH1F("jtGammaDRPassing_h", ";Jet-Gamma #DeltaR;Counts", 44, 0.0, 4.4);
  TH1F* jtPtPassing_p = new TH1F("jtPtPassing_h", ";Jet p_{T};Counts", 255, -5, 505);
  TH1F* jtDPhiPassing_p = new TH1F("jtDPhiPassing_h", ";#Delta#phi_{#gamma,jet};Counts", 100, -0.1, TMath::Pi()+0.1);

  mixMachine::mixMode inclusiveFlag = mixMachine::mixMode::INCLUSIVE;
  mixMachine::mixMode multiFlag = mixMachine::mixMode::MULTI;
  if(isPP){
    inclusiveFlag = mixMachine::mixMode::NONE;
    multiFlag = mixMachine::mixMode::NONE;
  }
  
  
  //Type and origin defined here
  //https://gitlab.cern.ch/atlas/athena/-/blob/21.2/PhysicsAnalysis/MCTruthClassifier/MCTruthClassifier/MCTruthClassifierDefs.h
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

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaVCentPt_p[cI][pI] = new TH1F(("photonEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #eta;Counts", nGammaEtaBins, gammaEtaBins);
      photonPhiVCentPt_p[cI][pI] = new TH1F(("photonPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #phi;Counts", nPhiBins, phiBins);

      photonJtDPhiVCentPt_p[cI][pI] = new TH1F(("photonJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";#Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);


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
  

    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      //Since this needs to be fine grain we will do 1.0 GeV bins over the range
      photonPtJtPtVCent_Config.SetValue("ISMC", isMC);
      photonPtJtXJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiJJGVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDRJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtXJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtAJJVCent_Config.SetValue("ISMC", isMC);
      leadingJtPtJtDPhiJJGVCent_Config.SetValue("ISMC", isMC);
      leadingJtPtJtDPhiJJVCent_Config.SetValue("ISMC", isMC);

      photonPtJtPtVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtPtVCent_Config);
      photonPtJtXJVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtXJVCent_Config);
      photonPtJtDPhiVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtDPhiVCent_Config);
      photonPtJtXJJVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtXJJVCent_Config);
      photonPtJtAJJVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtAJJVCent_Config);
      photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDPhiJJGVCent_Config);
      photonPtJtDPhiJJVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDPhiJJVCent_Config);
      photonPtJtDRJJVCent_MixMachine_p[cI][eI] = new mixMachine("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDRJJVCent_Config);


      leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI] = new mixMachine("leadingJtPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &leadingJtPtJtDPhiJJGVCent_Config);
      leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI] = new mixMachine("leadingJtPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr, multiFlag, &leadingJtPtJtDPhiJJVCent_Config);

      photonPtJtPtVCent_Config.SetValue("ISMC", 0);
      photonPtJtXJVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiJJGVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtDRJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtXJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtAJJVCent_Config.SetValue("ISMC", 0);

      leadingJtPtJtDPhiJJGVCent_Config.SetValue("ISMC", 0);
      leadingJtPtJtDPhiJJVCent_Config.SetValue("ISMC", 0);

      photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtPtVCent_Config);
      photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtXJVCent_Config);
      photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtDPhiVCent_Config);
      photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtXJJVCent_Config);
      photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtAJJVCent_Config);
      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDPhiJJGVCent_Config);
      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDPhiJJVCent_Config);
      photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDRJJVCent_Config);

      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("leadingJtPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &leadingJtPtJtDPhiJJGVCent_Config);
      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI] = new mixMachine("leadingJtPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &leadingJtPtJtDPhiJJVCent_Config);

      
      photonPtVCent_RAW_p[cI][eI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAW_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaFinePtBins, gammaFinePtBins);
      
      photonPtJtPtVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtDPhiVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. #Delta#phi_{#gamma,J};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtAJJVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
      photonPtJtDPhiJJGVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
       photonPtJtDPhiJJVCent_RAW_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAW_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
      
      photonPtVCent_RAWSideband_p[cI][eI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAWSideband_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaFinePtBins, gammaFinePtBins);
      
      
      photonPtJtPtVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtDPhiVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. #Delta#phi_{#gamma,J};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
      photonPtJtXJJVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
      photonPtJtAJJVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
      photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
      photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_RAWSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
      

      setSumW2({photonPtVCent_RAW_p[cI][eI], photonPtVCent_RAWSideband_p[cI][eI]});
      setSumW2({photonPtJtPtVCent_RAW_p[cI][eI], photonPtJtXJVCent_RAW_p[cI][eI], photonPtJtDPhiVCent_RAW_p[cI][eI], photonPtJtXJJVCent_RAW_p[cI][eI], photonPtJtAJJVCent_RAW_p[cI][eI], photonPtJtDPhiJJGVCent_RAW_p[cI][eI], photonPtJtDPhiJJVCent_RAW_p[cI][eI], photonPtJtPtVCent_RAWSideband_p[cI][eI], photonPtJtXJVCent_RAWSideband_p[cI][eI], photonPtJtDPhiVCent_RAWSideband_p[cI][eI], photonPtJtXJJVCent_RAWSideband_p[cI][eI], photonPtJtAJJVCent_RAWSideband_p[cI][eI], photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI], photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI]});
    
      if(doMix){
	photonPtJtPtVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	photonPtJtXJJVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXCorrection_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrection_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXCorrected_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrected_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

	photonPtJtXJJVCent_MIXPURE_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXPURE_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXMIX_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIX_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXMIXCorrected_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIXCorrected_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);


	photonPtJtAJJVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXCorrection_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrection_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXCorrected_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrected_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXPURE_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXPURE_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXMIX_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIX_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXMIXCorrected_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIXCorrected_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);


	photonPtJtDPhiJJGVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrection_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXPURE_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXPURE_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXMIX_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIX_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIXCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	photonPtJtDPhiJJVCent_MIX_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIX_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrection_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXPURE_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXPURE_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXMIX_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIX_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXMIXCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXMIXCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	

	photonPtJtPtVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	photonPtJtXJJVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXSidebandCorrection_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrection_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrected_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);


	photonPtJtAJJVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXSidebandCorrection_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrection_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrected_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	
	photonPtJtDPhiJJGVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrection_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	photonPtJtDPhiJJVCent_MIXSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrection_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXSidebandCorrected_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	

	setSumW2({photonPtJtPtVCent_MIX_p[cI][eI], photonPtJtXJVCent_MIX_p[cI][eI], photonPtJtDPhiVCent_MIX_p[cI][eI], photonPtJtXJJVCent_MIX_p[cI][eI], photonPtJtXJJVCent_MIXCorrection_p[cI][eI], photonPtJtXJJVCent_MIXCorrected_p[cI][eI], photonPtJtXJJVCent_MIXPURE_p[cI][eI], photonPtJtXJJVCent_MIXMIX_p[cI][eI], photonPtJtXJJVCent_MIXMIXCorrected_p[cI][eI], photonPtJtAJJVCent_MIX_p[cI][eI], photonPtJtAJJVCent_MIXCorrection_p[cI][eI], photonPtJtAJJVCent_MIXCorrected_p[cI][eI], photonPtJtAJJVCent_MIXPURE_p[cI][eI], photonPtJtAJJVCent_MIXMIX_p[cI][eI], photonPtJtAJJVCent_MIXMIXCorrected_p[cI][eI], photonPtJtDPhiJJGVCent_MIX_p[cI][eI], photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI], photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI], photonPtJtDPhiJJGVCent_MIXPURE_p[cI][eI], photonPtJtDPhiJJGVCent_MIXMIX_p[cI][eI], photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[cI][eI], photonPtJtDPhiJJVCent_MIX_p[cI][eI], photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI], photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI], photonPtJtDPhiJJVCent_MIXPURE_p[cI][eI], photonPtJtDPhiJJVCent_MIXMIX_p[cI][eI], photonPtJtDPhiJJVCent_MIXMIXCorrected_p[cI][eI], photonPtJtPtVCent_MIXSideband_p[cI][eI], photonPtJtXJVCent_MIXSideband_p[cI][eI], photonPtJtDPhiVCent_MIXSideband_p[cI][eI], photonPtJtXJJVCent_MIXSideband_p[cI][eI], photonPtJtXJJVCent_MIXSidebandCorrection_p[cI][eI], photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI], photonPtJtAJJVCent_MIXSideband_p[cI][eI], photonPtJtAJJVCent_MIXSidebandCorrection_p[cI][eI], photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI], photonPtJtDPhiJJGVCent_MIXSideband_p[cI][eI], photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[cI][eI], photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI], photonPtJtDPhiJJVCent_MIXSideband_p[cI][eI], photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[cI][eI], photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI]});
      }
      
      if(doMix || isPP){
	photonPtJtPtVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_SUB_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUB_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	
	photonPtJtPtVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. #Delta#phi_{J#gamam};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

	photonPtJtAJJVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);

	photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_SUBSideband_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	photonPtVCent_PURCORR_p[cI][eI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_PURCORR_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaFinePtBins, gammaFinePtBins);

	photonPtJtPtVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtXJJVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	photonPtJtAJJVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDPhiJJVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	photonPtJtDRJJVCent_PURCORR_p[cI][eI] = new TH2F(("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #DeltaR_{JJ};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][eI] = new TH2F(("leadingJtPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][eI] = new TH2F(("leadingJtPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	setSumW2({photonPtJtPtVCent_SUB_p[cI][eI], photonPtJtXJVCent_SUB_p[cI][eI], photonPtJtDPhiVCent_SUB_p[cI][eI], photonPtJtXJJVCent_SUB_p[cI][eI], photonPtJtAJJVCent_SUB_p[cI][eI], photonPtJtDPhiJJGVCent_SUB_p[cI][eI], photonPtJtDPhiJJVCent_SUB_p[cI][eI], photonPtJtPtVCent_SUBSideband_p[cI][eI], photonPtJtXJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiVCent_SUBSideband_p[cI][eI], photonPtJtXJJVCent_SUBSideband_p[cI][eI], photonPtJtAJJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI], photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI], photonPtJtDRJJVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][eI]});

	if(isMC){
	  photonPtVCent_TRUTH_p[cI][eI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_TRUTH_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaFinePtBins, gammaFinePtBins);
	  photonPtJtPtVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);


	  photonPtVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_TRUTHMATCHEDRECO_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaFinePtBins, gammaFinePtBins);
	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	  photonPtJtPtVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_MIXBKGD_p[cI][eI] = new TH2F(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);

	  photonPtJtAJJVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_MIXBKGD_p[cI][eI] = new TH2F(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);

	  photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI] = new TH2F(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	  photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI] = new TH2F(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dphiBins, nGammaPtBins, gammaPtBins);

	  setSumW2({photonPtJtPtVCent_TRUTH_p[cI][eI], photonPtJtXJVCent_TRUTH_p[cI][eI], photonPtJtDPhiVCent_TRUTH_p[cI][eI], photonPtJtXJJVCent_TRUTH_p[cI][eI], photonPtJtAJJVCent_TRUTH_p[cI][eI], photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI], photonPtJtDPhiJJVCent_TRUTH_p[cI][eI], photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtPtVCent_PUREBKGD_p[cI][eI], photonPtJtXJVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiVCent_PUREBKGD_p[cI][eI], photonPtJtXJJVCent_PUREBKGD_p[cI][eI], photonPtJtXJJVCent_MIXBKGD_p[cI][eI], photonPtJtAJJVCent_PUREBKGD_p[cI][eI], photonPtJtAJJVCent_MIXBKGD_p[cI][eI], photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI], photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI]});
	}
      }
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

  Int_t sampleTag;
  std::map<Int_t, Int_t> sampleTagCounter;
  std::vector<Float_t> uniqueMinPthats;

  sampleHandler sHandler;
  sHandler.Init("mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217");//we wont actually use this sample but i need to initialize w/ something

  
  if(isMC){
    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("sampleTag", 1);
    
    inTree_p->SetBranchAddress("sampleTag", &sampleTag);
    
    for(Int_t entry = 0; entry < inTree_p->GetEntries(); ++entry){
      inTree_p->GetEntry(entry);
      
      ++(sampleTagCounter[sampleTag]);
    }
    
    
    std::cout << "SAMPLE TAG, COUNTS:" << std::endl;
    for(auto const& tempTag : sampleTagCounter){
      std::cout << " " << tempTag.first << ", " << tempTag.second << std::endl;
      
      Float_t pthatMin = sHandler.GetMinPthat(tempTag.first);
      uniqueMinPthats.push_back(pthatMin);
    }
    std::sort(std::begin(uniqueMinPthats), std::end(uniqueMinPthats));
    
    std::cout << "UNIQUE PTHATMINS: " << std::endl;
    for(unsigned int uI = 0; uI < uniqueMinPthats.size(); ++uI){
      std::cout << " " << uniqueMinPthats[uI] << std::endl;
    }
    uniqueMinPthats.push_back(10000.);
  }

  outFile_p->cd();


  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  runNumber_p = new TH1F(("runNumber_" + systStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);
  if(!isMC) lumiFractionPerRun_p = new TH1F(("lumiFractionPerRun_" + systStr + "_h").c_str(), ";Run;Fraction of Lumiblocks", nRunBins+1, runMinF, runMaxF);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  

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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
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

  /*  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;
  std::vector<int>* truth_type_p=nullptr;
  std::vector<int>* truth_origin_p=nullptr;
  */

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
  
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;  
  std::vector<float>* photon_etcone_p=nullptr;

  std::vector<float>* aktRhi_etajes_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_etajes_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_etajes_jet_phi_p=nullptr;
  std::vector<std::vector<float>* > aktRhi_etajes_jet_pt_sysJES_p;
  std::vector<std::vector<float>* > aktRhi_etajes_jet_pt_sysJER_p;
  std::vector<float>* aktRhi_insitu_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_insitu_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_insitu_jet_phi_p=nullptr;
  std::vector<int>* aktRhi_truthpos_p=nullptr;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t jI = 0; jI < nJESSys; ++jI){aktRhi_etajes_jet_pt_sysJES_p.push_back(nullptr);}
  for(Int_t jI = 0; jI < nJERSys; ++jI){aktRhi_etajes_jet_pt_sysJER_p.push_back(nullptr);}  

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

    //Fix to prev. bug - used etajes instead of fully corrected jets
    //since for the mixing we are working with overlay, always take insitu corrections
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_pt").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_eta").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_phi").c_str(), 1);

    mixTree_p->SetBranchAddress("vert_z", &vert_z_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_pt").c_str(), &aktRhi_insitu_jet_pt_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_eta").c_str(), &aktRhi_insitu_jet_eta_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_phi").c_str(), &aktRhi_insitu_jet_phi_p);

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
    //    if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, (ULong64_t)nMaxEvt*100);
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
      for(unsigned int jI = 0; jI < aktRhi_insitu_jet_pt_p->size(); ++jI){
	if(aktRhi_insitu_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	if(aktRhi_insitu_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;

	Float_t jtEtaToCut = aktRhi_insitu_jet_eta_p->at(jI);
	if(jtEtaBinsDoAbs) jtEtaToCut = TMath::Abs(jtEtaToCut);

	if(jtEtaToCut <= jtEtaBinsLow) continue;
	if(jtEtaToCut >= jtEtaBinsHigh) continue;

	TLorentzVector temp;
	temp.SetPtEtaPhiM(aktRhi_insitu_jet_pt_p->at(jI), aktRhi_insitu_jet_eta_p->at(jI), aktRhi_insitu_jet_phi_p->at(jI), 0.0);

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
  inTree_p->SetBranchStatus("sampleTag", 1);

  if(isMC){
    inTree_p->SetBranchStatus("pthat", 1);
    inTree_p->SetBranchStatus("sampleWeight", 1);
    if(!isPP) inTree_p->SetBranchStatus("ncollWeight", 1);
    inTree_p->SetBranchStatus("fullWeight", 1);

    /*
    inTree_p->SetBranchStatus("truth_pt", 1);
    inTree_p->SetBranchStatus("truth_eta", 1);
    inTree_p->SetBranchStatus("truth_phi", 1);
    inTree_p->SetBranchStatus("truth_pdg", 1);
    inTree_p->SetBranchStatus("truth_origin", 1);
    inTree_p->SetBranchStatus("truth_type", 1);
    */

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
  inTree_p->SetBranchStatus(("photon_etcone" + std::to_string(isolationR) + "0").c_str(), 1);
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(isMC){
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_etajes_jet_eta").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_etajes_jet_phi").c_str(), 1);

    for(Int_t jI = 0; jI < nJESSys; ++jI){
      inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt_sys_JES_" + std::to_string(jI)).c_str(), 1);
    }

    for(Int_t jI = 0; jI < nJERSys; ++jI){
      inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt_sys_JER_" + std::to_string(jI)).c_str(), 1);
    }

  }
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_insitu_jet_phi").c_str(), 1);
  
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
  inTree_p->SetBranchAddress("sampleTag", &sampleTag);

  if(isMC){
    inTree_p->SetBranchAddress("pthat", &pthat);
    inTree_p->SetBranchAddress("sampleWeight", &sampleWeight);
    if(!isPP) inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
    inTree_p->SetBranchAddress("fullWeight", &fullWeight);

    /*
    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
    inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);
    inTree_p->SetBranchAddress("truth_origin", &truth_origin_p);
    inTree_p->SetBranchAddress("truth_type", &truth_type_p);
    */

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
  inTree_p->SetBranchAddress(("photon_etcone" + std::to_string(isolationR) + "0").c_str(), &photon_etcone_p);

  if(isMC){
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt").c_str(), &aktRhi_etajes_jet_pt_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_etajes_jet_eta").c_str(), &aktRhi_etajes_jet_eta_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_etajes_jet_phi").c_str(), &aktRhi_etajes_jet_phi_p);

    for(Int_t jI = 0; jI < nJESSys; ++jI){
      inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt_sys_JES_" + std::to_string(jI)).c_str(), &(aktRhi_etajes_jet_pt_sysJES_p[jI]));
    }

    for(Int_t jI = 0; jI < nJERSys; ++jI){
      inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_etajes_jet_pt_sys_JER_" + std::to_string(jI)).c_str(), &(aktRhi_etajes_jet_pt_sysJER_p[jI]));
    }
  }
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_pt").c_str(), &aktRhi_insitu_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_eta").c_str(), &aktRhi_insitu_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_insitu_jet_phi").c_str(), &aktRhi_insitu_jet_phi_p);

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

  std::map<int, int> eventCounter;
  for(int i = 0; i < 10; ++i){
    eventCounter[i] = 0;
  }


  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    double vert_z = vert_z_p->at(0);
    vert_z /= mmToCMDivFactor;

    //Cut 1: z vertex
    if(vert_z <= -15. || vert_z >= 15.) continue;      

    vzPassing_p->Fill(vert_z);

    ++(eventCounter[0]);

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
    
    Float_t minPthat = -1.0;
    Float_t maxPthat = -1.0;
   
    if(isMC){
      minPthat = sHandler.GetMinPthat(sampleTag);
      int minPtHatPos = -1;    
      for(unsigned int pI = 0; pI < uniqueMinPthats.size()-1; ++pI){
	if(TMath::Abs(minPthat - uniqueMinPthats[pI]) < 1.0){
	  minPtHatPos = pI;
	  break;
	}
      }
      
      if(minPtHatPos < 0){
	std::cout << "Min pthat of event, " << minPthat << ", not found in vector of unique minimum pthats. Please fix. return 1" << std::endl;
	return 1;
      }
      maxPthat = uniqueMinPthats[minPtHatPos+1];
    }

    //Check prescale is 1 in case I made a mistake on first unprescaled
    for(unsigned int hI = 0; hI < hltPrescaleVect.size(); ++hI){
      if(TMath::Abs((*(hltPrescaleVect[hI])) - 1.0) > hltPrescaleDelta){
	std::cout << "WARNING - prescale for \'" << hltList[hI] << "\' has non-unity value, \'" << (*(hltPrescaleVect[hI])) << "\'." << std::endl;
      }
    }

    //Grab the centrality position in your array of histograms for this event
    //in the case of p+p, this is always zero
    Int_t centPos = -1;
    Double_t cent = -1;
    if(!isPP){
      cent = centTable.GetCent(fcalA_et + fcalC_et);
      centPos = ghostPos(centBins, cent, true, doGlobalDebug);
    }
    else centPos = 0;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    //Cut 2: Centrality within selected range
    if(centPos < 0){
      bool vectContainsCent = vectContainsInt((Int_t)cent, &skippedCent);

      if(!vectContainsCent){
	std::cout << "gdjNTupleToHist Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
	skippedCent.push_back((Int_t)cent);
      }
      
      continue;
    }

    centPassing_p->Fill(cent);

    ++(eventCounter[1]);

    //In data, since we deal in unprescaled data, the event-by-event weights should always be 1
    if(!isMC) fullWeight = 1.0;

    //2021.08.12: 
    //Cut 3: MC only cut on is there a good truth photon
    // CFM NOTE YOU NEED TO FIX THIS W/ YEONJU'S CODE:
    // https://github.com/YeonjuGo/GDJ/blob/master/src/gdjNtuplePreProc_phoTaggedJetRaa.C#L1362-L1375
    // We need to veto on events w/ undressed background photons sneaking thru - we cannot handle them in MC w/ purity correction as in data
    if(isMC){
      /*
      Int_t leadingTruthPhotonPos = -1;
      Int_t leadingTruthPhotonPt = -1;
      for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	if(truth_pdg_p->at(tI) != 22) continue;
	
	if(truth_pt_p->at(tI) > leadingTruthPhotonPt){
	  leadingTruthPhotonPos = tI;
	  leadingTruthPhotonPt = truth_pt_p->at(tI);
	}
      }
      
      if(leadingTruthPhotonPos < 0) continue; //You need at least 1 leading photon
      if(truth_type_p->at(leadingTruthPhotonPos) != 14) continue; // Must be prompt iso
      if(truth_origin_p->at(leadingTruthPhotonPos) != 37) continue; // Must be origin prompt photon
      */

      
      //Make sure the photon is w/in the DP range
      if(truthPhotonPt < minPthat) continue;
      if(truthPhotonPt >= maxPthat) continue;

      truPhoPtPassing_p->Fill(truthPhotonPt);
    }


    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    ++(eventCounter[2]);
  
    Float_t mixWeight = fullWeight/(double)nMixEvents;
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
  
    //Check for good photon
    Float_t leadingPhoPt = -1.0;
    Int_t leadingPhoPos = -1;

    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      if(leadingPhoPt < photon_pt_p->at(pI)){
	leadingPhoPt = photon_pt_p->at(pI);
	leadingPhoPos = pI;
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << entry << std::endl; 

    //Cut 4: Is there at least one photon candidate within acceptance 
    if(leadingPhoPos < 0) continue; // is there at least one possible photon

    ++(eventCounter[3]);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << entry << std::endl; 

    //Get rid of this - Messing w/ under and overflow handling

    //Cut 5: Does leading photon pass pT cuts
    if(photon_pt_p->at(leadingPhoPos) < gammaPtBins[0]) continue;
    if(photon_pt_p->at(leadingPhoPos) >= gammaPtBins[nGammaPtBins]) continue;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << entry << std::endl; 

    //We need to correct the photon isolation according to functions in photonUtil.h
    Float_t correctedIso = photon_etcone_p->at(leadingPhoPos);
    if(doPtIsoCorrection){
      if(doCentIsoCorrection) correctedIso = getCorrectedPhotonIsolation(isPP, correctedIso, photon_pt_p->at(leadingPhoPos), photon_eta_p->at(leadingPhoPos), cent);
      else correctedIso = getPtCorrectedPhotonIsolation(correctedIso, photon_pt_p->at(leadingPhoPos), photon_eta_p->at(leadingPhoPos));
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << entry << std::endl; 

    
    //Cut 6: Leading reco photon passes full photon selection
    //Isolation (corrected), either in sideband or 'good' region, in good eta
    if(!isGoodPhoton(isPP, doPtIsoCorrection, sidebandType, photon_tight_p->at(leadingPhoPos), correctedIso, photon_eta_p->at(leadingPhoPos))) continue;

    ++(eventCounter[4]);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    //Cut 7: isMC we must make sure the leading reco photon matches a valid truth
    if(isMC){
      if(getDR(truthPhotonEta, truthPhotonPhi, photon_eta_p->at(leadingPhoPos), photon_phi_p->at(leadingPhoPos)) > 0.2) continue;
    }
    

    leadingPhoPtPassing_p->Fill(photon_pt_p->at(leadingPhoPos));
    leadingPhoEtaPassing_p->Fill(photon_eta_p->at(leadingPhoPos));
    leadingPhoIsoPassing_p->Fill(correctedIso);
    leadingPhoTightPassing_p->Fill(photon_tight_p->at(leadingPhoPos));

    bool isSideband = isSidebandPhoton(isPP, doPtIsoCorrection, sidebandType, photon_tight_p->at(leadingPhoPos), correctedIso);
    const Int_t phoPos = leadingPhoPos;
     
    Float_t etaValMain = photon_eta_p->at(phoPos);
    Float_t etaValSub = etaValMain;
    if(gammaEtaBinsDoAbs) etaValMain = TMath::Abs(etaValMain);
    if(gammaEtaBinsSubDoAbs) etaValSub = TMath::Abs(etaValSub);

    //Cut 8: Redundant cut on gamma eta position
    if(etaValMain <= gammaEtaBins[0]) continue;
    if(etaValMain >= gammaEtaBins[nGammaEtaBins]) continue;

    ++(eventCounter[5]);  //Reject outside of eta bounds

    //Position of the photon in barrel or endcap - some corrections were originally formulated with these two regions separate so they are handled separately at this stage
    Int_t barrelOrECPos = 0;
    if(TMath::Abs(etaValMain) > 1.52) barrelOrECPos = 1;
      
    //Position of the photon in histogram pT,Eta arrays
    Int_t ptPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(phoPos), true, doGlobalDebug);
    Int_t etaPos = ghostPos(nGammaEtaBinsSub, gammaEtaBinsSub, etaValSub, true, doGlobalDebug);

  
    if(etaPos >= 0){
      if(!isSideband){
	fillTH1(photonPtVCentEta_p[centPos][etaPos], photon_pt_p->at(phoPos), fullWeight);
	fillTH1(photonPtVCentEta_p[centPos][nGammaEtaBinsSub], photon_pt_p->at(phoPos), fullWeight);
      }

      if(isMC && truthPhotonPt > 0){
	if(!isSideband){
	  fillTH2(photonGenResVCentEta_p[centPos][etaPos], photon_pt_p->at(phoPos), truthPhotonPt, fullWeight);
	  fillTH2(photonGenResVCentEta_p[centPos][nGammaEtaBinsSub], photon_pt_p->at(phoPos), truthPhotonPt, fullWeight);	  
	}
      }
    }
    
    if(!isMC){
      ++(gammaCountsPerPtCent[ptPos][centPos]);
      ++(gammaCountsPerPtCent[nGammaPtBins][centPos]);
    }
    else{
      gammaCountsPerPtCent[ptPos][centPos] += fullWeight;
      gammaCountsPerPtCent[nGammaPtBins][centPos] += fullWeight;
    }
      
    if(!isSideband){
      fillTH1(photonEtaVCentPt_p[centPos][ptPos], etaValMain, fullWeight);
      fillTH1(photonPhiVCentPt_p[centPos][ptPos], photon_phi_p->at(phoPos), fullWeight);
      fillTH2(photonEtaPhiVCentPt_p[centPos][ptPos], etaValMain, photon_phi_p->at(phoPos), fullWeight);
	
      fillTH1(photonEtaVCentPt_p[centPos][nGammaPtBins], etaValMain, fullWeight);
      fillTH1(photonPhiVCentPt_p[centPos][nGammaPtBins], photon_phi_p->at(phoPos), fullWeight);
      fillTH2(photonEtaPhiVCentPt_p[centPos][nGammaPtBins], etaValMain, photon_phi_p->at(phoPos), fullWeight);
    }

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
    
    //Start w/ unpartnered truth
    if(isMC && !isSideband){
      for(unsigned int jI = 0; jI < aktR_truth_jet_pt_p->size(); ++jI){
	if(aktR_truth_jet_recopos_p->at(jI) >= 0) continue;
	
	recoGammaPt_ = photon_pt_p->at(phoPos);
	recoGammaPhi_ = photon_phi_p->at(phoPos);
	truthGammaPt_ = truthPhotonPt;	
	truthGammaPhi_ = truthPhotonPhi;

	for(Int_t jsI = 0; jsI < nJetSysAndNom; ++jsI){
	  recoXJ_[jsI] = -1;
	  recoJtPt_[jsI] = -1;
	}
	recoJtPhi_ = -100;
	
	truthXJ_ = aktR_truth_jet_pt_p->at(jI)/truthPhotonPt;
	truthJtPt_ = aktR_truth_jet_pt_p->at(jI);
	truthJtPhi_ = aktR_truth_jet_phi_p->at(jI);
	
	unfoldWeight_ = fullWeight;
	unfoldCent_ = cent;
	
	unfoldXJTree_p->Fill();
      }
    }
  
    //Fill for photons
    if(!isSideband){
      if(isMC){
	recoGammaPt_ = photon_pt_p->at(phoPos);
	truthGammaPt_ = truthPhotonPt;

	unfoldWeight_ = fullWeight;
	unfoldCent_ = cent;

	unfoldPhotonPtTree_p->Fill();
      }

      fillTH1(photonPtVCent_RAW_p[centPos][barrelOrECPos], photon_pt_p->at(phoPos), fullWeight);
      if(isMC){
	fillTH1(photonPtVCent_TRUTH_p[centPos][barrelOrECPos], truthPhotonPt, fullWeight);
	fillTH1(photonPtVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], photon_pt_p->at(phoPos), fullWeight);
      }
    }
    else fillTH1(photonPtVCent_RAWSideband_p[centPos][barrelOrECPos], photon_pt_p->at(phoPos), fullWeight);

    //Now go thru reco jets 
    for(unsigned int jI = 0; jI < aktRhi_insitu_jet_pt_p->size(); ++jI){
      Float_t ptToUse = aktRhi_insitu_jet_pt_p->at(jI);
      Float_t etaToUse = aktRhi_insitu_jet_eta_p->at(jI);
      Float_t phiToUse = aktRhi_insitu_jet_phi_p->at(jI);
      
      if(isMC){
	//If MC w/ a truth match we use etajes level of the correction
	if(aktRhi_truthpos_p->at(jI) >= 0){
	  ptToUse = aktRhi_etajes_jet_pt_p->at(jI);
	  etaToUse = aktRhi_etajes_jet_eta_p->at(jI);
	  phiToUse = aktRhi_etajes_jet_phi_p->at(jI);
	}	  
      }
      
      Float_t jtEtaToCut = etaToUse;
      if(jtEtaBinsDoAbs) jtEtaToCut = TMath::Abs(jtEtaToCut);
      
      //continue if eta is outside of cut bounds
      //Jet Cut 1: Jet must be within eta acceptance      
      if(jtEtaToCut <= jtEtaBinsLow) continue;
      if(jtEtaToCut >= jtEtaBinsHigh) continue;

      //Fill to check that the cut is applied correctly
      jtEtaPassing_p->Fill(etaToUse);
      
      //Jet Cut 2: continue if jet is on top of the photon (usual choice is dR)
      Float_t dR = getDR(etaToUse, phiToUse, photon_eta_p->at(phoPos), photon_phi_p->at(phoPos));
      if(dR < gammaExclusionDR) continue;
      
      //Fill to check that the cut is applied correctly
      jtGammaDRPassing_p->Fill(dR);

      if(recoJtPtMin > ptToUse) recoJtPtMin = ptToUse;//This is just to check what the input minimum recojtpt is
      
      if(isMC){
	int pos = aktRhi_truthpos_p->at(jI);
	if(pos >= 0){
	  if(aktR_truth_jet_pt_p->at(pos) >= jtPtBinsLow && aktR_truth_jet_pt_p->at(pos) < jtPtBinsHigh){
	    Int_t genJtPtPos = ghostPos(nJtPtBins, jtPtBins, aktR_truth_jet_pt_p->at(pos), true, doGlobalDebug);		
	    if(!isSideband) fillTH1(photonJtRecoOverGenVCentJtPt_p[centPos][genJtPtPos], ptToUse/aktR_truth_jet_pt_p->at(pos), fullWeight);
	  }
	}
      }
     
      //Jet Cut 3: Continue if reco is outside fiducial region
      if(ptToUse < jtPtBinsLow) continue;
      if(ptToUse >= jtPtBinsHigh) continue;
      
      jtPtPassing_p->Fill(ptToUse);

      TLorentzVector tempJet;
      tempJet.SetPtEtaPhiM(ptToUse, etaToUse, phiToUse, 0.0);
      goodJets.push_back(tempJet);
      goodJetsTPos.push_back(-1);
      if(isMC){
	if(aktRhi_truthpos_p->at(jI) >= 0){  
	  if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
	    goodJetsTPos[goodJetsTPos.size()-1] = aktRhi_truthpos_p->at(jI);
	  }
	}
      }

      Float_t dPhi = TMath::Abs(getDPHI(phiToUse, photon_phi_p->at(phoPos)));
      if(!isSideband){
	fillTH1(photonJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	fillTH1(photonJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);
      }
      
      if(isMC){
	if(aktRhi_truthpos_p->at(jI) >= 0){
	  if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
	    if(!isSideband){
	      fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	      fillTH1(photonGenMatchedJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);
	    }
	  }
	}
      }

            
      //Fill for the inclusive DPhi
      if(!isSideband){
	fillTH2(photonPtJtDPhiVCent_RAW_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), fullWeight);
	photonPtJtDPhiVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(dPhi, photon_pt_p->at(phoPos), fullWeight);	
      }
      else{
	fillTH2(photonPtJtDPhiVCent_RAWSideband_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), fullWeight);
	photonPtJtDPhiVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(dPhi, photon_pt_p->at(phoPos), fullWeight);	
      }
    
      if(!isSideband){
	if(isMC){	  
	  if(truthPhotonPt > 0){
	    int truthPos = aktRhi_truthpos_p->at(jI);
	      
	    if(truthPos >= 0){
	      Float_t truthPt = aktR_truth_jet_pt_p->at(truthPos);
	      Float_t truthPhi = aktR_truth_jet_phi_p->at(truthPos);
	    
	      if(truthPt > assocGenMinPt){
		Float_t dPhiTruth = TMath::Abs(getDPHI(truthPhi, truthPhotonPhi));

		fillTH2(photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), fullWeight);
		fillTH2(photonPtJtDPhiVCent_TRUTH_p[centPos][barrelOrECPos], dPhiTruth, truthPhotonPt, fullWeight);

		photonPtJtDPhiVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(dPhi, photon_pt_p->at(phoPos), fullWeight);
		photonPtJtDPhiVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(dPhiTruth, truthPhotonPt, fullWeight);
	      }
	      else fillTH2(photonPtJtDPhiVCent_PUREBKGD_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), fullWeight);
	    }
	    else fillTH2(photonPtJtDPhiVCent_PUREBKGD_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), fullWeight);
	  }
	}
      }
           
      //Require jets be 'back-to-back' w/ the photon     
      if(dPhi >= gammaJtDPhiCut){
	jtDPhiPassing_p->Fill(dPhi);
	
	goodJetsDPhi.push_back(tempJet);
	goodJetsDPhiTPos.push_back(-1);
	if(isMC){
	  if(aktRhi_truthpos_p->at(jI) >= 0){  
	    if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
	      goodJetsDPhiTPos[goodJetsDPhiTPos.size()-1] = aktRhi_truthpos_p->at(jI);
	    }
	  }
	}

	  
	if(!isSideband){
	  fillTH1(photonJtPtVCentPt_p[centPos][ptPos], ptToUse, fullWeight);
	  fillTH1(photonJtPtVCentPt_p[centPos][nGammaPtBins], ptToUse, fullWeight);
	  fillTH1(photonJtEtaVCentPt_p[centPos][ptPos], etaToUse, fullWeight);
	  fillTH1(photonJtEtaVCentPt_p[centPos][nGammaPtBins], etaToUse, fullWeight);
	  fillTH1(photonJtXJVCentPt_p[centPos][ptPos], ptToUse/photon_pt_p->at(phoPos), fullWeight);
	  fillTH1(photonJtXJVCentPt_p[centPos][nGammaPtBins], ptToUse/photon_pt_p->at(phoPos), fullWeight);
	}	
      
	if(!isSideband){       	  
	  photonPtJtPtVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(ptToUse, photon_pt_p->at(phoPos), fullWeight);
	  photonPtJtXJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

	  fillTH2(photonPtJtPtVCent_RAW_p[centPos][barrelOrECPos], ptToUse, photon_pt_p->at(phoPos), fullWeight);
	  fillTH2(photonPtJtXJVCent_RAW_p[centPos][barrelOrECPos], ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	  
	  if(isMC){	  
	    if(truthPhotonPt > 0){
	      int truthPos = aktRhi_truthpos_p->at(jI);
	      
	      if(truthPos >= 0){
		Float_t truthPt = aktR_truth_jet_pt_p->at(truthPos);
		Float_t truthEta = aktR_truth_jet_eta_p->at(truthPos);
		Float_t truthPhi = aktR_truth_jet_phi_p->at(truthPos);
				
		if(centPos == 0){
		  if(ptToUse >= jtPtBins[0] && ptToUse < jtPtBins[nJtPtBins]){
		    if(photon_pt_p->at(phoPos) >= gammaPtBins[3] && photon_pt_p->at(phoPos) < gammaPtBins[5]){
		      Int_t recoPtPos = ghostPos(nJtPtBins, jtPtBins, ptToUse, true, doGlobalDebug);		      
		      Float_t deltaRTruthReco = getDR(truthEta, truthPhi, etaToUse, phiToUse);
		      
		      truthPtOfMatches_Unweighted_p[recoPtPos]->Fill(truthPt);
		      truthPtOfMatches_Weighted_p[recoPtPos]->Fill(truthPt, fullWeight);
		      
		      truthPtDeltaROfMatches_Unweighted_p[recoPtPos]->Fill(truthPt, deltaRTruthReco);
		      truthPtDeltaROfMatches_Weighted_p[recoPtPos]->Fill(truthPt, deltaRTruthReco, fullWeight);
		      
		      if(ptToUse >= jtPtBins[0] && ptToUse < jtPtBins[1] && false){
			std::cout << "EVENT OF INTEREST: " << std::endl;
			std::cout << " Entry: " << entry << std::endl;
			//			std::cout << " RUN, LUMI, EVENT: " << runNumber << ", " << lumiBlock << std::endl;
			std::cout << " Photon pt, eta, phi: " << photon_pt_p->at(phoPos) << ", " << photon_eta_p->at(phoPos) << ", " << photon_phi_p->at(phoPos) << std::endl;
			std::cout << " Reco. jet pt, eta, phi: " << ptToUse << ", " << etaToUse << ", " << phiToUse << std::endl;
			std::cout << " Truth jet pt, eta, phi: " << truthPt << ", " << truthEta << ", " << truthPhi << std::endl;
			std::cout << std::endl;
		      }
		    }
		  }
		}
	   		
		if(truthPt > assocGenMinPt){
		  /*
		  if(truthPt < jtPtBins[0]) truthPt = jtPtBins[0] + 0.01;
		  else if(truthPt > jtPtBins[nJtPtBins]) truthPt = jtPtBins[nJtPtBins] - 0.01;
		  */

		  fillTH2(photonPtJtPtVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], ptToUse, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtXJVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

		  fillTH2(photonPtJtPtVCent_TRUTH_p[centPos][barrelOrECPos], truthPt, truthPhotonPt, fullWeight);
		  fillTH2(photonPtJtXJVCent_TRUTH_p[centPos][barrelOrECPos], truthPt/truthPhotonPt, truthPhotonPt, fullWeight);

		  photonPtJtPtVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(ptToUse, photon_pt_p->at(phoPos), fullWeight);
		  photonPtJtXJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

		  photonPtJtPtVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthPt, truthPhotonPt, fullWeight);
		  photonPtJtXJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthPt/truthPhotonPt, truthPhotonPt, fullWeight);		  
		}
		else{
		  fillTH2(photonPtJtPtVCent_PUREBKGD_p[centPos][barrelOrECPos], ptToUse, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtXJVCent_PUREBKGD_p[centPos][barrelOrECPos], ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		}
	      }
	      else{
		fillTH2(photonPtJtPtVCent_PUREBKGD_p[centPos][barrelOrECPos], ptToUse, photon_pt_p->at(phoPos), fullWeight);
		fillTH2(photonPtJtXJVCent_PUREBKGD_p[centPos][barrelOrECPos], ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	      }
	    }
	  }	    
	}
	else{	  
	  fillTH2(photonPtJtPtVCent_RAWSideband_p[centPos][barrelOrECPos], ptToUse, photon_pt_p->at(phoPos), fullWeight);
	  fillTH2(photonPtJtXJVCent_RAWSideband_p[centPos][barrelOrECPos], ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

	  photonPtJtPtVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(ptToUse, photon_pt_p->at(phoPos), fullWeight);
	  photonPtJtXJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(ptToUse/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	}
	
      
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
      
	if(isMC && truthPhotonPt > 0){
	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

	  recoGammaPt_ = photon_pt_p->at(phoPos);
	  recoXJ_[0] = ptToUse/recoGammaPt_;
	  recoJtPt_[0] = ptToUse;
	  recoJtPhi_ = phiToUse;

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

	  if(aktRhi_truthpos_p->at(jI) >= 0){	      
	    Float_t tempTruthJtPt = aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI));
	    Float_t tempTruthJtPhi= aktR_truth_jet_phi_p->at(aktRhi_truthpos_p->at(jI));
	    Float_t tempTruthJtEta = aktR_truth_jet_eta_p->at(aktRhi_truthpos_p->at(jI));
	      
	    TLorentzVector tempRecoJet1, tempRecoJet2;
	    TLorentzVector tempTruthJet1, tempTruthJet2;
	    tempTruthJet1.SetPtEtaPhiM(tempTruthJtPt, tempTruthJtEta, tempTruthJtPhi, 0.0);
	    tempRecoJet1.SetPtEtaPhiM(ptToUse, etaToUse, phiToUse, 0.0);
		
	    //		    if(tempTruthJtPt >= jtPtBins[0] && tempTruthJtPt < jtPtBins[nJtPtBins]){
	    truthGammaPt_ = truthPhotonPt;
	    truthXJ_ = tempTruthJtPt/truthGammaPt_;
	    truthJtPt_ = tempTruthJtPt;
	    truthJtPhi_ = tempTruthJtPhi;

	    unfoldWeight_ = fullWeight;
	    if(isPP) unfoldCent_ = 0;
	    else unfoldCent_ = cent;

	    for(Int_t jsI = 0; jsI < nJESSys; ++jsI){
	      recoXJ_[1 + jsI] = aktRhi_etajes_jet_pt_sysJES_p[jsI]->at(jI)/recoGammaPt_;
	      recoJtPt_[1 + jsI] = aktRhi_etajes_jet_pt_sysJES_p[jsI]->at(jI);
	    }

	    for(Int_t jsI = 0; jsI < nJERSys; ++jsI){
	      recoXJ_[1 + nJESSys + jsI] = aktRhi_etajes_jet_pt_sysJER_p[jsI]->at(jI)/recoGammaPt_;
	      recoJtPt_[1 + nJESSys + jsI] = aktRhi_etajes_jet_pt_sysJER_p[jsI]->at(jI);
	    }

	    
	    if(!isSideband) unfoldXJTree_p->Fill();

	    //Now do double XJJ
	    for(unsigned int jI2 = jI+1; jI2 < aktRhi_etajes_jet_pt_p->size(); ++jI2){
	      if(aktRhi_truthpos_p->at(jI2) < 0) continue;
	      
	      Float_t ptToUse2 = aktRhi_insitu_jet_pt_p->at(jI2);
	      if(isMC){
		if(aktRhi_truthpos_p->at(jI2) >= 0) ptToUse = aktRhi_etajes_jet_pt_p->at(jI2);
	      }
	      
	      Float_t dPhi = TMath::Abs(getDPHI(aktRhi_etajes_jet_phi_p->at(jI2), photon_phi_p->at(phoPos)));
	      if(dPhi < gammaJtDPhiCut) continue;			
	      //			if(aktRhi_etajes_jet_pt_p->at(jI2) < jtPtBins[0] && aktRhi_etajes_jet_pt_p->at(jI2) >= jtPtBins[nJtPtBins]) continue;
		
	      Float_t tempTruthJtPt2 = aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI2));
	      Float_t tempTruthJtPhi2 = aktR_truth_jet_phi_p->at(aktRhi_truthpos_p->at(jI2));
	      Float_t tempTruthJtEta2 = aktR_truth_jet_eta_p->at(aktRhi_truthpos_p->at(jI2));
	      
	      tempTruthJet2.SetPtEtaPhiM(tempTruthJtPt2, tempTruthJtEta2, tempTruthJtPhi2, 0.0);
	      tempRecoJet2.SetPtEtaPhiM(ptToUse2, aktRhi_etajes_jet_eta_p->at(jI2), aktRhi_etajes_jet_phi_p->at(jI2), 0.0);
		
	      recoJt1Pt_[0] = tempRecoJet1.Pt();
	      recoJt2Pt_[0] = tempRecoJet2.Pt();
	      
	      recoJt1Phi_ = tempRecoJet1.Phi();
	      recoJt2Phi_ = tempRecoJet2.Phi();

	      truthJt1Pt_ = tempTruthJet1.Pt();
	      truthJt2Pt_ = tempTruthJet2.Pt();

	      truthJt1Phi_ = tempTruthJet1.Phi();
	      truthJt2Phi_ = tempTruthJet2.Phi();
	      
	      tempTruthJet2 += tempTruthJet1;
	      tempRecoJet2 += tempRecoJet1;
	      truthXJJ_ = tempTruthJet2.Pt()/truthGammaPt_;
	      recoXJJ_[0] = tempRecoJet2.Pt()/recoGammaPt_;

	      for(Int_t jsI = 0; jsI < nJESSys; ++jsI){
		tempRecoJet1.SetPtEtaPhiM(aktRhi_etajes_jet_pt_sysJES_p[jsI]->at(jI), aktRhi_etajes_jet_eta_p->at(jI), aktRhi_etajes_jet_phi_p->at(jI), 0.0);
		tempRecoJet2.SetPtEtaPhiM(aktRhi_etajes_jet_pt_sysJES_p[jsI]->at(jI2), aktRhi_etajes_jet_eta_p->at(jI2), aktRhi_etajes_jet_phi_p->at(jI2), 0.0);

		recoJt1Pt_[1 + jsI] = tempRecoJet1.Pt();
		recoJt2Pt_[1 + jsI] = tempRecoJet2.Pt();

		tempRecoJet2 += tempRecoJet1;
		
		recoXJJ_[1 + jsI] = tempRecoJet2.Pt()/recoGammaPt_;
	      }

	      for(Int_t jsI = 0; jsI < nJERSys; ++jsI){
		tempRecoJet1.SetPtEtaPhiM(aktRhi_etajes_jet_pt_sysJER_p[jsI]->at(jI), aktRhi_etajes_jet_eta_p->at(jI), aktRhi_etajes_jet_phi_p->at(jI), 0.0);
		tempRecoJet2.SetPtEtaPhiM(aktRhi_etajes_jet_pt_sysJER_p[jsI]->at(jI2), aktRhi_etajes_jet_eta_p->at(jI2), aktRhi_etajes_jet_phi_p->at(jI2), 0.0);

		recoJt1Pt_[1 + nJESSys + jsI] = tempRecoJet1.Pt();
		recoJt2Pt_[1 + nJESSys + jsI] = tempRecoJet2.Pt();

		tempRecoJet2 += tempRecoJet1;
		
		recoXJJ_[1 + nJESSys + jsI] = tempRecoJet2.Pt()/recoGammaPt_;
	      }
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
	      if(!isSideband) unfoldXJJTree_p->Fill();
	    }
	  }	
	}

	++multCounter;

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
      	  
	if(isMC){
	  if(aktRhi_truthpos_p->at(jI) < 0){
	    if(!isSideband) fillTH1(photonJtFakeVCentPt_p[centPos][ptPos], ptToUse, fullWeight);
	  }
	  else{
	    if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) >= assocGenMinPt){
	      if(!isSideband){
		fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][ptPos], ptToUse, fullWeight);
		fillTH1(photonGenMatchedJtPtVCentPt_p[centPos][nGammaPtBins], ptToUse, fullWeight);
		
		fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][ptPos], aktRhi_etajes_jet_eta_p->at(jI), fullWeight);
		fillTH1(photonGenMatchedJtEtaVCentPt_p[centPos][nGammaPtBins], aktRhi_etajes_jet_eta_p->at(jI), fullWeight);
		
		fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][ptPos], ptToUse/photon_pt_p->at(phoPos), fullWeight);
		fillTH1(photonGenMatchedJtXJVCentPt_p[centPos][nGammaPtBins], ptToUse/photon_pt_p->at(phoPos), fullWeight);
		
		++multCounterGenMatched;
	      }
	      
	      if(truthPhotonPt > 0){
		if(truthPhotonPt >= gammaPtBins[0] && truthPhotonPt < gammaPtBins[nGammaPtBins]){
		  Int_t genPtPos = ghostPos(nGammaPtBins, gammaPtBins, truthPhotonPt, true, doGlobalDebug);
		  Int_t recoPtPos = ghostPos(nGammaPtBins, gammaPtBins, photon_pt_p->at(phoPos), true, doGlobalDebug);		
		  
		  if(!isSideband) fillTH2(photonJtGenResVCentGenPtRecoPt_p[centPos][genPtPos][recoPtPos], ptToUse, aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)), fullWeight); 
		}
	      }
	    }
	  }
	}
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
    
    if(isMC){
      std::vector<TLorentzVector> goodTruthJets, goodTruthJetsDPhi;
      
      for(unsigned int tI = 0; tI < aktR_truth_jet_pt_p->size(); ++tI){
	if(aktR_truth_jet_pt_p->at(tI) < jtPtBinsLow) continue;
	if(aktR_truth_jet_pt_p->at(tI) >= jtPtBinsHigh) continue;
	if(aktR_truth_jet_eta_p->at(tI) <= jtEtaBinsLow) continue;
	if(aktR_truth_jet_eta_p->at(tI) >= jtEtaBinsHigh) continue;
	
	Float_t dR = getDR(aktR_truth_jet_eta_p->at(tI), aktR_truth_jet_phi_p->at(tI), photon_eta_p->at(phoPos), photon_phi_p->at(phoPos));
	if(dR < gammaExclusionDR) continue;
	  
	Float_t dPhi = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(tI), photon_phi_p->at(phoPos)));
	  
	if(!isSideband){
	  fillTH1(photonGenJtDPhiVCentPt_p[centPos][ptPos], dPhi, fullWeight);
	  fillTH1(photonGenJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, fullWeight);
	}

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
	
	TLorentzVector tempJet;
	tempJet.SetPtEtaPhiM(aktR_truth_jet_pt_p->at(tI), aktR_truth_jet_eta_p->at(tI), aktR_truth_jet_phi_p->at(tI), 0.0);
	goodTruthJets.push_back(tempJet);
	
	if(dPhi >= gammaJtDPhiCut){
	  if(!isSideband){
	    fillTH1(photonGenJtPtVCentPt_p[centPos][ptPos], aktR_truth_jet_pt_p->at(tI), fullWeight);
	    fillTH1(photonGenJtPtVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_pt_p->at(tI), fullWeight);
	    fillTH1(photonGenJtEtaVCentPt_p[centPos][ptPos], aktR_truth_jet_eta_p->at(tI), fullWeight);
	    fillTH1(photonGenJtEtaVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_eta_p->at(tI), fullWeight);
	    fillTH1(photonGenJtXJVCentPt_p[centPos][ptPos], aktR_truth_jet_pt_p->at(tI)/photon_pt_p->at(phoPos), fullWeight);
	    fillTH1(photonGenJtXJVCentPt_p[centPos][nGammaPtBins], aktR_truth_jet_pt_p->at(tI)/photon_pt_p->at(phoPos), fullWeight);
	    ++multCounterGen;
		
	    goodTruthJetsDPhi.push_back(tempJet);
	  }
	}
      }
	
      for(unsigned int gI = 0; gI < goodTruthJetsDPhi.size(); ++gI){
	TLorentzVector jet1 = goodTruthJetsDPhi[gI];
	  
	for(unsigned int gI2 = gI+1; gI2 < goodTruthJetsDPhi.size(); ++gI2){
	  TLorentzVector jet2 = goodTruthJetsDPhi[gI2];
	  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	  jet2 += jet1;
	  
	  if(!isSideband){
	    fillTH1(photonGenMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
	    fillTH1(photonGenMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
	      
	    fillTH1(photonGenMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
	    fillTH1(photonGenMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);
	  }	    	    
	}
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
    
    if(!isSideband){
      fillTH1(photonJtMultVCentPt_p[centPos][ptPos], multCounter, fullWeight);
      fillTH1(photonJtMultVCentPt_p[centPos][nGammaPtBins], multCounter, fullWeight);	
    }
  
  
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

    for(unsigned int jI = 0; jI < goodJets.size(); ++jI){
      TLorentzVector jet1 = goodJets[jI];
      
      for(unsigned int jI2 = jI+1; jI2 < goodJets.size(); ++jI2){
	TLorentzVector jet2 = goodJets[jI2];
	

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

	Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
	if(!isSideband) photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(dR, photon_pt_p->at(phoPos), fullWeight);	  
	else photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(dR, photon_pt_p->at(phoPos), fullWeight);	  

	//	if(dR < mixJetExclusionDR) continue;	    

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 
	
	Float_t interJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), jet1.Phi()));    

	//	  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	jet2 += jet1;
	
	Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));    
	if(dR >= mixJetExclusionDR){
	  if(!isSideband){
	    photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    
	    
	    leadingJtPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(multiJtDPhi, goodJets[jI].Pt(), fullWeight);
	    leadingJtPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(interJtDPhi, goodJets[jI].Pt(), fullWeight);
	    
	    fillTH2(photonPtJtDPhiJJGVCent_RAW_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    fillTH2(photonPtJtDPhiJJVCent_RAW_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	  }
	  else{
	    photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    
	    leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(multiJtDPhi, goodJets[jI].Pt(), fullWeight);
	    leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(interJtDPhi, goodJets[jI].Pt(), fullWeight);

	    
	    fillTH2(photonPtJtDPhiJJGVCent_RAWSideband_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	    fillTH2(photonPtJtDPhiJJVCent_RAWSideband_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	  }
	}
	
	if(isMC){		
	  if(goodJetsTPos[jI] >= 0){
	    Float_t truthJt1Pt = aktR_truth_jet_pt_p->at(goodJetsTPos[jI]);
	    Float_t truthJt1Eta = aktR_truth_jet_eta_p->at(goodJetsTPos[jI]);
	    Float_t truthJt1Phi = aktR_truth_jet_phi_p->at(goodJetsTPos[jI]);
	    
	    TLorentzVector truthJt1;
	    truthJt1.SetPtEtaPhiM(truthJt1Pt, truthJt1Eta, truthJt1Phi, 0.0);
	    
	    if(goodJetsTPos[jI2] >= 0){
	      Float_t truthJt2Pt = aktR_truth_jet_pt_p->at(goodJetsTPos[jI2]);
	      Float_t truthJt2Eta = aktR_truth_jet_eta_p->at(goodJetsTPos[jI2]);
	      Float_t truthJt2Phi = aktR_truth_jet_phi_p->at(goodJetsTPos[jI2]);
	      
	      TLorentzVector truthJt2;
	      truthJt2.SetPtEtaPhiM(truthJt2Pt, truthJt2Eta, truthJt2Phi, 0.0);

	      Float_t truthDR = getDR(truthJt1Eta, truthJt1Phi, truthJt2Eta, truthJt2Phi);

	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

	      if(!isSideband){
		photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(dR, photon_pt_p->at(phoPos), fullWeight);
		photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthDR, truthPhotonPt, fullWeight);
	      }

	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

	      if(dR >= mixJetExclusionDR){
		if(!isSideband){
		  if(truthPhotonPt > 0){
		    if(truthJt1Pt > assocGenMinPt && truthJt2Pt > assocGenMinPt){
		      fillTH2(photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		      fillTH2(photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);

		      photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		      photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(interJtDPhi, photon_pt_p->at(phoPos), fullWeight);

		      leadingJtPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(multiJtDPhi, truthJt1Pt, fullWeight);
		      leadingJtPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(interJtDPhi, truthJt1Pt, fullWeight);
		    
		      Float_t truthAJJ = TMath::Abs(truthJt1Pt - truthJt2Pt);
		      truthJt1 += truthJt2;
		      
		      fillTH2(photonPtJtXJJVCent_TRUTH_p[centPos][barrelOrECPos], truthJt1.Pt()/truthPhotonPt, truthPhotonPt, fullWeight);
		      fillTH2(photonPtJtAJJVCent_TRUTH_p[centPos][barrelOrECPos], truthAJJ/truthPhotonPt, truthPhotonPt, fullWeight);
		      
		      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << std::endl; 

		      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthJt1.Pt()/truthPhotonPt, truthPhotonPt, fullWeight);
		      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthAJJ/truthPhotonPt, truthPhotonPt, fullWeight);
		    }
		  }
		}
	      }


	      if(!isSideband){
		if(truthJt2Pt > assocGenMinPt){
		  fillTH2(photonPtJtDPhiJJGVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		   fillTH2(photonPtJtDPhiJJVCent_MIXBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
		else{
		  fillTH2(photonPtJtDPhiJJGVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtDPhiJJVCent_PUREBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
	      }
	    }
	    else{
	      if(!isSideband){
		if(truthJt1Pt > assocGenMinPt){
		  fillTH2(photonPtJtDPhiJJGVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtDPhiJJVCent_MIXBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
		else{
		  fillTH2(photonPtJtDPhiJJGVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtDPhiJJVCent_PUREBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
	      }
	    }
	  }
	  else{
	    if(goodJetsTPos[jI2] >= 0){
	      Float_t truthJt2Pt = aktR_truth_jet_pt_p->at(goodJetsTPos[jI2]);

	      if(!isSideband){
		if(truthJt2Pt > assocGenMinPt){
		  fillTH2(photonPtJtDPhiJJGVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtDPhiJJVCent_MIXBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
		else{
		  fillTH2(photonPtJtDPhiJJGVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtDPhiJJVCent_PUREBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		}
	      }
	    }
	    else{
	      if(!isSideband){
		fillTH2(photonPtJtDPhiJJGVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), fullWeight);
		fillTH2(photonPtJtDPhiJJVCent_PUREBKGD_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), fullWeight);
	      }
	    }
	  }
	  
	}
      }
    }
  
    if(goodJetsDPhi.size() >= 2){
      for(auto const & jet : goodJetsDPhi){
	if(!isSideband){
	  fillTH1(photonMultiJtPtVCentPt_p[centPos][ptPos], jet.Pt(), fullWeight);
	  fillTH1(photonMultiJtPtVCentPt_p[centPos][nGammaPtBins], jet.Pt(), fullWeight);
	  fillTH1(photonMultiJtXJVCentPt_p[centPos][ptPos], jet.Pt()/photon_pt_p->at(phoPos), fullWeight);
	  fillTH1(photonMultiJtXJVCentPt_p[centPos][nGammaPtBins], jet.Pt()/photon_pt_p->at(phoPos), fullWeight);
	}
      }     

      for(unsigned int jI = 0; jI < goodJetsDPhi.size(); ++jI){
	TLorentzVector jet1 = goodJetsDPhi[jI];
	
	for(unsigned int jI2 = jI+1; jI2 < goodJetsDPhi.size(); ++jI2){
	  TLorentzVector jet2 = goodJetsDPhi[jI2];

	  Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
	  if(dR < mixJetExclusionDR) continue;	    
	  
	  Float_t multiJtAJ = TMath::Abs(jet1.Pt() - jet2.Pt());
	
	  double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	  jet2 += jet1;
	  
	  Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	  if(multiJtDPhi < gammaMultiJtDPhiCut) continue;

	  if(!isSideband){
	    fillTH1(photonMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
	    fillTH1(photonMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
	    
	    fillTH1(photonMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
	    fillTH1(photonMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);

	    photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	    photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYRaw(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

	    fillTH2(photonPtJtXJJVCent_RAW_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);	      
	    fillTH2(photonPtJtAJJVCent_RAW_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);	      
	  }
	  else{
	    photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	    photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYRaw(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

	    fillTH2(photonPtJtXJJVCent_RAWSideband_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	    fillTH2(photonPtJtAJJVCent_RAWSideband_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
	  }
	  
	  if(isMC){		
	    if(goodJetsDPhiTPos[jI] >= 0){
	      Float_t truthJt1Pt = aktR_truth_jet_pt_p->at(goodJetsDPhiTPos[jI]);
	      Float_t truthJt1Eta = aktR_truth_jet_eta_p->at(goodJetsDPhiTPos[jI]);
	      Float_t truthJt1Phi = aktR_truth_jet_phi_p->at(goodJetsDPhiTPos[jI]);
	      
	      TLorentzVector truthJt1;
	      truthJt1.SetPtEtaPhiM(truthJt1Pt, truthJt1Eta, truthJt1Phi, 0.0);
	      
	      if(goodJetsDPhiTPos[jI2] >= 0){
		Float_t truthJt2Pt = aktR_truth_jet_pt_p->at(goodJetsDPhiTPos[jI2]);
		Float_t truthJt2Eta = aktR_truth_jet_eta_p->at(goodJetsDPhiTPos[jI2]);
		Float_t truthJt2Phi = aktR_truth_jet_phi_p->at(goodJetsDPhiTPos[jI2]);
		
		TLorentzVector truthJt2;
		truthJt2.SetPtEtaPhiM(truthJt2Pt, truthJt2Eta, truthJt2Phi, 0.0);
		
		if(!isSideband){
		  fillTH1(photonGenMatchedMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
		  fillTH1(photonGenMatchedMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(phoPos), fullWeight);
		  
		  fillTH1(photonGenMatchedMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, fullWeight);
		  fillTH1(photonGenMatchedMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, fullWeight);
		  
		  if(truthPhotonPt > 0){
		    if(truthJt1Pt > assocGenMinPt && truthJt2Pt > assocGenMinPt){
		      fillTH2(photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		      fillTH2(photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);

		      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruthMatchedReco(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		     
		      Float_t truthJtAJ = TMath::Abs(truthJt1.Pt() - truthJt2.Pt());
		      truthJt1 += truthJt2;			
		      fillTH2(photonPtJtXJJVCent_TRUTH_p[centPos][barrelOrECPos], truthJt1.Pt()/truthPhotonPt, truthPhotonPt, fullWeight);
		      fillTH2(photonPtJtAJJVCent_TRUTH_p[centPos][barrelOrECPos], truthJtAJ/truthPhotonPt, truthPhotonPt, fullWeight);

		      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthJt1.Pt()/truthPhotonPt, truthPhotonPt, fullWeight);
		      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYTruth(truthJtAJ/truthPhotonPt, truthPhotonPt, fullWeight);
		    }
		  }
		}
		
		if(!isSideband){
		  if(truthJt2Pt > assocGenMinPt){
		    fillTH2(photonPtJtXJJVCent_MIXBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		  else{
		    fillTH2(photonPtJtXJJVCent_PUREBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		}
	      }
	      else{
		if(!isSideband){
		  if(truthJt1Pt > assocGenMinPt){
		    fillTH2(photonPtJtXJJVCent_MIXBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		  else{
		    fillTH2(photonPtJtXJJVCent_PUREBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		}
	      }
	    }
	    else{
	      if(!isSideband){
		if(goodJetsDPhiTPos[jI2] >= 0){
		  Float_t truthJt2Pt = aktR_truth_jet_pt_p->at(goodJetsDPhiTPos[jI2]);		
		  
		  if(truthJt2Pt > assocGenMinPt){
		    fillTH2(photonPtJtXJJVCent_MIXBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_MIXBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		  else{
		    fillTH2(photonPtJtXJJVCent_PUREBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		    fillTH2(photonPtJtAJJVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  }
		}
		else{
		  fillTH2(photonPtJtXJJVCent_PUREBKGD_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		  fillTH2(photonPtJtAJJVCent_PUREBKGD_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), fullWeight);
		}
	      }
	    }	  
	  }
	}
      }
    }
    
    if(isMC){
      if(!isSideband){
	fillTH1(photonGenJtMultVCentPt_p[centPos][ptPos], multCounterGen, fullWeight);
	fillTH1(photonGenJtMultVCentPt_p[centPos][nGammaPtBins], multCounterGen, fullWeight);
	
	fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][ptPos], multCounterGenMatched, fullWeight);
	fillTH1(photonGenMatchedJtMultVCentPt_p[centPos][nGammaPtBins], multCounterGenMatched, fullWeight);
      }
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
	  		  
	unsigned long long jetPos2 = maxPos;
	while(jetPos2 == jetPos || jetPos2 == maxPos){jetPos2 = randGen_p->Uniform(0, maxPos-1);}
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	
	std::vector<TLorentzVector> jets2 = mixingMap[key][jetPos2];
	
	int multCounterMix = 0;
	goodJetsMix[0].clear();
	goodJetsDPhiMix[0].clear();
	goodJetsMix[1].clear();
	goodJetsDPhiMix[1].clear();
	
	for(unsigned int jI = 0; jI < jets.size(); ++jI){
	  if(jets[jI].Pt() < jtPtBinsLow) continue;
	  if(jets[jI].Pt() >= jtPtBinsHigh) continue;

	  Float_t jtEtaToCut = jets[jI].Eta();
	  if(jtEtaBinsDoAbs) jtEtaToCut = TMath::Abs(jtEtaToCut);
	  
	  if(jtEtaToCut <= jtEtaBinsLow) continue;
	  if(jtEtaToCut >= jtEtaBinsHigh) continue;
	  
	  Float_t dR = getDR(jets[jI].Eta(), jets[jI].Phi(), photon_eta_p->at(phoPos), photon_phi_p->at(phoPos));
	  if(dR < gammaExclusionDR) continue;
	  	  
	  goodJetsMix[0].push_back(jets[jI]);
	  
	  Float_t dPhi = TMath::Abs(getDPHI(jets[jI].Phi(), photon_phi_p->at(phoPos)));	  
	    
	  if(!isSideband){
	    fillTH1(photonMixJtDPhiVCentPt_p[centPos][ptPos], dPhi, mixWeight);
	    fillTH1(photonMixJtDPhiVCentPt_p[centPos][nGammaPtBins], dPhi, mixWeight);
	  }

	  if(!isSideband){
	    photonPtJtDPhiVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(dPhi, photon_pt_p->at(phoPos), mixWeight);	    

	    fillTH2(photonPtJtDPhiVCent_MIX_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), mixWeight);
	  }
	  else{
	    photonPtJtDPhiVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(dPhi, photon_pt_p->at(phoPos), mixWeight);	    

	    fillTH2(photonPtJtDPhiVCent_MIXSideband_p[centPos][barrelOrECPos], dPhi, photon_pt_p->at(phoPos), mixWeight);
	  }
	  
	  if(dPhi >= gammaJtDPhiCut){
	    goodJetsDPhiMix[0].push_back(jets[jI]);
	      
	    if(!isSideband){
	      fillTH1(photonMixJtPtVCentPt_p[centPos][ptPos], jets[jI].Pt(), mixWeight);
	      fillTH1(photonMixJtPtVCentPt_p[centPos][nGammaPtBins], jets[jI].Pt(), mixWeight);
	      fillTH1(photonMixJtEtaVCentPt_p[centPos][ptPos], jets[jI].Eta(), mixWeight);
	      fillTH1(photonMixJtEtaVCentPt_p[centPos][nGammaPtBins], jets[jI].Eta(), mixWeight);
	      fillTH1(photonMixJtXJVCentPt_p[centPos][ptPos], jets[jI].Pt()/photon_pt_p->at(phoPos), mixWeight);
	      fillTH1(photonMixJtXJVCentPt_p[centPos][nGammaPtBins], jets[jI].Pt()/photon_pt_p->at(phoPos), mixWeight);
	      
	      photonPtJtPtVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(jets[jI].Pt(), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtXJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(jets[jI].Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);

	      fillTH2(photonPtJtPtVCent_MIX_p[centPos][barrelOrECPos], jets[jI].Pt(), photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtXJVCent_MIX_p[centPos][barrelOrECPos], jets[jI].Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	    }
	    else{
	      photonPtJtPtVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(jets[jI].Pt(), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtXJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(jets[jI].Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      

	      fillTH2(photonPtJtPtVCent_MIXSideband_p[centPos][barrelOrECPos], jets[jI].Pt(), photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtXJVCent_MIXSideband_p[centPos][barrelOrECPos], jets[jI].Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	    }
	    
	    ++multCounterMix;
	  }	    
	}
		
	for(unsigned int jI = 0; jI < jets2.size(); ++jI){
	  if(jets2[jI].Pt() < jtPtBinsLow) continue;
	  if(jets2[jI].Pt() >= jtPtBinsHigh) continue;

	  Float_t jtEtaToCut = jets2[jI].Eta();
	  if(jtEtaBinsDoAbs) jtEtaToCut = TMath::Abs(jtEtaToCut);

	  if(jtEtaToCut <= jtEtaBinsLow) continue;
	  if(jtEtaToCut >= jtEtaBinsHigh) continue;
	  
	  Float_t dR = getDR(jets2[jI].Eta(), jets2[jI].Phi(), photon_eta_p->at(phoPos), photon_phi_p->at(phoPos));
	  if(dR < gammaExclusionDR) continue;
	  
	  goodJetsMix[1].push_back(jets2[jI]);
	    
	  Float_t dPhi = TMath::Abs(getDPHI(jets2[jI].Phi(), photon_phi_p->at(phoPos)));
	  if(dPhi >= gammaJtDPhiCut){
	    goodJetsDPhiMix[1].push_back(jets2[jI]);
	  }	    
	}
	
      
	//	if(goodJetsDPhiMix[0].size() >= 2){
	for(auto const & jet : goodJetsDPhiMix[0]){
	  if(!isSideband){
	    fillTH1(photonMixMultiJtPtVCentPt_p[centPos][ptPos], jet.Pt(), mixWeight);
	    fillTH1(photonMixMultiJtPtVCentPt_p[centPos][nGammaPtBins], jet.Pt(), mixWeight);
	    fillTH1(photonMixMultiJtXJVCentPt_p[centPos][ptPos], jet.Pt()/photon_pt_p->at(phoPos), mixWeight);
	    fillTH1(photonMixMultiJtXJVCentPt_p[centPos][nGammaPtBins], jet.Pt()/photon_pt_p->at(phoPos), mixWeight);
	  }
	}
     	
	//First do pure background
	for(unsigned int bI = 0; bI < goodJetsMix[0].size(); ++bI){
	  TLorentzVector jet1 = goodJetsMix[0][bI];
	    
	  for(unsigned int bI2 = bI+1; bI2 < goodJetsMix[0].size(); ++bI2){
	    TLorentzVector jet2 = goodJetsMix[0][bI2];

	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());	    
	    if(!isSideband) photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(dR, photon_pt_p->at(phoPos), mixWeight);
	    else photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(dR, photon_pt_p->at(phoPos), mixWeight);


	    if(dR < mixJetExclusionDR) continue;	    	    

	    //	    double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	    Float_t interJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), jet1.Phi()));

	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	      
	    if(!isSideband){
	      fillTH2(photonPtJtDPhiJJGVCent_MIX_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtDPhiJJGVCent_MIXPURE_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtDPhiJJVCent_MIX_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtDPhiJJVCent_MIXPURE_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      leadingJtPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, jet1.Pt(), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, jet1.Pt(), mixWeight);
	    }
	    else{
	      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, jet1.Pt(), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, jet1.Pt(), mixWeight);


	      fillTH2(photonPtJtDPhiJJGVCent_MIXSideband_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtDPhiJJVCent_MIXSideband_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);	
	    }
	  }	
	}

	for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
	  TLorentzVector jet1 = goodJetsDPhiMix[0][bI];
	    
	  for(unsigned int bI2 = bI+1; bI2 < goodJetsDPhiMix[0].size(); ++bI2){
	    TLorentzVector jet2 = goodJetsDPhiMix[0][bI2];

	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());

	    if(dR < mixJetExclusionDR) continue;	    

	    Float_t multiJtAJ = TMath::Abs(jet1.Pt() - jet2.Pt());
	    
	    double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
	      
	    if(!isSideband){
	      fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
	      fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
		
	      fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
	      fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);

	      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      
	      fillTH2(photonPtJtXJJVCent_MIX_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtXJJVCent_MIXPURE_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	

	      fillTH2(photonPtJtAJJVCent_MIX_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtAJJVCent_MIXPURE_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	

	    }	  
	    else{
	      photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);


	      fillTH2(photonPtJtXJJVCent_MIXSideband_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtAJJVCent_MIXSideband_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	    }
	  }
	}
        
	//Now do mixed background (signal jet + associated w/ fake jets)
	for(unsigned int sI = 0; sI < goodJets.size(); ++sI){
	  TLorentzVector jet1 = goodJets[sI];
	  
	  for(unsigned int bI = 0; bI < goodJetsMix[0].size(); ++bI){
	    TLorentzVector jet2 = goodJetsMix[0][bI];
	    
	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
	    if(!isSideband) photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(dR, photon_pt_p->at(phoPos), mixWeight);
	    else photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(dR, photon_pt_p->at(phoPos), mixWeight);

	    if(dR < mixJetExclusionDR) continue;
	      
	    //	    double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	    Float_t interJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), jet1.Phi()));
	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	    
	    if(!isSideband){
	      photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      leadingJtPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, jet1.Pt(), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, jet1.Pt(), mixWeight);

	      fillTH2(photonPtJtDPhiJJGVCent_MIX_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight); 	   
	      fillTH2(photonPtJtDPhiJJGVCent_MIXMIX_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight); 	   

	      fillTH2(photonPtJtDPhiJJVCent_MIX_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight); 	   
	      fillTH2(photonPtJtDPhiJJVCent_MIXMIX_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight); 	   
	    }
	    else{
	      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtDPhi, jet1.Pt(), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(interJtDPhi, jet1.Pt(), mixWeight);

	      fillTH2(photonPtJtDPhiJJGVCent_MIXSideband_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);		      
	      fillTH2(photonPtJtDPhiJJVCent_MIXSideband_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);		      
	    }
	  }
	}

	for(unsigned int sI = 0; sI < goodJetsDPhi.size(); ++sI){
	  TLorentzVector jet1 = goodJetsDPhi[sI];
	  
	  for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
	    TLorentzVector jet2 = goodJetsDPhiMix[0][bI];
	    
	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
	    if(dR < mixJetExclusionDR) continue;
	      
	    double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
	    Float_t multiJtAJ = TMath::Abs(jet1.Pt() - jet2.Pt());

	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
	    
	    if(!isSideband){
	      fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
	      fillTH1(photonMixMultiJtXJJVCentPt_p[centPos][nGammaPtBins],jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
	      
	      fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
	      fillTH1(photonMixMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);

	      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMix(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);

	      fillTH2(photonPtJtXJJVCent_MIX_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtXJJVCent_MIXMIX_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	

	      fillTH2(photonPtJtAJJVCent_MIX_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	      fillTH2(photonPtJtAJJVCent_MIXMIX_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);	
	    }	  
	    else{
	      photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMix(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);

	      fillTH2(photonPtJtXJJVCent_MIXSideband_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);		      
	      fillTH2(photonPtJtAJJVCent_MIXSideband_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);		      
	    }
	  }
	}	  	

	//Now calculate the mixed event correction, for cases where your gamma + single jet embed accidentally picked a fake jet
	for(unsigned int bI = 0; bI < goodJetsMix[0].size(); ++bI){
	  TLorentzVector jet1 = goodJetsMix[0][bI];
	  
	  for(unsigned int bI2 = 0; bI2 < goodJetsMix[1].size(); ++bI2){
	    TLorentzVector jet2 = goodJetsMix[1][bI2];
	    
	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());

	    if(!isSideband) photonPtJtDRJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(dR, photon_pt_p->at(phoPos), mixWeight);
	    else photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(dR, photon_pt_p->at(phoPos), mixWeight);

	    if(dR < mixJetExclusionDR) continue;

	    Float_t interJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), jet1.Phi()));
	    
	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	    
	    if(!isSideband){	     
	      photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      leadingJtPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtDPhi, jet1.Pt(), mixWeight);

	      fillTH2(photonPtJtDPhiJJGVCent_MIXCorrection_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);		  
	      fillTH2(photonPtJtDPhiJJGVCent_MIXCorrection_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);		  

	      photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(interJtDPhi, jet1.Pt(), mixWeight);

	      fillTH2(photonPtJtDPhiJJVCent_MIXCorrection_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);		  
	      fillTH2(photonPtJtDPhiJJVCent_MIXCorrection_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);		  
	    }
	    else{
	      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(interJtDPhi, photon_pt_p->at(phoPos), mixWeight);

	      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtDPhi, jet1.Pt(), mixWeight);
	      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(interJtDPhi, jet1.Pt(), mixWeight);


	      fillTH2(photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[centPos][barrelOrECPos], multiJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[centPos][barrelOrECPos], interJtDPhi, photon_pt_p->at(phoPos), mixWeight);
	    }
	  }
	}

	for(unsigned int bI = 0; bI < goodJetsDPhiMix[0].size(); ++bI){
	  TLorentzVector jet1 = goodJetsDPhiMix[0][bI];
	  
	  for(unsigned int bI2 = 0; bI2 < goodJetsDPhiMix[1].size(); ++bI2){
	    TLorentzVector jet2 = goodJetsDPhiMix[1][bI2];
	    
	    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
	    if(dR < mixJetExclusionDR) continue;
	    
	    Float_t multiJtAJ = TMath::Abs(jet1.Pt() - jet2.Pt());
	    double dPhiJJ = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));

	    jet2 += jet1;

	    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(phoPos)));
	    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
	    
	    if(!isSideband){
 	      fillTH1(photonMixCorrectionMultiJtXJJVCentPt_p[centPos][ptPos], jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
	      fillTH1(photonMixCorrectionMultiJtXJJVCentPt_p[centPos][nGammaPtBins], jet2.Pt()/photon_pt_p->at(phoPos), mixWeight);
	      
	      fillTH1(photonMixCorrectionMultiJtDPhiJJVCentPt_p[centPos][ptPos], dPhiJJ, mixWeight);
	      fillTH1(photonMixCorrectionMultiJtDPhiJJVCentPt_p[centPos][nGammaPtBins], dPhiJJ, mixWeight);
	      
	      photonPtJtXJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);

	      fillTH2(photonPtJtXJJVCent_MIXCorrection_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);		  
	      fillTH2(photonPtJtAJJVCent_MIXCorrection_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);		  
	    }
	    else{
	      photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelOrECPos]->FillXYMixCorrection(multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);


	      fillTH2(photonPtJtXJJVCent_MIXSidebandCorrection_p[centPos][barrelOrECPos], jet2.Pt()/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	      fillTH2(photonPtJtAJJVCent_MIXSidebandCorrection_p[centPos][barrelOrECPos], multiJtAJ/photon_pt_p->at(phoPos), photon_pt_p->at(phoPos), mixWeight);
	    }
	  }
	}
	  
	if(!isSideband){
	  fillTH1(photonMixJtMultVCentPt_p[centPos][ptPos], multCounterMix, mixWeight);
	  fillTH1(photonMixJtMultVCentPt_p[centPos][nGammaPtBins], multCounterMix, mixWeight);	  
	  
	  fillTH1(photonSubJtMultModVCentPt_p[centPos][ptPos], TMath::Max(0, multCounter - multCounterMix), mixWeight);
	  fillTH1(photonSubJtMultModVCentPt_p[centPos][nGammaPtBins], TMath::Max(0, multCounter - multCounterMix), mixWeight);
	}
	
	++nCurrentMixEvents;
      }
    }
    else if(isPP){
      
      if(!isSideband){
	fillTH1(photonSubJtMultModVCentPt_p[centPos][ptPos], TMath::Max(0, multCounter), fullWeight);
	fillTH1(photonSubJtMultModVCentPt_p[centPos][nGammaPtBins], TMath::Max(0, multCounter), fullWeight);
      }    
    }
  
    if(!isSideband) fillTH2(photonEtaPt_p[centPos], etaValMain, photon_pt_p->at(phoPos), fullWeight);
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

  doGlobalDebug = false;
  if(doGlobalDebug) std::cout << "DOGLOBALDEBUG MANUALLY TURNED ON AT LINE: " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  std::cout << inFile_p << std::endl;
  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(isMC){
    for(Int_t i = 0; i < nJtPtBins; ++i){
      truthPtOfMatches_Unweighted_p[i]->Write("", TObject::kOverwrite);
      truthPtOfMatches_Weighted_p[i]->Write("", TObject::kOverwrite);
      
      truthPtDeltaROfMatches_Unweighted_p[i]->Write("", TObject::kOverwrite);
      truthPtDeltaROfMatches_Weighted_p[i]->Write("", TObject::kOverwrite);
      
      delete truthPtOfMatches_Unweighted_p[i];
      delete truthPtOfMatches_Weighted_p[i];
      
      delete truthPtDeltaROfMatches_Unweighted_p[i];
      delete truthPtDeltaROfMatches_Weighted_p[i];
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(isMC){
    unfoldPhotonPtTree_p->Write("", TObject::kOverwrite);
    delete unfoldPhotonPtTree_p;

    unfoldXJTree_p->Write("", TObject::kOverwrite);
    delete unfoldXJTree_p;

    unfoldXJJTree_p->Write("", TObject::kOverwrite);
    delete unfoldXJJTree_p;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* purityFile_p = new TFile(inPurityFileName.c_str(), "READ");

  //Purity correction definitions via Yeonju Go:

  // purity correction on jet observables 
  // O_sig = ( O_raw - (1-p)*O_bkg ) / p , where p = purity, O_x = photon normalized yield after event mixing subtraction
  // O_sig*nPhoton_sigPho*p = ( O_raw*nPhoton_sigPho - (1-p)*O_bkg*(nPhoton_sigPho) ) , here, nPhoton_sigPho is actually the number of RAW photon!

  //We will use the second definition as we are not yet normalizing to per photon

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
  

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
    //Combine the separated barrel and endcap
    for(Int_t bI = 0; bI < nBarrelAndEC; ++bI){
      if(barrelAndECComboPos == bI) continue;

      photonPtVCent_RAW_p[cI][barrelAndECComboPos]->Add(photonPtVCent_RAW_p[cI][bI]);
      photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->Add(photonPtVCent_RAWSideband_p[cI][bI]);


      photonPtJtPtVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_MixMachine_p[cI][bI]);
      photonPtJtXJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_MixMachine_p[cI][bI]);
      photonPtJtDPhiVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_MixMachine_p[cI][bI]);
      photonPtJtXJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_MixMachine_p[cI][bI]);
      photonPtJtAJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_MixMachine_p[cI][bI]);
      photonPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_MixMachine_p[cI][bI]);
      photonPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_MixMachine_p[cI][bI]);
      photonPtJtDRJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDRJJVCent_MixMachine_p[cI][bI]);

      leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][bI]);
      leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][bI]);
      
      photonPtJtPtVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtXJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtXJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtAJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][bI]);
      photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][bI]);

      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][bI]);
      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][bI]);

      
      if(isMC){
	photonPtVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtVCent_TRUTH_p[cI][bI]);
	photonPtJtPtVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_TRUTH_p[cI][bI]);
	photonPtJtXJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_TRUTH_p[cI][bI]);
	photonPtJtXJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_TRUTH_p[cI][bI]);
	photonPtJtAJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiJJGVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_TRUTH_p[cI][bI]);

	photonPtVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][bI]);
	photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][bI]);
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::cout << "WE ARE HERE " << __LINE__ << std::endl;
    //Compute the subtracted via the mix-machine
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      photonPtJtPtVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtXJVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtXJJVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtAJJVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiJJVCent_MixMachine_p[cI][eI]->ComputeSub();       
      photonPtJtDRJJVCent_MixMachine_p[cI][eI]->ComputeSub();       

      leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->ComputeSub();       
      leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI]->ComputeSub();       

      photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       

      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->ComputeSub();       
    }
    


    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::cout << "LINE: " << __LINE__ << std::endl;
  
    if(!isMC){//Purity correction; Data only
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      if(doUnifiedPurity){	
	//The unified purity cannot do barrel and endcap corrections separately 
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_Eta0p00to2p37";
	if(isStrSame(centBinsStr[cI], "Cent30to50")) purityFitStr = "fit_purity_Cent30to80_Eta0p00to2p37";
	else if(isStrSame(centBinsStr[cI], "Cent50to80")) purityFitStr = "fit_purity_Cent30to80_Eta0p00to2p37";
	
	TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());
	TH2F* tempRawHist_p = photonPtJtPtVCent_MixMachine_p[cI][barrelAndECComboPos]->GetTH2FPtr("RAW");

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//    	for(Int_t bIY = 0; bIY < tempRawHist_p->GetYaxis()->GetNbins(); ++bIY){
	  for(Int_t bIY = 0; bIY < tempRawHist_p->GetYaxis()->GetNbins(); ++bIY){
            Float_t lowVal = tempRawHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
            Float_t hiVal = tempRawHist_p->GetYaxis()->GetBinLowEdge(bIY+2);

            Float_t totalNum = 0.0;
            Float_t totalDenom = 0.0;//Total denom is also the total number of photons

            for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetNbins(); ++bIX){
              Float_t center = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
              if(center >= hiVal) continue;

              Float_t lowVal2 = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinLowEdge(bIX+1);
              Float_t hiVal2 = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinLowEdge(bIX+2);

              if(lowVal2 < lowVal - 0.0001 || hiVal2 >= hiVal + 0.0001){
		std::cout << "BIG ERROR BINS ARE NOT SELFCONSISTENT SUBSETS DURING PURITY CALC. return \
1" << std::endl;
                return 1;
              }

              Float_t content = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
              totalNum += content*center;
              totalDenom += content;
            }

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    //Get number of sideband photons
            Float_t totalSideband = 0.0;
            for(Int_t bIX = 0; bIX < photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetXaxis()->GetNbins(); ++bIX){
              Float_t center = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
              if(center >= hiVal) continue;
	      
	      Float_t lowVal2 = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinLowEdge(bIX+1);
              Float_t hiVal2 = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinLowEdge(bIX+2);
	      
	      if(lowVal2 < lowVal - 0.0001 || hiVal2 >= hiVal + 0.0001){
		std::cout << "BIG ERROR BINS ARE NOT SELFCONSISTENT SUBSETS DURING PURITY CALC. return 1" << std::endl;
                return 1;
              }
	      
	      Float_t content = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
	      totalSideband += content;
	    }
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


	    Float_t meanVal = totalNum/totalDenom;
            Float_t purity = purityFit_p->Eval(meanVal);

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	    std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtXJVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtDPhiVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtXJJVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtAJJVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos], photonPtJtDRJJVCent_MixMachine_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos]};
	    std::vector<mixMachine*> inMixMachines_Sideband_p = {photonPtJtPtVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtXJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtXJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtAJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]};
	    std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtXJVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtDPhiVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtXJJVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtAJJVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos], photonPtJtDRJJVCent_PURCORR_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos], leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos]};	  
	    
	    for(unsigned int hI = 0; hI < inMixMachines_p.size(); ++hI){
	      TH2F* subHist_p = inMixMachines_p[hI]->GetTH2FPtr("SUB");
	      TH2F* subHist_Sideband_p = inMixMachines_Sideband_p[hI]->GetTH2FPtr("SUB");	      
	      
	      for(Int_t bIX = 0; bIX < subHist_p->GetXaxis()->GetNbins(); ++bIX){
		Float_t binContent = subHist_p->GetBinContent(bIX+1, bIY+1);
		Float_t binError = subHist_p->GetBinError(bIX+1, bIY+1);
		
		Float_t binSidebandContent = subHist_Sideband_p->GetBinContent(bIX+1, bIY+1);
                Float_t binSidebandError = subHist_Sideband_p->GetBinError(bIX+1, bIY+1);

		//		std::cout << "CONTENT: " << binContent << std::endl;

		binSidebandContent = (1.0 - purity)*binSidebandContent*totalDenom/totalSideband;
                binSidebandError = (1.0 - purity)*binSidebandError*totalDenom/totalSideband;
		
                binContent -= binSidebandContent;
		
                binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
                outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
	      }
	    }	 	   
	  }       
      
	
	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t center = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetBinCenter(bIX+1);
            Float_t purity = purityFit_p->Eval(center);

            Float_t binContent = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
            Float_t binError = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinError(bIX+1);

	    Float_t binSidebandContent = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
	    binSidebandContent = (1.0 - purity)*binSidebandContent*binContent/binSidebandContent;
	    
	    binContent -= binSidebandContent;

	    photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinContent(bIX+1, binContent);
            photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinError(bIX+1, binError);
	  }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    //	}
      }
      else{
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  if(isStrSame(barrelAndECStr[eI], "BarrelAndEC")) continue;

	  std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_" + barrelAndECFitStr[eI];
	  TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());
	  
	  TH2F* tempRawHist_p = photonPtJtPtVCent_MixMachine_p[cI][eI]->GetTH2FPtr("RAW");

	  for(Int_t bIY = 0; bIY < tempRawHist_p->GetYaxis()->GetNbins(); ++bIY){
	    Float_t lowVal = tempRawHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
            Float_t hiVal = tempRawHist_p->GetYaxis()->GetBinLowEdge(bIY+2);

	    Float_t totalNum = 0.0;
	    Float_t totalDenom = 0.0;//Total denom is also the total number of photons
	    for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
	      if(center >= hiVal) continue;

	      Float_t lowVal2 = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinLowEdge(bIX+1);
	      Float_t hiVal2 = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinLowEdge(bIX+2);
	      
	      if(lowVal2 < lowVal - 0.0001 || hiVal2 >= hiVal + 0.0001){
		std::cout << "BIG ERROR BINS ARE NOT SELFCONSISTENT SUBSETS DURING PURITY CALC. return 1" << std::endl;
		return 1;
	      }

	      Float_t content = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	      totalNum += content*center;
              totalDenom += content;
	    }

	    //Get number of sideband photons
            Float_t totalSideband = 0.0;
            for(Int_t bIX = 0; bIX < photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
              Float_t center = photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
              if(center >= hiVal) continue;

	      Float_t lowVal2 = photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetBinLowEdge(bIX+1);
              Float_t hiVal2 = photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetBinLowEdge(bIX+2);
	      
	      if(lowVal2 < lowVal - 0.0001 || hiVal2 >= hiVal + 0.0001){
		std::cout << "BIG ERROR BINS ARE NOT SELFCONSISTENT SUBSETS DURING PURITY CALC. return 1" << std::endl;
                return 1;
              }

	      Float_t content = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	      totalSideband += content;
	    }
	    Float_t meanVal = totalNum/totalDenom;
            Float_t purity = purityFit_p->Eval(meanVal);

	    std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][eI], photonPtJtXJVCent_MixMachine_p[cI][eI], photonPtJtDPhiVCent_MixMachine_p[cI][eI], photonPtJtXJJVCent_MixMachine_p[cI][eI], photonPtJtAJJVCent_MixMachine_p[cI][eI], photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI], photonPtJtDPhiJJVCent_MixMachine_p[cI][eI], photonPtJtDRJJVCent_MixMachine_p[cI][eI], leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI], leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI]};
	    std::vector<mixMachine*> inMixMachines_Sideband_p = {photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI], photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI], photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI], photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI], photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI], photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI], photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI], photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI], leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI], leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]};
	    std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI], photonPtJtDRJJVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][eI]};	  

	    for(unsigned int hI = 0; hI < inMixMachines_p.size(); ++hI){
	      TH2F* subHist_p = inMixMachines_p[hI]->GetTH2FPtr("SUB");
	      TH2F* subHist_Sideband_p = inMixMachines_Sideband_p[hI]->GetTH2FPtr("SUB");	      
	      
	      std::cout << "PURITY CORRECTION HISTS: " << std::endl;
	      std::cout << " " << subHist_p->GetName() << std::endl;
	      std::cout << " " << subHist_Sideband_p->GetName() << std::endl;

	      for(Int_t bIX = 0; bIX < subHist_p->GetXaxis()->GetNbins(); ++bIX){
		Float_t binContent = subHist_p->GetBinContent(bIX+1, bIY+1);
		Float_t binError = subHist_p->GetBinError(bIX+1, bIY+1);

		Float_t binSidebandContent = subHist_Sideband_p->GetBinContent(bIX+1, bIY+1);
                Float_t binSidebandError = subHist_Sideband_p->GetBinError(bIX+1, bIY+1);

		binSidebandContent = (1.0 - purity)*binSidebandContent*totalDenom/totalSideband;
                binSidebandError = (1.0 - purity)*binSidebandError*totalDenom/totalSideband;

                binContent -= binSidebandContent;

                binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);

		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
                outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
	      }
	    }
	  }	  
	
	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	    if(isStrSame(barrelAndECStr[eI], "BarrelAndEC")) continue;

	    Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
            Float_t purity = purityFit_p->Eval(center);

            Float_t binContent = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
            Float_t binError = photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1);

	    Float_t binSidebandContent = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	    binSidebandContent = (1.0 - purity)*binSidebandContent*binContent/binSidebandContent;
	    
	    binContent -= binSidebandContent;

	    photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, binContent);
            photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, binError);
	  }
	}

      
	//Gotta manually sum to get the Barrel and EC
	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  if(eI == barrelAndECComboPos) continue;
	  
	  photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtVCent_PURCORR_p[cI][eI]);
      
	  photonPtJtPtVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_PURCORR_p[cI][eI]);
	  photonPtJtXJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_PURCORR_p[cI][eI]);
	  photonPtJtDPhiVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_PURCORR_p[cI][eI]);
	  photonPtJtXJJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_PURCORR_p[cI][eI]);
	  photonPtJtAJJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_PURCORR_p[cI][eI]);
	  photonPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI]);
	  photonPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]);
	  photonPtJtDRJJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(photonPtJtDRJJVCent_PURCORR_p[cI][eI]);

	  leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][eI]);
	  leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos]->Add(leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][eI]);
	}
      }
    }
    else{
      //If its MC we just pass thru; Purity correction goes on data only
      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){	
	std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][eI], photonPtJtXJVCent_MixMachine_p[cI][eI], photonPtJtDPhiVCent_MixMachine_p[cI][eI], photonPtJtXJJVCent_MixMachine_p[cI][eI], photonPtJtAJJVCent_MixMachine_p[cI][eI], photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI], photonPtJtDPhiJJVCent_MixMachine_p[cI][eI], photonPtJtDRJJVCent_MixMachine_p[cI][eI], leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI], leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI]};
	std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI], photonPtJtDRJJVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJGVCent_PURCORR_p[cI][eI], leadingJtPtJtDPhiJJVCent_PURCORR_p[cI][eI]};
	
	for(unsigned int mI = 0; mI < inMixMachines_p.size(); ++mI){
	  TH2F* subHist_p = inMixMachines_p[mI]->GetTH2FPtr("SUB");
	  
	  for(unsigned int bIX = 0; bIX < outHists_p[mI]->GetXaxis()->GetNbins(); ++bIX){
	    outHists_p[mI]->SetBinContent(bIX+1, subHist_p->GetBinContent(bIX+1));
	    outHists_p[mI]->SetBinError(bIX+1, subHist_p->GetBinError(bIX+1));
	  }
	}
      }
    }
  
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
    if(doMix){
      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){         
	//Correct DPhiJJ and XJJ Mix Event
	photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_MIX_p[cI][eI], photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_MIXMIX_p[cI][eI], photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_MIXSideband_p[cI][eI], photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[cI][eI], 1.0, -1.0);

	photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI]->Add(photonPtJtDPhiJJVCent_MIX_p[cI][eI], photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJVCent_MIXMIXCorrected_p[cI][eI]->Add(photonPtJtDPhiJJVCent_MIXMIX_p[cI][eI], photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI]->Add(photonPtJtDPhiJJVCent_MIXSideband_p[cI][eI], photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[cI][eI], 1.0, -1.0);

	photonPtJtXJJVCent_MIXCorrected_p[cI][eI]->Add(photonPtJtXJJVCent_MIX_p[cI][eI], photonPtJtXJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtXJJVCent_MIXMIXCorrected_p[cI][eI]->Add(photonPtJtXJJVCent_MIXMIX_p[cI][eI], photonPtJtXJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI]->Add(photonPtJtXJJVCent_MIXSideband_p[cI][eI], photonPtJtXJJVCent_MIXSidebandCorrection_p[cI][eI], 1.0, -1.0);

	photonPtJtAJJVCent_MIXCorrected_p[cI][eI]->Add(photonPtJtAJJVCent_MIX_p[cI][eI], photonPtJtAJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtAJJVCent_MIXMIXCorrected_p[cI][eI]->Add(photonPtJtAJJVCent_MIXMIX_p[cI][eI], photonPtJtAJJVCent_MIXCorrection_p[cI][eI], 1.0, -1.0);
	photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI]->Add(photonPtJtAJJVCent_MIXSideband_p[cI][eI], photonPtJtAJJVCent_MIXSidebandCorrection_p[cI][eI], 1.0, -1.0);

	//Subtract off the mixed event
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	
	photonPtJtPtVCent_SUB_p[cI][eI]->Add(photonPtJtPtVCent_RAW_p[cI][eI], photonPtJtPtVCent_MIX_p[cI][eI], 1.0, -1.0);
	photonPtJtXJVCent_SUB_p[cI][eI]->Add(photonPtJtXJVCent_RAW_p[cI][eI], photonPtJtXJVCent_MIX_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiVCent_RAW_p[cI][eI], photonPtJtDPhiVCent_MIX_p[cI][eI], 1.0, -1.0);
	photonPtJtXJJVCent_SUB_p[cI][eI]->Add(photonPtJtXJJVCent_RAW_p[cI][eI], photonPtJtXJJVCent_MIXCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtAJJVCent_SUB_p[cI][eI]->Add(photonPtJtAJJVCent_RAW_p[cI][eI], photonPtJtAJJVCent_MIXCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJGVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_RAW_p[cI][eI], photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiJJVCent_RAW_p[cI][eI], photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI], 1.0, -1.0);
	
	
	photonPtJtPtVCent_SUBSideband_p[cI][eI]->Add(photonPtJtPtVCent_RAWSideband_p[cI][eI], photonPtJtPtVCent_MIXSideband_p[cI][eI], 1.0, -1.0);
	photonPtJtXJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtXJVCent_RAWSideband_p[cI][eI], photonPtJtXJVCent_MIXSideband_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiVCent_RAWSideband_p[cI][eI], photonPtJtDPhiVCent_MIXSideband_p[cI][eI], 1.0, -1.0);
	photonPtJtXJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtXJJVCent_RAWSideband_p[cI][eI], photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtAJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtAJJVCent_RAWSideband_p[cI][eI], photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI], photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI], 1.0, -1.0);
	photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI], photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI], 1.0, -1.0);
       	           	
	if(!isMC && eI < nBarrelAndEC-1 && false){
	  //Purity correction
	  std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_" + barrelAndECFitStr[eI];
	  TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());
	  
	  for(Int_t bIY = 0; bIY < photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetNbins(); ++bIY){
	    Float_t lowVal = photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetBinLowEdge(bIY+1);
	    Float_t hiVal = photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetBinLowEdge(bIY+2);
	    
	    Float_t totalNum = 0.0;
	    Float_t totalDenom = 0.0;//Total denom is also the total number of photons
	    for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
	      if(center >= hiVal) continue;
	      
	      Float_t content = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	      
	      totalNum += content*center;
	      totalDenom += content;
	    }
	    
	    //Get number of sideband photons
	    Float_t totalSideband = 0.0;
	    for(Int_t bIX = 0; bIX < photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
	      if(center >= hiVal) continue;
	      
	      Float_t content = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	      
	      totalSideband += content;
	    }	 
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    Float_t meanVal = totalNum/totalDenom;
	    Float_t purity = purityFit_p->Eval(meanVal);

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    std::vector<TH2F*> inSignalHists_p = {photonPtJtPtVCent_SUB_p[cI][eI], photonPtJtXJVCent_SUB_p[cI][eI], photonPtJtDPhiVCent_SUB_p[cI][eI], photonPtJtXJJVCent_SUB_p[cI][eI], photonPtJtAJJVCent_SUB_p[cI][eI], photonPtJtDPhiJJGVCent_SUB_p[cI][eI], photonPtJtDPhiJJVCent_SUB_p[cI][eI]};
	    std::vector<TH2F*> inSidebandHists_p = {photonPtJtPtVCent_SUBSideband_p[cI][eI], photonPtJtXJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiVCent_SUBSideband_p[cI][eI], photonPtJtXJJVCent_SUBSideband_p[cI][eI], photonPtJtAJJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]};
	    std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]};

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	   
	    //We will correct the values on a loop
	    for(unsigned int hI = 0; hI < inSignalHists_p.size(); ++hI){
	      for(Int_t bIX = 0; bIX < inSignalHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
		Float_t binContent = inSignalHists_p[hI]->GetBinContent(bIX+1, bIY+1);	 
		Float_t binError = inSignalHists_p[hI]->GetBinError(bIX+1, bIY+1);	
		
		Float_t binSidebandContent = inSidebandHists_p[hI]->GetBinContent(bIX+1, bIY+1);
		Float_t binSidebandError = inSidebandHists_p[hI]->GetBinError(bIX+1, bIY+1);
		
		binSidebandContent = (1.0 - purity)*binSidebandContent*totalDenom/totalSideband;
		binSidebandError = (1.0 - purity)*binSidebandError*totalDenom/totalSideband;
		
		binContent -= binSidebandContent;

		binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
		//		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
		//		outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
	      }
	    }	  	    
	  }

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	    Float_t purity = purityFit_p->Eval(center);	   

	    Float_t binContent = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	    Float_t binError = photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1);

	    if(binContent <= TMath::Power(10, -50)) continue;

 	    Float_t binSidebandContent = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	    //	    Float_t binSidebandError = photonPtVCent_RAWSideband_p[cI][eI]->GetBinError(bIX+1);

	    binSidebandContent = (1.0 - purity)*binSidebandContent*binContent/binSidebandContent;
	    //	    binSidebandError = (1.0 - purity)*binSidebandError*binContent/binSidebandContent;
		
	    binContent -= binSidebandContent;	   
	    //	    binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, binContent);
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, binError);
	  }
	}
	else{
	  std::vector<TH2F*> inSignalHists_p = {photonPtJtPtVCent_SUB_p[cI][eI], photonPtJtXJVCent_SUB_p[cI][eI], photonPtJtDPhiVCent_SUB_p[cI][eI], photonPtJtXJJVCent_SUB_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_SUB_p[cI][eI], photonPtJtDPhiJJVCent_SUB_p[cI][eI]};
	  std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]};

	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t content = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	    Float_t error = photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1);
	    
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, content);
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, error);
	  }

	  for(unsigned int hI = 0; hI < inSignalHists_p.size(); ++hI){
	    for(Int_t bIX = 0; bIX < inSignalHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
	      for(Int_t bIY = 0; bIY < inSignalHists_p[hI]->GetYaxis()->GetNbins(); ++bIY){
		Float_t content = inSignalHists_p[hI]->GetBinContent(bIX+1, bIY+1);
		Float_t error = inSignalHists_p[hI]->GetBinError(bIX+1, bIY+1);

		//		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, content);
		//		outHists_p[hI]->SetBinError(bIX+1, bIY+1, error);
	      }
	    }
	  }	 
	}
      }
    }
    else if(isPP && false){

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	/*
	photonPtJtPtVCent_SUB_p[cI][eI]->Add(photonPtJtPtVCent_RAW_p[cI][eI]);
	photonPtJtXJVCent_SUB_p[cI][eI]->Add(photonPtJtXJVCent_RAW_p[cI][eI]);
	photonPtJtDPhiVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiVCent_RAW_p[cI][eI]);
	photonPtJtXJJVCent_SUB_p[cI][eI]->Add(photonPtJtXJJVCent_RAW_p[cI][eI]);
	photonPtJtAJJVCent_SUB_p[cI][eI]->Add(photonPtJtAJJVCent_RAW_p[cI][eI]);
	photonPtJtDPhiJJGVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_RAW_p[cI][eI]);
	photonPtJtDPhiJJVCent_SUB_p[cI][eI]->Add(photonPtJtDPhiJJVCent_RAW_p[cI][eI]);
		
	photonPtJtPtVCent_SUBSideband_p[cI][eI]->Add(photonPtJtPtVCent_RAWSideband_p[cI][eI]);
	photonPtJtXJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtXJVCent_RAWSideband_p[cI][eI]);
	photonPtJtDPhiVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiVCent_RAWSideband_p[cI][eI]);
	photonPtJtXJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtXJJVCent_RAWSideband_p[cI][eI]);
	photonPtJtAJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtAJJVCent_RAWSideband_p[cI][eI]);
	photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI]);
	photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]->Add(photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI]);
	*/

	if(!isMC){
	  std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_" + barrelAndECFitStr[eI];
	  TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());
	  
	  for(Int_t bIY = 0; bIY < photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetNbins(); ++bIY){
	    Float_t lowVal = photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetBinLowEdge(bIY+1);
	    Float_t hiVal = photonPtJtPtVCent_RAW_p[cI][eI]->GetYaxis()->GetBinLowEdge(bIY+2);
	    
	    Float_t totalNum = 0.0;
	    Float_t totalDenom = 0.0;
	    for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
	      if(center >= hiVal) continue;
	      
	      Float_t content = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	      
	      totalNum += content*center;
	      totalDenom += content;
	    }
	    
	    //Get number of sidband photons
	    Float_t totalSideband = 0.0;
	    for(Int_t bIX = 0; bIX < photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = photonPtVCent_RAWSideband_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	      
	      if(center < lowVal) continue;
	      if(center >= hiVal) continue;
	      
	      Float_t content = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	      
	      totalSideband += content;
	    }	 
	    
	    
	    Float_t meanVal = totalNum/totalDenom;
	    Float_t purity = purityFit_p->Eval(meanVal);
	    
	    std::vector<TH2F*> inSignalHists_p = {photonPtJtPtVCent_SUB_p[cI][eI], photonPtJtXJVCent_SUB_p[cI][eI], photonPtJtDPhiVCent_SUB_p[cI][eI], photonPtJtXJJVCent_SUB_p[cI][eI], photonPtJtAJJVCent_SUB_p[cI][eI], photonPtJtDPhiJJGVCent_SUB_p[cI][eI], photonPtJtDPhiJJVCent_SUB_p[cI][eI]};
	    std::vector<TH2F*> inSidebandHists_p = {photonPtJtPtVCent_SUBSideband_p[cI][eI], photonPtJtXJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiVCent_SUBSideband_p[cI][eI], photonPtJtXJJVCent_SUBSideband_p[cI][eI], photonPtJtAJJVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI], photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]};
	    std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]};
	      
	    //We will correct the values on a loop
	    for(unsigned int hI = 0; hI < inSignalHists_p.size(); ++hI){
	      for(Int_t bIX = 0; bIX < inSignalHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
		Float_t binContent = inSignalHists_p[hI]->GetBinContent(bIX+1, bIY+1);	 
		Float_t binError = inSignalHists_p[hI]->GetBinError(bIX+1, bIY+1);	

		Float_t binSidebandContent = inSidebandHists_p[hI]->GetBinContent(bIX+1, bIY+1);
		Float_t binSidebandError = inSidebandHists_p[hI]->GetBinError(bIX+1, bIY+1);

		//		std::cout << " CALC VAL (1.0 - " << purity << ")*" << binSidebandContent << "*" << totalDenom << "/" << totalSideband << std::endl;
		//		std::cout << " CALC ERR (1.0 - " << purity << ")*" << binSidebandError << "*" << totalDenom << "/" << totalSideband << std::endl;
		
		binSidebandContent = (1.0 - purity)*binSidebandContent*totalDenom/totalSideband;
		binSidebandError = (1.0 - purity)*binSidebandError*totalDenom/totalSideband;
		

		binContent -= binSidebandContent;
		binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
		//		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
		//		outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
	      }
	    }	  
	  }
       
	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t center = photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetBinCenter(bIX+1);
	    Float_t purity = purityFit_p->Eval(center);

	    Float_t binContent = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	    Float_t binError = photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1);

	    if(binContent <= TMath::Power(10, -50)) continue;

	    Float_t binSidebandContent = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(bIX+1);
	    //	    Float_t binSidebandError = photonPtVCent_RAWSideband_p[cI][eI]->GetBinError(bIX+1);

	    binSidebandContent = (1.0 - purity)*binSidebandContent*binContent/binSidebandContent;
	    //	    binSidebandError = (1.0 - purity)*binSidebandError*binContent/binSidebandContent;
		
	    binContent -= binSidebandContent;	   
	    //	    binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, binContent);
	    //	    photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, binError);
	  }	  
	}	
	else{
	  std::vector<TH2F*> inSignalHists_p = {photonPtJtPtVCent_SUB_p[cI][eI], photonPtJtXJVCent_SUB_p[cI][eI], photonPtJtDPhiVCent_SUB_p[cI][eI], photonPtJtXJJVCent_SUB_p[cI][eI], photonPtJtAJJVCent_SUB_p[cI][eI], photonPtJtDPhiJJGVCent_SUB_p[cI][eI], photonPtJtDPhiJJVCent_SUB_p[cI][eI]};
	  std::vector<TH2F*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI], photonPtJtXJVCent_PURCORR_p[cI][eI], photonPtJtDPhiVCent_PURCORR_p[cI][eI], photonPtJtXJJVCent_PURCORR_p[cI][eI], photonPtJtAJJVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]};

	  for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t content = photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1);
	    Float_t error = photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1);
	    
	    photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, content);
	    photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, error);
	  }


	  for(unsigned int hI = 0; hI < inSignalHists_p.size(); ++hI){
	    for(Int_t bIX = 0; bIX < inSignalHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
	      for(Int_t bIY = 0; bIY < inSignalHists_p[hI]->GetYaxis()->GetNbins(); ++bIY){
		Float_t content = inSignalHists_p[hI]->GetBinContent(bIX+1, bIY+1);
		Float_t error = inSignalHists_p[hI]->GetBinError(bIX+1, bIY+1);

		//		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, content);
		//		outHists_p[hI]->SetBinError(bIX+1, bIY+1, error);
	      }
	    }
	  }
	}

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  purityFile_p->Close();
  delete purityFile_p;

  outFile_p->cd();
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
  
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      photonPtVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);

      photonPtJtPtVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtPtVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtXJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtXJVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtDPhiVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtXJJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtXJJVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtAJJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtAJJVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtDPhiJJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiJJVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtDRJJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDRJJVCent_MixMachine_p[cI][eI]->Clean();

      leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      leadingJtPtJtDPhiJJGVCent_MixMachine_p[cI][eI]->Clean();

      leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI]->WriteToDirectory(centDir_p);
      leadingJtPtJtDPhiJJVCent_MixMachine_p[cI][eI]->Clean();

      photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      leadingJtPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI]->Clean();

      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->WriteToDirectory(centDir_p);
      leadingJtPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI]->Clean();


      photonPtJtPtVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtXJVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtAJJVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiJJGVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiJJVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);

      photonPtVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);      
      photonPtJtPtVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtXJVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtXJJVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtAJJVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);

      photonPtVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
      if(isMC){
	photonPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
      }

      if(doMix){
	photonPtJtPtVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtXJJVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXPURE_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXMIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXMIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtAJJVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXPURE_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXMIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXMIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtDPhiJJGVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXPURE_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXMIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtDPhiJJVCent_MIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXPURE_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXMIX_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXMIXCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtPtVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtPtVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtXJJVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXSidebandCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtAJJVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXSidebandCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtDPhiJJGVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtDPhiJJVCent_MIXSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtPtVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
      	
	photonPtJtPtVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);

	if(isMC){
	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);

	  photonPtJtPtVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);

	  
	  photonPtJtPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	}
      }
      else if(isPP){
	photonPtJtPtVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_SUB_p[cI][eI]->Write("", TObject::kOverwrite);

	photonPtJtPtVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI]->Write("", TObject::kOverwrite);
	
	photonPtJtPtVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtXJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtAJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtJtDPhiJJVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);

	if(isMC){
	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);

	  photonPtJtPtVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);

	  photonPtJtPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	}
      }
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
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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

    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      delete photonPtVCent_RAW_p[cI][eI];

      delete photonPtJtPtVCent_RAW_p[cI][eI];
      delete photonPtJtDPhiVCent_RAW_p[cI][eI];
      delete photonPtJtXJVCent_RAW_p[cI][eI];
      delete photonPtJtXJJVCent_RAW_p[cI][eI];
      delete photonPtJtAJJVCent_RAW_p[cI][eI];
      delete photonPtJtDPhiJJGVCent_RAW_p[cI][eI];
      delete photonPtJtDPhiJJVCent_RAW_p[cI][eI];

      delete photonPtVCent_RAWSideband_p[cI][eI];
      
      delete photonPtJtPtVCent_RAWSideband_p[cI][eI];
      delete photonPtJtXJVCent_RAWSideband_p[cI][eI];
      delete photonPtJtDPhiVCent_RAWSideband_p[cI][eI];
      delete photonPtJtXJJVCent_RAWSideband_p[cI][eI];
      delete photonPtJtAJJVCent_RAWSideband_p[cI][eI];
      delete photonPtJtDPhiJJGVCent_RAWSideband_p[cI][eI];
      delete photonPtJtDPhiJJVCent_RAWSideband_p[cI][eI];

      delete photonPtVCent_PURCORR_p[cI][eI];
      if(isMC){
	delete photonPtVCent_TRUTH_p[cI][eI];
	delete photonPtVCent_TRUTHMATCHEDRECO_p[cI][eI];
      }
      
      if(doMix){
	delete photonPtJtPtVCent_MIX_p[cI][eI];
	delete photonPtJtXJVCent_MIX_p[cI][eI];
	delete photonPtJtDPhiVCent_MIX_p[cI][eI];

	delete photonPtJtXJJVCent_MIX_p[cI][eI];
	delete photonPtJtXJJVCent_MIXCorrection_p[cI][eI];
	delete photonPtJtXJJVCent_MIXCorrected_p[cI][eI];
	delete photonPtJtXJJVCent_MIXPURE_p[cI][eI];
	delete photonPtJtXJJVCent_MIXMIX_p[cI][eI];
	delete photonPtJtXJJVCent_MIXMIXCorrected_p[cI][eI];

	delete photonPtJtAJJVCent_MIX_p[cI][eI];
	delete photonPtJtAJJVCent_MIXCorrection_p[cI][eI];
	delete photonPtJtAJJVCent_MIXCorrected_p[cI][eI];
	delete photonPtJtAJJVCent_MIXPURE_p[cI][eI];
	delete photonPtJtAJJVCent_MIXMIX_p[cI][eI];
	delete photonPtJtAJJVCent_MIXMIXCorrected_p[cI][eI];

	delete photonPtJtDPhiJJGVCent_MIX_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXCorrection_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXCorrected_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXPURE_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXMIX_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXMIXCorrected_p[cI][eI];

	delete photonPtJtDPhiJJVCent_MIX_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXCorrection_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXCorrected_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXPURE_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXMIX_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXMIXCorrected_p[cI][eI];
	
	delete photonPtJtPtVCent_SUB_p[cI][eI];
	delete photonPtJtXJVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiVCent_SUB_p[cI][eI];
	delete photonPtJtXJJVCent_SUB_p[cI][eI];
	delete photonPtJtAJJVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiJJVCent_SUB_p[cI][eI];
	
	delete photonPtJtPtVCent_MIXSideband_p[cI][eI];
	delete photonPtJtXJVCent_MIXSideband_p[cI][eI];
	delete photonPtJtDPhiVCent_MIXSideband_p[cI][eI];
	delete photonPtJtXJJVCent_MIXSideband_p[cI][eI];
	delete photonPtJtXJJVCent_MIXSidebandCorrection_p[cI][eI];
	delete photonPtJtXJJVCent_MIXSidebandCorrected_p[cI][eI];
	delete photonPtJtAJJVCent_MIXSideband_p[cI][eI];
	delete photonPtJtAJJVCent_MIXSidebandCorrection_p[cI][eI];
	delete photonPtJtAJJVCent_MIXSidebandCorrected_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXSideband_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXSidebandCorrection_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_MIXSidebandCorrected_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXSideband_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXSidebandCorrection_p[cI][eI];
	delete photonPtJtDPhiJJVCent_MIXSidebandCorrected_p[cI][eI];
	
	delete photonPtJtPtVCent_SUBSideband_p[cI][eI];
	delete photonPtJtXJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiVCent_SUBSideband_p[cI][eI];
	delete photonPtJtXJJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtAJJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI];

	delete photonPtJtPtVCent_PURCORR_p[cI][eI];
	delete photonPtJtXJVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiVCent_PURCORR_p[cI][eI];
	delete photonPtJtXJJVCent_PURCORR_p[cI][eI];
	delete photonPtJtAJJVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiJJVCent_PURCORR_p[cI][eI];

	if(isMC){
	  delete photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI];

	  delete photonPtJtPtVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtXJVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtDPhiVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtXJJVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtXJJVCent_MIXBKGD_p[cI][eI];
	  delete photonPtJtAJJVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtAJJVCent_MIXBKGD_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI];

	  delete photonPtJtPtVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTH_p[cI][eI];
	}
      }
      else if(isPP){
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	delete photonPtJtPtVCent_SUB_p[cI][eI];
	delete photonPtJtXJVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiVCent_SUB_p[cI][eI];
	delete photonPtJtXJJVCent_SUB_p[cI][eI];
	delete photonPtJtAJJVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_SUB_p[cI][eI];
	delete photonPtJtDPhiJJVCent_SUB_p[cI][eI];

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	delete photonPtJtPtVCent_SUBSideband_p[cI][eI];
	delete photonPtJtXJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiVCent_SUBSideband_p[cI][eI];
	delete photonPtJtXJJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtAJJVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_SUBSideband_p[cI][eI];
	delete photonPtJtDPhiJJVCent_SUBSideband_p[cI][eI];
	
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	delete photonPtJtPtVCent_PURCORR_p[cI][eI];
	delete photonPtJtXJVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiVCent_PURCORR_p[cI][eI];
	delete photonPtJtXJJVCent_PURCORR_p[cI][eI];
	delete photonPtJtAJJVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI];
	delete photonPtJtDPhiJJVCent_PURCORR_p[cI][eI];

	if(isMC){
	  delete photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI];

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  delete photonPtJtPtVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTH_p[cI][eI];
	}
      }
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

  outFile_p->cd();
  config_p->SetValue("RECOJTPTMIN", prettyString(recoJtPtMin, 1, false).c_str());
  config_p->SetValue("NGAMMAFINEPTBINS", nGammaFinePtBins);
  config_p->SetValue("GAMMAFINEPTBINS", gammaFinePtBinsStr.c_str());


  config_p->SetValue("NJETSYSANDNOM", nJetSysAndNom);
  std::string jetSysAndNom = "";
  for(Int_t sI = 0; sI < nJetSysAndNom; ++sI){
    jetSysAndNom = jetSysAndNom + jetSysAndNomStr[sI] + ",";
  }
  if(jetSysAndNom.find(",") != std::string::npos) jetSysAndNom.replace(jetSysAndNom.size()-1, 1, "");
  config_p->SetValue("JETSYSANDNOM", jetSysAndNom.c_str());
  
  config_p->Write("config", TObject::kOverwrite);

  TEnv labelEnv;
  for(auto const & lab : binsToLabelStr){
    labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
  }
  labelEnv.Write("label", TObject::kOverwrite);

  vzPassing_p->Write("", TObject::kOverwrite);  
  centPassing_p->Write("", TObject::kOverwrite);
  truPhoPtPassing_p->Write("", TObject::kOverwrite);
  leadingPhoPtPassing_p->Write("", TObject::kOverwrite);
  leadingPhoEtaPassing_p->Write("", TObject::kOverwrite);
  leadingPhoIsoPassing_p->Write("", TObject::kOverwrite);
  leadingPhoTightPassing_p->Write("", TObject::kOverwrite);

  jtEtaPassing_p->Write("", TObject::kOverwrite);
  jtPtPassing_p->Write("", TObject::kOverwrite);
  jtDPhiPassing_p->Write("", TObject::kOverwrite);
  jtGammaDRPassing_p->Write("", TObject::kOverwrite);

  delete vzPassing_p;  
  delete centPassing_p;
  delete truPhoPtPassing_p;
  delete leadingPhoPtPassing_p;
  delete leadingPhoEtaPassing_p;
  delete leadingPhoIsoPassing_p;
  delete leadingPhoTightPassing_p;

  delete jtEtaPassing_p;
  delete jtPtPassing_p;
  delete jtDPhiPassing_p;
  delete jtGammaDRPassing_p;

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
