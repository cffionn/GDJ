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
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TTree.h"

//Local
#include "include/binFlattener.h"
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/fillUtils.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/mixMachine.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/returnFileList.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"

bool sortJetVectAndTruthPos(std::vector<TLorentzVector>* jetVect, std::vector<std::vector<int>* > jetTruthPosVect)
{  
  if(jetTruthPosVect.size() == 0){
    std::cout << "SORTJETVECTANDTRUTHPOS: NO PAIRED VECTORS GIVEN. RETURN FALSE" << std::endl;
    return false;
  }
  
  for(unsigned int jI = 0; jI < jetTruthPosVect.size(); ++jI){
    if(jetVect->size() != jetTruthPosVect[jI]->size()){
      std::cout << "SORTJETVECTANDTRUTHPOS: SIZE MISMATCH. RETURN FALSE" << std::endl;
      return false;
    }
  }

  unsigned int maxPos = jetVect->size();
  unsigned int jI = 0;
  while(jI < maxPos){
    bool passes = true;
    for(unsigned int jI2 = jI+1; jI2 < maxPos; ++jI2){
      if((*jetVect)[jI].Pt() < (*jetVect)[jI2].Pt()){
	passes = false;

	TLorentzVector tempJet = (*jetVect)[jI];
	std::vector<int> tempTruthPos;
	for(unsigned int pI = 0; pI < jetTruthPosVect.size(); ++pI){
	  tempTruthPos.push_back((*(jetTruthPosVect[pI]))[jI]);
	}

	(*jetVect)[jI] = (*jetVect)[jI2];
	for(unsigned int pI = 0; pI < jetTruthPosVect.size(); ++pI){
	  (*(jetTruthPosVect[pI]))[jI] = (*(jetTruthPosVect[pI]))[jI2];
	}
	
	(*jetVect)[jI2] = tempJet;
	for(unsigned int pI = 0; pI < jetTruthPosVect.size(); ++pI){
	  (*(jetTruthPosVect[pI]))[jI2] = tempTruthPos[pI];
	}
      }
    }

    if(passes) ++jI;
  }

  return true;
}


void createMixMachineTEnv(TEnv* inEnv_p, bool is2DUnfoldValue, bool isMCValue, int nBinsXValue, std::string binsXValue, std::string titleXValue, int nBinsYValue, std::string binsYValue, std::string titleYValue)
{
  inEnv_p->SetValue("IS2DUNFOLD", is2DUnfoldValue);
  inEnv_p->SetValue("ISMC", isMCValue);

  inEnv_p->SetValue("NBINSX", nBinsXValue);
  inEnv_p->SetValue("BINSX", binsXValue.c_str());
  inEnv_p->SetValue("TITLEX", titleXValue.c_str());

  inEnv_p->SetValue("NBINSY", nBinsYValue);
  inEnv_p->SetValue("BINSY", binsYValue.c_str());
  inEnv_p->SetValue("TITLEY", titleYValue.c_str());

  return;
}


bool readBins(TEnv* config_p, std::string varStr, Int_t nMaxBins, Int_t* nBins, Float_t* low, Float_t* high, Bool_t* doLog, Bool_t* doCustom, Bool_t* doAbs, Double_t bins[], std::string* binsStrForConfig, Float_t* lowReco = nullptr, Float_t* highReco = nullptr) //Last two args are optional
{
  std::string inConfigFileName = config_p->GetValue("CONFIGFILENAME", "");
  (*nBins) = config_p->GetValue(("N" + varStr + "BINS").c_str(), nMaxBins+1);
  if(!goodBinning(inConfigFileName, nMaxBins, *nBins, "N" + varStr + "BINS")) return false;

  (*low) = config_p->GetValue((varStr + "BINSLOW").c_str(), 80.0);
  (*high) = config_p->GetValue((varStr + "BINSHIGH").c_str(), 280.0);
  if(doLog != nullptr) (*doLog) = (bool)config_p->GetValue((varStr + "BINSDOLOG").c_str(), 0);
  if(doCustom != nullptr) (*doCustom) = (bool)config_p->GetValue((varStr + "BINSDOCUSTOM").c_str(), 0);
  if(doAbs != nullptr) (*doAbs) = (bool)config_p->GetValue((varStr + "BINSDOABS").c_str(), 0);

  if(lowReco != nullptr && highReco != nullptr){
    const float defaultVal = -1000.0;
    const float defaultValTest = -999.0;

    (*lowReco) = config_p->GetValue((varStr + "BINSLOWRECO").c_str(), defaultVal);
    (*highReco) = config_p->GetValue((varStr + "BINSHIGHRECO").c_str(), defaultVal);

    if(*lowReco < defaultValTest) (*lowReco) = (*low);
    if(*highReco < defaultValTest) (*highReco) = (*high);
  }

  Float_t tempLow = (*low);
  //Handle absolute bins
  if(doAbs != nullptr){
    if(doAbs){
      if(tempLow < 0.0) tempLow = 0.0;      
    }
  }

  bool binsFilled = false;
  if(doLog != nullptr){
    if(*doLog){
      getLogBins(tempLow, *high, *nBins, bins);
      binsFilled = true;
    }
  }
  
  if(!binsFilled && doCustom != nullptr){
    if(*doCustom){
      std::string binsStr = config_p->GetValue((varStr + "BINSCUSTOM").c_str(), "");
      if(binsStr.size() == 0){
	std::cout << "CUSTOM " << varStr << " BINS REQUESTED, BUT \'" << varStr << "BINSDOCUSTOM\' not defined in config \'" << inConfigFileName << "\'. return false" << std::endl;
	return false;
      }
      std::vector<float> binsCustom = strToVectF(binsStr);
      Int_t nBinsTemp = ((Int_t)binsCustom.size()) - 1;
      if(nBinsTemp != *nBins){
	std::cout << "REQUESTED CUSTON BINS \'" << binsStr << "\' HAS SIZE \'" << nBinsTemp << "\' does not match given n" << varStr << "Bins \'" << *nBins << "\'. return false" << std::endl;
	return false;
      }
      
      for(unsigned int gI = 0; gI < binsCustom.size(); ++gI){
	bins[gI] = binsCustom[gI];
      }   
      binsFilled = true;
    }
  }
  if(!binsFilled) getLinBins(tempLow, *high, *nBins, bins);    

  (*binsStrForConfig) = "";
  for(int bI = 0; bI < (*nBins)+1; ++bI){
    std::string newStr = std::to_string(bins[bI]);
    if(newStr.find(".") != std::string::npos){
      while(newStr.rfind("0") == newStr.size()-1){
	if(newStr.rfind(".") == newStr.size()-2) break;
	newStr = newStr.substr(0, newStr.size()-1);
      }
    }	
    
    (*binsStrForConfig) = (*binsStrForConfig) + newStr + ",";    
  }
  if(binsStrForConfig->size() > 0) binsStrForConfig->replace(binsStrForConfig->size()-1, 1, "");  
  return binsFilled;
}

//Taken from Run2 dijet asymmetry, ATL-COM-PHY-2020-138
Double_t getRTrkJESSysPt(float cent, Double_t jtPt)
{
  Double_t factor = 1 - 0.000153*(80.0 - cent);
  return jtPt*factor;
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
  config_p->SetValue("CONFIGFILENAME", inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INFILENAME",
                                              "OUTFILENAME",
                                              "CENTFILENAME",
					      "MIXFILENAME",
					      "DOUNIFIEDPURITY",
					      "PURITYFILENAME",
					      "JETR",
					      "GAMMAEXCLUSIONDR",
					      "MIXJETEXCLUSIONDR",
					      "DOMIX",
					      "MIXCAP",
					      "NMIXEVENTS",
					      "DOSTRICTMIX",
                                              "ISPP",
					      "ISMC",
					      "KEEPRESPONSETREE",
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
					      "GAMMAPTBINSLOWRECO",
					      "GAMMAPTBINSHIGHRECO",
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
					      "JTPTBINSLOWRECO",
					      "JTPTBINSLOWRECOSYST",
					      "JTPTBINSHIGHRECO",
					      "NSUBJTPTBINS",
					      "SUBJTPTBINSLOW",
					      "SUBJTPTBINSHIGH",
					      "SUBJTPTBINSDOLOG",
					      "SUBJTPTBINSDOCUSTOM",
					      "SUBJTPTBINSCUSTOM",
					      "SUBJTPTBINSLOWRECO",
					      "SUBJTPTBINSLOWRECOSYST",
					      "SUBJTPTBINSHIGHRECO",
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
					      "XJBINSLOWRECO",
					      "XJBINSHIGHRECO",
					      "NXJJBINS",
					      "XJJBINSLOW",
					      "XJJBINSHIGH",
					      "XJJBINSDOCUSTOM",
					      "XJJBINSCUSTOM",
					      "XJJBINSLOWRECO",
					      "XJJBINSHIGHRECO",
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
  const std::string nStartEvtStr = config_p->GetValue("NEVTSTART", "");

  ULong64_t nMaxEvt = 0;
  if(nMaxEvtStr.size() != 0){nMaxEvt = std::stol(nMaxEvtStr);} 

  ULong64_t nStartEvt = 0;
  if(nStartEvtStr.size() != 0){nStartEvt = std::stol(nStartEvtStr);} 

  const int jetR = config_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);
  //  const double assocGenMinPt = config_p->GetValue("ASSOCGENMINPT", 15.0);
  const double gammaExclusionDR = config_p->GetValue("GAMMAEXCLUSIONDR", 0.5);
  const double mixJetExclusionDR = config_p->GetValue("MIXJETEXCLUSIONDR", 0.5);

  const bool doMix = config_p->GetValue("DOMIX", 0);

  if(doMix){
    if(!checkEnvForParams(config_p, mixParams)) return 1;
  } 
  
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
  const bool keepResponseTree = config_p->GetValue("KEEPRESPONSETREE", 0);
  
  //multiple files are possibly input; check if input is dir or file, and create vector of filenames
  std::vector<std::string> inROOTFileNames;
  if(check.checkFileExt(inROOTFileName, "root")){
    inROOTFileNames.push_back(inROOTFileName);
  }
  else if(check.checkDir(inROOTFileName)){
    inROOTFileNames = returnFileList(inROOTFileName, ".root");
    
    //check we have at least 1 valid file
    if(inROOTFileNames.size() == 0){
      std::cout << "No valid files found for INFILENAME \'" << inROOTFileName << "\'. return 1" << std::endl;
      return 1;
    }
  }
  else{
    std::cout << "Given INFILENAME: \'" << inROOTFileName << "\' is neither file nor directory. return 1" << std::endl;
    return 1;
  }

  //  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
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
  std::string systemStr = "PP";
  if(!isPP) systemStr = "PbPb";
  
  std::vector<int> centBins;
  std::vector<std::string> centBinsStr = {systemStr};
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

  //Block of variable declarations
  Int_t nGammaPtBins, nJtPtBins, nSubJtPtBins, nXJBins, nXJJBins, nDPhiBins, nAJBins, nDRBins;
  Float_t gammaPtBinsLow, gammaPtBinsHigh, jtPtBinsLow, jtPtBinsHigh, subJtPtBinsLow, subJtPtBinsHigh, xjBinsLow, xjBinsHigh, xjjBinsLow, xjjBinsHigh, dPhiBinsLow, dPhiBinsHigh, ajBinsLow, ajBinsHigh, drBinsLow, drBinsHigh;
  Bool_t gammaPtBinsDoLog, gammaPtBinsDoCustom, jtPtBinsDoLog, jtPtBinsDoCustom, subJtPtBinsDoLog, subJtPtBinsDoCustom, xjBinsDoLog, xjBinsDoCustom, xjjBinsDoLog, xjjBinsDoCustom, dPhiBinsDoLog, dPhiBinsDoCustom, ajBinsDoLog, ajBinsDoCustom, drBinsDoLog, drBinsDoCustom;
  Double_t gammaPtBins[nMaxPtBins+1], jtPtBins[nMaxPtBins+1], subJtPtBins[nMaxPtBins+1], xjBins[nMaxPtBins+1], xjjBins[nMaxPtBins+1], dPhiBins[nMaxPtBins+1], ajBins[nMaxPtBins+1], drBins[nMaxPtBins+1];
  std::string gammaPtBinsStrForConfig, jtPtBinsStrForConfig, subJtPtBinsStrForConfig, xjBinsStrForConfig, xjjBinsStrForConfig, dPhiBinsStrForConfig, ajBinsStrForConfig, drBinsStrForConfig;
  Float_t gammaPtBinsLowReco, gammaPtBinsHighReco, jtPtBinsLowReco, jtPtBinsHighReco, subJtPtBinsLowReco, subJtPtBinsHighReco, xjBinsLowReco, xjBinsHighReco, xjjBinsLowReco, xjjBinsHighReco, ajBinsLowReco, ajBinsHighReco;

  readBins(config_p, "GAMMAPT", nMaxPtBins, &nGammaPtBins, &gammaPtBinsLow, &gammaPtBinsHigh, &gammaPtBinsDoLog, &gammaPtBinsDoCustom, nullptr, gammaPtBins, &gammaPtBinsStrForConfig, &gammaPtBinsLowReco, &gammaPtBinsHighReco);

  //Grabbing jtpt bins is more involved because of the pt systematic
  readBins(config_p, "JTPT", nMaxPtBins, &nJtPtBins, &jtPtBinsLow, &jtPtBinsHigh, &jtPtBinsDoLog, &jtPtBinsDoCustom, nullptr, jtPtBins, &jtPtBinsStrForConfig, &jtPtBinsLowReco, &jtPtBinsHighReco);
  const Float_t jtPtBinsLowRecoSyst = config_p->GetValue("JTPTBINSLOWRECOSYST", -1000.0);
  if(jtPtBinsLowRecoSyst < 0 || jtPtBinsLowRecoSyst < jtPtBinsLowReco){
    std::cout << "ERROR: Given value for jtPtBinsLowRecoSyst (systematic variation of pt minimum for unfolding sensitivity) \'" << jtPtBinsLowRecoSyst << "\' is not valid. Must be > 0 and > jtPtBinsLowReco, \'" << jtPtBinsLowReco << "\'. return 1" << std::endl;
    return 1;
  }

  //Repeat for the subjtpt bins
  readBins(config_p, "SUBJTPT", nMaxPtBins, &nSubJtPtBins, &subJtPtBinsLow, &subJtPtBinsHigh, &subJtPtBinsDoLog, &subJtPtBinsDoCustom, nullptr, subJtPtBins, &subJtPtBinsStrForConfig, &subJtPtBinsLowReco, &subJtPtBinsHighReco);
  const Float_t subJtPtBinsLowRecoSyst = config_p->GetValue("SUBJTPTBINSLOWRECOSYST", -1000.0);
  if(subJtPtBinsLowRecoSyst < 0 || subJtPtBinsLowRecoSyst < subJtPtBinsLowReco){
    std::cout << "ERROR: Given value for subJtPtBinsLowRecoSyst (systematic variation of pt minimum for unfolding sensitivity) \'" << subJtPtBinsLowRecoSyst << "\' is not valid. Must be > 0 and > subJtPtBinsLowReco, \'" << subJtPtBinsLowReco << "\'. return 1" << std::endl;
    return 1;
  }

  //Before continuing, institute a bunch of checks to make sure everything maps correctly between jtptbins and subjtptbins
  const double deltaPt = 0.001; //should be sufficient for integer level pt bins
  if(!isFloatSame(jtPtBinsLow, subJtPtBinsLow, deltaPt)){
    std::cout << "ERROR gdjNTupleToHist: jtPtBinsLow \'" << jtPtBinsLow << "\' doesn't match subJtPtBinsLow \'" << subJtPtBinsLow << "\' to precision \'" << deltaPt << "\'. return 1" << std::endl;
    return 1;
  }
  if(!isFloatSame(jtPtBinsLowReco, subJtPtBinsLowReco, deltaPt)){
    std::cout << "ERROR gdjNTupleToHist: jtPtBinsLowReco \'" << jtPtBinsLowReco << "\' doesn't match subJtPtBinsLowReco \'" << subJtPtBinsLowReco << "\' to precision \'" << deltaPt << "\'. return 1" << std::endl;
    return 1;
  }
  if(!isFloatSame(jtPtBinsLowRecoSyst, subJtPtBinsLowRecoSyst, deltaPt)){
    std::cout << "ERROR gdjNTupleToHist: jtPtBinsLowRecoSyst \'" << jtPtBinsLowRecoSyst << "\' doesn't match subJtPtBinsLowRecoSyst \'" << subJtPtBinsLowRecoSyst << "\' to precision \'" << deltaPt << "\'. return 1" << std::endl;
    return 1;
  }
  
  readBins(config_p, "XJ", nMaxPtBins, &nXJBins, &xjBinsLow, &xjBinsHigh, &xjBinsDoLog, &xjBinsDoCustom, nullptr, xjBins, &xjBinsStrForConfig, &xjBinsLowReco, &xjBinsHighReco);
  readBins(config_p, "XJJ", nMaxPtBins, &nXJJBins, &xjjBinsLow, &xjjBinsHigh, &xjjBinsDoLog, &xjjBinsDoCustom, nullptr, xjjBins, &xjjBinsStrForConfig, &xjjBinsLowReco, &xjjBinsHighReco);  
   readBins(config_p, "DPHI", nMaxPtBins, &nDPhiBins, &dPhiBinsLow, &dPhiBinsHigh, &dPhiBinsDoLog, &dPhiBinsDoCustom, nullptr, dPhiBins, &dPhiBinsStrForConfig);
  readBins(config_p, "AJ", nMaxPtBins, &nAJBins, &ajBinsLow, &ajBinsHigh, &ajBinsDoLog, &ajBinsDoCustom, nullptr, ajBins, &ajBinsStrForConfig, &ajBinsLowReco, &ajBinsHighReco);
  readBins(config_p, "DR", nMaxPtBins, &nDRBins, &drBinsLow, &drBinsHigh, &drBinsDoLog, &drBinsDoCustom, nullptr, drBins, &drBinsStrForConfig);
 
 
  std::vector<std::string> genGammaPtBinsStr, recoGammaPtBinsStr, gammaPtBinsStr;
  for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
    genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(pI));
    recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(pI));
    gammaPtBinsStr.push_back("GammaPt" + std::to_string(pI));

    //gammaPtBinsStrForConfig = gammaPtBinsStrForConfig + std::to_string(gammaPtBins[pI]) + ", ";

    binsToLabelStr[genGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[recoGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    binsToLabelStr[gammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
  }
  //  gammaPtBinsStrForConfig = gammaPtBinsStrForConfig + std::to_string(gammaPtBins[nGammaPtBins]);

  genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(nGammaPtBins));
  recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(nGammaPtBins));
  gammaPtBinsStr.push_back("GammaPt" + std::to_string(nGammaPtBins));

  binsToLabelStr[genGammaPtBinsStr[genGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  binsToLabelStr[recoGammaPtBinsStr[recoGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
   binsToLabelStr[gammaPtBinsStr[gammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
  //Also set the fine grain photon bins right now
  
 
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


  std::string jtPtBinsGlobalStr = "GlobalJtPt0";
  std::string jtPtBinsGlobalLabel = prettyString(jtPtBinsLow,1,false) + " < p_{T,jet} < " + prettyString(jtPtBinsHigh,1,false);
  binsToLabelStr[jtPtBinsGlobalStr] = jtPtBinsGlobalLabel;   

  std::string multiJtCutGlobalStr = "MultiJt0";
  std::string multiJtCutGlobalLabel = "N_{Jet,Reco.} >= 2";
  binsToLabelStr[multiJtCutGlobalStr] = multiJtCutGlobalLabel;   

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //  std::string jtPtBinsStrForConfig = "";
  std::vector<std::string> jtPtBinsStr;
  for(Int_t pI = 0; pI < nJtPtBins; ++pI){
    if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
    else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));
    //    jtPtBinsStrForConfig = jtPtBinsStrForConfig + std::to_string(jtPtBins[pI]) + ", ";

    binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
  }
  //  jtPtBinsStrForConfig = jtPtBinsStrForConfig + std::to_string(jtPtBins[nJtPtBins]);

  jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
  binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);


  const Double_t gammaJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);  
  const Double_t gammaMultiJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAMULTIJTDPHI", ""), doGlobalDebug);  

  const std::string gammaJtDPhiStr = "DPhi0";
  std::string gammaJtDPhiLabel = returnAllCapsString(config_p->GetValue("GAMMAJTDPHI", ""));
  if(gammaJtDPhiLabel.find("PI") != std::string::npos) gammaJtDPhiLabel.replace(gammaJtDPhiLabel.find("PI"), 2, "#pi");
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  gammaJtDPhiLabel = "|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiLabel;
  binsToLabelStr[gammaJtDPhiStr] = gammaJtDPhiLabel;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  //Initial bin grabbing is done, now make the multi-dim flattened bins
  //Block of variable declarations
  Int_t nSubJtGammaPtBins = -1;
  Float_t subJtGammaPtBinsLow, subJtGammaPtBinsHigh;
  Double_t subJtGammaPtBins[nMaxPtBins+1];
  std::string subJtGammaPtBinsStrForConfig = "";

  binFlattener subJtGammaPtBinFlattener;
  //Now create your flattened bins using the binFlattener class
  if(!subJtGammaPtBinFlattener.Init("subJtGammaPtBinFlattener", nGammaPtBins, gammaPtBins, nSubJtPtBins, subJtPtBins)){
    std::cout << "Error gdjNTupleToHist: Failure to init binFlattener, L" << __LINE__ << ". return 1" << std::endl;
    return 1;
  }

  double subJtGammaPtMax = 1.0;
  while(subJtGammaPtMax/(nGammaPtBins*nSubJtPtBins) < 1){subJtGammaPtMax *= 10.0;}

  std::vector<double> subJtGammaPtBinsV = subJtGammaPtBinFlattener.GetFlattenedBins(0.0, subJtGammaPtMax);

 
  for(unsigned int sgI = 0; sgI < subJtGammaPtBinsV.size(); ++sgI){
    ++nSubJtGammaPtBins;
    subJtGammaPtBins[nSubJtGammaPtBins] = subJtGammaPtBinsV[sgI];

    std::string binStr = std::to_string(subJtGammaPtBinsV[sgI]);

    if(binStr.find(".") != std::string::npos){
      while(binStr.rfind("0") == binStr.size()-1){
        if(binStr.rfind(".") == binStr.size()-2) break;
        binStr = binStr.substr(0, binStr.size()-1);
      }
    }

    subJtGammaPtBinsStrForConfig = subJtGammaPtBinsStrForConfig + binStr + ",";
  }

  /*
  std::cout << "subJtGammaPtBins: " << std ::endl;
  for(Int_t sgI = 0; sgI < nSubJtGammaPtBins; ++sgI){
    std::cout << " " << sgI << ": " << subJtGammaPtBins[sgI] << "-" << subJtGammaPtBins[sgI+1] << " (gammaBin=" << subJtGammaPtBinFlattener.GetBin1PosFromGlobal(sgI) << ", subJtPtBins=" << subJtGammaPtBinFlattener.GetBin2PosFromGlobal(sgI) << ")" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Gamma,subjtptbins: " << std::endl;
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    std::cout << " gI=" << gI << ", " << gammaPtBins[gI] << "-" << gammaPtBins[gI+1] << ": " << std::endl;
    for(Int_t sI = 0; sI < nSubJtPtBins; ++sI){
      std::cout << "  sI=" << sI << ", " << subJtPtBins[sI] << "-" << subJtPtBins[sI+1] << ", globalBin=" << subJtGammaPtBinFlattener.GetGlobalFromBin12Pos(gI, sI) << std::endl;
    }  
  }

  return 1;
  */
  

  const unsigned long long mixCap = config_p->GetValue("MIXCAP", 100000000);
  const unsigned long long nMixEvents = config_p->GetValue("NMIXEVENTS", 1);
  const bool doStrictMix = config_p->GetValue("DOSTRICTMIX", true);

  std::string outFileNameForTxt = outFileName;
  outFileNameForTxt.replace(outFileNameForTxt.find(".root"), 5, ".txt");

  int outFileFills = 0;
  std::ofstream outFileTxt(outFileNameForTxt.c_str());

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  //Define some numbers and names for photon es/er systematics
  const Int_t nGESSys = 2;
  const Int_t nGERSys = 2;

  const Int_t nGammaSysAndNom = 1 + nGESSys + nGERSys;
  std::string gammaSysAndNomStr[nGammaSysAndNom] = {"Nominal", "GES0", "GES1", "GER0", "GER1"};

  //Define some numbers and names for jet es/er systematics
  const Int_t nJESSys = 18;
  const Int_t nJERSys = 9;

  //1 (nominal) + nJESSys + 1 (centrality dependent rtrk JESSys manual) + nJERSys
  const Int_t nJetSysAndNom = 1 + nJESSys + 1 + nJERSys;

  std::string jetSysAndNomStr[nJetSysAndNom] = {"Nominal"}; 
  for(Int_t sI = 0; sI < nJESSys; ++sI){
    jetSysAndNomStr[1 + sI] = "JES" + std::to_string(sI);
  }
  jetSysAndNomStr[1+nJESSys] = "JESRTRK";
  for(Int_t sI = 0; sI < nJERSys; ++sI){
    jetSysAndNomStr[1 + nJESSys + 1 + sI] = "JER" + std::to_string(sI);
  }

  //We need to declare some systematics (like jet pt cut unfolding variation)
  std::vector<std::string> systStrVect = {"NOMINAL", "JTPTCUT"};
  std::vector<std::string> systTypeVect = {"NOMINAL", "UNFOLDING"};

  if(systStrVect.size() != systTypeVect.size()){
    std::cout << "ERROR: systStrVect size \'" << systStrVect.size() << "\' does not match size of systTypeVect, \'" << systTypeVect.size() << "\'. All declared syst must have corresponding type. return 1" << std::endl;
    return 1; 
  }

  Int_t sampleTag;
  UInt_t runNumber;
  UInt_t lumiBlock;
  ULong64_t eventNumber;

  const Int_t nMaxJets = 100;

  Float_t recoGammaPt_[nGammaSysAndNom];
  Float_t recoGammaPhi_;
  Float_t recoGammaEta_;

  Float_t truthGammaPt_;
  Float_t truthGammaPhi_;
  Float_t truthGammaEta_;

  Int_t nRecoJt_;
  Float_t recoJtPt_[nMaxJets][nJetSysAndNom];
  Float_t recoJtPhi_[nMaxJets];
  Float_t recoJtEta_[nMaxJets];

  Float_t truthJtPt_[nMaxJets];
  Float_t truthJtPhi_[nMaxJets];
  Float_t truthJtEta_[nMaxJets];
  Int_t truthJtFlavor_[nMaxJets];

  Int_t nTruthJtUnmatched_;
  Float_t truthJtUnmatchedPt_[nMaxJets];
  Float_t truthJtUnmatchedPhi_[nMaxJets];
  Float_t truthJtUnmatchedEta_[nMaxJets];
  Int_t truthJtUnmatchedFlavor_[nMaxJets];

  Double_t unfoldWeight_;
  Float_t unfoldCent_;
  
  TTree* unfoldTree_p = nullptr;

  if(isMC && keepResponseTree){
    unfoldTree_p = new TTree("unfoldTree_p", "");

    unfoldTree_p->Branch("sampleTag", &sampleTag, "sampleTag/I");

    unfoldTree_p->Branch("runNumber", &runNumber, "runNumber/i");
    unfoldTree_p->Branch("lumiBlock", &lumiBlock, "lumiBlock/i");
    unfoldTree_p->Branch("eventNumber", &eventNumber, "eventNumber/l");

    unfoldTree_p->Branch("recoGammaPt", recoGammaPt_, ("recoGammaPt[" + std::to_string(nGammaSysAndNom) + "]/F").c_str()); 
    unfoldTree_p->Branch("recoGammaPhi", &recoGammaPhi_, "recoGammaPhi/F"); 
    unfoldTree_p->Branch("recoGammaEta", &recoGammaEta_, "recoGammaEta/F"); 
    unfoldTree_p->Branch("truthGammaPt", &truthGammaPt_, "truthGammaPt/F"); 
    unfoldTree_p->Branch("truthGammaPhi", &truthGammaPhi_, "truthGammaPhi/F"); 
    unfoldTree_p->Branch("truthGammaEta", &truthGammaEta_, "truthGammaEta/F"); 

    unfoldTree_p->Branch("nRecoJt", &nRecoJt_, "nRecoJt/I");     
    unfoldTree_p->Branch("recoJtPt", recoJtPt_, ("recoJtPt[nRecoJt][" + std::to_string(nJetSysAndNom) + "]/F").c_str()); 
    unfoldTree_p->Branch("recoJtPhi", recoJtPhi_, "recoJtPhi[nRecoJt]/F"); 
    unfoldTree_p->Branch("recoJtEta", recoJtEta_, "recoJtEta[nRecoJt]/F"); 

    unfoldTree_p->Branch("truthJtPt", truthJtPt_, "truthJtPt[nRecoJt]/F"); 
    unfoldTree_p->Branch("truthJtPhi", truthJtPhi_, "truthJtPhi[nRecoJt]/F"); 
    unfoldTree_p->Branch("truthJtEta", truthJtEta_, "truthJtEta[nRecoJt]/F"); 
    unfoldTree_p->Branch("truthJtFlavor", truthJtFlavor_, "truthJtFlavor[nRecoJt]/I"); 

    unfoldTree_p->Branch("nTruthJtUnmatched", &nTruthJtUnmatched_, "nTruthJtUnmatched/I");
    unfoldTree_p->Branch("truthJtUnmatchedPt", truthJtUnmatchedPt_, "truthJtUnmatchedPt[nTruthJtUnmatched]/F"); 
    unfoldTree_p->Branch("truthJtUnmatchedPhi", truthJtUnmatchedPhi_, "truthJtUnmatchedPhi[nTruthJtUnmatched]/F"); 
    unfoldTree_p->Branch("truthJtUnmatchedEta", truthJtUnmatchedEta_, "truthJtUnmatchedEta[nTruthJtUnmatched]/F"); 
    unfoldTree_p->Branch("truthJtUnmatchedFlavor", truthJtUnmatchedFlavor_, "truthJtUnmatchedFlavor[nTruthJtUnmatched]/I"); 

    unfoldTree_p->Branch("unfoldWeight", &unfoldWeight_, "unfoldWeight/D"); 
    unfoldTree_p->Branch("unfoldCent", &unfoldCent_, "unfoldCent/F"); 
  }
  
  TH1D* mixingCentrality_p = nullptr;
  TH1D* mixingPsi2_p = nullptr;
  TH1D* mixingVz_p = nullptr;
  TH2D* mixingCentralityPsi2_p = nullptr;
  TH2D* mixingCentralityVz_p = nullptr;
  TH2D* mixingPsi2Vz_p = nullptr;
  TH2D* mixingCentralityPsi2_VzBinned_p[nMaxMixBins];

  if(doMix){
    if(doMixCent){
      mixingCentrality_p = new TH1D("mixingCentrality_h", ";Centrality;Counts", nMixCentBins, mixCentBins);
      if(doMixPsi2){
	mixingCentralityPsi2_p = new TH2D("mixingCentralityPsi2_h", ";Centrality;#Psi_{2}", nMixCentBins, mixCentBins, nMixPsi2Bins, mixPsi2Bins);

	if(doMixVz){
	  for(Int_t vI = 0; vI < nMixVzBins; ++vI){
	    mixingCentralityPsi2_VzBinned_p[vI] = new TH2D(("mixingCentralityPsi2_Vz" + std::to_string(vI) + "_h").c_str(), ";;", nMixCentBins, mixCentBins, nMixPsi2Bins, mixPsi2Bins);
	  }
	}
      }

      if(doMixVz) mixingCentralityVz_p = new TH2D("mixingCentralityVz_h", ";Centrality;v_{z}", nMixCentBins, mixCentBins, nMixVzBins, mixVzBins);
    }
    if(doMixPsi2){
      mixingPsi2_p = new TH1D("mixingPsi2_h", ";#Psi_{2};Counts", nMixPsi2Bins, mixPsi2Bins);
      if(doMixVz) mixingPsi2Vz_p = new TH2D("mixingPsi2Vz_h", ";#Psi_{2};v_{z}", nMixPsi2Bins, mixPsi2Bins, nMixVzBins, mixVzBins);
    }
    if(doMixVz) mixingVz_p = new TH1D("mixingVz_h", ";v_{z};Counts", nMixVzBins, mixVzBins);
  }

  

  TH1D* runNumber_p = nullptr;
  TH1D* lumiFractionPerRun_p = nullptr;
  TH1D* pthat_p = nullptr;
  TH1D* pthat_Unweighted_p = nullptr;
  TH1D* centrality_p = nullptr;
  TH1D* centrality_Unweighted_p = nullptr;
  TH1D* photonPtVCentEta_p[nMaxCentBins][nMaxSubBins+1];
  TH2D* photonGenResVCentEta_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH2D* photonEtaPt_p[nMaxCentBins];
  TH2D* photonEtaPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  
  TH1D* photonJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded
  //Make this a class, one for multijet one for single jet, so you dont have to carry stuff around
  //Inclusive jet mixmachine TEnv
  TEnv photonPtJtPtVCent_Config;
  createMixMachineTEnv(&photonPtJtPtVCent_Config, true, isMC, nJtPtBins, jtPtBinsStrForConfig, "Jet p_{T}", nGammaPtBins, gammaPtBinsStrForConfig, "Photon p_{T}");



  TEnv photonPtJtXJVCent_Config;
  createMixMachineTEnv(&photonPtJtXJVCent_Config, true, isMC, nXJBins, xjBinsStrForConfig, "Jet x_{J#gamma}", nGammaPtBins, gammaPtBinsStrForConfig, "Photon p_{T}");

  TEnv photonPtJtDPhiVCent_Config;
  createMixMachineTEnv(&photonPtJtDPhiVCent_Config, true, isMC, nDPhiBins, dPhiBinsStrForConfig, "Jet #Delta#phi_{J#gamma}", nGammaPtBins, gammaPtBinsStrForConfig, "Photon p_{T}");

  //Multi-jet mixmachine TEnv
  TEnv photonPtJtXJJVCent_Config;
  createMixMachineTEnv(&photonPtJtXJJVCent_Config, true, isMC, nXJJBins, xjjBinsStrForConfig, "Jet #vec{x}_{JJ#gamma}", nSubJtGammaPtBins, subJtGammaPtBinsStrForConfig, "Sub. Jet x Photon p_{T}");

  TEnv photonPtJtAJJVCent_Config;
  createMixMachineTEnv(&photonPtJtAJJVCent_Config, true, isMC, nAJBins, ajBinsStrForConfig, "Jet A_{JJ#gamma}", nSubJtGammaPtBins, subJtGammaPtBinsStrForConfig, "Sub. Jet x Photon p_{T}");
 
  TEnv photonPtJtDPhiJJGVCent_Config;
  createMixMachineTEnv(&photonPtJtDPhiJJGVCent_Config, true, isMC, nDPhiBins, dPhiBinsStrForConfig, "Jet #Delta#phi_{JJ#gamma}", nSubJtGammaPtBins, subJtGammaPtBinsStrForConfig, "Sub. Jet x Photon p_{T}");

  TEnv photonPtJtDPhiJJVCent_Config;
  createMixMachineTEnv(&photonPtJtDPhiJJVCent_Config, true, isMC, nDPhiBins, dPhiBinsStrForConfig, "Jet #Delta#phi_{JJ}", nSubJtGammaPtBins, subJtGammaPtBinsStrForConfig, "Sub. Jet x Photon p_{T}");

  TEnv photonPtJtDRJJVCent_Config;
 createMixMachineTEnv(&photonPtJtDRJJVCent_Config, true, isMC, nDRBins, drBinsStrForConfig, "Jet #DeltaR_{JJ}", nSubJtGammaPtBins, subJtGammaPtBinsStrForConfig, "Sub. Jet x Photon p_{T}");

  const Int_t nMaxSyst = 100;
  mixMachine* photonPtJtPtVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtXJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtXJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtAJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiJJGVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDRJJVCent_MixMachine_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];

  mixMachine* photonPtJtPtVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtXJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtXJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtAJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDPhiJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  mixMachine* photonPtJtDRJJVCent_MixMachine_Sideband_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];

  TH2D* drJJ_OneJetNoTruth_p[nMaxCentBins];
  TH2D* singleJetResponse_OneJetNoTruth_p[nMaxCentBins];
  TH2D* doubleJetResponse_OneJetNoTruth_p[nMaxCentBins];

  TH1D* photonPtVCent_RAW_p[nMaxCentBins][nBarrelAndEC];//This is needed for the purity correction - it is a pure partner of JtPt and JtXJ to be filled on every one of their fills. unclear yet whether it also works for XJJ
  TH1D* photonPtVCent_ValXWeightSum_p[nMaxCentBins][nBarrelAndEC];  

  TH1D* photonPtVCent_RAWNoTruthMatch_p[nMaxCentBins][nBarrelAndEC];
  TH1D* photonPtVCent_RAWWithTruthMatch_p[nMaxCentBins][nBarrelAndEC];
  TH1D* photonPtVCent_RAWSideband_p[nMaxCentBins][nBarrelAndEC];

  TH1D* photonPtVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH1D* photonPtVCent_TRUTHWithRecoMatch_p[nMaxCentBins][nBarrelAndEC];
  TH1D* photonPtVCent_TRUTHNoRecoMatch_p[nMaxCentBins][nBarrelAndEC];  
      
  TH1D* photonMixJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixCorrectionMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixCorrectedMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixCorrectionMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixCorrectedMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonMixJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded
  TH1D* photonSubJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubMultiJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubMultiJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonSubJtMultModVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  //2-D histograms to be unfolded  
  TH1D* photonPtVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtPtVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtXJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtDPhiVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtXJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtAJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtDPhiJJGVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtDPhiJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  TH2D* photonPtJtDRJJVCent_PURCORR_p[nMaxCentBins][nBarrelAndEC][nMaxSyst];
  
  TH2D* photonPtJtPtVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtAJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJGVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJVCent_TRUTH_p[nMaxCentBins][nBarrelAndEC];


  TH2D* photonPtJtPtVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[nMaxCentBins][nBarrelAndEC];
  
  TH2D* photonPtJtPtVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtXJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtAJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtAJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJGVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJGVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJVCent_PUREBKGD_p[nMaxCentBins][nBarrelAndEC];
  TH2D* photonPtJtDPhiJJVCent_MIXBKGD_p[nMaxCentBins][nBarrelAndEC];
  

  TH1D* photonGenJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1D* photonGenMatchedJtDPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedJtPtVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedJtEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedJtXJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedMultiJtXJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedMultiJtDPhiJJVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1D* photonGenMatchedJtMultVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  
  TH1D* photonJtFakeVCentPt_p[nMaxCentBins][nMaxSubBins+1];

  TH1D* vzPassing_p = new TH1D("vzPassing_h", ";vz;Counts", 80, -20.0, 20.0);
  TH1D* centPassing_p = new TH1D("centPassing_h", ";Centrality;Counts", 80, -0.5, 79.5);
  TH1D* truPhoPtPassing_p = new TH1D("truPhoPt_h", ";True Photon p_{T};Counts", 255, -5.0, 505.0);
  TH1D* leadingPhoPtPassing_p = new TH1D("leadingPhoPt_h", ";Leading Photon p_{T};Counts", 255, -5.0, 505.0);
  TH1D* leadingPhoEtaPassing_p = new TH1D("leadingPhoEta_h", ";Leading Photon #eta;Counts", 84, -2.1, 2.1);
  TH1D* leadingPhoIsoPassing_p = new TH1D("leadingPhoIso_h", ";Leading Photon Iso;Counts", 80, -10, 30);
  TH1D* leadingPhoTightPassing_p = new TH1D("leadingPhoTight_h", ";Leading Photon Tight;Counts", 2, -0.5, 1.5);

  TH1D* jtEtaPassing_p = new TH1D("jtEtaPassing_h", ";Jet #eta;Counts", 124, -3.1, 3.1);
  TH1D* jtGammaDRPassing_p = new TH1D("jtGammaDRPassing_h", ";Jet-Gamma #DeltaR;Counts", 44, 0.0, 4.4);
  TH1D* jtPtPassing_p = new TH1D("jtPtPassing_h", ";Jet p_{T};Counts", 255, -5, 505);
  TH1D* jtDPhiPassing_p = new TH1D("jtDPhiPassing_h", ";#Delta#phi_{#gamma,jet};Counts", 100, -0.1, TMath::Pi()+0.1);

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
    pthat_p = new TH1D(("pthat_" + systemStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    pthat_Unweighted_p = new TH1D(("pthat_Unweighted_" + systemStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    centerTitles({pthat_p, pthat_Unweighted_p});
  }
  
  if(!isPP){
    centrality_p = new TH1D(("centrality_" + systemStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
    centerTitles(centrality_p);

    if(isMC){
      centrality_Unweighted_p = new TH1D(("centrality_Unweighted_" + systemStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
      centerTitles(centrality_Unweighted_p);
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t eI = 0; eI < nGammaEtaBinsSub+1; ++eI){
      photonPtVCentEta_p[cI][eI] = new TH1D(("photonPtVCentEta_" + centBinsStr[cI] + "_" + gammaEtaBinsSubStr[eI]+ "_h").c_str(), ";#gamma p_{T} [GeV];Counts", nGammaPtBins, gammaPtBins);
      centerTitles(photonPtVCentEta_p[cI][eI]);

      if(isMC){
	photonGenResVCentEta_p[cI][eI] = new TH2D(("photonGenResVCentEta_" + centBinsStr[cI] + "_" + gammaEtaBinsSubStr[eI] + "_h").c_str(), ";Reco. Photon p_{T};Gen. Photon p_{T}", nGammaPtBins, gammaPtBins, nGammaPtBins, gammaPtBins);
      }
    }

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaVCentPt_p[cI][pI] = new TH1D(("photonEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #eta;Counts", nGammaEtaBins, gammaEtaBins);
      photonPhiVCentPt_p[cI][pI] = new TH1D(("photonPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #phi;Counts", nPhiBins, phiBins);

      photonJtDPhiVCentPt_p[cI][pI] = new TH1D(("photonJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";#Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);


      photonJtPtVCentPt_p[cI][pI] = new TH1D(("photonJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
      photonMultiJtPtVCentPt_p[cI][pI] = new TH1D(("photonMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
      photonJtEtaVCentPt_p[cI][pI] = new TH1D(("photonJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
      photonJtXJVCentPt_p[cI][pI] = new TH1D(("photonJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
      photonMultiJtXJVCentPt_p[cI][pI] = new TH1D(("photonMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

      photonMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);

      photonMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";#Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

      photonJtMultVCentPt_p[cI][pI] = new TH1D(("photonJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

      setSumW2({photonJtDPhiVCentPt_p[cI][pI], photonJtPtVCentPt_p[cI][pI], photonMultiJtPtVCentPt_p[cI][pI], photonJtEtaVCentPt_p[cI][pI], photonJtXJVCentPt_p[cI][pI], photonMultiJtXJVCentPt_p[cI][pI], photonMultiJtXJJVCentPt_p[cI][pI], photonMultiJtDPhiJJVCentPt_p[cI][pI]});
  
      if(doMix){
	photonMixJtDPhiVCentPt_p[cI][pI] = new TH1D(("photonMixJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Mixed event #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonMixJtPtVCentPt_p[cI][pI] = new TH1D(("photonMixJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonMixMultiJtPtVCentPt_p[cI][pI] = new TH1D(("photonMixMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonMixJtEtaVCentPt_p[cI][pI] = new TH1D(("photonMixJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonMixJtXJVCentPt_p[cI][pI] = new TH1D(("photonMixJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Mixed event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonMixMultiJtXJVCentPt_p[cI][pI] = new TH1D(("photonMixMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);

	photonMixMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonMixMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);

	photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonMixCorrectionMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event correction #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);

	photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonMixCorrectedMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event corrected #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);


	photonMixMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonMixMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonMixCorrectionMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event correction #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonMixCorrectedMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Mixed event corrected #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonMixJtMultVCentPt_p[cI][pI] = new TH1D(("photonMixJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Mixed Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	setSumW2({photonMixJtDPhiVCentPt_p[cI][pI], photonMixJtPtVCentPt_p[cI][pI], photonMixMultiJtPtVCentPt_p[cI][pI], photonMixJtEtaVCentPt_p[cI][pI], photonMixJtXJVCentPt_p[cI][pI], photonMixMultiJtXJVCentPt_p[cI][pI], photonMixMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtXJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtXJJVCentPt_p[cI][pI], photonMixMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectionMultiJtDPhiJJVCentPt_p[cI][pI], photonMixCorrectedMultiJtDPhiJJVCentPt_p[cI][pI], photonMixJtMultVCentPt_p[cI][pI]});
      }

      if(doMix || isPP){
	photonSubJtDPhiVCentPt_p[cI][pI] = new TH1D(("photonSubJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Subtracted event #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonSubJtPtVCentPt_p[cI][pI] = new TH1D(("photonSubJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonSubMultiJtPtVCentPt_p[cI][pI] = new TH1D(("photonSubMultiJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);

	photonSubJtEtaVCentPt_p[cI][pI] = new TH1D(("photonSubJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonSubJtXJVCentPt_p[cI][pI] = new TH1D(("photonSubJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Subtracted event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonSubMultiJtXJVCentPt_p[cI][pI] = new TH1D(("photonSubMultiJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonSubMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonSubMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);
	photonSubMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonSubMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Subtracted event #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonSubJtMultVCentPt_p[cI][pI] = new TH1D(("photonSubJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Subtracted Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	photonSubJtMultModVCentPt_p[cI][pI] = new TH1D(("photonSubJtMultModVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Subtracted Event Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);
      
	setSumW2({photonSubJtDPhiVCentPt_p[cI][pI], photonSubJtPtVCentPt_p[cI][pI], photonSubMultiJtPtVCentPt_p[cI][pI], photonSubJtEtaVCentPt_p[cI][pI], photonSubJtXJVCentPt_p[cI][pI], photonSubMultiJtXJVCentPt_p[cI][pI], photonSubMultiJtXJJVCentPt_p[cI][pI], photonSubMultiJtDPhiJJVCentPt_p[cI][pI], photonSubJtMultVCentPt_p[cI][pI], photonSubJtMultModVCentPt_p[cI][pI]});
      }
  
      if(isMC){
	photonGenJtDPhiVCentPt_p[cI][pI] = new TH1D(("photonGenJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Gen. #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonGenJtPtVCentPt_p[cI][pI] = new TH1D(("photonGenJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonGenJtEtaVCentPt_p[cI][pI] = new TH1D(("photonGenJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonGenJtXJVCentPt_p[cI][pI] = new TH1D(("photonGenJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen. x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonGenMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen. #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);
	photonGenMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonGenMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen. #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonGenJtMultVCentPt_p[cI][pI] = new TH1D(("photonGenJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen. Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);

	photonGenMatchedJtDPhiVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedJtDPhiVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_h").c_str(), ";Gen.-matched #Delta#phi_{#gamma,jet};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);	
	photonGenMatchedJtPtVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedJtPtVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}/N_{#gamma}", nJtPtBins, jtPtBins);
	photonGenMatchedJtEtaVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedJtEtaVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched #gamma-tagged Jet #eta;N_{#gamma,jet}/N_{#gamma}", nJtEtaBins, jtEtaBins);
	photonGenMatchedJtXJVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedJtXJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";Gen.-matched x_{J,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJBins, xjBins);
	photonGenMatchedMultiJtXJJVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedMultiJtXJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen.-matched #vec{x}_{JJ,#gamma};N_{#gamma,jet}/N_{#gamma}", nXJJBins, xjjBins);
	photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedMultiJtDPhiJJVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_" + multiJtCutGlobalStr + "_h").c_str(), ";Gen.-matched #Delta#phi_{JJ};N_{#gamma,jet}/N_{#gamma}", nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

	photonGenMatchedJtMultVCentPt_p[cI][pI] = new TH1D(("photonGenMatchedJtMultVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), (";Gen.-matched Jet multiplicity (" + prettyString(jtPtBins[0], 1, false) + " < p_{T} < " + prettyString(jtPtBins[nJtPtBins], 1, false) + ");Counts").c_str(), 7, -0.5, 6.5);
      
             

	photonJtFakeVCentPt_p[cI][pI] = new TH1D(("photonJtFakeVCentPt_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_" + jtPtBinsGlobalStr + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} ;#frac{N_{Fake jets}}{N_{All jets}}", nJtPtBins, jtPtBins);
	
      	setSumW2({photonGenJtDPhiVCentPt_p[cI][pI], photonGenJtPtVCentPt_p[cI][pI], photonGenJtEtaVCentPt_p[cI][pI], photonGenJtXJVCentPt_p[cI][pI], photonGenMultiJtXJJVCentPt_p[cI][pI], photonGenMultiJtDPhiJJVCentPt_p[cI][pI], photonGenJtMultVCentPt_p[cI][pI], photonGenMatchedJtDPhiVCentPt_p[cI][pI], photonGenMatchedJtPtVCentPt_p[cI][pI], photonGenMatchedJtEtaVCentPt_p[cI][pI], photonGenMatchedJtXJVCentPt_p[cI][pI], photonGenMatchedMultiJtXJJVCentPt_p[cI][pI], photonGenMatchedMultiJtDPhiJJVCentPt_p[cI][pI], photonGenMatchedJtMultVCentPt_p[cI][pI], photonJtFakeVCentPt_p[cI][pI]});
      }
    }
	    
    if(isMC){
      drJJ_OneJetNoTruth_p[cI] = new TH2D(("drJJ_OneJetNoTruth_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Jet p_{T};#DeltaR_{JJ,Lead}", 70, 30.0, 100.0, 50, 0.0, 1.0);      

      singleJetResponse_OneJetNoTruth_p[cI] = new TH2D(("singleJetResponse_OneJetNoTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Jet p_{T};Response (Reco./Gen.)", nJtPtBins, jtPtBins, 100, 0.0, 2.0);
      doubleJetResponse_OneJetNoTruth_p[cI] = new TH2D(("doubleJetResponse_OneJetNoTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Jet p_{T};Response (Reco./Gen.)", nJtPtBins, jtPtBins, 100, 0.0, 2.0);
    }
    
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      photonPtJtPtVCent_Config.SetValue("ISMC", isMC);
      photonPtJtXJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiJJGVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDPhiJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtDRJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtXJJVCent_Config.SetValue("ISMC", isMC);
      photonPtJtAJJVCent_Config.SetValue("ISMC", isMC);
      

      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	photonPtJtPtVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtPtVCent_Config);
	photonPtJtXJVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtXJVCent_Config);
	photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, inclusiveFlag, &photonPtJtDPhiVCent_Config);
	photonPtJtXJJVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtXJJVCent_Config);
	photonPtJtAJJVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtAJJVCent_Config);
	photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDPhiJJGVCent_Config);
	photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDPhiJJVCent_Config);
	photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI] = new mixMachine("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr, multiFlag, &photonPtJtDRJJVCent_Config);
      }
    

      //Data or MC, Sideband should always be treated as data
      photonPtJtPtVCent_Config.SetValue("ISMC", 0);
      photonPtJtXJVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiJJGVCent_Config.SetValue("ISMC", 0);
      photonPtJtDPhiJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtDRJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtXJJVCent_Config.SetValue("ISMC", 0);
      photonPtJtAJJVCent_Config.SetValue("ISMC", 0);
      
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtPtVCent_Config);
	photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtXJVCent_Config);
	photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", inclusiveFlag, &photonPtJtDPhiVCent_Config);
	photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtXJJVCent_Config);
	photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtAJJVCent_Config);
	photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDPhiJJGVCent_Config);
	photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDPhiJJVCent_Config);
	photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI] = new mixMachine("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_Sideband", multiFlag, &photonPtJtDRJJVCent_Config);       
      }               
      photonPtVCent_RAW_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAW_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
      photonPtVCent_ValXWeightSum_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_ValXWeightSum_h").c_str(), ";Photon p_{T};Weight x Value", nGammaPtBins, gammaPtBins);
      
      if(isMC){
	photonPtVCent_RAWWithTruthMatch_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAWWithTruthMatch_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
	photonPtVCent_RAWNoTruthMatch_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAWNoTruthMatch_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
      }

      photonPtVCent_RAWSideband_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_RAWSideband_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);


      if(isMC){
	photonPtVCent_TRUTH_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_TRUTH_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
	photonPtVCent_TRUTHWithRecoMatch_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_TRUTHWithRecoMatch_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
	photonPtVCent_TRUTHNoRecoMatch_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_TRUTHNoRecoMatch_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);
      }
      
      setSumW2({photonPtVCent_RAW_p[cI][eI], photonPtVCent_ValXWeightSum_p[cI][eI], photonPtVCent_RAWSideband_p[cI][eI]});
      if(isMC) setSumW2({photonPtVCent_RAWWithTruthMatch_p[cI][eI], photonPtVCent_RAWNoTruthMatch_p[cI][eI], photonPtVCent_TRUTH_p[cI][eI], photonPtVCent_TRUTHWithRecoMatch_p[cI][eI], photonPtVCent_TRUTHNoRecoMatch_p[cI][eI]});	       
      
      if(/*doMix*/true || isPP){
	photonPtVCent_PURCORR_p[cI][eI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_PURCORR_h").c_str(), ";Photon p_{T};Counts (Weighted)", nGammaPtBins, gammaPtBins);

	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  photonPtJtPtVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);

	  photonPtJtXJVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJJBins, xjjBins, nSubJtGammaPtBins, subJtGammaPtBins);
	  photonPtJtAJJVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nSubJtGammaPtBins, subJtGammaPtBins);
	  photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nSubJtGammaPtBins, subJtGammaPtBins);
	  photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #Delta#phi_{JJ};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nSubJtGammaPtBins, subJtGammaPtBins);
	  photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI] = new TH2D(("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + systStrVect[systI] + "_" + gammaJtDPhiStr + "_PURCORR_h").c_str(), ";Reco. #DeltaR_{JJ};Reco. Photon p_{T}", nDRBins, drBins, nSubJtGammaPtBins, subJtGammaPtBins);
	  
	  setSumW2({photonPtJtPtVCent_PURCORR_p[cI][eI][systI], photonPtJtXJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI], photonPtJtXJJVCent_PURCORR_p[cI][eI][systI], photonPtJtAJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]});
	}
	
	if(isMC){
	  photonPtJtPtVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJJBins, xjjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI] = new TH2D(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTH_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  

	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJJBins, xjjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI] = new TH2D(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_TRUTHMATCHEDRECO_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);

	  photonPtJtPtVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtPtVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. Jet p_{T};Reco. Photon p_{T}", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtXJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. x_{J};Reco. Photon p_{T}", nXJBins, xjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtDPhiVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{J#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJJBins, xjjBins, nGammaPtBins, gammaPtBins);
	  photonPtJtXJJVCent_MIXBKGD_p[cI][eI] = new TH2D(("photonPtJtXJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #vec{x}_{JJ#gamma};Reco. Photon p_{T}", nXJJBins, xjjBins, nGammaPtBins, gammaPtBins);

	  photonPtJtAJJVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
	  photonPtJtAJJVCent_MIXBKGD_p[cI][eI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. A_{JJ#gamma};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);

	  photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI] = new TH2D(("photonPtJtDPhiJJGVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);

	  photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI] = new TH2D(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_PUREBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);
	  photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI] = new TH2D(("photonPtJtDPhiJJVCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI] + "_" + gammaJtDPhiStr + "_MIXBKGD_h").c_str(), ";Reco. #Delta#phi_{JJ#gamma};Reco. Photon p_{T}", nDPhiBins, dPhiBins, nGammaPtBins, gammaPtBins);

	  setSumW2({photonPtJtPtVCent_TRUTH_p[cI][eI], photonPtJtXJVCent_TRUTH_p[cI][eI], photonPtJtDPhiVCent_TRUTH_p[cI][eI], photonPtJtXJJVCent_TRUTH_p[cI][eI], photonPtJtAJJVCent_TRUTH_p[cI][eI], photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI], photonPtJtDPhiJJVCent_TRUTH_p[cI][eI], photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI], photonPtJtPtVCent_PUREBKGD_p[cI][eI], photonPtJtXJVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiVCent_PUREBKGD_p[cI][eI], photonPtJtXJJVCent_PUREBKGD_p[cI][eI], photonPtJtXJJVCent_MIXBKGD_p[cI][eI], photonPtJtAJJVCent_PUREBKGD_p[cI][eI], photonPtJtAJJVCent_MIXBKGD_p[cI][eI], photonPtJtDPhiJJGVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiJJGVCent_MIXBKGD_p[cI][eI], photonPtJtDPhiJJVCent_PUREBKGD_p[cI][eI], photonPtJtDPhiJJVCent_MIXBKGD_p[cI][eI]});
	}
      }
    }

    photonEtaPt_p[cI] = new TH2D(("photonEtaPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;#gamma p_{T} [GeV]", nGammaEtaBins, gammaEtaBins, nGammaPtBins, gammaPtBins);

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI] = new TH2D(("photonEtaPhi_" + centBinsStr[cI] + "_" + gammaPtBinsStr[pI] + "_h").c_str(), ";#gamma #eta;#gamma #phi", nGammaEtaBins, gammaEtaBins, nPhiBins, phiBins);
    }
  }

  TFile* inFile_p = nullptr;
  TTree* inTree_p = nullptr;
  std::map<Int_t, Int_t> sampleTagCounter;
  std::vector<Float_t> uniqueMinPthats;
  sampleHandler sHandler;
  sHandler.Init("mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217");//we wont actually use this sample but i need to initialize w/ something

  Int_t runMin = TMath::Power(10,10);
  Int_t runMax = -1;

  ULong64_t nEntriesAllFiles = 0;

  for(unsigned int fileI = 0; fileI < inROOTFileNames.size(); ++fileI){
    inFile_p = new TFile(inROOTFileNames[fileI].c_str(), "READ");
    inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

    nEntriesAllFiles += (ULong64_t)inTree_p->GetEntries();
  
    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("runNumber", 1);
    runMin = TMath::Min((Int_t)inTree_p->GetMinimum("runNumber"), runMin);
    runMax = TMath::Max((Int_t)inTree_p->GetMaximum("runNumber"), runMax);    
  
    if(isMC){
      inTree_p->SetBranchStatus("*", 0);
      inTree_p->SetBranchStatus("sampleTag", 1);
      
      inTree_p->SetBranchAddress("sampleTag", &sampleTag);
      for(Int_t entry = 0; entry < inTree_p->GetEntries(); ++entry){
	inTree_p->GetEntry(entry);
	
	++(sampleTagCounter[sampleTag]);
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }  

  Int_t nRunBins = runMax - runMin;
  Float_t runMinF = ((Float_t)runMin) - 0.5;
  Float_t runMaxF = ((Float_t)runMax) + 0.5;
  
  if(isMC){    
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
  outFile_p->cd();

  runNumber_p = new TH1D(("runNumber_" + systemStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);
  if(!isMC) lumiFractionPerRun_p = new TH1D(("lumiFractionPerRun_" + systemStr + "_h").c_str(), ";Run;Fraction of Lumiblocks", nRunBins+1, runMinF, runMaxF);

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
    
  //define the trigger names we will restrict data samples to
  const std::string hltNamePbPb = "HLT_g20_loose_ion";
  const std::string hltNamePP = "HLT_g35_loose_L1EM15";

  std::string hltName = "";
  int hltPos = -1;
  if(!isMC){
    if(isPP) hltName = hltNamePP;
    else hltName = hltNamePbPb;
  }

  //Grab the hltbranches for some basic prescale checks
  //Now that we do this as a loop over many files, re-open to get branches
  inFile_p = new TFile(inROOTFileNames[0].c_str(), "READ");
  inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  std::vector<std::string> listOfBranches = getVectBranchList(inTree_p);
  std::vector<std::string> hltList;
  std::vector<std::string> hltListPres;
  std::string hltStr = "HLT_";
  std::string prescaleStr = "_prescale";

  //Close file and delete now that we have our branch list
  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //HLTLists are built
  for(auto const & branchStr : listOfBranches){
    if(branchStr.size() < hltStr.size()) continue;
    if(!isStrSame(branchStr.substr(0, hltStr.size()), hltStr)) continue;

    if(branchStr.find(prescaleStr) != std::string::npos) hltListPres.push_back(branchStr);
    else hltList.push_back(branchStr);

    if(isStrSame(branchStr, hltName)){
      hltPos = hltList.size()-1;
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  bool allHLTPrescalesFound = true;
  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    std::string modStr = hltList[hI] + prescaleStr;
    if(!vectContainsStr(modStr, &hltListPres)){
      allHLTPrescalesFound = false;
      std::cout << "HLT " << hltList[hI] << " has no prescale. return 1" << std::endl;
    }
  }
  if(!allHLTPrescalesFound) return 1;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  float hltPrescaleDelta = 0.01;
  std::vector<bool*> hltVect;
  std::vector<float*> hltPrescaleVect;

  Float_t pthat;
  Float_t sampleWeight;
  Float_t ncollWeight;
  Float_t fullWeight;
  Float_t fcalA_et, fcalC_et;
  Float_t evtPlane2Phi;
  std::vector<float>* vert_z_p=nullptr;
  const Double_t mmToCMDivFactor = 10.;

  std::vector<float> fullWeightVals;

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta, truthPhotonIso4;
  
  std::vector<float>* photon_pt_p=nullptr;
  //Sys1 and 2 correspond to Photon ES variations up and down
  std::vector<float>* photon_pt_sys1_p=nullptr;
  std::vector<float>* photon_pt_sys2_p=nullptr;
  //Sys3 and 4 correspond to Photon ER variations up and down
  std::vector<float>* photon_pt_sys3_p=nullptr;
  std::vector<float>* photon_pt_sys4_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;  
  std::vector<float>* photon_etcone_p=nullptr;
  std::vector<float>* photon_correctedIso_p=nullptr;
  
  std::vector<float>* aktRhi_etajes_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_etajes_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_etajes_jet_phi_p=nullptr;
  std::vector<std::vector<float>* > aktRhi_etajes_jet_pt_sysJES_p;
  std::vector<std::vector<float>* > aktRhi_etajes_jet_pt_sysJER_p;
  std::vector<float>* aktRhi_insitu_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_insitu_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_insitu_jet_phi_p=nullptr;
  std::vector<int>* aktRhi_truthpos_p=nullptr;

  std::vector<std::vector<float> >* aktRhi_jetconstit_pt_p=nullptr;
  std::vector<std::vector<float> >* aktRhi_jetconstit_eta_p=nullptr;
  std::vector<std::vector<float> >* aktRhi_jetconstit_phi_p=nullptr;

  for(Int_t jI = 0; jI < nJESSys; ++jI){aktRhi_etajes_jet_pt_sysJES_p.push_back(nullptr);}
  for(Int_t jI = 0; jI < nJERSys; ++jI){aktRhi_etajes_jet_pt_sysJER_p.push_back(nullptr);}  

  std::vector<float>* aktR_truth_jet_pt_p=nullptr;
  std::vector<float>* aktR_truth_jet_eta_p=nullptr;
  std::vector<float>* aktR_truth_jet_phi_p=nullptr;
  std::vector<int>* aktR_truth_jet_recopos_p=nullptr;
  std::vector<int>* aktR_truth_jet_partonid_p=nullptr;

  TFile* mixFile_p = nullptr;
  TTree* mixTree_p = nullptr;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doMix){
    //First create a map of signal events for different categories
    for(unsigned int fileI = 0; fileI < inROOTFileNames.size(); ++fileI){    
      //If we are targetting events for studies we should not do a full multi-file processing
      if((nStartEvtStr.size() != 0 || nMaxEvtStr.size() != 0) && fileI != 0) continue;

      inFile_p = new TFile(inROOTFileNames[fileI].c_str(), "READ");
      inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
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
      
      //TTree processing; Range of entry handling - default is full range but you can set a start and an endpoint - this just handles that in a somewhat smart way
      ULong64_t nEntriesTemp = inTree_p->GetEntries();
      ULong64_t nEntriesStart = 0;     
      
      if(nStartEvtStr.size() != 0) nEntriesStart = TMath::Min(nEntriesTemp-1, (ULong64_t)nStartEvt);
      if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, (ULong64_t)nMaxEvt);
      const ULong64_t nSigEntries = nEntriesTemp;
      
      for(ULong64_t entry = nEntriesStart; entry < nEntriesStart + nSigEntries; ++entry){
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

      inFile_p->Close();
      delete inFile_p;
    }
  
    //Now create the corresponding map of mixed events for different mixing categories + save the relevant jet collections      
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
    
    ULong64_t nEntriesTemp = mixTree_p->GetEntries();
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
	if(aktRhi_insitu_jet_pt_p->at(jI) < jtPtBinsLowReco) continue;
	if(aktRhi_insitu_jet_pt_p->at(jI) >= jtPtBinsHighReco) continue;
	
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
    if(tempCounts.size() != 0) std::cout << "MEDIAN NUMBER TO MIX: " << tempCounts[tempCounts.size()/2] << std::endl;
    else std::cout << "NO MEDIAN COUNTS" << std::endl;
  }//end if(doMix){

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //We need to initialize some things before getting to the primary analysis loop
  //Init. hlt branch vectors
  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    hltVect.push_back(new bool(false));
    hltPrescaleVect.push_back(new float(0.0));
  }

  //Initialize maps of photon counts
  std::map<int, int> perEventPhotonCounts; //how many events w/ 0, 1, 2 photons, etc.
  for(unsigned int i = 0; i < 10; ++i){
    perEventPhotonCounts[i] = 0;
  }

  //How many events w/ photon fall into cent, pt category
  std::vector<std::vector<Double_t> > gammaCountsPerPtCent;
  for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
    gammaCountsPerPtCent.push_back({});

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      gammaCountsPerPtCent[pI].push_back(0.0);
    }
  }

  //Init recoJtPtMin
  Double_t recoJtPtMin = 100000.;

  //Define nDiv outputs + needed iterator across files
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntriesAllFiles/20);
  ULong64_t currEntry = 0;

  //bool for checking if you are missing trigger fires in data
  bool didOneFireMiss = false;

  std::vector<int> skippedCent;
  std::vector<std::string> eventCutStr = {"v_{z} Cut",
					  "HLT req. (Data Only)",
					  "Centrality in Range (Pb+Pb Only)",
					  "Truth photon p_{T} in MC range"};  
  
  std::map<int, int> eventCounter;
  for(int i = 0; i < 100; ++i){
    eventCounter[i] = 0;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::cout << "Begin processing files..." << std::endl;
  for(unsigned int fileI = 0; fileI < inROOTFileNames.size(); ++fileI){
    //If we are targetting events for studies we should not do a full multi-file processing           
    if((nStartEvtStr.size() != 0 || nMaxEvtStr.size() != 0) && fileI != 0) continue;

    std::cout << "Begin processing file " << fileI << "/" << inROOTFileNames.size() << ":" << std::endl;
    std::cout << "  \'" << inROOTFileNames[fileI] << "\'" << std::endl;
    
    inFile_p = new TFile(inROOTFileNames[fileI].c_str(), "READ");
    inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    inTree_p->SetBranchStatus("*", 0);
    
    for(unsigned int hI = 0; hI < hltList.size(); ++hI){
      inTree_p->SetBranchStatus(hltList[hI].c_str(), 1);
      inTree_p->SetBranchStatus(hltListPres[hI].c_str(), 1);
    }  
    
    inTree_p->SetBranchStatus("runNumber", 1);
    inTree_p->SetBranchStatus("eventNumber", 1);
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
      inTree_p->SetBranchStatus("truthPhotonIso4", 1);
    }
    
    if(!isPP){
      inTree_p->SetBranchStatus("fcalA_et", 1);
      inTree_p->SetBranchStatus("fcalC_et", 1);
      if(doMixPsi2) inTree_p->SetBranchStatus("evtPlane2Phi", 1);
    }
    
    inTree_p->SetBranchStatus("vert_z", 1);
    
    inTree_p->SetBranchStatus("photon_pt", 1);
    inTree_p->SetBranchStatus("photon_pt_sys1", 1);
    inTree_p->SetBranchStatus("photon_pt_sys2", 1);
    inTree_p->SetBranchStatus("photon_pt_sys3", 1);
    inTree_p->SetBranchStatus("photon_pt_sys4", 1);
    inTree_p->SetBranchStatus("photon_eta", 1);
    inTree_p->SetBranchStatus("photon_phi", 1);
    inTree_p->SetBranchStatus("photon_tight", 1);
    inTree_p->SetBranchStatus(("photon_etcone" + std::to_string(isolationR) + "0").c_str(), 1);
    
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
      inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_partonid").c_str(), 1);
    }
    
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(unsigned int hI = 0; hI < hltList.size(); ++hI){
      inTree_p->SetBranchAddress(hltList[hI].c_str(), hltVect[hI]);
      inTree_p->SetBranchAddress(hltListPres[hI].c_str(), hltPrescaleVect[hI]);
    }
    
    inTree_p->SetBranchAddress("runNumber", &runNumber);
    inTree_p->SetBranchAddress("eventNumber", &eventNumber);
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
      inTree_p->SetBranchAddress("truthPhotonIso4", &truthPhotonIso4);
    }
    
    if(!isPP){
      inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
      inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
      if(doMixPsi2) inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
    }
    
    inTree_p->SetBranchAddress("vert_z", &vert_z_p);
    
    inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
    inTree_p->SetBranchAddress("photon_pt_sys1", &photon_pt_sys1_p);
    inTree_p->SetBranchAddress("photon_pt_sys2", &photon_pt_sys2_p);
    inTree_p->SetBranchAddress("photon_pt_sys3", &photon_pt_sys3_p);
    inTree_p->SetBranchAddress("photon_pt_sys4", &photon_pt_sys4_p);
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
      inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_partonid").c_str(), &aktR_truth_jet_partonid_p);
    }
    
    ULong64_t nEntriesTemp = inTree_p->GetEntries();
    ULong64_t nEntriesStart = 0;
    if(nStartEvtStr.size() != 0) nEntriesStart = TMath::Min(nEntriesTemp-1, (ULong64_t)nStartEvt);
    if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, nMaxEvt);

    const ULong64_t nEntries = nEntriesTemp;
    std::cout << "  Processing " << nEntries << " events in file..." << std::endl;
    
    //Main signal processing loop
    for(ULong64_t entry = nEntriesStart; entry < nEntries+nEntriesStart; ++entry){
      if(currEntry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries+nEntriesStart << "..." << std::endl;
      ++currEntry;
      inTree_p->GetEntry(entry);

      //Cut 1: z vertex
      double vert_z = vert_z_p->at(0);
      vert_z /= mmToCMDivFactor;
      if(vert_z <= -15. || vert_z >= 15.) continue;      

      vzPassing_p->Fill(vert_z);
      ++(eventCounter[0]);
      
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
      
      //Trigger processing
      if(!isMC){
	if(!didOneFireMiss){//only check this once per input
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
	
	//Cut 2: data only require single designated trigger fire
	if(!(*(hltVect[hltPos]))) continue;
      }
      ++(eventCounter[1]);
      
      //Grab the pthat of the sample if applicable (i.e. MC)
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

      //Cut 3: Centrality within selected range
      if(centPos < 0){
	bool vectContainsCent = vectContainsInt((Int_t)cent, &skippedCent);
	if(!vectContainsCent){
	  std::cout << "gdjNTupleToHist Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
	  skippedCent.push_back((Int_t)cent);
	}     
	continue;
      }
      centPassing_p->Fill(cent);
      ++(eventCounter[2]);
      
      //In data, since we deal in unprescaled data, the event-by-event weights should always be 1
      if(!isMC) fullWeight = 1.0;
      
      bool fullWeightFound = false;
      for(unsigned int fI = 0; fI < fullWeightVals.size(); ++fI){
	if(TMath::Abs(fullWeight - fullWeightVals[fI]) < TMath::Power(10,-50)){
	  fullWeightFound = true;
	  break;
	}
      }
      if(!fullWeightFound) fullWeightVals.push_back(fullWeight);
      
      if(isMC){
	//Cut 4: Make sure the photon is w/in the designated MC sample range
	if(truthPhotonPt < minPthat) continue;
	if(truthPhotonPt >= maxPthat) continue;
	
	truPhoPtPassing_p->Fill(truthPhotonPt);
      }
      ++(eventCounter[3]);
      
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
    
      //Apply pt range cut and isolation now
      //Values hardcoded, find definitions here ATL-COM-PHYS-2021-215
      Bool_t truthPhoIsGood = false;
      const Float_t truthIsoCutInGeV = 5.0;
      if(isMC){
	truthPhoIsGood = photonEtaIsGood(truthPhotonEta);
	truthPhoIsGood = truthPhoIsGood && truthPhotonPt >= gammaPtBinsLow && truthPhotonPt < gammaPtBinsHigh;
	truthPhoIsGood = truthPhoIsGood && truthPhotonIso4 < truthIsoCutInGeV;
	
	//Cut 5: InMC only require the truth level photon fall in your bins and is isolated
	if(!truthPhoIsGood){
	  //	std::cout << "photon  no good, here is iso: " << truthPhotonIso4 << std::endl;
	  continue;
	}
      }
      ++(eventCounter[4]);
    
      //0.2 or less is the defined matching dR found in ATL-COM-PHYS-2021-215
      const Float_t truthPhoRecoDR = 0.2;
      Int_t truthPhoRecoPos = -1;
      Bool_t truthPhoHasGoodReco = false;
          
      //Clean and re-new the vector for corrected isolated
      if(photon_correctedIso_p != nullptr){
	photon_correctedIso_p->clear();
	delete photon_correctedIso_p;
      }
      photon_correctedIso_p = new std::vector<float>;

      //Since truth is independent of syst., define truth pho. position here
      std::vector<int> barrelECFillTruth = {2}; //Always fill at the inclusive position
      if(isMC){
	if(photonEtaIsBarrel(truthPhotonEta)) barrelECFillTruth.push_back(0);
	else if(photonEtaIsEC(truthPhotonEta)) barrelECFillTruth.push_back(1);
	else{
	  std::cout << "Photon w/ eta=\'" << truthPhotonEta << "\' is not registering as barrel or endcap check for bugs" << std::endl;
	  return 1;
	}
      }      
    
      //We will have to construct some alt histograms for systematics so we will do this in a loop
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	//Go thru the photons to find the reco. match of truthPhoton + calc corrected iso
	for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
	  //Should only be done once - get the correctedIso and populate the vector
	  if(systI == 0){
	    Float_t correctedIso = photon_etcone_p->at(pI);
	    
	    if(doPtIsoCorrection){
	      if(doCentIsoCorrection) correctedIso = getCorrectedPhotonIsolation(isPP, correctedIso, photon_pt_p->at(pI), photon_eta_p->at(pI), cent);
	      else correctedIso = getPtCorrectedPhotonIsolation(correctedIso, photon_pt_p->at(pI), photon_eta_p->at(pI));
	    }       
	    
	    photon_correctedIso_p->push_back(correctedIso);
	  }

	  if(isMC){
	  //Continue on anything < 15
	    if(photon_pt_p->at(pI) < 15.0) continue;
	    
	    Float_t tempDR = getDR(truthPhotonEta, truthPhotonPhi, photon_eta_p->at(pI), photon_phi_p->at(pI), "Entry: " + std::to_string(entry));
	    
	    Float_t truthPhoRecoMatchPt = -999;
	    if(truthPhoRecoPos >= 0) truthPhoRecoMatchPt = photon_pt_p->at(truthPhoRecoPos);
	    
	    if(tempDR < truthPhoRecoDR && photon_pt_p->at(pI) > truthPhoRecoMatchPt){
	      if(truthPhoRecoPos >= 0 && false){ //Convinced that 
		std::cout << "UHOH MULTIPLE PHOTON RECO MATCH TRUTH (ENTRY: " << entry << ")" << std::endl;
		std::cout << " 1st Match (pt, eta, phi): " << photon_pt_p->at(truthPhoRecoPos) << ", " << photon_eta_p->at(truthPhoRecoPos) << ", " << photon_phi_p->at(truthPhoRecoPos) << std::endl;
		std::cout << " 2nd Match (pt, eta, phi): " << photon_pt_p->at(pI) << ", " << photon_eta_p->at(pI) << ", " << photon_phi_p->at(pI) << std::endl;
	      }
	      
	      truthPhoRecoPos = pI;
	      truthPhoHasGoodReco = true;
	    }
	  }
	}//End first photon loop after finding reco match for truth and doing good iso calc
	
	//Check the reco is good
	if(truthPhoRecoPos >= 0){
	  recoGammaPt_[0] = photon_pt_p->at(truthPhoRecoPos);
	  recoGammaPt_[1] = photon_pt_sys1_p->at(truthPhoRecoPos);
	  recoGammaPt_[2] = photon_pt_sys2_p->at(truthPhoRecoPos);
	  recoGammaPt_[3] = photon_pt_sys3_p->at(truthPhoRecoPos);
	  recoGammaPt_[4] = photon_pt_sys4_p->at(truthPhoRecoPos);

	  recoGammaEta_ = photon_eta_p->at(truthPhoRecoPos);
	  Float_t tempCorrectedIso = photon_correctedIso_p->at(truthPhoRecoPos);
	  
	  if(!photonEtaIsGood(recoGammaEta_)) truthPhoHasGoodReco = false;
	  else if(!photon_tight_p->at(truthPhoRecoPos)) truthPhoHasGoodReco = false;
	  else if(!isIsolatedPhoton(isPP, doPtIsoCorrection, tempCorrectedIso)) truthPhoHasGoodReco = false;    
	  else if(recoGammaPt_[0] < gammaPtBinsLowReco) truthPhoHasGoodReco = false;    
	  else if(recoGammaPt_[0] >= gammaPtBinsHighReco) truthPhoHasGoodReco = false;          
	}
	
	//Response TTree filling
	if(isMC && systI == 0){
	  bool isSignal = false;

	  if(truthPhoRecoPos >= 0){
	    isSignal = photon_tight_p->at(truthPhoRecoPos);
	    isSignal = isSignal && isIsolatedPhoton(isPP, doPtIsoCorrection, photon_correctedIso_p->at(truthPhoRecoPos));
	  }
      	  
	  if(keepResponseTree){
	    for(Int_t gI = 0; gI < nGammaSysAndNom; ++gI){
	      recoGammaPt_[gI] = -999;
	    }	    
	    recoGammaPhi_ = -999;
	    recoGammaEta_ = -999;//photon_phi_p->at(phoPos);
	    truthGammaPt_ = truthPhotonPt;	
	    truthGammaPhi_ = truthPhotonPhi;
	    truthGammaEta_ = truthPhotonEta;
	    
	    if(truthPhoRecoPos >= 0 && truthPhoHasGoodReco && isSignal){	  	  
	      recoGammaPt_[0] = photon_pt_p->at(truthPhoRecoPos);
	      recoGammaPt_[1] = photon_pt_sys1_p->at(truthPhoRecoPos);
	      recoGammaPt_[2] = photon_pt_sys2_p->at(truthPhoRecoPos);
	      recoGammaPt_[3] = photon_pt_sys3_p->at(truthPhoRecoPos);
	      recoGammaPt_[4] = photon_pt_sys4_p->at(truthPhoRecoPos);
	      recoGammaPhi_ = photon_phi_p->at(truthPhoRecoPos);
	      recoGammaEta_ = photon_eta_p->at(truthPhoRecoPos);
	    }
	    
	    unfoldWeight_ = fullWeight;
	    unfoldCent_ = cent;
	    
	    nTruthJtUnmatched_ = 0;      
	    nRecoJt_ = 0;      
	    
	    for(unsigned int jI = 0; jI < aktR_truth_jet_pt_p->size(); ++jI){
	      int recoPos = aktR_truth_jet_recopos_p->at(jI);
	      
	      if(recoPos >= 0){
		recoJtPt_[nRecoJt_][0] = aktRhi_etajes_jet_pt_p->at(recoPos);
		
		for(Int_t jsI = 0; jsI < nJESSys; ++jsI){
		  recoJtPt_[nRecoJt_][1 + jsI] = aktRhi_etajes_jet_pt_sysJES_p[jsI]->at(recoPos);
		}
		//Manually insert the rtrk centrality dependent Pb+Pb specific uncertainty
		//in pp, just set it to zero for simplicity
		if(!isPP) recoJtPt_[nRecoJt_][1+nJESSys] = getRTrkJESSysPt(cent, aktRhi_etajes_jet_pt_p->at(recoPos));
		else recoJtPt_[nRecoJt_][1+nJESSys] = aktRhi_etajes_jet_pt_p->at(recoPos);
		
		for(Int_t jsI = 0; jsI < nJERSys; ++jsI){
		  recoJtPt_[nRecoJt_][1 + nJESSys + 1 + jsI] = aktRhi_etajes_jet_pt_sysJER_p[jsI]->at(recoPos);
		}
		
		recoJtPhi_[nRecoJt_] = aktRhi_etajes_jet_phi_p->at(recoPos);
		recoJtEta_[nRecoJt_] = aktRhi_etajes_jet_eta_p->at(recoPos);
		
		truthJtPt_[nRecoJt_] = aktR_truth_jet_pt_p->at(jI);
		truthJtPhi_[nRecoJt_] = aktR_truth_jet_phi_p->at(jI);
		truthJtEta_[nRecoJt_] = aktR_truth_jet_eta_p->at(jI);
		truthJtFlavor_[nRecoJt_] = aktR_truth_jet_partonid_p->at(jI);
		
		++nRecoJt_;
	      }
	      else{
		truthJtUnmatchedPt_[nTruthJtUnmatched_] = aktR_truth_jet_pt_p->at(jI);
		truthJtUnmatchedPhi_[nTruthJtUnmatched_] = aktR_truth_jet_phi_p->at(jI);
		truthJtUnmatchedEta_[nTruthJtUnmatched_] = aktR_truth_jet_eta_p->at(jI);
		truthJtUnmatchedFlavor_[nTruthJtUnmatched_] = aktR_truth_jet_partonid_p->at(jI);
		
		++nTruthJtUnmatched_;
	      }
	    }
	    unfoldTree_p->Fill();
	  }
	
	  for(auto const barrelEC : barrelECFillTruth){
	    fillTH1(photonPtVCent_TRUTH_p[centPos][barrelEC], truthPhotonPt, fullWeight);
	    if(truthPhoHasGoodReco) fillTH1(photonPtVCent_TRUTHWithRecoMatch_p[centPos][barrelEC], truthPhotonPt, fullWeight);
	    else fillTH1(photonPtVCent_TRUTHNoRecoMatch_p[centPos][barrelEC], truthPhotonPt, fullWeight);
	  }
	}//End response filling	

	//Now we populate the histograms
	//All photons passing requirements must be included
	//If no reco match exists just fill
	for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
	  //Passes basic fiducial cuts
	  bool isGoodReco = photonEtaIsGood(photon_eta_p->at(pI));
	  isGoodReco = isGoodReco && photon_pt_p->at(pI) >= gammaPtBinsLowReco;
	  isGoodReco = isGoodReco && photon_pt_p->at(pI) < gammaPtBinsHighReco;
	  
	  //Signal is always defined as tightid + isolation so hardcode
	  bool isSignal = photon_tight_p->at(pI);
	  isSignal = isSignal && isIsolatedPhoton(isPP, doPtIsoCorrection, photon_correctedIso_p->at(pI));
	  //Sideband is not orthogonal to isGoodRecoSignal and definition can change - use function
	  bool isSideband = isSidebandPhoton(isPP, doPtIsoCorrection, sidebandType, photon_tight_p->at(pI), photon_correctedIso_p->at(pI));
	  
	  //goodreco signal
	  bool isGoodRecoSignal = isGoodReco && isSignal;
	  //goodreco sideband
	  bool isGoodRecoSideband = isGoodReco && isSideband;
	  	  
	  //Check if this reco is the truth match
	  bool isTruthPhotonMatched = isMC && pI == truthPhoRecoPos && truthPhoHasGoodReco;
	  if(isGoodReco){
	    std::vector<int> barrelECFill = {2}; //Always fill at the inclusive position
	    //Is this barrel or ec for the other fill
	    
	    if(photonEtaIsBarrel(photon_eta_p->at(pI))) barrelECFill.push_back(0);
	    else if(photonEtaIsEC(photon_eta_p->at(pI))) barrelECFill.push_back(1);
	    else{
	      std::cout << "Photon w/ eta=\'" << photon_eta_p->at(pI) << "\' is not registering as barrel or endcap check for bugs" << std::endl;
	      return 1;
	    }

	    if(systI == 0){	    
	      for(auto const barrelEC : barrelECFill){
		if(isGoodRecoSignal){
		  fillTH1(photonPtVCent_RAW_p[centPos][barrelEC], photon_pt_p->at(pI), fullWeight);
		  fillTH1(photonPtVCent_ValXWeightSum_p[centPos][barrelEC], photon_pt_p->at(pI), fullWeight*photon_pt_p->at(pI));
		  
		  if(isMC){
		    if(isTruthPhotonMatched) fillTH1(photonPtVCent_RAWWithTruthMatch_p[centPos][barrelEC], photon_pt_p->at(pI), fullWeight);
		    else fillTH1(photonPtVCent_RAWNoTruthMatch_p[centPos][barrelEC], photon_pt_p->at(pI), fullWeight);
		  }
		}
		else if(isGoodRecoSideband) fillTH1(photonPtVCent_RAWSideband_p[centPos][barrelEC], photon_pt_p->at(pI), fullWeight);
	      }
	    }

	    //We will keep collections of jets passing cuts for multijet observable construction
	    std::vector<TLorentzVector> goodRecoJets;
	    std::vector<int> goodRecoJetsPos, goodRecoJetsTruthPos;
	    
	    //Now produce photon+jet observables
	    for(unsigned int jI = 0; jI < aktRhi_insitu_jet_pt_p->size(); ++jI){
	      Float_t jtPtToUse = aktRhi_insitu_jet_pt_p->at(jI);
	      Float_t jtPhiToUse = aktRhi_insitu_jet_phi_p->at(jI);
	      Float_t jtEtaToUse = aktRhi_insitu_jet_eta_p->at(jI);	  
	      
	      bool isGoodTruthJet = false;
	      int truthPos = -1;

	      if(isMC){
		//MC corrections only apply to PYTHIA jets inserted in overlay
		if(aktRhi_truthpos_p->at(jI) >= 0){
		  jtPtToUse = aktRhi_etajes_jet_pt_p->at(jI);
		  jtPhiToUse = aktRhi_etajes_jet_phi_p->at(jI);
		  jtEtaToUse = aktRhi_etajes_jet_eta_p->at(jI);	    
		  
		  truthPos = aktRhi_truthpos_p->at(jI);
		  
		  if(aktR_truth_jet_pt_p->at(truthPos) < jtPtBinsLow || aktR_truth_jet_pt_p->at(truthPos) > jtPtBinsHigh) truthPos = -1;
		  else isGoodTruthJet = true;
		  
		//Add eta cut
		  if(isGoodTruthJet){
		    if(aktR_truth_jet_eta_p->at(truthPos) < jtEtaBinsLow || aktR_truth_jet_eta_p->at(truthPos) > jtEtaBinsHigh){
		      truthPos = -1;
		      isGoodTruthJet = false;
		    }
		  }
		}
		
		
		if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		//Add dR Cut
		if(isGoodTruthJet){
		  Float_t dRTruthGammaJet = getDR(aktR_truth_jet_eta_p->at(truthPos), aktR_truth_jet_phi_p->at(truthPos), truthPhotonEta, truthPhotonPhi);
		  
		  if(dRTruthGammaJet < gammaExclusionDR){
		    isGoodTruthJet = false;
		    truthPos = -1;
		  }
		}
	      }
	      
	      //Good reco jet requirements - pass eta cuts, pass pt cuts, pass dr cut and (later) pass dphi w/ gamma cut
	      bool isGoodRecoJet = jtPtToUse >= jtPtBinsLowReco && jtPtToUse < jtPtBinsHighReco;
	      if(isStrSame(systStrVect[systI], "JTPTCUT")) isGoodRecoJet = jtPtToUse >= jtPtBinsLowRecoSyst  && jtPtToUse < jtPtBinsHighReco;
	      isGoodRecoJet = isGoodRecoJet && jtEtaToUse >= jtEtaBinsLow && jtEtaToUse <= jtEtaBinsHigh;
	      Float_t dPhiRecoGammaJet = TMath::Abs(getDPHI(jtPhiToUse, photon_phi_p->at(pI)));
	      Float_t dRRecoGammaJet = getDR(jtEtaToUse, jtPhiToUse, photon_eta_p->at(pI), photon_phi_p->at(pI));
	      
	      //First do dphi - just copy the shit below
	      //Exclude any jet around the photon candidate
	      if(dRRecoGammaJet < gammaExclusionDR) isGoodRecoJet = false;
	      
	      if(isGoodRecoJet){	    
		for(auto const barrelEC : barrelECFill){
		  if(isGoodRecoSignal){
		    photonPtJtDPhiVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(dPhiRecoGammaJet, photon_pt_p->at(pI), fullWeight);
		    if(isMC){
		      if(isTruthPhotonMatched && isGoodTruthJet){
			photonPtJtDPhiVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(dPhiRecoGammaJet, photon_pt_p->at(pI), fullWeight);
			photonPtJtDPhiVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap(truthPos);
		      }
		      else{
			photonPtJtDPhiVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(dPhiRecoGammaJet, photon_pt_p->at(pI), fullWeight);		  
		      }		
		    }//end if(isGoodRecoSignal)
		  }
		  else if(isGoodRecoSideband){
		    photonPtJtDPhiVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(dPhiRecoGammaJet, photon_pt_p->at(pI), fullWeight);
		  }
		}//End reco barrelECFill
	      }//end if(isGoodRecoJet)
	      
	      if(isMC && isGoodTruthJet){
		//DPhi truth Cut only
		Float_t dPhiTruthGammaJet = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(truthPos), truthPhotonPhi));
		
		if(dPhiTruthGammaJet < gammaJtDPhiCut){
		  isGoodTruthJet = false;
		  truthPos = -1;
		}	  
	      }
	      
	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
	      //Now construct non-dphi based observables by enforcing the dPhi gamma-jet cut
	      if(dPhiRecoGammaJet < gammaJtDPhiCut) isGoodRecoJet = false;
	      
	      if(isGoodRecoJet){
		Float_t xJValue = jtPtToUse / photon_pt_p->at(pI);
		Bool_t xJValueGood = xJValue >= xjBinsLowReco && xJValue < xjBinsHighReco;
		Float_t xJValueTruth = -999;
		Bool_t xJValueTruthGood = false;
	      
		if(isGoodTruthJet && isTruthPhotonMatched){
		  xJValueTruth = aktR_truth_jet_pt_p->at(truthPos) / truthPhotonPt;
		  xJValueTruthGood = xJValueTruth > xjBinsLow && xJValueTruth < xjBinsHigh;
		}
		
		TLorentzVector tL;
		tL.SetPtEtaPhiM(jtPtToUse, jtEtaToUse, jtPhiToUse, 0.0);
		goodRecoJets.push_back(tL);
		goodRecoJetsPos.push_back(jI);
		goodRecoJetsTruthPos.push_back(truthPos);
		
		if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	      
		
		for(auto const barrelEC : barrelECFill){
		  if(isGoodRecoSignal){
		    photonPtJtPtVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(jtPtToUse, photon_pt_p->at(pI), fullWeight);
		    if(xJValueGood) photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(xJValue, photon_pt_p->at(pI), fullWeight);
		    
		    if(isMC){
		      if(isTruthPhotonMatched && isGoodTruthJet){
			photonPtJtPtVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(jtPtToUse, photon_pt_p->at(pI), fullWeight);
			photonPtJtPtVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap(truthPos);
			
			if(xJValueGood){
			  if(xJValueTruthGood){
			    photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(xJValue, photon_pt_p->at(pI), fullWeight);
			    photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap(truthPos);
			  }
			  else photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(xJValue, photon_pt_p->at(pI), fullWeight);
			}		  
		      }
		      else{
			photonPtJtPtVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(jtPtToUse, photon_pt_p->at(pI), fullWeight);
			if(xJValueGood) photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(xJValue, photon_pt_p->at(pI), fullWeight);
		      }
		    }
		  }
		  else if(isGoodRecoSideband){
		    photonPtJtPtVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(jtPtToUse, photon_pt_p->at(pI), fullWeight);
		    if(xJValueGood) photonPtJtXJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(xJValue, photon_pt_p->at(pI), fullWeight);
		  }
		}//End barrelECFill loop for reco
	      }
	    }//End reco jet loop
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	  
	    //Now we start another loop but one for multijet events
	    for(unsigned int gI = 0; gI < goodRecoJets.size(); ++gI){
	      //Get first  reco jet + corresponding truth info
	      TLorentzVector goodRecoJet1 = goodRecoJets[gI];
	      int truthPos1 = goodRecoJetsTruthPos[gI];
	      bool isTruthMatched1 = truthPos1 >= 0;	  
	    
	      for(unsigned int gI2 = gI+1; gI2 < goodRecoJets.size(); ++gI2){
		//Get second reco jet + corresponding truth info
		TLorentzVector goodRecoJet2 = goodRecoJets[gI2];
		int truthPos2 = goodRecoJetsTruthPos[gI2];
		bool isTruthMatched2 = truthPos2 >= 0;

		//combine truth info
		bool isTruthMatched = isTruthMatched1 && isTruthMatched2;       	
		int truthID1 = truthPos1*1000 + truthPos2;
		int truthID2 = truthPos2*1000 + truthPos1;

		//Construct variables before 4-vector sum of jets
		Float_t aJJValue = TMath::Abs(goodRecoJet1.Pt() - goodRecoJet2.Pt())/photon_pt_p->at(pI);
		Float_t dPhiJJValue = TMath::Abs(getDPHI(goodRecoJet1.Phi(), goodRecoJet2.Phi()));
		Float_t dRJJValue = getDR(goodRecoJet1.Eta(), goodRecoJet1.Phi(), goodRecoJet2.Eta(), goodRecoJet2.Phi());

		//Construct Booleans from existing variables
		Bool_t aJJValueGood = aJJValue >= ajBinsLowReco && aJJValue < ajBinsHighReco;
		Bool_t dRJJPasses = dRJJValue >= mixJetExclusionDR;	     	      

		//Use the binflattener for mixmachines before 4-vector sum of jets
		Float_t subJtGammaPtValReco = -999.0;
		if(isGoodReco) subJtGammaPtValReco = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(photon_pt_p->at(pI), goodRecoJet2.Pt(), __LINE__);
		Float_t subJtGammaPtValTruth = -999.0;

		//4-vector sum of the two jets
		goodRecoJet2 += goodRecoJet1;
		
		if(!dRJJPasses) continue;
		
		//Does the multijet pass multijtdphi cut at truth level?	    
		Float_t multiJtDPhiTruth = -999.0;
		Float_t aJJValueTruth = -999.0;
		Float_t xJJValueTruth = -999.0;
		
		bool isTruthMatchedDPhi = isTruthMatched;
		Bool_t xJJValueTruthGood = false;
		
		if(isTruthMatchedDPhi){
		  TLorentzVector goodTruthJet1;
		  goodTruthJet1.SetPtEtaPhiM(aktR_truth_jet_pt_p->at(truthPos1), aktR_truth_jet_eta_p->at(truthPos1), aktR_truth_jet_phi_p->at(truthPos1), 0.0);
		  
		  TLorentzVector goodTruthJet2;
		  goodTruthJet2.SetPtEtaPhiM(aktR_truth_jet_pt_p->at(truthPos2), aktR_truth_jet_eta_p->at(truthPos2), aktR_truth_jet_phi_p->at(truthPos2), 0.0);
		  
		  aJJValueTruth = TMath::Abs(goodTruthJet1.Pt() - goodTruthJet2.Pt())/truthPhotonPt;

		  //Swapping effects may seriously matter here
		  subJtGammaPtValTruth = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(truthPhotonPt, goodTruthJet2.Pt(), __LINE__);
		  
		  goodTruthJet2 += goodTruthJet1;
		  
		  xJJValueTruth = goodTruthJet2.Pt()/truthPhotonPt;
		  xJJValueTruthGood = xJJValueTruth >= xjjBinsLow && xJJValueTruth < xjjBinsHigh;
		  
		  multiJtDPhiTruth = TMath::Abs(getDPHI(truthPhotonPhi, goodTruthJet2.Phi()));
		  if(multiJtDPhiTruth < gammaMultiJtDPhiCut) isTruthMatchedDPhi = false;	     
		}
		
		Float_t multiJtDPhiReco = TMath::Abs(getDPHI(photon_phi_p->at(pI), goodRecoJet2.Phi()));
		Float_t xJJValue = goodRecoJet2.Pt() / photon_pt_p->at(pI);
		Bool_t xJJValueGood = xJJValue >= xjjBinsLowReco && xJJValue < xjjBinsHighReco;	       
		
		for(auto const barrelEC : barrelECFill){	    
		  if(isGoodRecoSignal){
		    photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(multiJtDPhiReco, subJtGammaPtValReco, fullWeight);
		    
		    if(isMC){
		      if(isTruthPhotonMatched && isTruthMatched){
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(multiJtDPhiReco, subJtGammaPtValReco, fullWeight);
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap({truthID1, truthID2});
		      }
		      else photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(multiJtDPhiReco, subJtGammaPtValReco, fullWeight);
		    }
		  }
		  else if(isGoodRecoSideband) photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(multiJtDPhiReco, subJtGammaPtValReco, fullWeight);
		  
		  //If it fails the reco cut we do not fill
		  if(multiJtDPhiReco < gammaMultiJtDPhiCut) continue;	    

		  
		  if(isGoodRecoSignal){
		    if(xJJValueGood) photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(xJJValue, subJtGammaPtValReco, fullWeight);
		    if(aJJValueGood) photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(aJJValue, subJtGammaPtValReco, fullWeight);
		  
		    photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(dPhiJJValue, subJtGammaPtValReco, fullWeight);
		    photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRaw(dRJJValue, subJtGammaPtValReco, fullWeight);
		    
		    if(isMC){
		      //xJJ Handling
		      if(xJJValueGood){
			bool goodTruth = isTruthPhotonMatched && isTruthMatchedDPhi && xJJValueTruthGood;
			if(goodTruth){
			  photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(xJJValue, subJtGammaPtValReco, fullWeight);		  		      
			  photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap({truthID1, truthID2}); 
			}
			else{
			  photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(xJJValue, subJtGammaPtValReco, fullWeight);		  		      
			}
		      }
		      
		      //aJJ Handling
		      if(aJJValueGood){
			bool goodTruth = isTruthPhotonMatched && isTruthMatchedDPhi;
			if(goodTruth){
			  photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(aJJValue, subJtGammaPtValReco, fullWeight);		  		      
			  photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap({truthID1, truthID2}); 
			}
			else{
			  photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(aJJValue, subJtGammaPtValReco, fullWeight);		  		      
			}
		      }
		      
		      //dphiJJ Handling
		      {
			bool goodTruth = isTruthPhotonMatched && isTruthMatchedDPhi;
			if(goodTruth){
			  photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(dPhiJJValue, subJtGammaPtValReco, fullWeight);		  		      
			  photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap({truthID1, truthID2}); 
			}
			else{
			  photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(dPhiJJValue, subJtGammaPtValReco, fullWeight);		  		      
			}
		      }
		      
		      //drJJ Handling
		      {
			bool goodTruth = isTruthPhotonMatched && isTruthMatchedDPhi;
			if(goodTruth){
			  photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawWithTruthMatch(dRJJValue, subJtGammaPtValReco, fullWeight);		  		      
			  photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->PushTrackingMap({truthID1, truthID2}); 
			}
			else{
			  photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYRawNoTruthMatch(dRJJValue, subJtGammaPtValReco, fullWeight);		  		      
			}
		      }
		    }
		  }//isGoodRecoSignal if statement ends
		  else if(isGoodRecoSideband){
		    if(xJJValueGood) photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(xJJValue, subJtGammaPtValReco, fullWeight);
		    if(aJJValueGood) photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(aJJValue, subJtGammaPtValReco, fullWeight);
		    
		    photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(dPhiJJValue, subJtGammaPtValReco, fullWeight);
		    photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYRaw(dRJJValue, subJtGammaPtValReco, fullWeight);
		  }
		}//End barrelECFill	      
	      }//for(goodRecoJets2)
	    }//End multijet loop, for(goodRecoJets)
	  
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
	    
	    //If doMix, begin running the mixing
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
	      
	      //Create the key to grab the mixed event jets
	      unsigned long long key = keyBoy.GetKey(eventKeyVect);
	      unsigned long long maxPos = mixingMap[key].size();
	      if(maxPos == 0){
		std::cout << "WHOOPS NO AVAILABLE MIXED EVENT. bailing" << std::endl;
		std::cout << key << ", " << mixCentPos << ", " << cent << std::endl;
		return 1;
	      }
	      unsigned long long nCurrentMixEvents = 0;
	      
	      std::vector<unsigned long long> jetPos1s, jetPos2s;
	      
	      while(nCurrentMixEvents < nMixEvents){
		unsigned long long jetPos = maxPos;	 
		bool goodJetPos = false;
		while(!goodJetPos){
		  jetPos = randGen_p->Uniform(0, maxPos-1);
		  goodJetPos = jetPos != maxPos && !vectContainsULL(jetPos, &jetPos1s);
		}
		++(signalMapCounterPost[key]);
		std::vector<TLorentzVector> jets = mixingMap[key][jetPos];
		
		unsigned long long jetPos2 = maxPos;
		bool goodJetPos2 = false;
		while(!goodJetPos2){
		  jetPos2 = randGen_p->Uniform(0, maxPos-1);
		  goodJetPos2 = jetPos2 != maxPos && jetPos2 != jetPos;
		  
		  if(goodJetPos2){
		    //Make sure we haven't used this combo yet -- only need to test swapped case as jetPos != jetPos1s[jpI] by earlier test !vectContainsULL
		    for(unsigned int jpI = 0; jpI < jetPos1s.size(); ++jpI){
		      if(jetPos2s[jpI] == jetPos && jetPos1s[jpI] == jetPos2){
			goodJetPos2 = false;
			break;
		      }
		    }
		  }//End if(goodJetPos2){
		}//End while(!goodJetPos2){
		
		jetPos1s.push_back(jetPos);
		jetPos2s.push_back(jetPos2);
		
		std::vector<TLorentzVector> jets2 = mixingMap[key][jetPos2];
		
		//Go thru and select jets passing cuts from first mixed event
		std::vector<TLorentzVector> passingJets1, passingJets2;
		for(unsigned int jI = 0; jI < jets.size(); ++jI){
		  Float_t dRRecoGammaJet = getDR(jets[jI].Eta(), jets[jI].Phi(), photon_eta_p->at(pI), photon_phi_p->at(pI));
		  Float_t dPhiRecoGammaJet = TMath::Abs(getDPHI(jets[jI].Phi(), photon_phi_p->at(pI)));
		  bool isGoodRecoJet = jets[jI].Pt() >= jtPtBinsLowReco && jets[jI].Pt() < jtPtBinsHighReco;
		  if(isStrSame(systStrVect[systI], "JTPTCUT")) isGoodRecoJet = jets[jI].Pt() >= jtPtBinsLowRecoSyst && jets[jI].Pt() < jtPtBinsHighReco;

		  isGoodRecoJet = isGoodRecoJet && jets[jI].Eta() >= jtEtaBinsLow && jets[jI].Eta() <= jtEtaBinsHigh;
		  isGoodRecoJet = isGoodRecoJet && dRRecoGammaJet >= gammaExclusionDR;
		  
		  if(isGoodRecoJet){
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			photonPtJtDPhiVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(dPhiRecoGammaJet, photon_pt_p->at(pI), mixWeight);
		      }
		      else if(isGoodRecoSideband){
			photonPtJtDPhiVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(dPhiRecoGammaJet, photon_pt_p->at(pI), mixWeight);
		      }
		    }//End for(auto const barrelEC...
		  }//End if(isGoodRecoJet)
		  
		  //Add in the dphi cut
		  isGoodRecoJet = isGoodRecoJet && dPhiRecoGammaJet >= gammaJtDPhiCut;
		  if(isGoodRecoJet){
		    //Since jet passes fill passingJets1
		    passingJets1.push_back(jets[jI]);
		    
		    Float_t xJValue = jets[jI].Pt() / photon_pt_p->at(pI);
		    Bool_t xJValueGood = xJValue >= xjBinsLowReco && xJValue < xjBinsHighReco;
		    for(auto const barrelEC : barrelECFill){		
		      if(isGoodRecoSignal){
			photonPtJtPtVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(jets[jI].Pt(), photon_pt_p->at(pI), mixWeight);
			if(xJValueGood) photonPtJtXJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(jets[jI].Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);		  
		      }
		      else if(isGoodRecoSideband){
			photonPtJtPtVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(jets[jI].Pt(), photon_pt_p->at(pI), mixWeight);
			if(xJValueGood) photonPtJtXJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(jets[jI].Pt()/photon_pt_p->at(pI), photon_pt_p->at(pI), mixWeight);		  
		      }
		    }//End for(auto const barrelEC : barrelECFill){
		  }//End if(isGoodRecoJet){
		}//End for(unsigned int jI = 0; jI < jets.size...	    
	      		
		//For multijet mixing we must process a second event
		for(unsigned int jI = 0; jI < jets2.size(); ++jI){
		  Float_t dRRecoGammaJet = getDR(jets2[jI].Eta(), jets2[jI].Phi(), photon_eta_p->at(pI), photon_phi_p->at(pI));
		  Float_t dPhiRecoGammaJet = TMath::Abs(getDPHI(jets2[jI].Phi(), photon_phi_p->at(pI)));
		  bool isGoodRecoJet = jets2[jI].Pt() >= jtPtBinsLowReco && jets2[jI].Pt() < jtPtBinsHighReco;
		  if(isStrSame(systStrVect[systI], "JTPTCUT")) isGoodRecoJet = jets2[jI].Pt() >= jtPtBinsLowRecoSyst && jets2[jI].Pt() < jtPtBinsHighReco;

		  isGoodRecoJet = isGoodRecoJet && jets2[jI].Eta() >= jtEtaBinsLow && jets2[jI].Eta() <= jtEtaBinsHigh;
		  isGoodRecoJet = isGoodRecoJet && dRRecoGammaJet >= gammaExclusionDR;
		  //Add in the dphi cut
		  isGoodRecoJet = isGoodRecoJet && dPhiRecoGammaJet >= gammaJtDPhiCut;
		  
		  if(isGoodRecoJet) passingJets2.push_back(jets2[jI]);
		}//End for(unsigned int jI = 0; jI < jets2.size...	    
		
		//We have 2 valid jet collections now for this photon - do multijet mixing
		//First pure background, single mixed event w/ itself
		for(unsigned int jI = 0; jI < passingJets1.size(); ++jI){
		  TLorentzVector jet1 = passingJets1[jI];
		  
		  for(unsigned int jI2 = jI+1; jI2 < passingJets1.size(); ++jI2){
		    TLorentzVector jet2 = passingJets1[jI2];
		    
		    //enforce dR exclusion region
		    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		    if(dR < mixJetExclusionDR) continue;
		    
		    Float_t aJJValue = TMath::Abs(jet1.Pt() - jet2.Pt()) / photon_pt_p->at(pI);
		    Bool_t aJJValueGood = aJJValue >= ajBinsLowReco && aJJValue < ajBinsHighReco;
		    
		    Float_t dPhiJJValue = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
		    Float_t dRJJValue = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		    
		    Float_t subJtGammaPtVal = -999.0;
		    if(isGoodReco) subJtGammaPtVal = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(photon_pt_p->at(pI), jet2.Pt(), __LINE__);

		    jet2 += jet1;
		    //enforce multijetdphi cut
		    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(pI)));
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		    }
		    
		    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
		    
		    Float_t xJJValue = jet2.Pt() / photon_pt_p->at(pI);
		    Bool_t xJJValueGood = xJJValue >= xjjBinsLowReco && xJJValue < xjjBinsHighReco;
		    
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(xJJValue, subJtGammaPtVal, mixWeight);		  
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(aJJValue, subJtGammaPtVal, mixWeight);
			
			photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(xJJValue, subJtGammaPtVal, mixWeight);
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(aJJValue, subJtGammaPtVal, mixWeight);
		      
			photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		    }//end for(barrelECFill){
		  }//end for(jI2 < passingJets2.size()){
		}//end for(jI < passingJets1.size()){

		if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
		
		//Now do mixed event crossed w/ current event i.e. one real one fake jet	    
		for(unsigned int jI = 0; jI < goodRecoJets.size(); ++jI){
		  TLorentzVector signalJet = goodRecoJets[jI];
		
		  for(unsigned int jI2 = 0; jI2 < passingJets1.size(); ++jI2){
		    TLorentzVector mixJet = passingJets1[jI2];
		  
		    //enforce dR exclusion region
		    Float_t dR = getDR(signalJet.Eta(), signalJet.Phi(), mixJet.Eta(), mixJet.Phi());
		    if(dR < mixJetExclusionDR) continue;
		    
		    Float_t aJJValue = TMath::Abs(signalJet.Pt() - mixJet.Pt()) / photon_pt_p->at(pI);
		    Bool_t aJJValueGood = aJJValue >= ajBinsLowReco && aJJValue < ajBinsHighReco;
		    
		    Float_t dPhiJJValue = TMath::Abs(getDPHI(signalJet.Phi(), mixJet.Phi()));
		    Float_t dRJJValue = getDR(signalJet.Eta(), signalJet.Phi(), mixJet.Eta(), mixJet.Phi());

		    Float_t subJtGammaPtVal = -999.0;
		    if(isGoodReco) subJtGammaPtVal = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(photon_pt_p->at(pI), mixJet.Pt(), __LINE__);
		    
		    mixJet += signalJet;
		    //enforce multijetdphi cut
		    Float_t multiJtDPhi = TMath::Abs(getDPHI(mixJet.Phi(), photon_phi_p->at(pI)));

		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		    }		
		    
		    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
		    
		    Float_t xJJValue = mixJet.Pt() / photon_pt_p->at(pI);
		    Bool_t xJJValueGood = xJJValue >= xjjBinsLowReco && xJJValue < xjjBinsHighReco;
		    
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(xJJValue, subJtGammaPtVal, mixWeight);		  
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(aJJValue, subJtGammaPtVal, mixWeight);
			
			photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMix(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(xJJValue, subJtGammaPtVal, mixWeight);
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(aJJValue, subJtGammaPtVal, mixWeight);
			
			photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMix(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		    }//end for(barrelECFill){
		  }//end for(jI2<passingJets1.size){
		}//end for(jI<goodRecoJets.size()){
		
		//Finally we need to correct the mixed event for instances where we took a fake jet from the signal event and mixed it with a fake jet from the mixed event, by using 2 mixed events
		for(unsigned int jI = 0; jI < passingJets1.size(); ++jI){
		  TLorentzVector jet1 = passingJets1[jI];
		  
		  for(unsigned int jI2 = 0; jI2 < passingJets2.size(); ++jI2){
		    TLorentzVector jet2 = passingJets2[jI2];	
	    
		    //enforce dR exclusion region
		    Float_t dR = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		    if(dR < mixJetExclusionDR) continue;
		    
		    Float_t aJJValue = TMath::Abs(jet1.Pt() - jet2.Pt()) / photon_pt_p->at(pI);
		    Bool_t aJJValueGood = aJJValue >= ajBinsLowReco && aJJValue < ajBinsHighReco;
		    
		    Float_t dPhiJJValue = TMath::Abs(getDPHI(jet1.Phi(), jet2.Phi()));
		    Float_t dRJJValue = getDR(jet1.Eta(), jet1.Phi(), jet2.Eta(), jet2.Phi());
		    
		    Float_t subJtGammaPtVal = -999.0;
		    if(isGoodReco) subJtGammaPtVal = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(photon_pt_p->at(pI), jet2.Pt(), __LINE__);

		    jet2 += jet1;
		    //enforce multijetdphi cut
		    Float_t multiJtDPhi = TMath::Abs(getDPHI(jet2.Phi(), photon_phi_p->at(pI)));		    
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(multiJtDPhi, subJtGammaPtVal, mixWeight);
		      }
		    }//end for(barrelECFill){
		    
		    if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
		    
		    Float_t xJJValue = jet2.Pt() / photon_pt_p->at(pI);
		    Bool_t xJJValueGood = xJJValue >= xjjBinsLowReco && xJJValue < xjjBinsHighReco;
		    
		    for(auto const barrelEC : barrelECFill){
		      if(isGoodRecoSignal){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(xJJValue, subJtGammaPtVal, mixWeight);
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(aJJValue, subJtGammaPtVal, mixWeight);
			
			photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(multiJtDPhi, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_p[centPos][barrelEC][systI]->FillXYMixCorrection(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		      else if(isGoodRecoSideband){
			if(xJJValueGood) photonPtJtXJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(xJJValue, subJtGammaPtVal, mixWeight);
			if(aJJValueGood) photonPtJtAJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(aJJValue, subJtGammaPtVal, mixWeight);
			
			photonPtJtDPhiJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(dPhiJJValue, subJtGammaPtVal, mixWeight);
			photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(multiJtDPhi, subJtGammaPtVal, mixWeight);
			photonPtJtDRJJVCent_MixMachine_Sideband_p[centPos][barrelEC][systI]->FillXYMixCorrection(dRJJValue, subJtGammaPtVal, mixWeight);
		      }
		    }//end for(barrelECFill){
		  }
		}
	      
		++nCurrentMixEvents; // increment to next mixed event
	      }//End while(nCurrentMixEvents < nMixEvents){		  
	    }//End if(doMix){
	  }//End fills for pure photon good reco   
	}//End Photon loop
            
	if(isMC){
	  std::vector<TLorentzVector> goodTruthJets;
	  std::vector<int> goodTruthJetsPos, goodTruthJetsRecoPos;
	  
	  for(unsigned int jI = 0; jI < aktR_truth_jet_pt_p->size(); ++jI){	  
	    //Continue if the Truth jet is no good
	    Float_t dRTruthGammaJet = TMath::Abs(getDR(aktR_truth_jet_eta_p->at(jI), aktR_truth_jet_phi_p->at(jI), truthPhotonEta, truthPhotonPhi));
	    bool isGoodTruthJet = aktR_truth_jet_pt_p->at(jI) >= jtPtBinsLow && aktR_truth_jet_pt_p->at(jI) < jtPtBinsHigh;
	    isGoodTruthJet = isGoodTruthJet && aktR_truth_jet_eta_p->at(jI) >= jtEtaBinsLow && aktR_truth_jet_eta_p->at(jI) < jtEtaBinsHigh;
	    isGoodTruthJet = isGoodTruthJet && dRTruthGammaJet >= gammaExclusionDR;
	    Float_t dPhiTruthGammaJet = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(jI), truthPhotonPhi));
	    
	    if(!isGoodTruthJet) continue;  
	    
	    for(auto const barrelECTruth : barrelECFillTruth){
	      photonPtJtDPhiVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(dPhiTruthGammaJet, truthPhotonPt, fullWeight);	   	  
	      
	      if(photonPtJtDPhiVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(jI)) photonPtJtDPhiVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(dPhiTruthGammaJet, truthPhotonPt, fullWeight);	    
	      else photonPtJtDPhiVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(dPhiTruthGammaJet, truthPhotonPt, fullWeight);	    
	    }//end barrelECTruth for loop
	    
	    //Now add dphi cutting
	    if(isGoodTruthJet){
	      if(dPhiTruthGammaJet < gammaJtDPhiCut) isGoodTruthJet = false;
	    }
	    if(!isGoodTruthJet) continue;  
	    
	    Float_t xJValueTruth = -999.0;
	    Bool_t xJValueTruthGood = false;
	    
	    //Populate the good Truth jets vector
	    TLorentzVector tL;
	    tL.SetPtEtaPhiM(aktR_truth_jet_pt_p->at(jI), aktR_truth_jet_eta_p->at(jI), aktR_truth_jet_phi_p->at(jI), 0.0);
	    goodTruthJets.push_back(tL);
	    goodTruthJetsPos.push_back(jI);
	    goodTruthJetsRecoPos.push_back(aktR_truth_jet_recopos_p->at(jI));
	    
	    for(auto const barrelECTruth : barrelECFillTruth){
	      photonPtJtPtVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(aktR_truth_jet_pt_p->at(jI), truthPhotonPt, fullWeight);	   
	      if(photonPtJtPtVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(jI)) photonPtJtPtVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(aktR_truth_jet_pt_p->at(jI), truthPhotonPt, fullWeight);
	      else photonPtJtPtVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(aktR_truth_jet_pt_p->at(jI), truthPhotonPt, fullWeight);
	      
	      photonPtJtXJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(aktR_truth_jet_pt_p->at(jI)/truthPhotonPt, truthPhotonPt, fullWeight);	   
	      if(photonPtJtXJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(jI)) photonPtJtXJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(aktR_truth_jet_pt_p->at(jI)/truthPhotonPt, truthPhotonPt, fullWeight);	    
	      else photonPtJtXJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(aktR_truth_jet_pt_p->at(jI)/truthPhotonPt, truthPhotonPt, fullWeight);	    
	    }//end barrelECTruth for loop
	  }
	  
	  //Construct multijet truth w/o reco fills
	  for(unsigned int jI = 0; jI < goodTruthJets.size(); ++jI){
	    TLorentzVector goodTruthJet1 = goodTruthJets[jI];
	    int truthID1 = goodTruthJetsPos[jI];
	    
	    for(unsigned int jI2 = jI+1; jI2 < goodTruthJets.size(); ++jI2){
	      TLorentzVector goodTruthJet2 = goodTruthJets[jI2];
	      int truthID2 = goodTruthJetsPos[jI2];
	      
	      int truthCompID1 = 1000*truthID1 + truthID2;
	      int truthCompID2 = 1000*truthID2 + truthID1;
	      
	      Float_t aJJValue = TMath::Abs(goodTruthJet1.Pt() - goodTruthJet2.Pt()) / truthPhotonPt;
	      Bool_t aJJValueGood = aJJValue >= ajBinsLow && aJJValue < ajBinsHigh;
	      
	      Float_t dPhiJJValue = TMath::Abs(getDPHI(goodTruthJet1.Phi(), goodTruthJet2.Phi()));
	      Float_t dRJJValue = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodTruthJet2.Eta(), goodTruthJet2.Phi());
	      
	      //	      if(dRJJValue < mixJetExclusionDR) continue;
	      
	      
	      Float_t subLeadingJetPt = goodTruthJet1.Pt();
	      if(goodTruthJet2.Pt() < subLeadingJetPt) subLeadingJetPt = goodTruthJet2.Pt();
	      Float_t subJtGammaPtValTruth = subJtGammaPtBinFlattener.GetGlobalBinCenterFromBin12Val(truthPhotonPt, subLeadingJetPt, __LINE__);
	      
	      goodTruthJet2 += goodTruthJet1;
	      
	      Float_t multiJtDPhiValue = TMath::Abs(getDPHI(truthPhotonPhi, goodTruthJet2.Phi()));
	      Float_t xJJValue = goodTruthJet2.Pt() / truthPhotonPt;
	      Bool_t xJJValueGood = xJJValue >= xjjBinsLow && xJJValue < xjjBinsHigh;
	      
	      for(auto const barrelECTruth : barrelECFillTruth){
		photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(multiJtDPhiValue, subJtGammaPtValTruth, fullWeight);
		
		bool truthFillWithRecoPos = photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID1) || photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID2);
			      
		if(truthFillWithRecoPos) photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(multiJtDPhiValue, subJtGammaPtValTruth, fullWeight);
		else photonPtJtDPhiJJGVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(multiJtDPhiValue, subJtGammaPtValTruth, fullWeight);
	      }	  
	      
	      if(multiJtDPhiValue < gammaMultiJtDPhiCut) continue; //We dont fill truth that fails this cut
	      
	      for(auto const barrelECTruth : barrelECFillTruth){
	      
		bool truthFillWithRecoPos = photonPtJtXJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID1) || photonPtJtXJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID2);
		photonPtJtXJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(xJJValue, subJtGammaPtValTruth, fullWeight);
		if(truthFillWithRecoPos) photonPtJtXJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(xJJValue, subJtGammaPtValTruth, fullWeight);
		else{
		  photonPtJtXJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(xJJValue, subJtGammaPtValTruth, fullWeight);
		}
		
		truthFillWithRecoPos = photonPtJtAJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID1) || photonPtJtAJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID2);
		photonPtJtAJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(aJJValue, subJtGammaPtValTruth, fullWeight);
		if(truthFillWithRecoPos) photonPtJtAJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(aJJValue, subJtGammaPtValTruth, fullWeight);
		else photonPtJtAJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(aJJValue, subJtGammaPtValTruth, fullWeight);
	      
		
		truthFillWithRecoPos = photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID1) || photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID2);
		photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(dPhiJJValue, subJtGammaPtValTruth, fullWeight);
		if(truthFillWithRecoPos) photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(dPhiJJValue, subJtGammaPtValTruth, fullWeight);
		else photonPtJtDPhiJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(dPhiJJValue, subJtGammaPtValTruth, fullWeight);
		
		
		truthFillWithRecoPos = photonPtJtDRJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID1) || photonPtJtDRJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->IsInTrackingMap(truthCompID2);
		photonPtJtDRJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruth(dRJJValue, subJtGammaPtValTruth, fullWeight);
		if(truthFillWithRecoPos) photonPtJtDRJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthWithRecoMatch(dRJJValue, subJtGammaPtValTruth, fullWeight);
		else photonPtJtDRJJVCent_MixMachine_p[centPos][barrelECTruth][systI]->FillXYTruthNoRecoMatch(dRJJValue, subJtGammaPtValTruth, fullWeight);
	      }	  
	    }
	  }      
	}//end if(isMC)
      }//End systStrVect loop

      //Clean the tracking maps of all mixmachines
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	    photonPtJtPtVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtXJVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtXJJVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtAJJVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	    photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI]->CleanTrackingMap();
	  }
	}
      }
    }
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

    if(tempRatioVals.size() != 0) std::cout << "MEDIAN RATIO MIX-PER-SIGNAL: " << tempRatioVals[tempRatioVals.size()/2] << std::endl;
    else std::cout << "NO MEDIAN RATIO MIX-PER-SIGNAL." << std::endl;
  }
  
  std::cout << std::endl;
  std::cout << "Number of events w/ N photons" << std::endl;
  for(auto const & val: perEventPhotonCounts){
    std::cout << " N=" << val.first << ": " << val.second << std::endl;
  }
  std::cout << std::endl;

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "FULL WEIGHTS USED: " << std::endl;
  for(unsigned int i = 0; i < fullWeightVals.size(); ++i){
    std::cout << " " << i << ": " << fullWeightVals[i] << std::endl;
  }
  std::cout << std::endl;

  //  doGlobalDebug = false;
  if(doGlobalDebug) std::cout << "DOGLOBALDEBUG MANUALLY TURNED ON AT LINE: " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  std::cout << inFile_p << std::endl;
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  
  if(isMC && keepResponseTree){
    unfoldTree_p->Write("", TObject::kOverwrite);
    delete unfoldTree_p;
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

      //      photonPtVCent_RAW_p[cI][barrelAndECComboPos]->Add(photonPtVCent_RAW_p[cI][bI]);
      //      photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->Add(photonPtVCent_RAWSideband_p[cI][bI]);
      
      bool mixMachineAddsAllGood = true;
      /*
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtPtVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtXJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtXJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtAJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_MixMachine_p[cI][bI]);

      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDRJJVCent_MixMachine_p[cI][barrelAndECComboPos]->Add(photonPtJtDRJJVCent_MixMachine_p[cI][bI]);


      //Sidebands
      //JtPt
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtPtVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_MixMachine_Sideband_p[cI][bI]);

      //JtXJ
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtXJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_MixMachine_Sideband_p[cI][bI]);

      //JtXJJ
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtXJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_MixMachine_Sideband_p[cI][bI]);

      //JtDPhi
     mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][bI]);

     //JtAJJ
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtAJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_MixMachine_Sideband_p[cI][bI]);

      //JtDPhiJJG
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][bI]);

      //JtDPhiJJ
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][bI]);

      //JtDRJJ
      mixMachineAddsAllGood = mixMachineAddsAllGood && photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos]->Add(photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][bI]);
      */
      std::cout << "MixMachineAddsAllGood: " << mixMachineAddsAllGood << std::endl;
      
      if(isMC){
	//	photonPtVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtVCent_TRUTH_p[cI][bI]);

	photonPtJtPtVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtPtVCent_TRUTH_p[cI][bI]);
	photonPtJtXJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtXJVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiVCent_TRUTH_p[cI][bI]);
	photonPtJtXJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtXJJVCent_TRUTH_p[cI][bI]);
	photonPtJtAJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtAJJVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiJJGVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJGVCent_TRUTH_p[cI][bI]);
	photonPtJtDPhiJJVCent_TRUTH_p[cI][barrelAndECComboPos]->Add(photonPtJtDPhiJJVCent_TRUTH_p[cI][bI]);


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

    //    std::cout << "WE ARE HERE " << __LINE__ << std::endl;
    //Compute the subtracted via the mix-machine
    //    /*
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	photonPtJtPtVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtXJVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtXJJVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtAJJVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI]->ComputeSub();            
	photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]->ComputeSub();       
	
	photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
	photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI]->ComputeSub();       
      }
    }

    if(!isMC){//Purity correction; Data only
      if(doUnifiedPurity){	
	//The unified purity cannot do barrel and endcap corrections separately 
	std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_Eta0p00to2p37";
	if(isStrSame(centBinsStr[cI], "Cent30to50")) purityFitStr = "fit_purity_Cent30to80_Eta0p00to2p37";
	else if(isStrSame(centBinsStr[cI], "Cent50to80")) purityFitStr = "fit_purity_Cent30to80_Eta0p00to2p37";
	
	TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());      
	std::vector<double> meanValPerGammaBin, rawPerGammaBin, sidebandPerGammaBin;
	
	for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
	  double numVal = photonPtVCent_ValXWeightSum_p[cI][barrelAndECComboPos]->GetBinContent(gI+1);
	  double denomVal = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinContent(gI+1);

	  if(numVal < TMath::Power(10,-100)) meanValPerGammaBin.push_back(0.0);
	  else meanValPerGammaBin.push_back(numVal/denomVal);

	  rawPerGammaBin.push_back(denomVal);

	  double sidebandVal = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetBinContent(gI+1);

	  if(sidebandVal < TMath::Power(10,-100)) sidebandPerGammaBin.push_back(0.0);
	  else sidebandPerGammaBin.push_back(sidebandVal);
	}
		  	  
	
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtXJVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtXJJVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtAJJVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJGVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJVCent_MixMachine_p[cI][barrelAndECComboPos][systI], photonPtJtDRJJVCent_MixMachine_p[cI][barrelAndECComboPos][systI]};

	  std::vector<mixMachine*> inMixMachines_Sideband_p = {photonPtJtPtVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtXJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtXJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtAJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI], photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][barrelAndECComboPos][systI]};
	    
	  std::vector<TH2D*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtXJVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtXJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtAJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI], photonPtJtDRJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]};	  
	  
	  //Need to know whether it is multijet or not + check the sizes make sense
	  std::vector<bool> isMultijet_p = {false, false, false, true, true, true, true, true};
	  if(isMultijet_p.size() != outHists_p.size()){
	    std::cout << "gdjNtupleToHist Error L" << __LINE__ << ": Mismatch between hist vector and bool isMultijet_p vector. fix. return 1" << std::endl;
	    return 1;
	  }
	  
	  for(unsigned int hI = 0; hI < inMixMachines_p.size(); ++hI){
	    TH2D* subHist_p = inMixMachines_p[hI]->GetTH2DPtr("SUB");
	    TH2D* subHist_Sideband_p = inMixMachines_Sideband_p[hI]->GetTH2DPtr("SUB");	      
	      
	    bool isMultijet = isMultijet_p[hI];

	    for(Int_t bIX = 0; bIX < subHist_p->GetXaxis()->GetNbins(); ++bIX){
	      for(Int_t bIY = 0; bIY < subHist_p->GetYaxis()->GetNbins(); ++bIY){
		Float_t binContent = subHist_p->GetBinContent(bIX+1, bIY+1);
		Float_t binError = subHist_p->GetBinError(bIX+1, bIY+1);
		
		Float_t binSidebandContent = subHist_Sideband_p->GetBinContent(bIX+1, bIY+1);
		Float_t binSidebandError = subHist_Sideband_p->GetBinError(bIX+1, bIY+1);

		//Continue on a zero bin
		if(binContent < TMath::Power(10,-50) && binSidebandContent < TMath::Power(10,-50)){
		  outHists_p[hI]->SetBinContent(bIX+1, bIY+1, 0.0);
		  outHists_p[hI]->SetBinError(bIX+1, bIY+1, 0.0);
		  
		  continue;
		}
	
		//Non-multijet, y-axis should just track gammaPtBins i.e.
		Int_t gammaPtBin = bIY;
		if(isMultijet){//Multijet need to find the appropriate gammaPtBin for purity
		  gammaPtBin = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		}
		
		Float_t purity = purityFit_p->Eval(meanValPerGammaBin[gammaPtBin]);
		
		binSidebandContent = (1.0 - purity)*binSidebandContent*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];
		binSidebandError = (1.0 - purity)*binSidebandError*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];
		
		binContent -= binSidebandContent;
		
		binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		
		outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
		outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
	      }
	    }	 	   
	  }
	}
	
	for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetXaxis()->GetNbins(); ++bIX){
	  Float_t purity = purityFit_p->Eval(meanValPerGammaBin[bIX]);
	  
	  Float_t binContent = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
	  Float_t binError = photonPtVCent_RAW_p[cI][barrelAndECComboPos]->GetBinError(bIX+1);
	  
	  Float_t binSidebandContent = photonPtVCent_RAWSideband_p[cI][barrelAndECComboPos]->GetBinContent(bIX+1);
	  
	  if(binContent < TMath::Power(10,-50) && binSidebandContent < TMath::Power(10,-50)){
	    photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinContent(bIX+1, 0.0);
	    photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinError(bIX+1, 0.0);
	    
	    continue;
	  }
	    
	  binSidebandContent = (1.0 - purity)*binSidebandContent*binContent/binSidebandContent;	  
	  binContent -= binSidebandContent;
	  
	  photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinContent(bIX+1, binContent);
	  photonPtVCent_PURCORR_p[cI][barrelAndECComboPos]->SetBinError(bIX+1, binError);
	}
      }
      else{	
	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  if(isStrSame(barrelAndECStr[eI], "BarrelAndEC")) continue;
	  
	  std::string purityFitStr = "fit_purity_" + centBinsStr[cI] + "_" + barrelAndECFitStr[eI];
	  TF1* purityFit_p = (TF1*)purityFile_p->Get(purityFitStr.c_str());
	  std::vector<double> meanValPerGammaBin, rawPerGammaBin, sidebandPerGammaBin;

	  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
	    double numVal = photonPtVCent_ValXWeightSum_p[cI][eI]->GetBinContent(gI+1);
	    double denomVal = photonPtVCent_RAW_p[cI][eI]->GetBinContent(gI+1);

	    if(numVal < TMath::Power(10,-100)) meanValPerGammaBin.push_back(0.0);
	    else meanValPerGammaBin.push_back(numVal/denomVal);

	    rawPerGammaBin.push_back(denomVal);

	    double sidebandVal = photonPtVCent_RAWSideband_p[cI][eI]->GetBinContent(gI+1);

	    if(sidebandVal < TMath::Power(10,-100)) sidebandPerGammaBin.push_back(0.0);
	    else sidebandPerGammaBin.push_back(sidebandVal);
	  }

	  //EDITING HERE 2022.09.19	  
	  for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	    std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][eI][systI], photonPtJtXJVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI], photonPtJtXJJVCent_MixMachine_p[cI][eI][systI], photonPtJtAJJVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI], photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]};
	    std::vector<mixMachine*> inMixMachines_Sideband_p = {photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI], photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI]};
	    std::vector<TH2D*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI][systI], photonPtJtXJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI], photonPtJtXJJVCent_PURCORR_p[cI][eI][systI], photonPtJtAJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]};	  
	    //Need to know whether it is multijet or not + check the sizes make sense
	    std::vector<bool> isMultijet_p = {false, false, false, true, true, true, true, true};
	    if(isMultijet_p.size() != outHists_p.size()){
	      std::cout << "gdjNtupleToHist Error L" << __LINE__ << ": Mismatch between hist vector and bool isMultijet_p vector. fix. return 1" << std::endl;
	      return 1;
	    }

	      
	    for(unsigned int hI = 0; hI < inMixMachines_p.size(); ++hI){
	      TH2D* subHist_p = inMixMachines_p[hI]->GetTH2DPtr("SUB");
	      TH2D* subHist_Sideband_p = inMixMachines_Sideband_p[hI]->GetTH2DPtr("SUB");	      
	      
	      bool isMultijet = isMultijet_p[hI];
	
	      for(Int_t bIX = 0; bIX < subHist_p->GetXaxis()->GetNbins(); ++bIX){
		for(Int_t bIY = 0; bIY < subHist_p->GetYaxis()->GetNbins(); ++bIY){
		  Float_t binContent = subHist_p->GetBinContent(bIX+1, bIY+1);
		  Float_t binError = subHist_p->GetBinError(bIX+1, bIY+1);
		  
		  Float_t binSidebandContent = subHist_Sideband_p->GetBinContent(bIX+1, bIY+1);
		  Float_t binSidebandError = subHist_Sideband_p->GetBinError(bIX+1, bIY+1);

		  //Continue on a zero bin
		  if(binContent < TMath::Power(10,-50) && binSidebandContent < TMath::Power(10,-50)){
		    outHists_p[hI]->SetBinContent(bIX+1, bIY+1, 0.0);
		    outHists_p[hI]->SetBinError(bIX+1, bIY+1, 0.0);
		    continue;
		  }

		  //Non-multijet, y-axis should just track gammaPtBins i.e.
		  Int_t gammaPtBin = bIY;
		  if(isMultijet){//Multijet need to find the appropriate gammaPtBin for purity
		    gammaPtBin = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		  }

		  Float_t purity = purityFit_p->Eval(meanValPerGammaBin[gammaPtBin]);		  
		  binSidebandContent = (1.0 - purity)*binSidebandContent*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];
		  binSidebandError = (1.0 - purity)*binSidebandError*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];
		  
		  binContent -= binSidebandContent;
		  
		  binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);
		  
		  outHists_p[hI]->SetBinContent(bIX+1, bIY+1, binContent);
		  outHists_p[hI]->SetBinError(bIX+1, bIY+1, binError);
		}
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
	    
	  for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	    photonPtJtPtVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtPtVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtXJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtXJVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtDPhiVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtXJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtXJJVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtAJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtAJJVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtDPhiJJGVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtDPhiJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI]);
	    photonPtJtDRJJVCent_PURCORR_p[cI][barrelAndECComboPos][systI]->Add(photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]);
	  }
	}
      }
    }
    else{
      //If its MC we just pass thru; Purity correction goes on data only
      std::cout << "PASS THRU" << std::endl;

      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){	
	
	//Do Photons, then follow up with mixmachines
	for(Int_t bIX = 0; bIX < photonPtVCent_RAW_p[cI][eI]->GetXaxis()->GetNbins(); ++bIX){
	  photonPtVCent_PURCORR_p[cI][eI]->SetBinContent(bIX+1, photonPtVCent_RAW_p[cI][eI]->GetBinContent(bIX+1));
	  photonPtVCent_PURCORR_p[cI][eI]->SetBinError(bIX+1, photonPtVCent_RAW_p[cI][eI]->GetBinError(bIX+1));
	}
	
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  std::vector<mixMachine*> inMixMachines_p = {photonPtJtPtVCent_MixMachine_p[cI][eI][systI], photonPtJtXJVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI], photonPtJtXJJVCent_MixMachine_p[cI][eI][systI], photonPtJtAJJVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI], photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI], photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]};
	  std::vector<TH2D*> outHists_p = {photonPtJtPtVCent_PURCORR_p[cI][eI][systI], photonPtJtXJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI], photonPtJtXJJVCent_PURCORR_p[cI][eI][systI], photonPtJtAJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI], photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI], photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]};
	  
	  for(unsigned int mI = 0; mI < inMixMachines_p.size(); ++mI){
	    TH2D* subHist_p = inMixMachines_p[mI]->GetTH2DPtr("SUB");	
	    for(unsigned int bIX = 0; bIX < outHists_p[mI]->GetXaxis()->GetNbins(); ++bIX){
	      for(unsigned int bIY = 0; bIY < outHists_p[mI]->GetYaxis()->GetNbins(); ++bIY){
		outHists_p[mI]->SetBinContent(bIX+1, bIY+1, subHist_p->GetBinContent(bIX+1, bIY+1));
		outHists_p[mI]->SetBinError(bIX+1, bIY+1, subHist_p->GetBinError(bIX+1, bIY+1));
	      }
	    }	  
	  }
	}
      }
    }
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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
  
    if(isMC){
      drJJ_OneJetNoTruth_p[cI]->Write("", TObject::kOverwrite);
    }
  
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      //Pure photon values
      photonPtVCent_RAW_p[cI][eI]->Write("", TObject::kOverwrite);
      if(isMC){
	photonPtVCent_RAWWithTruthMatch_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtVCent_RAWNoTruthMatch_p[cI][eI]->Write("", TObject::kOverwrite);
      }
      photonPtVCent_RAWSideband_p[cI][eI]->Write("", TObject::kOverwrite);      
      if(isMC){
	photonPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtVCent_TRUTHWithRecoMatch_p[cI][eI]->Write("", TObject::kOverwrite);
	photonPtVCent_TRUTHNoRecoMatch_p[cI][eI]->Write("", TObject::kOverwrite);
      }
      
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	//JtPt
	photonPtJtPtVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtPtVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtXJ
	photonPtJtXJVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtXJVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtDPhi
	photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtXJJ
	photonPtJtXJJVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtXJJVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtAJJ
	photonPtJtAJJVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtAJJVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtDPhiJJG
	photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtDPhiJJ
	photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtDRJJ
	photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI]->Clean();
	
	//JtPt
	photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtPtVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtXJ
	photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtXJVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtDPhi
	photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtXJJ
	photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtXJJVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtAJJ
	photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtAJJVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtDPhiJJG
	photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiJJGVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtDPhiJJ
	photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDPhiJJVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
	
	//JtDRJJ
	photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI]->WriteToDirectory(centDir_p);
	photonPtJtDRJJVCent_MixMachine_Sideband_p[cI][eI][systI]->Clean();
      }

      photonPtVCent_PURCORR_p[cI][eI]->Write("", TObject::kOverwrite);

      if(doMix){
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  photonPtJtPtVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	}

	if(isMC){
	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);

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
	  //	  photonPtJtDRJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);

	  
	  photonPtJtPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	}
      }
      else if(isPP){	
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  photonPtJtPtVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	  photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI]->Write("", TObject::kOverwrite);
	}

	if(isMC){
	  photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_TRUTHMATCHEDRECO_p[cI][eI]->Write("", TObject::kOverwrite);

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
	  // 	  photonPtJtDRJJVCent_PUREBKGD_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_MIXBKGD_p[cI][eI]->Write("", TObject::kOverwrite);

	  photonPtJtPtVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtXJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtAJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  photonPtJtDPhiJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	  //	  photonPtJtDRJJVCent_TRUTH_p[cI][eI]->Write("", TObject::kOverwrite);
	}
      }
    }    

    photonEtaPt_p[cI]->Write("", TObject::kOverwrite);

    for(Int_t pI = 0; pI < nGammaPtBins+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
    }
    centDir_p->Close();
    delete centDir_p;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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
      if(isMC){
	delete photonPtVCent_RAWWithTruthMatch_p[cI][eI];
	delete photonPtVCent_RAWNoTruthMatch_p[cI][eI];
      }
      delete photonPtVCent_RAWSideband_p[cI][eI];
      if(isMC){
	delete photonPtVCent_TRUTH_p[cI][eI];
	delete photonPtVCent_TRUTHWithRecoMatch_p[cI][eI];
	delete photonPtVCent_TRUTHNoRecoMatch_p[cI][eI];
      }
      delete photonPtVCent_PURCORR_p[cI][eI];
      
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	delete photonPtJtPtVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtXJVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtDPhiVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtXJJVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtAJJVCent_MixMachine_p[cI][eI][systI];
	delete photonPtJtDPhiJJGVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtDPhiJJVCent_MixMachine_p[cI][eI][systI];
        delete photonPtJtDRJJVCent_MixMachine_p[cI][eI][systI];
      }

      if(doMix){	
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  delete photonPtJtPtVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtXJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtXJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtAJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI];
	}

	if(isMC){
	  delete photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  //	  delete photonPtJtDRJJVCent_TRUTHMATCHEDRECO_p[cI][eI];

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
	  //	  delete photonPtJtDRJJVCent_PUREBKGD_p[cI][eI];
	  //	  delete photonPtJtDRJJVCent_MIXBKGD_p[cI][eI];

	  delete photonPtJtPtVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTH_p[cI][eI];
	  //	  delete photonPtJtDRJJVCent_TRUTH_p[cI][eI];
	}
      }
      else if(isPP){
	for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	  delete photonPtJtPtVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtXJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtXJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtAJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiJJGVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDPhiJJVCent_PURCORR_p[cI][eI][systI];
	  delete photonPtJtDRJJVCent_PURCORR_p[cI][eI][systI];
	}

	if(isMC){
	  delete photonPtJtPtVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTHMATCHEDRECO_p[cI][eI];
	  //	  delete photonPtJtDRJJVCent_TRUTHMATCHEDRECO_p[cI][eI];

	  delete photonPtJtPtVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiVCent_TRUTH_p[cI][eI];
	  delete photonPtJtXJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtAJJVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJGVCent_TRUTH_p[cI][eI];
	  delete photonPtJtDPhiJJVCent_TRUTH_p[cI][eI];
	  //	  delete photonPtJtDRJJVCent_TRUTH_p[cI][eI];
	}
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

  config_p->SetValue("NGAMMASYSANDNOM", nGammaSysAndNom);
  std::string gammaSysAndNom = "";
  for(Int_t sI = 0; sI < nGammaSysAndNom; ++sI){
    gammaSysAndNom = gammaSysAndNom + gammaSysAndNomStr[sI] + ",";
  }
  if(gammaSysAndNom.find(",") != std::string::npos) gammaSysAndNom.replace(gammaSysAndNom.size()-1, 1, "");
  config_p->SetValue("GAMMASYSANDNOM", gammaSysAndNom.c_str());

  config_p->SetValue("NJETSYSANDNOM", nJetSysAndNom);
  std::string jetSysAndNom = "";
  for(Int_t sI = 0; sI < nJetSysAndNom; ++sI){
    jetSysAndNom = jetSysAndNom + jetSysAndNomStr[sI] + ",";
  }
  if(jetSysAndNom.find(",") != std::string::npos) jetSysAndNom.replace(jetSysAndNom.size()-1, 1, "");
  config_p->SetValue("JETSYSANDNOM", jetSysAndNom.c_str());
  
  std::string systStrForConfig = vectToStrComma(systStrVect);
  std::string systTypeForConfig = vectToStrComma(systTypeVect);

  config_p->SetValue("SYSTNAMES", systStrForConfig.c_str());
  config_p->SetValue("SYSTTYPES", systTypeForConfig.c_str());

  config_p->SetValue("SUBJTGAMMAPTMIN", 0.0);
  config_p->SetValue("SUBJTGAMMAPTMAX", subJtGammaPtMax);
  
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
  delete truPhoPtPassing_p;  delete leadingPhoPtPassing_p;
  delete leadingPhoEtaPassing_p;
  delete leadingPhoIsoPassing_p;
  delete leadingPhoTightPassing_p;

  delete jtEtaPassing_p;
  delete jtPtPassing_p;
  delete jtDPhiPassing_p;
  delete jtGammaDRPassing_p;

  outFile_p->Close();
  delete outFile_p;

    outFileTxt.close();

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
