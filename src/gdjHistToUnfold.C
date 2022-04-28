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

//Primary function called in main (should be the only one - support functions go elsewhere)
int gdjHistToUnfold(std::string inConfigFileName)
{
  //Watch just for timing tests if things look low
  cppWatch globalTimer;
  globalTimer.start();

  //Random number generator seed is fixed but determined randomly for debuggable results
  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101             
  TRandom3* randGen_p = new TRandom3(randSeed);

  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  //Debugging is set as an environment variable - we use a simple class to pick up
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();


  //Global debug default spits out file and line as follows
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INRESPONSEFILENAME",
					      "INUNFOLDFILENAME",
					      "VARNAME",
                                              "OUTFILENAME",
					      "NITER",
					      "DOREWEIGHT",
					      "UNFOLDERRTYPE",
					      "NTOYS",
					      "DOREBIN"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string varName = config_p->GetValue("VARNAME", "");
  std::string varNameLower = returnAllLowercaseString(varName);

  //Multijet variables are a well defined list and we must handle differently, so define a bool
  const bool isMultijet = isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");
  
  
  const int nIter = config_p->GetValue("NITER", -1);
  const bool doReweight = config_p->GetValue("DOREWEIGHT", 0);

  //Rebinning requires additional params i.e. what is the rebin in X and Y
  const bool doRebin = (bool)config_p->GetValue("DOREBIN", 0);
  std::vector<std::string> rebinParams = {"REBINX",
					  "REBINY"};

  //Since arrays require predefined sizes, just define a large cap on them lie 1000 - we will limit ourselves to the subset of the array we need when we get to it
  const int nMaxBins = 1000;
  Int_t nRebinX, nRebinY;
  Double_t rebinX[nMaxBins+1], rebinY[nMaxBins+1];  
  std::vector<float> rebinXVect, rebinYVect;

  //We need some param delta defining when two bin edges are the same - here I pick something small
  const Float_t deltaValue = 0.001;

  //If we are doing rebin, grab the new arrays of bins 
  if(doRebin){
    if(!checkEnvForParams(config_p, rebinParams)) return 1;

    std::string rebinXStr = config_p->GetValue("REBINX", "");
    std::string rebinYStr = config_p->GetValue("REBINY", "");

    //Currently we require both rebin arrays specified even if we only wanna rebin one dimension
    if(rebinXStr.size() == 0){
      std::cout << "ERROR: 'DOREBIN' is " << doRebin << ", but rebinXStr is empty. return 1" << std::endl;
      return 1;
    }
    if(rebinYStr.size() == 0){
      std::cout << "ERROR: 'DOREBIN' is " << doRebin << ", but rebinYStr is empty. return 1" << std::endl;
      return 1;
    }

    //Convert the string to a vector of floats
    rebinXVect = strToVectF(rebinXStr);
    rebinYVect = strToVectF(rebinYStr);

    //Since we are making a choice w/ delta for bin uniqueness we need to check there is nothing below that level
    //Pre-check the bins requested are unique at the level of 2xDeltaValue
    if(!checkMinBinWidth(rebinXVect, 2*deltaValue)) return 1;
    if(!checkMinBinWidth(rebinYVect, 2*deltaValue)) return 1;

    //Finally, fill our final bin arrays
    nRebinX = rebinXVect.size() - 1;
    for(Int_t bIX = 0; bIX < nRebinX+1; ++bIX){
      rebinX[bIX] = rebinXVect[bIX];
    }
    nRebinY = rebinYVect.size() - 1;
    for(Int_t bIY = 0; bIY < nRebinY+1; ++bIY){
      rebinY[bIY] = rebinYVect[bIY];
    }   
  }

  //RooUnfold allows for many different methods for error calculation
  //My preference is toys - simple to understand, precision well defined by number of toys
  //Downside: Not exact, significant CPU time increase 
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

  //Define a cap on the number of iterations you can do to not waste time on absurd high amounts
  const int nIterMax = 100;
  if(nIter > nIterMax){
    std::cout << "Requested number of iterations \'" << nIter << "\' exceeds maximum allowed \'" << nIterMax << "\'" << std::endl;
  }

  //Modify the output file name for missing bits (.root, output directory, etc)
  if(outFileName.find(".root") != std::string::npos){
    outFileName.replace(outFileName.rfind(".root"), 5, "");
    outFileName = outFileName + "_" + dateStr + ".root";
  }
  if(outFileName.find("/") == std::string::npos){
    outFileName = outputStr + "/" + outFileName;
  }

  //Define some unfolding variables in arrays along with paired-by-position stylized versions
  //Previously i rebuilt everything from a stripped down unfolding TTree - ive re-worked this to be histogram based and temporarily comment out the *Tree variabls  
  std::vector<std::string> validUnfoldVar = {"Pt", "XJ", "XJJ"};
  //  std::vector<std::string> validUnfoldVarTree = {"XJ", "XJ", "XJJ"};
  std::vector<std::string> validUnfoldVarStyle = {"Jet p_{T}", "x_{J}", "#vec{x}_{JJ#gamma}"};

  //Check that our requested unfolding Var matches what is in the vectors + grab position
  const int vectPos = vectContainsStrPos(varName, &validUnfoldVar);
  if(vectPos < 0){
    std::cout << "GIVEN VARNAME \'" << varName << "\' is invalid. Please pick from:" << std::endl;
    for(unsigned int i = 0; i < validUnfoldVar.size(); ++i){
      std::cout << " \'" << validUnfoldVar[i] << "\'" << std::endl;
    }
    std::cout << "return 1" << std::endl;
    return 1;
  }
  //  std::string varNameTree = validUnfoldVarTree[vectPos];
  std::string varNameStyle = validUnfoldVarStyle[vectPos];
  
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;

  //We  need to define a cap here too - given we are extremely unlikely to ever drop below 10% level 10 is sufficient
  const Int_t nMaxCentBins = 10;

  //Grab the response file - this will also be the unfold file while working with MC - with data it will differ
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");  
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  
  //declare as nullptr; actual file / configs will depend on if its data or MC
  TFile* inUnfoldFile_p = nullptr;
  TEnv* inUnfoldFileConfig_p = nullptr;
  TEnv* inUnfoldFileLabels_p = nullptr;

  //If MC - same file; if data - different files
  const bool isSameInputFile = isStrSame(inResponseFileName, inUnfoldFileName);

  //Either have two var point to the same file and handle it that way or separate files (in case of data)
  if(isSameInputFile){
    inUnfoldFile_p = inResponseFile_p;
    inUnfoldFileConfig_p = inResponseFileConfig_p;
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }
  else{
    inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");  
    inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");    
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }

  //there is an additional set of parameters that we must make sure are defined from the file we are passed 
  std::vector<std::string> necessaryInFileParams = {"ISMC",
						    "DOUNIFIEDPURITY",
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
						    "GAMMAJTDPHI",
						    "ISPP"};

 
  std::vector<std::string> pbpbParams = {"CENTBINS"};
  std::vector<std::string> mcParams = {"ASSOCGENMINPT"};

  //Some of these params should be identical between data and MC else the unfold wont make any sense
  std::vector<std::string> paramsToCompare = {"DOUNIFIEDPURITY",
					       "JETR",
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

  //These are almost always true - throw them in the names and use the label map later to make sure we get the cuts on the plot
  //DPhi0 - the dphi cut
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string gammaJtDPhiStr = "DPhi0";
  //Multijtcut - the minimum pT cut on the jets to construct the multijet object
  const std::string multiJtCutGlobalStr = "MultiJt0";
  //dont remember
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";
 
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Return failed if we are missing required parameters
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  const bool isMCResponse = inResponseFileConfig_p->GetValue("ISMC", 1);
  //make sure response file given is MC
  if(!isMCResponse){
    std::cout << "Given responseFile \'" << inResponseFileName << "\' is not an MC file. return 1" << std::endl;
    return 1;
  }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Grab first set of parameters defined in necessaryInFileParams
  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  const bool isPP = inResponseFileConfig_p->GetValue("ISPP", 0);
  const bool doUnifiedPurity = inResponseFileConfig_p->GetValue("DOUNIFIEDPURITY", 0);

  //Unified purity refers to having purity defined barrel and endcap separately or together
  //Originally we were using a separately defined purity but have switched to a barrel + endcap combined definition
  std::cout << "UNIFIED PURITY?" << std::endl;
  std::cout << " " << doUnifiedPurity << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //MixMode - 0 is p+p (no mixing), 1 is normal inclusive jet mixing, 2 is multi-jet mixing
  std::string mixModeStr = "MIXMODE0";
  if(!isPP){
    if(isMultijet) mixModeStr = "MIXMODE2";
    else mixModeStr = "MIXMODE1";
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
  //Grab the in-file defined bins
  Int_t nGammaPtBins = inResponseFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxBins+1];
  const Float_t gammaPtBinsLowReco = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOWRECO", gammaPtBinsLow);
  const Float_t gammaPtBinsHighReco = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGHRECO", gammaPtBinsHigh);

  const Int_t nGammaFinePtBins = inResponseFileConfig_p->GetValue("NGAMMAFINEPTBINS", 0);
  std::string gammaFinePtBinsStr = inResponseFileConfig_p->GetValue("GAMMAFINEPTBINS", "");
  std::vector<float> gammaFinePtBinsVect = strToVectF(gammaFinePtBinsStr);
  Double_t gammaFinePtBins[nMaxBins+1];

  for(unsigned int i = 0; i < gammaFinePtBinsVect.size(); ++i){gammaFinePtBins[i] = gammaFinePtBinsVect[i];}
  
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nJtPtBins = inResponseFileConfig_p->GetValue("NJTPTBINS", 0);
  const Float_t jtPtBinsLow = inResponseFileConfig_p->GetValue("JTPTBINSLOW", -1.0);
  const Float_t jtPtBinsHigh = inResponseFileConfig_p->GetValue("JTPTBINSHIGH", -1.0);
  const Float_t jtPtBinsLowReco = inResponseFileConfig_p->GetValue("JTPTBINSLOWRECO", -1.0);
  const Float_t jtPtBinsHighReco = inResponseFileConfig_p->GetValue("JTPTBINSHIGHRECO", -1.0);
  const Bool_t jtPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOLOG", 1);
  const Bool_t jtPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOCUSTOM", 1);
  Double_t jtPtBins[nMaxBins+1];
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

  const Double_t gammaJtDPhiCut = mathStringToNum(inResponseFileConfig_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);

  std::cout << "GAMMAJTDPHICUT: " << gammaJtDPhiCut << std::endl;
  
  //  const Float_t assocGenMinPt = inResponseFileConfig_p->GetValue("ASSOCGENMINPT", -1.0);
  
  const Int_t nXJBins = inResponseFileConfig_p->GetValue("NXJBINS", 20);
  const Float_t xjBinsLow = inResponseFileConfig_p->GetValue("XJBINSLOW", 0.0);
  const Float_t xjBinsHigh = inResponseFileConfig_p->GetValue("XJBINSHIGH", 2.0);
  const Bool_t xjBinsDoCustom = inResponseFileConfig_p->GetValue("XJBINSDOCUSTOM", 0);
  Double_t xjBins[nMaxBins+1];
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
  if(isSameInputFile){
    if(!compEnvParams(inResponseFileConfig_p, inUnfoldFileConfig_p, paramsToCompare)) return 1;
  }

  //Gotta do this bit manually
  Int_t nVarBins = nXJBins;
  Double_t varBins[nMaxBins+1];
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
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  

  TH1F* photonPtReco_p[nMaxCentBins];
  TH1F* photonPtTruth_p[nMaxCentBins];  

  TH2F* photonPtJetVarReco_p[nMaxCentBins];
  TH2F* photonPtJetVarTruth_p[nMaxCentBins];  
  
  //  const Int_t nBarrelAndEC = 2;
  //  const std::string barrelAndECStr[nBarrelAndEC] = {"Barrel", "EC"};

  //Define some necessary variables for handling Barrel/EC separately or unified like labels
  const Int_t nBarrelAndECMax = 2;
  Int_t nBarrelAndEC = nBarrelAndECMax;
  std::string barrelAndECStr[nBarrelAndECMax] = {"Barrel", "EC"};
  if(doUnifiedPurity){
    nBarrelAndEC = 1;
    barrelAndECStr[0] = "BarrelAndEC";

    std::cout << "SWITCHING TO UNIFIED PURITY MODE" << std::endl;
  }
  
  TH1F* photonPtReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];  
  TH1F* photonPtReco_PURCORR_COMBINED_p[nMaxCentBins];  
  TH1F* photonPtReco_TRUTH_p[nMaxCentBins][nBarrelAndECMax];
  TH1F* photonPtReco_TRUTH_COMBINED_p[nMaxCentBins];

  TH2F* photonPtJetVarReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];
  TH2F* photonPtJetVarReco_PURCORR_COMBINED_p[nMaxCentBins];
  TH2F* photonPtJetVarReco_TRUTH_p[nMaxCentBins][nBarrelAndECMax];
  TH2F* photonPtJetVarReco_TRUTH_COMBINED_p[nMaxCentBins];
  
  RooUnfoldResponse* rooResGamma_p[nMaxCentBins];  
  TH2F* rooResGammaMatrix_p[nMaxCentBins];
  TH1F* rooResGammaMisses_p[nMaxCentBins];

  RooUnfoldResponse* rooResGammaJetVar_p[nMaxCentBins]; 
  TH2F* rooResGammaJetVarMatrixReco_p[nMaxCentBins];
  TH2F* rooResGammaJetVarMatrixTruth_p[nMaxCentBins];
  TH2F* rooResGammaJetVarMisses_p[nMaxCentBins];
  
  //Lets check if a rebin is requested that it is valid
  if(doRebin){
    std::string histName = centBinsStr[0] + "/photonPtJt" + varName + "VCent_" + centBinsStr[0] + "_" + barrelAndECStr[0] + "_" + gammaJtDPhiStr + "_" + mixModeStr + "_SUB_h";
    std::cout << "HISTNAME: " << histName  << std::endl;
    TH2F* tempTH2ForCheck_p = (TH2F*)inUnfoldFile_p->Get(histName.c_str());

    //Check that the macro histogram bins align with requested new binning                      
    if(!checkHistContainsBins(rebinXVect, tempTH2ForCheck_p, deltaValue, false)) return 1;
    if(!checkHistContainsBins(rebinYVect, tempTH2ForCheck_p, deltaValue, true)) return 1;    

    //Since we got out of previous check, reset the varbins!
    nVarBins = nRebinX;
    for(Int_t vI = 0; vI < nRebinX+1; ++vI){
      varBins[vI] = rebinXVect[vI];
    }

    nGammaPtBins = nRebinY;
    for(Int_t vI = 0; vI < nRebinY+1; ++vI){
      gammaPtBins[vI] = rebinYVect[vI];
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++ cI){    
    photonPtReco_p[cI] = new TH1F(("photonPtReco_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Photon p_{T};Counts", nGammaPtBins, gammaPtBins);
    photonPtTruth_p[cI] = new TH1F(("photonPtTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Photon p_{T};Counts", nGammaPtBins, gammaPtBins);
    
    photonPtJetVarReco_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetVarTruth_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      std::string baseName = centBinsStr[cI] + "/photonPt";
      std::string jtName = "Jt" + varName;
      std::string midName =  "VCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI];
      std::string photonPtPurCorrName = baseName + midName + "_PURCORR_h";
      std::string photonPtJtVarPurCorrName = baseName + jtName + midName + "_" + gammaJtDPhiStr + "_PURCORR_h";
      std::string photonPtTruthName = baseName + midName + "_TRUTH_h";
      std::string photonPtJtVarTruthName = baseName + jtName + midName + "_" + gammaJtDPhiStr + "_" + mixModeStr + "_TRUTH_h";

      TFile* tempFilePointer_p = inUnfoldFile_p;
      if(!isMC) tempFilePointer_p = inResponseFile_p;
      
      if(doRebin){
	photonPtReco_PURCORR_p[cI][eI] = (TH1F*)inUnfoldFile_p->Get(photonPtPurCorrName.c_str());     

	TH2F* tempHistTH2F_p = (TH2F*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());
	strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_h");		  
	std::string xTitle = tempHistTH2F_p->GetXaxis()->GetTitle();
	std::string yTitle = tempHistTH2F_p->GetYaxis()->GetTitle();
	std::string titleStr = ";" + xTitle + ";" + yTitle;
	std::cout << "LINE, NAME: " << __LINE__ << ", " << photonPtJtVarTruthName << std::endl;

	photonPtJetVarReco_PURCORR_p[cI][eI] = new TH2F(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nGammaPtBins, gammaPtBins, nVarBins, varBins);
	
	fineHistToCoarseHist(tempHistTH2F_p, photonPtJetVarReco_PURCORR_p[cI][eI]);
	photonPtReco_TRUTH_p[cI][eI] = (TH1F*)tempFilePointer_p->Get(photonPtTruthName.c_str());
	
	tempHistTH2F_p = (TH2F*)tempFilePointer_p->Get(photonPtJtVarTruthName.c_str());
	strReplace(&photonPtJtVarTruthName, "TRUTH_h", "TRUTH_REBIN_h");		  
	xTitle = tempHistTH2F_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH2F_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;
	photonPtJetVarReco_TRUTH_p[cI][eI] = new TH2F(photonPtJtVarTruthName.c_str(), titleStr.c_str(), nGammaPtBins, gammaPtBins, nVarBins, varBins);

	fineHistToCoarseHist(tempHistTH2F_p, photonPtJetVarReco_TRUTH_p[cI][eI]);
      }
      else{
	std::cout << "CORRNAME: " << photonPtPurCorrName << std::endl;
	photonPtReco_PURCORR_p[cI][eI] = (TH1F*)inUnfoldFile_p->Get(photonPtPurCorrName.c_str());     
	std::cout << "LINE, NAMEPURCORR: " << __LINE__ << ", " << photonPtJtVarPurCorrName << std::endl;
	photonPtJetVarReco_PURCORR_p[cI][eI] = (TH2F*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());     
	
	photonPtReco_TRUTH_p[cI][eI] = (TH1F*)tempFilePointer_p->Get(photonPtTruthName.c_str());
	std::cout << "TRUTH HIST: " << photonPtJtVarTruthName << std::endl;
	std::cout << "LINE, NAMETRUTH: " << __LINE__ << ", " << photonPtJtVarTruthName << std::endl;
	photonPtJetVarReco_TRUTH_p[cI][eI] = (TH2F*)tempFilePointer_p->Get(photonPtJtVarTruthName.c_str());
      }
    }
    
    photonPtReco_PURCORR_COMBINED_p[cI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), ";Reco. Photon p_{T};Counts (Purity Corrected)", nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_PURCORR_COMBINED_p[cI] = new TH2F(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    
    photonPtReco_TRUTH_COMBINED_p[cI] = new TH1F(("photonPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_TRUTH_COMBINED_p[cI] = new TH2F(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);   

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
  
    for(Int_t bIX = 0; bIX < photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      for(Int_t bIY = 0; bIY < photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins(); ++bIY){

	Float_t content = 0.0;
	Float_t error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetVarReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetVarReco_PURCORR_p[cI][eI]->GetBinError(bIX+1, bIY+1);

	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetVarReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetVarReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);	

	content = 0.0;
	error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetVarReco_TRUTH_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetVarReco_TRUTH_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	  
	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetVarReco_TRUTH_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetVarReco_TRUTH_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);		
      }
    }

    std::vector<TH1*> hists_p = {photonPtJetVarReco_p[cI], photonPtJetVarTruth_p[cI]};
    hists_p.push_back(photonPtJetVarReco_PURCORR_COMBINED_p[cI]);
    hists_p.push_back(photonPtJetVarReco_TRUTH_COMBINED_p[cI]);
    
    setSumW2(hists_p);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  //Grab the unfolding ttree for matrix filling

  TTree* unfoldTree_p = (TTree*)inResponseFile_p->Get("unfoldTree_p");

  //Before declaring variables define array size for JES/JER syst. + check it matches input
  const Int_t nJESSys = 18;
  const Int_t nJERSys = 9;
  const Int_t nJetSysAndNom = 1 + nJESSys + nJERSys;
  const Int_t nJetSysAndNomInput = inResponseFileConfig_p->GetValue("NJETSYSANDNOM", -1);

  if(nJetSysAndNom != nJetSysAndNomInput){
    std::cout << "Mismatch between number of jet systematics expected and input from file \'" << inResponseFileName << "\'. return 1" << std::endl;
    return 1;
  }

  //Run,lumi,event check for duplication
  Int_t runNumber;
  UInt_t lumiBlock;
  Int_t eventNumber;
  
  //Declare unfolding TTree variables
  const Int_t nMaxJets = 100;

  Float_t recoGammaPt_;
  Float_t recoGammaPhi_;
  Float_t truthGammaPt_;
  Float_t truthGammaPhi_;

  Int_t nRecoJt_;
  Float_t recoJtPt_[nMaxJets][nJetSysAndNom];
  Float_t recoJtPhi_[nMaxJets];
  Float_t recoJtEta_[nMaxJets];

  Float_t truthJtPt_[nMaxJets];
  Float_t truthJtPhi_[nMaxJets];
  Float_t truthJtEta_[nMaxJets];

  Int_t nTruthJtUnmatched_;
  Float_t truthJtUnmatchedPt_[nMaxJets];
  Float_t truthJtUnmatchedPhi_[nMaxJets];
  Float_t truthJtUnmatchedEta_[nMaxJets];

  Float_t unfoldWeight_;
  Float_t unfoldCent_;

  unfoldTree_p->SetBranchAddress("runNumber", &runNumber);
  unfoldTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  unfoldTree_p->SetBranchAddress("eventNumber", &eventNumber);
  
  unfoldTree_p->SetBranchAddress("recoGammaPt", &recoGammaPt_);
  unfoldTree_p->SetBranchAddress("recoGammaPhi", &recoGammaPhi_);
  unfoldTree_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldTree_p->SetBranchAddress("truthGammaPhi", &truthGammaPhi_);
  
  unfoldTree_p->SetBranchAddress("nRecoJt", &nRecoJt_);
  unfoldTree_p->SetBranchAddress("recoJtPt", recoJtPt_);
  unfoldTree_p->SetBranchAddress("recoJtPhi", recoJtPhi_);
  unfoldTree_p->SetBranchAddress("recoJtEta", recoJtEta_);
  
  unfoldTree_p->SetBranchAddress("truthJtPt", truthJtPt_);
  unfoldTree_p->SetBranchAddress("truthJtPhi", truthJtPhi_);
  unfoldTree_p->SetBranchAddress("truthJtEta", truthJtEta_);
  
  unfoldTree_p->SetBranchAddress("nTruthJtUnmatched", &nTruthJtUnmatched_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPt", truthJtUnmatchedPt_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPhi", truthJtUnmatchedPhi_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedEta", truthJtUnmatchedEta_);
  
  unfoldTree_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  unfoldTree_p->SetBranchAddress("unfoldCent", &unfoldCent_);    

  bool varIsPhi = varName.find("DPhiJJ") != std::string::npos;
  
  //Starting out we need to calc our gamma normalization
  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    rooResGamma_p[cI] = new RooUnfoldResponse(photonPtReco_p[cI], photonPtTruth_p[cI], ("rooResGamma_" + centBinsStr[cI]).c_str(), "");
    rooResGammaMatrix_p[cI] = new TH2F(("rooResGammaMatrix_" + centBinsStr[cI] + "_h").c_str(), ";Reco #gamma p_{T} [GeV]; Truth #gamma p_{T} [GeV]", nGammaPtBins, gammaPtBins, nGammaPtBins, gammaPtBins);
    rooResGammaMisses_p[cI] = new TH1F(("rooResGammaMisses_" + centBinsStr[cI] + "_h").c_str(), ";Truth #gamma p_{T} [GeV]; Counts (weighted)", nGammaPtBins, gammaPtBins);


    rooResGammaJetVar_p[cI] = new RooUnfoldResponse(photonPtJetVarReco_p[cI], photonPtJetVarTruth_p[cI], ("rooResGammaJetVar_" + centBinsStr[cI]).c_str(), "");
    rooResGammaJetVarMatrixReco_p[cI] = new TH2F(("rooResJetVarGammaMatrixReco_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Jet p_{T} [GeV]; Reco. #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
    rooResGammaJetVarMatrixTruth_p[cI] = new TH2F(("rooResJetVarGammaMatrixTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Jet p_{T} [GeV]; Truth #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
    rooResGammaJetVarMisses_p[cI] = new TH2F(("rooResJetVarGammaMisses_" + centBinsStr[cI] + "_h").c_str(), ";Truth Misses Jet p_{T} [GeV]; Truth Misses #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
  }

  std::cout << "Constructing response matrices..." << std::endl;
  const ULong64_t nEntriesUnfold = (ULong64_t)unfoldTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntriesUnfold/20);

  //Check we don't have duplicate events
  std::vector<Int_t> runNumbers, eventNumbers;
  std::vector<UInt_t> lumiBlocks;
  
  for(ULong64_t entry = 0; entry < nEntriesUnfold; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntriesUnfold << "..." << std::endl;
    unfoldTree_p->GetEntry(entry);

    //pushbacks
    runNumbers.push_back(runNumber);
    lumiBlocks.push_back(lumiBlock);
    eventNumbers.push_back(eventNumber);
    
    Int_t centPos = -1;
    if(isPP) centPos = 0;
    else centPos = ghostPos(centBinsF, unfoldCent_, true, doGlobalDebug);
    if(centPos < 0) continue;

    bool truthGammaOutOfBounds = truthGammaPt_ < gammaPtBinsLow || truthGammaPt_ >= gammaPtBinsHigh;
    if(truthGammaOutOfBounds) continue;
    
    //Construct 1-D unfold response matrix for gamma-pt (required in normalization)
    bool recoGammaOutOfBounds = recoGammaPt_ < gammaPtBinsLowReco || recoGammaPt_ >= gammaPtBinsHighReco;
    if(recoGammaOutOfBounds){
      rooResGamma_p[centPos]->Miss(truthGammaPt_, unfoldWeight_);
      rooResGammaMisses_p[centPos]->Fill(truthGammaPt_, unfoldWeight_);
    }
    else{
      rooResGamma_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_);
      rooResGammaMatrix_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_);
    }

    //Now construct 2-D unfold response matrix for gammapt-jetvariable    
    //Start w/ truth jets, no reco

    
    for(Int_t tI = 0; tI < nTruthJtUnmatched_; ++tI){
      bool truthJetOutOfBounds = truthJtUnmatchedPt_[tI] < jtPtBinsLow || truthJtUnmatchedPt_[tI] >= jtPtBinsHigh;
      if(truthJetOutOfBounds) continue;

      Float_t gammaJtDPhiTruth = TMath::Abs(getDPHI(truthJtUnmatchedPhi_[tI], truthGammaPhi_));
      bool gammaJtPassesDPhiTruth = gammaJtDPhiTruth >= gammaJtDPhiCut;
      if(!gammaJtPassesDPhiTruth) continue;
      
      if(varName.find("Pt") != std::string::npos){
	rooResGammaJetVar_p[centPos]->Miss(truthJtUnmatchedPt_[tI], truthGammaPt_, unfoldWeight_);
	rooResGammaJetVarMisses_p[centPos]->Fill(truthJtUnmatchedPt_[tI], truthGammaPt_, unfoldWeight_);
      }
    }
    
    
    for(Int_t jI = 0; jI < nRecoJt_; ++jI){
      bool recoJtOutOfBounds = recoJtPt_[jI][0] < jtPtBinsLowReco || recoJtPt_[jI][0] >= jtPtBinsHighReco;
      bool truthJetOutOfBounds = truthJtPt_[jI] < jtPtBinsLow || truthJtPt_[jI] >= jtPtBinsHigh;
      if(truthJetOutOfBounds) continue;
      
      Float_t gammaJtDPhiTruth = TMath::Abs(getDPHI(truthJtPhi_[jI], truthGammaPhi_));
      bool gammaJtPassesDPhiTruth = gammaJtDPhiTruth >= gammaJtDPhiCut;
      if(!gammaJtPassesDPhiTruth) continue;

      Float_t gammaJtDPhiReco = -999;
      bool gammaJtPassesDPhiReco = false;
      if(!recoGammaOutOfBounds){
	gammaJtDPhiReco = TMath::Abs(getDPHI(recoJtPhi_[jI], recoGammaPhi_));
	gammaJtPassesDPhiReco = gammaJtDPhiReco >= gammaJtDPhiCut;
      }
      
      if(varName.find("Pt") != std::string::npos){
	if(recoJtOutOfBounds || recoGammaOutOfBounds || !gammaJtPassesDPhiReco){
	  rooResGammaJetVar_p[centPos]->Miss(truthJtPt_[jI], truthGammaPt_, unfoldWeight_);
	  rooResGammaJetVarMisses_p[centPos]->Fill(truthJtPt_[jI], truthGammaPt_, unfoldWeight_);
	}
	else{
	  rooResGammaJetVar_p[centPos]->Fill(recoJtPt_[jI][0], recoGammaPt_, truthJtPt_[jI], truthGammaPt_, unfoldWeight_);
	  rooResGammaJetVarMatrixReco_p[centPos]->Fill(recoJtPt_[jI][0], recoGammaPt_, unfoldWeight_);
	  rooResGammaJetVarMatrixTruth_p[centPos]->Fill(truthJtPt_[jI], truthGammaPt_, unfoldWeight_);
	}
      }
    }    
  }
  std::cout << "Construction of response matrices is complete."   << std::endl;

  //This doesnt work because MC repeats numbers...
  /* 
  std::cout << "Check for event duplication..." << std::endl;

  bool hasDuplicates = false;
  for(unsigned int rI = 0; rI < runNumbers.size(); ++rI){
    int runN = runNumbers[rI];
    unsigned int lumiB = lumiBlocks[rI];
    int eventN = eventNumbers[rI];

    for(unsigned int rI2 = rI+1; rI2 < runNumbers.size(); ++rI2){
      if(runN != runNumbers[rI2]) continue;
      if(lumiB != lumiBlocks[rI2]) continue;
      if(eventN != eventNumbers[rI2]) continue;

      std::cout << " DUPLICATED RUN LUMI EVENT: " << runN << ", " << lumiB << ", " << eventN << std::endl;
      
      hasDuplicates = true;
    }
  }

  if(hasDuplicates){
    std::cout << "Duplication found! RETURN AND FIX" << std::endl;
    return 1;
  }
  else std::cout << "No duplication found" << std::endl;

  //  return 1;
  */
  
  TH2F* reweightingHist_p[nMaxCentBins];
  if(!isSameInputFile && doReweight){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      reweightingHist_p[cI] = new TH2F(("reweightingHist_" + centBinsStr[cI] + "_h").c_str(), (";Reco. Jet " + varNameStyle + ";Gamma p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

      for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
	for(Int_t bIY = 0; bIY < nGammaPtBins; ++bIY){
	
	  Float_t numeratorData = photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinContent(bIX+1, bIY+1);
	  Float_t denominatorMC = photonPtJetVarReco_p[cI]->GetBinContent(bIX+1, bIY+1);

	  Float_t ratioVal = 0.0;
	  if(denominatorMC > TMath::Power(10, -50)) ratioVal = numeratorData/denominatorMC;
	  
	  reweightingHist_p[cI]->SetBinContent(bIX+1, bIY+1, ratioVal);
	}
      }
    }

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      delete photonPtJetVarReco_p[cI];
      delete photonPtJetVarTruth_p[cI];

      photonPtJetVarReco_p[cI] = new TH2F(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
      photonPtJetVarTruth_p[cI] = new TH2F(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    }
  }
  
  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);
    
    for(int gI = 0; gI < nGammaPtBins; ++gI){
      TH1D* xjRecoProjection_PURCORR_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
      TH1D* xjRecoProjection_TRUTH_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);

      for(int xI = 0; xI < nVarBins; ++xI){
	xjRecoProjection_PURCORR_COMBINED_p->SetBinContent(xI+1, photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_PURCORR_COMBINED_p->SetBinError(xI+1, photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinError(xI+1, gI+1));

	xjRecoProjection_TRUTH_COMBINED_p->SetBinContent(xI+1, photonPtJetVarReco_TRUTH_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_TRUTH_COMBINED_p->SetBinError(xI+1, photonPtJetVarReco_TRUTH_COMBINED_p[cI]->GetBinError(xI+1, gI+1));
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
      
      xjRecoProjection_PURCORR_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_PURCORR_COMBINED_p;

      xjRecoProjection_TRUTH_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_TRUTH_COMBINED_p;
    }
    
    //Write the full histograms at reco level to file
    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtTruth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    photonPtJetVarReco_PURCORR_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtJetVarReco_TRUTH_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    TH2D* tempTruthHist_p = (TH2D*)rooResGammaJetVar_p[cI]->Htruth();
    tempTruthHist_p->Write(("photonPtJet" + varName + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH2_COMBINED_h").c_str(), TObject::kOverwrite); 

    TH2D* tempRecoHist_p = (TH2D*)rooResGammaJetVar_p[cI]->Hmeasured();
    tempRecoHist_p->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR2_COMBINED_h").c_str(), TObject::kOverwrite); 
    
    
    rooResGammaMatrix_p[cI]->Write("", TObject::kOverwrite);
    rooResGammaMisses_p[cI]->Write("", TObject::kOverwrite);

    std::cout << "HIST NAMES:" << std::endl;
    std::cout << " HIST " << rooResGammaJetVarMatrixReco_p[cI]->GetName() << std::endl;
    std::cout << " HIST " << rooResGammaJetVarMatrixTruth_p[cI]->GetName() << std::endl;
    std::cout << " HIST " << rooResGammaJetVarMisses_p[cI]->GetName() << std::endl;
    //    return 1;
    
    rooResGammaJetVarMatrixReco_p[cI]->Write("", TObject::kOverwrite);
    rooResGammaJetVarMatrixTruth_p[cI]->Write("", TObject::kOverwrite);
    rooResGammaJetVarMisses_p[cI]->Write("", TObject::kOverwrite);
    
    //Full unfold gamma only
    for(int i = 1; i <= nIter; ++ i){
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

    //Gamma-jet unfold
    for(int i = 1; i <= nIter; ++ i){
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGammaJetVar_p[cI], photonPtJetVarReco_PURCORR_COMBINED_p[cI], i);
      
      TH2D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH2D*)rooBayes_p->Hreco()->Clone(("photonPtJetVarReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      else if(unfoldErrType == 1){
        rooBayes_p->SetNToys(nToys);
        unfolded_p = (TH2D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtJetVarReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      }

      
      TH2D* refolded_p = (TH2D*)rooResGammaJetVar_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str());

      unfolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str(), TObject::kOverwrite);

      delete unfolded_p;
      delete refolded_p;
      delete rooBayes_p;
    }
    

    //Full unfold gamma+jet, should be instant when the unfold file is also mc    
    /*
    for(int i = 1; i < nIter; ++ i){
      std::string saveNameUnfold = "photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h";
      std::string saveNameRefold = "photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h";

      //2021.11.09 debug - added rebinning intermediary but appears there are serious problems
      if(saveNameUnfold.find("photonPtJetPtReco_Iter2_Cent0to10_PURCORR_COMBINED_h") != std::string::npos){

	std::cout << "Histogram to Unfold Print: " << std::endl;
	photonPtJetVarReco_PURCORR_COMBINED_p[cI]->Print("ALL");

	std::cout << "Hist to unfold nX, nY: " << photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins() << ", " << photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins() << std::endl;
	return 1;
      }
      
    
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGammaJetVar_p[cI], photonPtJetVarReco_PURCORR_COMBINED_p[cI], i);  
      
      TH2D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH2D*)rooBayes_p->Hreco()->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      else if(unfoldErrType == 1){
	rooBayes_p->SetNToys(nToys);
	unfolded_p = (TH2D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());//Very annoyingly Hreco is the name of the function that returns unfolded distribution
      }
      TH2D* refolded_p = (TH2D*)rooResGammaJetVar_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str());

      unfolded_p->Write(saveNameUnfold.c_str(), TObject::kOverwrite);
      refolded_p->Write(saveNameRefold.c_str(), TObject::kOverwrite);

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

    photonPtJetVarReco_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetVarTruth_p[cI]->Write("", TObject::kOverwrite);

    photonPtJetVarReco_PURCORR_COMBINED_p[cI]->Write("", TObject::kOverwrite);
    photonPtJetVarReco_TRUTH_COMBINED_p[cI]->Write("", TObject::kOverwrite);
    
    if(!isSameInputFile && doReweight){
      reweightingHist_p[cI]->Write("", TObject::kOverwrite);
    }
    */
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;  

    delete photonPtReco_p[cI];
    delete photonPtTruth_p[cI];

    delete photonPtJetVarReco_p[cI];
    delete photonPtJetVarTruth_p[cI];
    
    //    delete photonPtJetVarReco_p[cI];
    //    delete photonPtJetVarTruth_p[cI];

    //    delete photonPtJetVarReco_PURCORR_COMBINED_p[cI];
    //    delete photonPtJetVarReco_TRUTH_COMBINED_p[cI];
    
    //    delete rooResGammaJetVar_p[cI];

    if(!isSameInputFile && doReweight){
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
  
  if(!isSameInputFile){
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
