//Author: Chris McGinn (2021.01.13)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/configParser.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/unfoldingUtil.h"
// #include "include/treeUtil.h"

int gdjPlotUnfoldDiagnostics(std::string inConfigFileName)
{
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.2;
  
  cppWatch globalTimer;
  globalTimer.start();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INUNFOLDFILENAME",
					      "OUTUNFOLDCONFIG",
					      "TERMTHRESH",
					      "DORELATIVETERM",
					      "SAVETAG",
					      "SAVEEXT",
					      "DOLOGX",
					      "DOLOGY",
					      "XJMIN.GLOBAL",
					      "XJMAX.GLOBAL",
					      "XJLEGX.GLOBAL.0",
					      "XJLEGY.GLOBAL.0",
					      "XJLABELX.GLOBAL",
					      "XJLABELY.GLOBAL",
					      "XJLABELALIGNRIGHT.GLOBAL",
					      "DIAGLEGX.GLOBAL.0",
					      "DIAGLEGY.GLOBAL.0",
					      "DIAGLABELX.GLOBAL",
					      "DIAGLABELY.GLOBAL",
					      "DIAGLABELALIGNRIGHT.GLOBAL"};

  
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outUnfoldConfigName = config_p->GetValue("OUTUNFOLDCONFIG", "");
  TEnv* outUnfoldConfig_p = new TEnv();

  const double termThresh = config_p->GetValue("TERMTHRESH", -1.0);
  outUnfoldConfig_p->SetValue("TERMTHRESH", config_p->GetValue("TERMTHRESH", termThresh));

  const std::string saveTag = config_p->GetValue("SAVETAG", "");
  const std::string saveExt = config_p->GetValue("SAVEEXT", "");
  const bool doRelativeTerm = config_p->GetValue("DORELATIVETERM",  0);
  std::string deltaStr = "Rel";
  if(!doRelativeTerm) deltaStr = "Abs";

  
  const bool doLogX = config_p->GetValue("DOLOGX", 0);
  const bool doLogY = config_p->GetValue("DOLOGY", 0);
  
  std::vector<std::string> validExt = {"pdf", "png"};

  if(!vectContainsStr(saveExt, &validExt)){
    std::cout << "Given SAVEEXT \'" << saveExt << "\' is not valid. return 1" << std::endl;
    return 1;
  }  
  
  const Int_t nIterCap = config_p->GetValue("NITERPLOTCAP", 100000);
  
  const Float_t xjMinGlobal = config_p->GetValue("XJMIN.GLOBAL", -1000.);
  const Float_t xjMaxGlobal = config_p->GetValue("XJMAX.GLOBAL", 1000.);

  const Float_t xjRatMinGlobal = config_p->GetValue("XJRATMIN.GLOBAL", -1000.);
  const Float_t xjRatMaxGlobal = config_p->GetValue("XJRATMAX.GLOBAL", 1000.);
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<Float_t> xjLegXGlobal, xjLegYGlobal, xjLegNGlobal;
  for(Int_t i = 0; i < 100; ++i){
    Float_t tempLegX = config_p->GetValue(("XJLEGX.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegY = config_p->GetValue(("XJLEGY.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegN = config_p->GetValue(("XJLEGN.GLOBAL." + std::to_string(i)).c_str(), -1000.);

    if(tempLegX >= 0){
      xjLegXGlobal.push_back(tempLegX);
      xjLegYGlobal.push_back(tempLegY);
      if(tempLegN >= 0) xjLegNGlobal.push_back(tempLegN);
    }
  }
  if(xjLegNGlobal.size() != xjLegXGlobal.size()){
    std::cout << "Number of histograms per legend not specified manually for all legends - reverting to equal split" << std::endl;
    xjLegNGlobal.clear();
  }
  
  const Float_t xjLabelXGlobal = config_p->GetValue("XJLABELX.GLOBAL", -1000.);
  const Float_t xjLabelYGlobal = config_p->GetValue("XJLABELY.GLOBAL", -1000.);
  const Bool_t xjLabelAlignRightGlobal = config_p->GetValue("XJLABELALIGNRIGHT.GLOBAL", 0);

  std::vector<Float_t> diagLegXGlobal, diagLegYGlobal, diagLegNGlobal;
  for(Int_t i = 0; i < 100; ++i){
    Float_t tempLegX = config_p->GetValue(("DIAGLEGX.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegY = config_p->GetValue(("DIAGLEGY.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegN = config_p->GetValue(("DIAGLEGN.GLOBAL." + std::to_string(i)).c_str(), -1000.);

    if(tempLegX >= 0){
      diagLegXGlobal.push_back(tempLegX);
      diagLegYGlobal.push_back(tempLegY);
      if(tempLegN >= 0) diagLegNGlobal.push_back(tempLegN);
    }
  }
  if(diagLegNGlobal.size() != diagLegXGlobal.size()){
    std::cout << "Number of histograms per legend not specified manually for all legends - reverting to equal split" << std::endl;
    diagLegNGlobal.clear();
  }
  
  const Float_t diagLabelXGlobal = config_p->GetValue("DIAGLABELX.GLOBAL", -1000.);
  const Float_t diagLabelYGlobal = config_p->GetValue("DIAGLABELY.GLOBAL", -1000.);
  const Bool_t diagLabelAlignRightGlobal = config_p->GetValue("DIAGLABELALIGNRIGHT.GLOBAL", 0);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Get global labels
  std::vector<std::string> globalLabels;
  std::vector<std::string> globalLabelsDataOverrides;
  std::vector<std::string> globalLabelsMCOverrides;
  for(unsigned int i = 0; i < 100; ++i){
    std::string tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
    if(tempStr.size() != 0){
      globalLabels.push_back(tempStr);        
      tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i) + ".DATAOVERRIDE").c_str(), "");
      if(tempStr.size() != 0) globalLabelsDataOverrides.push_back(tempStr);
      else globalLabelsDataOverrides.push_back(globalLabels[globalLabels.size()-1]);

      tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i) + ".MCOVERRIDE").c_str(), "");
      if(tempStr.size() != 0) globalLabelsMCOverrides.push_back(tempStr);
      else globalLabelsMCOverrides.push_back(globalLabels[globalLabels.size()-1]);
    }
  }

  std::cout << "DEFAULT: " << std::endl;
  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    std::cout << " " << gI << "/" << globalLabels.size() << ": " << globalLabels[gI] << std::endl;
  }

  std::cout << "DEFAULT: " << std::endl;
  for(unsigned int gI = 0; gI < globalLabelsMCOverrides.size(); ++gI){
    std::cout << " " << gI << "/" << globalLabelsMCOverrides.size() << ": " << globalLabelsMCOverrides[gI] << std::endl;
  }
  


  //check that our given input files exist
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;
  
  TFile* inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");
  TEnv* inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");    
  TEnv* inUnfoldFileLabel_p = (TEnv*)inUnfoldFile_p->Get("label");    

  std::vector<std::string> necessaryUnfoldFileParams = {"ISMC",
							"NGAMMAPTBINS",
							"GAMMAPTBINSLOW",
							"GAMMAPTBINSHIGH",
							"GAMMAPTBINSDOLOG",
							"GAMMAPTBINSDOCUSTOM",
							"ISPP",
							"NITER",
							"JETR",
							"VARNAME",
							"ASSOCGENMINPT"};
  
  const std::string gammaJtDPhiStr = "DPhi0";
  //  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> pbpbParams = {"CENTBINS"};
  if(!checkEnvForParams(inUnfoldFileConfig_p, necessaryUnfoldFileParams)) return 1;

  std::string varName = inUnfoldFileConfig_p->GetValue("VARNAME", "");
  outUnfoldConfig_p->SetValue("VARNAME", varName.c_str());
  
  std::string varNameUpper = strLowerToUpper(varName);
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string varNameLabel = varNameLower;

  std::map<std::string, std::string> varNameLowerToLabel;
  varNameLowerToLabel["xj"] = "x_{J#gamma}";
  varNameLowerToLabel["xjj"] = "#vec{x}_{JJ#gamma}";

  if(varNameLowerToLabel.count(varNameLabel) != 0) varNameLabel = varNameLowerToLabel[varNameLabel];
  
  //Add to global labels
  int jetR = inUnfoldFileConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  else{
    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  }
  globalLabels.push_back(jetRStr);
  globalLabels.push_back(jtPtBinsGlobalStr);
  globalLabels.push_back(gammaJtDPhiStr);

  globalLabelsDataOverrides.push_back(jetRStr);
  globalLabelsDataOverrides.push_back(jtPtBinsGlobalStr);
  globalLabelsDataOverrides.push_back(gammaJtDPhiStr);

  globalLabelsMCOverrides.push_back(jetRStr);
  globalLabelsMCOverrides.push_back(jtPtBinsGlobalStr);
  globalLabelsMCOverrides.push_back(gammaJtDPhiStr);
   
  configParser labels(inUnfoldFileLabel_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();

  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    if(labelMap.count(globalLabels[gI]) != 0){
      globalLabels[gI] = labelMap[globalLabels[gI]];
    }
  }

  for(unsigned int gI = 0; gI < globalLabelsDataOverrides.size(); ++gI){
    if(labelMap.count(globalLabelsDataOverrides[gI]) != 0){
      globalLabelsDataOverrides[gI] = labelMap[globalLabelsDataOverrides[gI]];
    }
  }

  for(unsigned int gI = 0; gI < globalLabelsMCOverrides.size(); ++gI){
    if(labelMap.count(globalLabelsMCOverrides[gI]) != 0){
      globalLabelsMCOverrides[gI] = labelMap[globalLabelsMCOverrides[gI]];
    }
  }
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Grab first set of parameters defined in necessaryInFileParams
  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  const bool isPP = inUnfoldFileConfig_p->GetValue("ISPP", 0);
  std::vector<std::string> centBinsStr = {"PP"};
  std::vector<std::string> centBinsLabel = {"#bf{p+p}"};

  if(!isPP){
    std::vector<std::string> tempCentBinsStr = strToVect(inUnfoldFileConfig_p->GetValue("CENTBINS", ""));
    centBinsStr.clear();
    centBinsLabel.clear();
    for(unsigned int cI = 0; cI < tempCentBinsStr.size()-1; ++cI){
      centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
      centBinsLabel.push_back("#bf{Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%}");
    }
  }  
  
  const Int_t nIterMax = 100;
  const Int_t nIter = inUnfoldFileConfig_p->GetValue("NITER", 0);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(nIter > nIterMax){
    std::cout << "nIter \'" << nIter << "\' exceeds maximum allowed nIterMax \'" << nIterMax << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nIterForLoop = TMath::Min(nIter, nIterCap);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Define a bunch of parameters for the TCanvas
  const Int_t nPad = 3;
  Double_t padSplit1 = 0.43;
  Double_t padSplit2 = 0.27976;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit2;
  const Double_t topMargin = 0.01;
  const Double_t rightMargin = 0.01;
  const Double_t noMargin = 0.001;
  Double_t height = 900.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit1) - bottomMargin*padSplit2)/(1.0 - leftMargin - rightMargin);

  const Double_t leftMargin2D = 0.12;
  const Double_t bottomMargin2D = leftMargin2D;
  const Double_t topMargin2D = 0.01;
  const Double_t rightMargin2D = 0.12;
  Double_t height2D = 900.0;
  Double_t width2D = height2D*(1.0 - topMargin2D* - bottomMargin2D)/(1.0 - leftMargin2D - rightMargin2D);
  
  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = inUnfoldFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inUnfoldFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inUnfoldFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  std::vector<Float_t> xJMax, xJMin;
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    std::string gammaPtBinsStr = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBins != nGammaPtBinsTemp){
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  gStyle->SetOptStat(0);

  std::vector<std::string> recoMatrixNames, truthMatrixNames;
  std::vector<std::vector<std::string> > recoHistNames, truthHistNames, statsDeltaNames, iterDeltaNames, totalDeltaNames;
  std::vector<std::vector<std::vector<std::string> > > unfoldNames;

  std::vector<std::string> statsDeltaNames2D, iterDeltaNames2D, totalDeltaNames2D;
  std::vector<std::vector<std::string> > unfoldNames2D;

  for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
    statsDeltaNames2D.push_back("statsDelta_GammaPtAll_" + centBinsStr[cI] + "_h");
    iterDeltaNames2D.push_back("iterDelta_GammaPtAll_" + centBinsStr[cI] + "_h");
    totalDeltaNames2D.push_back("totalDelta_GammaPtAll_" + centBinsStr[cI] + "_h");
    
    unfoldNames2D.push_back({});

    for(Int_t i = 0; i < nIter; ++i){
      unfoldNames2D[cI].push_back(("photonPtJet" + varName + "Reco_HalfToUnfold_Iter" + std::to_string(i) + "_" +centBinsStr[cI] + "_h").c_str());
    }
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
    statsDeltaNames2D.push_back("statsDelta_GammaPtAll_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
    iterDeltaNames2D.push_back("iterDelta_GammaPtAll_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
    totalDeltaNames2D.push_back("totalDelta_GammaPtAll_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");

    unfoldNames2D.push_back({});

    for(Int_t i = 0; i < nIter; ++i){
      unfoldNames2D[cI + centBinsStr.size()].push_back(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
    }
  }
  
  for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
    recoMatrixNames.push_back("photonPtJet" + varName + "Reco_HalfResponse_" + centBinsStr[cI] + "_h");
    truthMatrixNames.push_back("photonPtJet" + varName + "Truth_HalfResponse_" + centBinsStr[cI] + "_h");

    recoHistNames.push_back({});
    truthHistNames.push_back({});
    statsDeltaNames.push_back({});
    iterDeltaNames.push_back({});
    totalDeltaNames.push_back({});
    unfoldNames.push_back({});
    
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){    
      recoHistNames[cI].push_back(varNameLower + "Reco_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");
      truthHistNames[cI].push_back(varNameLower + "Truth_HalfToUnfold_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");

      statsDeltaNames[cI].push_back("statsDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");
      iterDeltaNames[cI].push_back("iterDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");
      totalDeltaNames[cI].push_back("totalDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");

      unfoldNames[cI].push_back({});
      for(Int_t uI = 0; uI < nIter; ++uI){
	unfoldNames[cI][gI].push_back(varNameLower + "Reco_HalfToUnfold_Iter" + std::to_string(uI) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_h");
      }
    }    
  }   

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
    Int_t cI2 = cI + centBinsStr.size();
    recoMatrixNames.push_back("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h");
    truthMatrixNames.push_back("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h");

    recoHistNames.push_back({});
    truthHistNames.push_back({});
    statsDeltaNames.push_back({});
    iterDeltaNames.push_back({});
    totalDeltaNames.push_back({});
    unfoldNames.push_back({});
    
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){    
      recoHistNames[cI2].push_back(varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
      truthHistNames[cI2].push_back(varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_TRUTH_COMBINED_h");

      statsDeltaNames[cI2].push_back("statsDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
      iterDeltaNames[cI2].push_back("iterDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
      totalDeltaNames[cI2].push_back("totalDelta_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");

      unfoldNames[cI2].push_back({});
      for(Int_t uI = 0; uI < nIter; ++uI){
	unfoldNames[cI2][gI].push_back(varNameLower + "Reco_Iter" + std::to_string(uI) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h");
      }
    }
  }

  for(unsigned int cI = 0; cI < recoMatrixNames.size(); ++cI){    
    TH2F* recoMatrix_p = (TH2F*)inUnfoldFile_p->Get(recoMatrixNames[cI].c_str());
    TH2F* truthMatrix_p = (TH2F*)inUnfoldFile_p->Get(truthMatrixNames[cI].c_str());

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    std::string saveNameTag = dateStr;
    if(recoMatrixNames[cI].find("_Half") != std::string::npos) saveNameTag = "HalfTest_" + saveTag + "_" + saveNameTag;
    else saveNameTag = "PURCORR_" + saveTag + "_" + saveNameTag;

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    std::vector<TH2F*> matrices_p = {recoMatrix_p, truthMatrix_p};
    for(unsigned int mI = 0; mI < matrices_p.size(); ++mI){
      TCanvas* canv_p = new TCanvas("canv_p", "", width2D, height2D);
      canv_p->SetTopMargin(topMargin2D);
      canv_p->SetBottomMargin(bottomMargin2D);
      canv_p->SetLeftMargin(leftMargin2D);
      canv_p->SetRightMargin(rightMargin2D);
      
      canv_p->cd();

      matrices_p[mI]->DrawCopy("COLZ");
      gPad->SetLogz();
      
      std::string saveName = matrices_p[mI]->GetName();
      saveName.replace(saveName.rfind("_h"), 2, "");
      saveName = "pdfDir/" + dateStr + "/" + saveName + "_Delta" + deltaStr + "_" + saveNameTag + "." + saveExt;
      quietSaveAs(canv_p, saveName);    
      delete canv_p;
    }

    
    
    //2-D analysis of the unfolding results - this should be the overriding result - the 1-D analysis is cross-check + possible systematic
    //    for(int yExclude = -1; yExclude < nGammaPtBins; ++yExclude){
    TH1F* statsDelta2D_p = new TH1F(statsDeltaNames2D[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIter-1, 1.5, nIter+0.5);
    TH1F* iterDelta2D_p = new TH1F(iterDeltaNames2D[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIter-1, 1.5, ((double)nIter)+0.5);
    TH1F* totalDelta2D_p = new TH1F(totalDeltaNames2D[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" +deltaStr + "}}").c_str(), nIter-1, 1.5, ((double)nIter)+0.5);
    TH2F* unfold2D_p[nIterMax];
    
    //If its sub data lets save a nice 2-D version of the input before unfolding
    if(unfoldNames2D[cI][0].find("PURCORR") != std::string::npos){// && yExclude < 0){
      std::string inputName = unfoldNames2D[cI][0];
      
      std::string preName = inputName.substr(0, inputName.find("_Iter"));
      std::string postName = inputName;
      while(postName.find("Iter") != std::string::npos){
	postName.replace(0, postName.find("_")+1, "");
      }
      
      inputName = preName + "_PreUnfold_" + postName;
      
      TH2F* preUnfold_p = (TH2F*)inUnfoldFile_p->Get(inputName.c_str());  
      
      TCanvas* canv_p = new TCanvas("canv_p", "", width2D, height2D);
      canv_p->SetTopMargin(topMargin2D);
      canv_p->SetBottomMargin(bottomMargin2D);
      canv_p->SetLeftMargin(leftMargin2D);
      canv_p->SetRightMargin(rightMargin2D);
      
      canv_p->cd();
      
      preUnfold_p->DrawCopy("COLZ");
      gPad->SetLogz();
      
      std::string saveName = inputName;
      saveName.replace(saveName.rfind("_h"), 2, "");
      saveName = "pdfDir/" + dateStr + "/" + saveName + "_Delta" + deltaStr + "_" + saveNameTag + "." + saveExt;
      quietSaveAs(canv_p, saveName);    
      delete canv_p;
    }
    
    std::vector<TH2F*> unfoldedHists2D_p;
    for(Int_t i = 1; i < nIter; ++i){
      unfold2D_p[i] = (TH2F*)inUnfoldFile_p->Get(unfoldNames2D[cI][i].c_str());      
      
      TCanvas* canv_p = new TCanvas("canv_p", "", width2D, height2D);
      canv_p->SetTopMargin(topMargin2D);
      canv_p->SetBottomMargin(bottomMargin2D);
      canv_p->SetLeftMargin(leftMargin2D);
      canv_p->SetRightMargin(rightMargin2D);
      
      canv_p->cd();
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      unfold2D_p[i]->DrawCopy("COLZ");
      gPad->SetLogz();
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      std::string saveName = unfoldNames2D[cI][i];
      saveName.replace(saveName.rfind("_h"), 2, "");
      saveName = "pdfDir/" + dateStr + "/" + saveName + "_Delta" + deltaStr + "_" + saveNameTag + "." + saveExt;
      quietSaveAs(canv_p, saveName);    
      delete canv_p;
      
      unfoldedHists2D_p.push_back(unfold2D_p[i]);
    }
    
    getIterativeHists2D(unfoldedHists2D_p, statsDelta2D_p, iterDelta2D_p, totalDelta2D_p, doRelativeTerm);
    
    std::vector<std::string> tempLabels = globalLabels;
    if(!isMC && truthHistNames[cI][0].find("PURCORR") != std::string::npos){
      tempLabels = globalLabelsDataOverrides;
    }
    else if(!isMC && truthHistNames[cI][0].find("HalfTo") != std::string::npos){
      tempLabels = globalLabelsMCOverrides;
    }
    tempLabels.push_back("#bf{" + prettyString(gammaPtBins[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false) + "}");
    if(!isPP) tempLabels.push_back(centBinsLabel[cI%centBinsLabel.size()]);
      
    TLatex* label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(titleFont);
    label_p->SetTextSize(titleSize);
    label_p->SetTextColor(1);
    if(diagLabelAlignRightGlobal) label_p->SetTextAlign(31);    
    
    Int_t maxLegX = -1;
    std::vector<TLegend*> leg_p;
    std::vector<TH1F*> histsForLeg_p = {iterDelta2D_p, statsDelta2D_p, totalDelta2D_p};
    std::vector<std::string> stringsForLeg_p = {"Iter. Change", "Statistical Unc.", "Quad. Sum"};
    for(unsigned int lI = 0; lI < stringsForLeg_p.size(); ++lI){
      if(maxLegX < (Int_t)stringsForLeg_p[lI].size()) maxLegX = (Int_t)stringsForLeg_p[lI].size();  
    }
    
    Int_t numberHists = histsForLeg_p.size();
    Int_t histsPerLeg = 0;
    
    Int_t totalCustom = 0;
    for(unsigned int hI = 0; hI < diagLegNGlobal.size(); ++hI){
      totalCustom += diagLegNGlobal[hI];
    }
    if(totalCustom < numberHists){
      if(totalCustom != 0) std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
      diagLegNGlobal.clear();
    }
    
    if(diagLegNGlobal.size() == 0){
      while(histsPerLeg*(Int_t)diagLegXGlobal.size() < numberHists){
	++histsPerLeg;
      }
    }
    
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(unsigned int lI = 0; lI < diagLegXGlobal.size(); ++lI){
      leg_p.push_back(nullptr);
      
      if(diagLegNGlobal.size() == 0){
	if(lI+1 < diagLegXGlobal.size()) leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - histsPerLeg*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	else{
	  if(numberHists%histsPerLeg == 0) leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - histsPerLeg*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	  else leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - (numberHists%histsPerLeg)*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	}	
      }
      else{
	leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - diagLegNGlobal[lI]*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
      }
	
      leg_p[lI]->SetTextFont(titleFont);
      leg_p[lI]->SetTextSize(titleSize);
      leg_p[lI]->SetBorderSize(0);
      leg_p[lI]->SetFillColor(0);
      leg_p[lI]->SetFillStyle(0);
    }
    
    TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
    canv_p->SetTopMargin(topMargin);
    canv_p->SetRightMargin(rightMargin);
    canv_p->SetLeftMargin(leftMargin);
    canv_p->SetBottomMargin(leftMargin); //Because of the padsplit adjustment
    canv_p->cd();
    
    HIJet::Style::EquipHistogram(iterDelta2D_p, 0);
    HIJet::Style::EquipHistogram(statsDelta2D_p, 1);
    HIJet::Style::EquipHistogram(totalDelta2D_p, 2);
    totalDelta2D_p->DrawCopy("HIST E1 P");
    statsDelta2D_p->DrawCopy("HIST E1 P SAME");
    iterDelta2D_p->DrawCopy("HIST E1 P SAME");
    totalDelta2D_p->DrawCopy("HIST E1 P SAME");    
    
    if(diagLegNGlobal.size() == 0){
      for(unsigned int lI = 0; lI < histsForLeg_p.size(); ++lI){      
	leg_p[lI/histsPerLeg]->AddEntry(histsForLeg_p[lI], stringsForLeg_p[lI].c_str(), "P");
      }
    }
    else{
      Int_t legCounter = 0;
      for(unsigned int lI = 0; lI < diagLegNGlobal.size(); ++lI){
	for(unsigned int dI = 0; dI < diagLegNGlobal[lI]; ++dI){	
	  leg_p[lI]->AddEntry(histsForLeg_p[legCounter], stringsForLeg_p[legCounter].c_str(), "P");
	  ++legCounter;
	}
      }
    }
    
    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){leg_p[lI]->Draw("SAME");}
    for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){label_p->DrawLatex(diagLabelXGlobal, diagLabelYGlobal - tI*0.06, tempLabels[tI].c_str());}
    quietSaveAs(canv_p, "pdfDir/" + dateStr + "/diag" + varName + "_GammaPtAll_" + centBinsStr[cI%centBinsStr.size()] + "_Delta" + deltaStr + "_" + saveNameTag +  "." + saveExt);
    delete canv_p;
    delete label_p;
    
    std::string configStr2D = varNameUpper + "_GAMMAPTALL_" + centBinsStr[cI%centBinsStr.size()];
    int termPos2D = -1;
    for(Int_t i = 1; i < nIter; ++i){
      Double_t totalDelta = totalDelta2D_p->GetBinContent(i);
      Double_t statsDelta = statsDelta2D_p->GetBinContent(i);      
      if(statsDelta/totalDelta  >= termThresh){
	termPos2D = i+1;
	break;
      }
    }
    outUnfoldConfig_p->SetValue(configStr2D.c_str(), termPos2D);
    
    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){     
      delete leg_p[lI];
    }
    leg_p.clear();
    
    delete statsDelta2D_p;
    delete iterDelta2D_p;
    delete totalDelta2D_p;
  
    
       
    //1-D analysis of the unfolding results
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){    
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      TH1F* truth_p = nullptr;
      if(truthHistNames[cI][gI].find("TRUTH_COMBINED") == std::string::npos || isMC){
	truth_p = (TH1F*)inUnfoldFile_p->Get(truthHistNames[cI][gI].c_str());
	std::cout << "NAME Passed: " << truthHistNames[cI][gI] << std::endl;
	std::cout << " " << recoHistNames[cI][gI] << std::endl;
      }
      else if(!isMC){
	std::cout << "NEW NAME: " << truthHistNames[cI][gI] << std::endl;						   
	truth_p = (TH1F*)inUnfoldFile_p->Get(truthHistNames[cI][gI].c_str());
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      TH1F* reco_p = (TH1F*)inUnfoldFile_p->Get(recoHistNames[cI][gI].c_str());
      TH1F* unfold_p[nIterMax];
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    

      TH1F* statsDelta_p = new TH1F(statsDeltaNames[cI][gI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIter-1, 1.5, nIter+0.5);
      TH1F* iterDelta_p = new TH1F(iterDeltaNames[cI][gI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIter-1, 1.5, ((double)nIter)+0.5);
      TH1F* totalDelta_p = new TH1F(totalDeltaNames[cI][gI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIter-1, 1.5, ((double)nIter)+0.5);

      std::vector<TH1F*> unfoldedHists_p;
      for(Int_t i = 1; i < nIter; ++i){
	unfold_p[i] = (TH1F*)inUnfoldFile_p->Get(unfoldNames[cI][gI][i].c_str());
	unfoldedHists_p.push_back(unfold_p[i]);
      }
    
      getIterativeHists(unfoldedHists_p, statsDelta_p, iterDelta_p, totalDelta_p, doRelativeTerm);
      
      tempLabels = globalLabels;
      if(!isMC && truthHistNames[cI][gI].find("PURCORR") != std::string::npos){
	tempLabels = globalLabelsDataOverrides;
      }
      else if(!isMC && truthHistNames[cI][0].find("HalfTo") != std::string::npos){
	tempLabels = globalLabelsMCOverrides;
      }
      
      tempLabels.push_back("#bf{" + prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false) + "}");
      if(!isPP) tempLabels.push_back(centBinsLabel[cI%centBinsLabel.size()]);
    
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(titleFont);
      label_p->SetTextSize(titleSize);
      label_p->SetTextColor(1);
      if(diagLabelAlignRightGlobal) label_p->SetTextAlign(31);
      
      Int_t maxLegX = -1;
      histsForLeg_p = {iterDelta_p, statsDelta_p, totalDelta_p};
      stringsForLeg_p = {"Iter. Change", "Statistical Unc.", "Quad. Sum"};
      for(unsigned int lI = 0; lI < stringsForLeg_p.size(); ++lI){
	if(maxLegX < (Int_t)stringsForLeg_p[lI].size()) maxLegX = (Int_t)stringsForLeg_p[lI].size();
      }
      numberHists = histsForLeg_p.size();
      histsPerLeg = 0;

      totalCustom = 0;
      for(unsigned int hI = 0; hI < diagLegNGlobal.size(); ++hI){
	totalCustom += diagLegNGlobal[hI];
      }
      if(totalCustom < numberHists){
	if(totalCustom != 0) std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
	diagLegNGlobal.clear();
      }
      
      if(diagLegNGlobal.size() == 0){
	while(histsPerLeg*(Int_t)diagLegXGlobal.size() < numberHists){
	  ++histsPerLeg;
	}
      }

      for(unsigned int lI = 0; lI < diagLegXGlobal.size(); ++lI){
	leg_p.push_back(nullptr);
	
	if(diagLegNGlobal.size() == 0){
	  if(lI+1 < diagLegXGlobal.size()) leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - histsPerLeg*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	  else{
	    if(numberHists%histsPerLeg == 0) leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - histsPerLeg*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	    else leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - (numberHists%histsPerLeg)*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	  }	
	}
	else{
	  leg_p[lI] = new TLegend(diagLegXGlobal[lI], diagLegYGlobal[lI] - diagLegNGlobal[lI]*0.04, diagLegXGlobal[lI] + maxLegX*0.02, diagLegYGlobal[lI]);
	}
      
	leg_p[lI]->SetTextFont(titleFont);
	leg_p[lI]->SetTextSize(titleSize);
	leg_p[lI]->SetBorderSize(0);
	leg_p[lI]->SetFillColor(0);
	leg_p[lI]->SetFillStyle(0);
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
      canv_p->SetTopMargin(topMargin);
      canv_p->SetRightMargin(rightMargin);
      canv_p->SetLeftMargin(leftMargin);
      canv_p->SetBottomMargin(leftMargin); //Because of the padsplit adjustment
      canv_p->cd();

      HIJet::Style::EquipHistogram(iterDelta_p, 0);
      HIJet::Style::EquipHistogram(statsDelta_p, 1);
      HIJet::Style::EquipHistogram(totalDelta_p, 2);
      totalDelta_p->DrawCopy("HIST E1 P");
      statsDelta_p->DrawCopy("HIST E1 P SAME");
      iterDelta_p->DrawCopy("HIST E1 P SAME");
      totalDelta_p->DrawCopy("HIST E1 P SAME");    
      
      if(diagLegNGlobal.size() == 0){
	for(unsigned int lI = 0; lI < histsForLeg_p.size(); ++lI){      
	  leg_p[lI/histsPerLeg]->AddEntry(histsForLeg_p[lI], stringsForLeg_p[lI].c_str(), "P");
	}
      }
      else{
	Int_t legCounter = 0;
	for(unsigned int lI = 0; lI < diagLegNGlobal.size(); ++lI){
	  for(unsigned int dI = 0; dI < diagLegNGlobal[lI]; ++dI){	
	    leg_p[lI]->AddEntry(histsForLeg_p[legCounter], stringsForLeg_p[legCounter].c_str(), "P");
	    ++legCounter;
	  }
	}
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	leg_p[lI]->Draw("SAME");
      }
      
      for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	label_p->DrawLatex(diagLabelXGlobal, diagLabelYGlobal - tI*0.06, tempLabels[tI].c_str());
      }

      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/diag" + varName + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI%centBinsStr.size()] + "_Delta" + deltaStr + "_" + saveNameTag +  "." + saveExt);
      delete canv_p;

      //Write termination to a config
      std::string configStr = varNameUpper + "_GAMMAPT" + std::to_string(gI) + "_" + centBinsStr[cI%centBinsStr.size()];
      int termPos = -1;
      for(Int_t i = 1; i < nIter; ++i){
	Double_t totalDelta = totalDelta_p->GetBinContent(i);
	Double_t statsDelta = statsDelta_p->GetBinContent(i);
	  
	if(statsDelta/totalDelta  >= termThresh){
	  termPos = i+1;
	  break;
	}
      }      
      outUnfoldConfig_p->SetValue(configStr.c_str(), termPos);
          
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	delete leg_p[lI];
      }
      leg_p.clear();
      
      Float_t xjMinTemp = config_p->GetValue(("XJMIN." + std::to_string(gI)).c_str(), -1000.);
      Float_t xjMaxTemp = config_p->GetValue(("XJMAX." + std::to_string(gI)).c_str(), -1000.);
      if(truth_p == nullptr){
	Float_t xjMinTemp2 = config_p->GetValue(("XJMIN." + std::to_string(gI) + ".DATAOVERRIDE").c_str(), -1000.);
	Float_t xjMaxTemp2 = config_p->GetValue(("XJMAX." + std::to_string(gI) + ".DATAOVERRIDE").c_str(), -1000.);
	if(xjMinTemp2 > 0) xjMinTemp = xjMinTemp2;
	if(xjMaxTemp2 > 0) xjMaxTemp = xjMaxTemp2;
      }

      if(xjMaxTemp < 0){
	xjMinTemp = xjMinGlobal;
	xjMaxTemp = xjMaxGlobal;	
      }
            
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      canv_p = new TCanvas("canv_p", "", width, height);
      TCanvas* canvBest_p = new TCanvas("canvBest_p", "", width, height);
      TPad* pads_p[nPad];
      TPad* padsBest_p[nPad];
      canv_p->SetTopMargin(noMargin);
      canv_p->SetRightMargin(noMargin);
      canv_p->SetLeftMargin(noMargin);
      canv_p->SetBottomMargin(noMargin);
      canv_p->cd();
      
      canvBest_p->SetTopMargin(noMargin);
      canvBest_p->SetRightMargin(noMargin);
      canvBest_p->SetLeftMargin(noMargin);
      canvBest_p->SetBottomMargin(noMargin);
      canvBest_p->cd();
      
      pads_p[0] = new TPad("pad0", "", 0.0, padSplit1, 1.0, 1.0);
      pads_p[0]->SetTopMargin(topMargin);
      pads_p[0]->SetRightMargin(rightMargin);
      pads_p[0]->SetLeftMargin(leftMargin);
      pads_p[0]->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      pads_p[0]->cd();
      canv_p->cd();

      padsBest_p[0] = new TPad("padBest0", "", 0.0, padSplit1, 1.0, 1.0);
      padsBest_p[0]->SetTopMargin(topMargin);
      padsBest_p[0]->SetRightMargin(rightMargin);
      padsBest_p[0]->SetLeftMargin(leftMargin);
      padsBest_p[0]->SetBottomMargin(noMargin);
      
      canvBest_p->cd();
      padsBest_p[0]->Draw("SAME");
      padsBest_p[0]->cd();
      canvBest_p->cd();

      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", 0.0, padSplit2, 1.0, padSplit1);
      pads_p[1]->SetTopMargin(noMargin);
      pads_p[1]->SetRightMargin(rightMargin);
      pads_p[1]->SetLeftMargin(leftMargin);
      pads_p[1]->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[1]->Draw("SAME");
      pads_p[1]->cd(); 
      canv_p->cd();

      canvBest_p->cd();
      padsBest_p[1] = new TPad("padBest1", "", 0.0, padSplit2, 1.0, padSplit1);
      padsBest_p[1]->SetTopMargin(noMargin);
      padsBest_p[1]->SetRightMargin(rightMargin);
      padsBest_p[1]->SetLeftMargin(leftMargin);
      padsBest_p[1]->SetBottomMargin(noMargin);
      
      canvBest_p->cd();
      padsBest_p[1]->Draw("SAME");
      padsBest_p[1]->cd();
      canvBest_p->cd();

      
      canv_p->cd();
      pads_p[2] = new TPad("pad2", "", 0.0, 0.0, 1.0, padSplit2);      
      pads_p[2]->SetTopMargin(noMargin);
      pads_p[2]->SetRightMargin(rightMargin);
      pads_p[2]->SetLeftMargin(leftMargin);
      pads_p[2]->SetBottomMargin(bottomMargin);
      
      canv_p->cd();
      pads_p[2]->Draw("SAME");
      pads_p[2]->cd(); 
      canv_p->cd();
            
      canvBest_p->cd();
      padsBest_p[2] = new TPad("padBest2", "", 0.0, 0.0, 1.0, padSplit2);      
      padsBest_p[2]->SetTopMargin(noMargin);
      padsBest_p[2]->SetRightMargin(rightMargin);
      padsBest_p[2]->SetLeftMargin(leftMargin);
      padsBest_p[2]->SetBottomMargin(bottomMargin);
      
      canvBest_p->cd();
      padsBest_p[2]->Draw("SAME");
      padsBest_p[2]->cd();
      canvBest_p->cd();
      
      canv_p->cd();
      pads_p[0]->cd();

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(truth_p != nullptr) HIJet::Style::EquipHistogram(truth_p, 2); //TRUTH IS kAzure-3 lines (atlas style picked up here)
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      HIJet::Style::EquipHistogram(reco_p, 1); //RECO IS kRed-4 lines (atlas style picked up here)
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      reco_p->SetMarkerSize(0.00001);
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      if(truth_p != nullptr) truth_p->SetMarkerSize(0.00001);    
      
      maxLegX = -1;
      //    Int_t nLeg = nIterForLoop + 2 - 1;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(titleFont);
      label_p->SetTextSize(titleSize/(1.0 - padSplit1));
      label_p->SetTextColor(1);
      if(xjLabelAlignRightGlobal) label_p->SetTextAlign(31);

      std::vector<std::string> legFillStr;
      if(truth_p != nullptr){
	stringsForLeg_p = {"Truth", "Reco."};
	histsForLeg_p = {truth_p, reco_p};    
	legFillStr = {"L", "L"};
      }
      else{
	stringsForLeg_p = {"Reco."};
	histsForLeg_p = {reco_p};    
	legFillStr = {"L"};
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      for(Int_t i = 1; i < nIterForLoop; ++i){
	stringsForLeg_p.push_back("Unfolded, " + std::to_string(i));
	histsForLeg_p.push_back(unfold_p[i]);
	legFillStr.push_back("P L");
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      for(unsigned int lI = 0; lI < stringsForLeg_p.size(); ++lI){
	if(maxLegX < (Int_t)stringsForLeg_p[lI].size()) maxLegX = (Int_t)stringsForLeg_p[lI].size();      
      }
      numberHists = histsForLeg_p.size();
      histsPerLeg = 0;
      
      totalCustom = 0;
      for(unsigned int hI = 0; hI < xjLegNGlobal.size(); ++hI){
	totalCustom += xjLegNGlobal[hI];
      }
      if(totalCustom < numberHists){
	if(totalCustom != 0) std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
	xjLegNGlobal.clear();
      }
      
      if(xjLegNGlobal.size() == 0){
	while(histsPerLeg*(Int_t)xjLegXGlobal.size() < numberHists){
	  ++histsPerLeg;
	}
      }
      
      for(unsigned int lI = 0; lI < xjLegXGlobal.size(); ++lI){
	leg_p.push_back(nullptr);
	if(xjLegNGlobal.size() == 0){
	  if(lI+1 < xjLegXGlobal.size()) leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - histsPerLeg*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);
	  else{
	    if(numberHists%histsPerLeg == 0) leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - histsPerLeg*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);
	    else leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - (numberHists%histsPerLeg)*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);
	  }
	}
	else{
	  double nVert = xjLegNGlobal[lI];
	  if(truth_p == nullptr && lI == xjLegXGlobal.size()-1) nVert -= 1.0;

	  leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - nVert*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);	
	}
	
	leg_p[lI]->SetTextFont(titleFont);
	leg_p[lI]->SetTextSize(titleSize/(1.0 - padSplit1));
	leg_p[lI]->SetBorderSize(0);
	leg_p[lI]->SetFillColor(0);
	leg_p[lI]->SetFillStyle(0);
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      binWidthAndSelfNorm(reco_p);

      reco_p->SetMinimum(xjMinTemp);
      reco_p->SetMaximum(xjMaxTemp);

      reco_p->GetXaxis()->SetTitleFont(titleFont);
      reco_p->GetYaxis()->SetTitleFont(titleFont);
      reco_p->GetXaxis()->SetTitleSize(0.00001);
      reco_p->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit1));
      
      reco_p->GetXaxis()->SetLabelFont(titleFont);
      reco_p->GetYaxis()->SetLabelFont(titleFont);
      reco_p->GetXaxis()->SetLabelSize(0.00001);
      reco_p->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit1));
      
      reco_p->GetYaxis()->SetTitleOffset(yOffset);
      
      reco_p->GetYaxis()->SetNdivisions(505);
      reco_p->GetXaxis()->SetNdivisions(505);
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      reco_p->GetYaxis()->SetTitle(("#frac{1}{N_{" + varNameLabel + "}} #frac{dN}{d" + varNameLabel + "}").c_str());
      reco_p->DrawCopy("HIST E1");
      if(doLogX) gPad->SetLogx();
      if(doLogY) gPad->SetLogy();

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(truth_p != nullptr){
	binWidthAndSelfNorm(truth_p);      
	truth_p->DrawCopy("HIST E1 SAME");
      }

      for(Int_t i = 1; i < nIter; ++i){
	HIJet::Style::EquipHistogram(unfold_p[i], i-1);
	binWidthAndSelfNorm(unfold_p[i]);
      }
      
      for(Int_t i = 1; i < nIterForLoop; ++i){
	unfold_p[i]->DrawCopy("HIST E1 P SAME");
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      canvBest_p->cd();
      padsBest_p[0]->cd();

      reco_p->DrawCopy("HIST E1");
      if(truth_p != nullptr) truth_p->DrawCopy("HIST E1 SAME");

      if(doLogX) gPad->SetLogx();
      if(doLogY) gPad->SetLogy();

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      std::vector<double> deltaVals;
      for(Int_t i = 1; i < nIter; ++i){
	HIJet::Style::EquipHistogram(unfold_p[i], i-1);
	if(termPos2D == i+1 || termPos == i+1){
	  unfold_p[i]->DrawCopy("HIST E1 P SAME");
	  if(termPos2D == i+1){
	    Int_t tempColor = HIJet::Style::GetColor(i-1);
	    Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);	                

	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << i << "/" << nIter << std::endl;
	    
	    TBox* box_p = new TBox();
	    box_p->SetFillColorAlpha(tempColor, tempOpacity);
	    for(Int_t bIX = 0; bIX < unfold_p[i]->GetNbinsX(); ++bIX){
	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << i << "/" << nIter << ", " << bIX << "/" << unfold_p[i]->GetNbinsX() << std::endl;
	    	      
	      Double_t x1 = unfold_p[i]->GetBinLowEdge(bIX+1);
	      Double_t x2 = unfold_p[i]->GetBinLowEdge(bIX+2);

	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	      
	      Double_t yDelta1 = TMath::Abs(unfold_p[i+1]->GetBinContent(bIX+1) - unfold_p[i]->GetBinContent(bIX+1));
	      Double_t yDelta2 = 0.0;
	      if(i > 1) yDelta2 = TMath::Abs(unfold_p[i-1]->GetBinContent(bIX+1) - unfold_p[i]->GetBinContent(bIX+1));

	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	      	      
	      Double_t yLow = unfold_p[i]->GetBinContent(bIX+1) - TMath::Max(yDelta1, yDelta2);
	      Double_t yHigh = unfold_p[i]->GetBinContent(bIX+1) + TMath::Max(yDelta1, yDelta2);

	      if(unfold_p[i]->GetBinContent(bIX+1) < TMath::Power(10,-50)) deltaVals.push_back(0.0);
	      else deltaVals.push_back(TMath::Max(yDelta1, yDelta2)/(double)unfold_p[i]->GetBinContent(bIX+1));

	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	      
	      box_p->DrawBox(x1, yLow, x2, yHigh);
	    }
	    
	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    delete box_p;
	  }
	}
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      canv_p->cd();
      pads_p[0]->cd();
      
      if(xjLegNGlobal.size() == 0){
	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	for(unsigned int lI = 0; lI < histsForLeg_p.size(); ++lI){
	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << lI << ", " << stringsForLeg_p[lI] << std::endl;

	  leg_p[lI/histsPerLeg]->AddEntry(histsForLeg_p[lI], stringsForLeg_p[lI].c_str(), legFillStr[lI].c_str());
	}
      }
      else{
	unsigned int legCounter = 0;
	for(unsigned int lI = 0; lI < xjLegNGlobal.size(); ++lI){
	  for(unsigned int dI = 0; dI < xjLegNGlobal[lI]; ++dI){	
	    if(histsForLeg_p.size() == legCounter) break;
	    leg_p[lI]->AddEntry(histsForLeg_p[legCounter], stringsForLeg_p[legCounter].c_str(), legFillStr[legCounter].c_str());
	    ++legCounter;
	  }
	}
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
      
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	leg_p[lI]->Draw("SAME");
      }
      
      for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	label_p->DrawLatex(xjLabelXGlobal, xjLabelYGlobal - tI*0.083, tempLabels[tI].c_str());
      }    

      canvBest_p->cd();
      padsBest_p[0]->cd();

      for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	label_p->DrawLatex(xjLabelXGlobal, xjLabelYGlobal - tI*0.083, tempLabels[tI].c_str());
      }    
      
      canv_p->cd();
      pads_p[1]->cd();
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      Float_t xjRatMinTemp = config_p->GetValue(("XJRATMIN." + std::to_string(gI)).c_str(), -1000.);
      Float_t xjRatMaxTemp = config_p->GetValue(("XJRATMAX." + std::to_string(gI)).c_str(), -1000.);
      if(xjRatMaxTemp < 0){
	xjRatMinTemp = xjRatMinGlobal;
	xjRatMaxTemp = xjRatMaxGlobal;	
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t i = 1; i < nIterForLoop; ++i){
	TH1F* unfoldClone_p = (TH1F*)unfold_p[i]->Clone("tempClone");	
	unfoldClone_p->Divide(reco_p);
	
	if(i == 1){
	  unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Reco.}");
	  unfoldClone_p->SetMinimum(xjRatMinTemp);
	  unfoldClone_p->SetMaximum(xjRatMaxTemp);
	  
	  unfoldClone_p->GetXaxis()->SetTitleOffset(1.2);
	  unfoldClone_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit1 - padSplit2)/(1.0-padSplit1));
	  
	  unfoldClone_p->GetXaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  unfoldClone_p->GetYaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  
	  unfoldClone_p->GetXaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	  unfoldClone_p->GetYaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	  
	  unfoldClone_p->GetYaxis()->SetNdivisions(404);	
	  
	  unfoldClone_p->DrawCopy("HIST E1 P");
	  if(doLogX) gPad->SetLogx();
	}
	else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	delete unfoldClone_p;
      }
    
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;
      canvBest_p->cd();
      padsBest_p[1]->cd();
    
      bool hasDrawn = false;
      for(Int_t i = 1; i < nIter; ++i){
	if(i+1 != termPos && i+1 != termPos2D) continue;

	TH1F* unfoldClone_p = (TH1F*)unfold_p[i]->Clone("tempClone");    
	unfoldClone_p->Divide(reco_p);
	
	if(!hasDrawn){
	  unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Reco.}");
	  unfoldClone_p->SetMinimum(xjRatMinTemp);
	  unfoldClone_p->SetMaximum(xjRatMaxTemp);
	  
	  unfoldClone_p->GetXaxis()->SetTitleOffset(1.2);
	  unfoldClone_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit1 - padSplit2)/(1.0-padSplit1));
	  
	  unfoldClone_p->GetXaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  unfoldClone_p->GetYaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  
	  unfoldClone_p->GetXaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	  unfoldClone_p->GetYaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	  
	  unfoldClone_p->GetYaxis()->SetNdivisions(404);	
	  
	  unfoldClone_p->DrawCopy("HIST E1 P");	 
	  if(doLogX) gPad->SetLogx();
	  hasDrawn = true;
	}
	else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	if(i+1 == termPos2D){
	  Int_t tempColor = HIJet::Style::GetColor(i-1);
	  Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);
	  
	  TBox* box_p = new TBox();
	  box_p->SetFillColorAlpha(tempColor, tempOpacity);

	  for(Int_t bIX = 0; bIX < unfoldClone_p->GetXaxis()->GetNbins(); ++bIX){
	    Double_t x1 = unfoldClone_p->GetBinLowEdge(bIX+1);
	    Double_t x2 = unfoldClone_p->GetBinLowEdge(bIX+2);

	    box_p->DrawBox(x1, 1.0 - deltaVals[bIX], x2, 1.0 + deltaVals[bIX]);
	  }

	  delete box_p;
	}

	delete unfoldClone_p;
      }

      //Third panel, divide by truth
      canvBest_p->cd();
      pads_p[2]->cd();
      
      for(Int_t i = 1; i < nIterForLoop; ++i){
	TH1F* unfoldClone_p = (TH1F*)unfold_p[i]->Clone("tempClone");	
	unfoldClone_p->Divide(truth_p);
	
	if(i == 1){
	  if(isMC) unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Truth}");
	  else unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded Data}{Truth PYTHIA}");
	  unfoldClone_p->SetMinimum(xjRatMinTemp);
	  unfoldClone_p->SetMaximum(xjRatMaxTemp);
	  
	  unfoldClone_p->GetXaxis()->SetTitleOffset(1.2);
	  unfoldClone_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit2)/(1.0-padSplit1));
	  
	  unfoldClone_p->GetXaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetTitleSize(titleSize/(padSplit2));
	  unfoldClone_p->GetYaxis()->SetTitleSize(titleSize/(padSplit2));
	  
	  unfoldClone_p->GetXaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetLabelSize(labelSize/(padSplit2));
	  unfoldClone_p->GetYaxis()->SetLabelSize(labelSize/(padSplit2));
	  
	  unfoldClone_p->GetYaxis()->SetNdivisions(404);	
	  
	  unfoldClone_p->DrawCopy("HIST E1 P");
	  if(doLogX) gPad->SetLogx();
	}
	else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	delete unfoldClone_p;
      }
    
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;
      canvBest_p->cd();
      padsBest_p[2]->cd();
    
      hasDrawn = false;
      for(Int_t i = 1; i < nIter; ++i){
	if(i+1 != termPos && i+1 != termPos2D) continue;

	TH1F* unfoldClone_p = (TH1F*)unfold_p[i]->Clone("tempClone");    
	unfoldClone_p->Divide(truth_p);
	
	if(!hasDrawn){
	  if(isMC) unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Truth}");
	  else unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded Data}{Truth PYTHIA}");

	  unfoldClone_p->SetMinimum(xjRatMinTemp);
	  unfoldClone_p->SetMaximum(xjRatMaxTemp);
	  
	  unfoldClone_p->GetXaxis()->SetTitleOffset(1.2);
	  unfoldClone_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit2)/(1.0-padSplit1));
	  
	  unfoldClone_p->GetXaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetTitleFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetTitleSize(titleSize/(padSplit2));
	  unfoldClone_p->GetYaxis()->SetTitleSize(titleSize/(padSplit2));
	  
	  unfoldClone_p->GetXaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetYaxis()->SetLabelFont(titleFont);
	  unfoldClone_p->GetXaxis()->SetLabelSize(labelSize/(padSplit2));
	  unfoldClone_p->GetYaxis()->SetLabelSize(labelSize/(padSplit2));
	  
	  unfoldClone_p->GetYaxis()->SetNdivisions(404);	
	  
	  unfoldClone_p->DrawCopy("HIST E1 P");	 
	  if(doLogX) gPad->SetLogx();
	  hasDrawn = true;
	}
	else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	if(i+1 == termPos2D){
	  Int_t tempColor = HIJet::Style::GetColor(i-1);
	  Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);
	  
	  TBox* box_p = new TBox();
	  box_p->SetFillColorAlpha(tempColor, tempOpacity);

	  for(Int_t bIX = 0; bIX < unfoldClone_p->GetXaxis()->GetNbins(); ++bIX){
	    Double_t x1 = unfoldClone_p->GetBinLowEdge(bIX+1);
	    Double_t x2 = unfoldClone_p->GetBinLowEdge(bIX+2);

	    box_p->DrawBox(x1, 1.0 - deltaVals[bIX], x2, 1.0 + deltaVals[bIX]);
	  }

	  delete box_p;
	}

	delete unfoldClone_p;
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;      
          
      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(reco_p->GetXaxis()->GetBinLowEdge(1), 1.0, reco_p->GetXaxis()->GetBinLowEdge(reco_p->GetXaxis()->GetNbins()+1), 1.0);

      canv_p->cd();
      pads_p[2]->cd();
      line_p->DrawLine(reco_p->GetXaxis()->GetBinLowEdge(1), 1.0, reco_p->GetXaxis()->GetBinLowEdge(reco_p->GetXaxis()->GetNbins()+1), 1.0);      
      pads_p[1]->cd();
      line_p->DrawLine(reco_p->GetXaxis()->GetBinLowEdge(1), 1.0, reco_p->GetXaxis()->GetBinLowEdge(reco_p->GetXaxis()->GetNbins()+1), 1.0);      
      padsBest_p[1]->cd();
      line_p->DrawLine(reco_p->GetXaxis()->GetBinLowEdge(1), 1.0, reco_p->GetXaxis()->GetBinLowEdge(reco_p->GetXaxis()->GetNbins()+1), 1.0);      


      delete line_p;
      delete label_p;
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;      
    
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/unfold" + varName + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI%centBinsStr.size()] + "_Delta" + deltaStr + "_" + saveNameTag + "." + saveExt);

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;
      
      for(Int_t i = 0; i < nPad; ++i){
	delete pads_p[i];
      }
      delete canv_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", gI: " << gI << "/" << nGammaPtBins << std::endl;
      
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	delete leg_p[lI];      
      }
      leg_p.clear();      

      Int_t nBestHist = 4;
      if(truth_p == nullptr) nBestHist -= 1.0;
      if(termPos2D < 0 && termPos < 0) nBestHist -= 2.0;
      else if(termPos2D < 0 || termPos < 0) nBestHist -= 1.0;
      else if(termPos2D == termPos) nBestHist -= 1.0;      

      Int_t lastLegPos = -1;
      Int_t totalHist = 0;
      for(unsigned int lI = 0; lI < xjLegXGlobal.size(); ++lI){
	leg_p.push_back(nullptr);
	if(xjLegNGlobal.size() == 0){
	  int newHistsPerLeg = histsPerLeg;
	  if(newHistsPerLeg > nBestHist - totalHist) newHistsPerLeg = nBestHist - totalHist;

	  leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - newHistsPerLeg*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);

	  totalHist += newHistsPerLeg;
	}
	else{ 
	  double nVert = xjLegNGlobal[lI];
	  if(nVert > nBestHist - totalHist) nVert = nBestHist - totalHist;
	  
	  leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - nVert*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);

	  totalHist += nVert;
	}
	
	leg_p[lI]->SetTextFont(titleFont);
	leg_p[lI]->SetTextSize(titleSize/(1.0 - padSplit1));
	leg_p[lI]->SetBorderSize(0);
	leg_p[lI]->SetFillColor(0);
	leg_p[lI]->SetFillStyle(0);	

	lastLegPos = lI;
	if(nBestHist == totalHist) break;
      }

      histsForLeg_p.clear();
      stringsForLeg_p.clear();
      legFillStr.clear();
      if(truth_p != nullptr){
	histsForLeg_p.push_back(truth_p);
	stringsForLeg_p.push_back("Truth");
	legFillStr.push_back("L");
      }
      histsForLeg_p.push_back(reco_p);
      stringsForLeg_p.push_back("Reco.");
      legFillStr.push_back("L");

      for(Int_t i = 1; i < nIter; ++i){
	if(i+1 == termPos2D){
	  histsForLeg_p.push_back(unfold_p[i]);
	  if(i+1 == termPos) stringsForLeg_p.push_back("Iter. " + std::to_string(termPos2D) + ", Best 1 and 2D");
	  else stringsForLeg_p.push_back("Iter. " + std::to_string(termPos2D) + ", Best 2D");

	  legFillStr.push_back("L P");
	}
	else if(i+1 == termPos){
	  histsForLeg_p.push_back(unfold_p[i]);
	  stringsForLeg_p.push_back("Iter. " + std::to_string(termPos) + ", Best 1D");
	  legFillStr.push_back("L P");
	}
      }

      if(xjLegNGlobal.size() == 0){
	for(unsigned int lI = 0; lI < histsForLeg_p.size(); ++lI){
	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << lI << ", " << stringsForLeg_p[lI] << std::endl;

	  leg_p[lI/histsPerLeg]->AddEntry(histsForLeg_p[lI], stringsForLeg_p[lI].c_str(), legFillStr[lI].c_str());
	}
      }
      else{
	unsigned int legCounter = 0;
	for(unsigned int lI = 0; lI < xjLegNGlobal.size(); ++lI){
	  for(unsigned int dI = 0; dI < xjLegNGlobal[lI]; ++dI){	
	    if(histsForLeg_p.size() == legCounter) break;
	    leg_p[lI]->AddEntry(histsForLeg_p[legCounter], stringsForLeg_p[legCounter].c_str(), legFillStr[legCounter].c_str());
	    ++legCounter;
	  }

	  if(((int)lI) == lastLegPos) break;
	}
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    

      canvBest_p->cd();
      padsBest_p[0]->cd();
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	leg_p[lI]->Draw("SAME");
      }      
      
      quietSaveAs(canvBest_p, "pdfDir/" + dateStr + "/unfold" + varName + "_BESTIterations_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI%centBinsStr.size()] + "_Delta" + deltaStr + "_" + saveNameTag + "." + saveExt);

      for(Int_t i = 0; i < nPad; ++i){
	delete padsBest_p[i];
      }
      delete canvBest_p;
      
      for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	delete leg_p[lI];      
      }
      leg_p.clear();            
    }
  }  

  inUnfoldFile_p->Close();
  delete inUnfoldFile_p;

  if(outUnfoldConfigName.find(".config") != std::string::npos){
    outUnfoldConfigName.replace(outUnfoldConfigName.rfind(".config"), 7, "");
  }

  outUnfoldConfigName = outUnfoldConfigName + "_" + varName + ".config";

  if(outUnfoldConfigName.find("/") == std::string::npos){//means no path is supplied
    outUnfoldConfigName = "output/" + dateStr + "/" + outUnfoldConfigName;    
  }
  
  outUnfoldConfig_p->SetValue("UNFOLDFILENAME", inUnfoldFileName.c_str());
  
  outUnfoldConfig_p->WriteFile(outUnfoldConfigName.c_str());
  
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::cout << "GDJPLOTUNFOLDDIAGNOSTICS COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotUnfoldDiagnostics.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotUnfoldDiagnostics(argv[1]);
  return retVal;
}
