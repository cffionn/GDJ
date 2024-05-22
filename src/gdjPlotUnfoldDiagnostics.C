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
#include "TH1D.h"
#include "TH2D.h"
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
#include "include/binFlattener.h"
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
#include "include/varUtil.h"

template <typename T>
bool axisFix(T* input)
{
  if(input == nullptr){
    std::cout << "AXISFIX ERROR: Pad is nullptr, return false" << std::endl;
    return false;
  }

  input->RedrawAxis();
  input->SetTicks();

  return true;
}

bool tpadsAxisFix(std::vector<TPad*> pads_p)
{
  for(unsigned int pI = 0; pI < pads_p.size(); ++pI){
    if(!axisFix(pads_p[pI])) return false;
  }

  return true;
}

//Switch from void to bool because there are failure conditions to setup
//Needs to be a pointer to a vector of pointers
template <typename T>
bool setupLegends(int font, float size, std::vector<TLegend*>* leg_p, std::vector<T*> histsForLeg_p, std::vector<std::string> strForLeg, std::vector<std::string> legFillStr, std::vector<Int_t>* legN, std::vector<Float_t> legX, std::vector<Float_t> legY)
{
  Int_t nHistsForLeg = histsForLeg_p.size();
  if(nHistsForLeg != strForLeg.size()){
    std::cout << "Number of histrograms \'" << nHistsForLeg << "\' does not match number of stringsForLegend provided, \'" << strForLeg.size() << "\'. setupLegends() return false" << std::endl;
    return false;
  }

  if(legFillStr.size() != histsForLeg_p.size()){
    legFillStr.clear();
    for(unsigned int hI = 0; hI < histsForLeg_p.size(); ++hI){
      legFillStr.push_back("P L");
    }
  }

  //Float just for future multiplication
  Float_t maxLegXChar = -1;
  for(unsigned int i = 0; i < strForLeg.size(); ++i){
    if(strForLeg[i].size() > maxLegXChar) maxLegXChar = strForLeg[i].size();
  }

  //Define the number of histograms that will go to each legend
  Int_t numberHists = histsForLeg_p.size();
  Int_t histsPerLeg = 0;

  Int_t totalCustom = 0;
  for(unsigned int hI = 0; hI < legN->size(); ++hI){
    totalCustom += legN->at(hI);
  }

  if(totalCustom < numberHists){
    if(totalCustom != 0) std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
    legN->clear();
  }
  if(legN->size() == 0){
    while(histsPerLeg*(Int_t)legX.size() < numberHists){
      ++histsPerLeg;
    }
    for(Int_t dI = 0; dI < legX.size(); ++dI){
      legN->push_back(histsPerLeg);
    }
  }

 for(unsigned int lI = 0; lI < legX.size(); ++lI){
    leg_p->push_back(nullptr);
    leg_p->at(lI) = new TLegend(legX[lI], legY[lI] - ((Float_t)legN->at(lI))*0.0455, legX[lI] + maxLegXChar*0.02, legY[lI]);

    leg_p->at(lI)->SetTextFont(font);
    leg_p->at(lI)->SetTextSize(size);
    leg_p->at(lI)->SetBorderSize(0);
    leg_p->at(lI)->SetFillColor(0);
    leg_p->at(lI)->SetFillStyle(0);
  }

  unsigned int legCounter = 0;
  for(unsigned int lI = 0; lI < legN->size(); ++lI){
    for(unsigned int dI = 0; dI < legN->at(lI); ++dI){
      leg_p->at(lI)->AddEntry(histsForLeg_p[legCounter], strForLeg[legCounter].c_str(), legFillStr[legCounter].c_str());
      ++legCounter;

      if(legCounter == histsForLeg_p.size()) break;
    }
    if(legCounter == histsForLeg_p.size()) break;
  }

  return true;
}


//Middle hist in vector should be nominal - others are adjacent iterations
bool drawBoxAdjacentIter(TPad* pad_p, std::vector<TH1*> hists_p, Float_t opacity, std::vector<double>* deltaVals)
{
  if(pad_p == nullptr){
    std::cout << "drawBoxAdjacentIter ERROR: TPad ptr given is nullptr. return false" << std::endl;
    return false;
  }

  pad_p->cd();

  //Has to be 3, nominal + 2 adjacent iterations
  if(hists_p.size() != 3){
    std::cout << "drawBoxAdjacentIter ERROR: Number of hists must be three; given \'" << hists_p.size() << "\'. return false" << std::endl;
    return false;
  }

  //Color of the middle (nominal!) hist - others are adjacent iterations
  Int_t color = hists_p[1]->GetMarkerColor();

  TBox* box_p = new TBox();
  box_p->SetFillColorAlpha(color, opacity);
  for(Int_t bIX = 0; bIX < hists_p[0]->GetNbinsX(); ++bIX){
    Double_t x1 = hists_p[0]->GetBinLowEdge(bIX+1);
    Double_t x2 = hists_p[0]->GetBinLowEdge(bIX+2);

    Double_t yDelta1 = TMath::Abs(hists_p[0]->GetBinContent(bIX+1) - hists_p[1]->GetBinContent(bIX+1));
    Double_t yDelta2 = TMath::Abs(hists_p[1]->GetBinContent(bIX+1) - hists_p[2]->GetBinContent(bIX+1));
    Double_t yDelta = TMath::Max(yDelta1, yDelta2);

    Double_t yLow = TMath::Max(hists_p[1]->GetBinContent(bIX+1) - yDelta, hists_p[1]->GetMinimum());
    Double_t yHigh = TMath::Min(hists_p[1]->GetBinContent(bIX+1) + yDelta, hists_p[1]->GetMaximum());

   if(deltaVals != nullptr){
      if(hists_p[1]->GetBinContent(bIX+1) < TMath::Power(10,-50)) deltaVals->push_back(0.0);
      else deltaVals->push_back(yDelta/(double)hists_p[1]->GetBinContent(bIX+1));
    }

    box_p->DrawBox(x1, yLow, x2, yHigh);
  }

  delete box_p;

  return true;
}

bool drawBoxAdjacentIterRel(TPad* pad_p, TH1* hist_p, Float_t opacity, std::vector<double>* deltaVals)
{
  if(pad_p == nullptr){
    std::cout << "drawBoxAdjacentIterRel ERROR: TPad ptr given is nullptr. return false" << std::endl;
    return false;
  }

  if(hist_p == nullptr){
    std::cout << "drawBoxAdjacentIterRel ERROR: THist ptr given is nullptr. return false" << std::endl;
    return false;
  }

  pad_p->cd();

  Int_t color = hist_p->GetMarkerColor();

  TBox* box_p = new TBox();
  box_p->SetFillColorAlpha(color, opacity);

  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    Double_t x1 = hist_p->GetBinLowEdge(bIX+1);
    Double_t x2 = hist_p->GetBinLowEdge(bIX+2);

    Double_t yLow = TMath::Max(1.0 - deltaVals->at(bIX), hist_p->GetMinimum());
    Double_t yHigh = TMath::Min(1.0 + deltaVals->at(bIX), hist_p->GetMaximum());

    box_p->DrawBox(x1, yLow, x2, yHigh);
  }
  delete box_p;

  return true;
}

int gdjPlotUnfoldDiagnostics(std::string inConfigFileName)
{
  //Hard-code some latex/legend parameters
  const int titleFont = 42;
  const double titleSize = 0.03;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;

  //Internal stopwatch for code timing tests
  cppWatch globalTimer;
  globalTimer.start();

  //Check config exists
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Create date-stamped output and pdf dir
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  //Global debug handler - do 'export DOGLOBALDEBUGROOT=1' to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  //Grab our input config, define necessary internal params
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INUNFOLDFILENAME",
					      "OUTUNFOLDCONFIG",
					      "TERMTHRESH",
					      "TERM1D",
					      "TERM2D",
					      "TERMBOTH",
					      "DORELATIVETERM",
					      "DOBINWIDTHANDSELFNORM",
					      "SAVETAG",
					      "SAVEEXT",
					      "DOLOGX",
					      "DOLOGY",
					      "DOJTPTCUT",
					      "JETVARMIN.GLOBAL",
					      "JETVARMAX.GLOBAL",
					      "JETVARLEGX.GLOBAL.0",
					      "JETVARLEGY.GLOBAL.0",
					      "JETVARLABELX.GLOBAL",
					      "JETVARLABELY.GLOBAL",
					      "JETVARLABELALIGNRIGHT.GLOBAL",
					      "DIAGLEGX.GLOBAL.0",
					      "DIAGLEGY.GLOBAL.0",
					      "DIAGLABELX.GLOBAL",
					      "DIAGLABELY.GLOBAL",
					      "DIAGLABELALIGNRIGHT.GLOBAL"};

  //Check that the required parameters exist
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Start grabbing parameters from config file
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outUnfoldConfigName = config_p->GetValue("OUTUNFOLDCONFIG", "");

  TEnv* outUnfoldConfig_p = new TEnv();

  //Threshold by which we terminate the unfolding - deltaStat contribution to quad sum exceeds this value, term
  const double termThresh = config_p->GetValue("TERMTHRESH", -1.0);
  outUnfoldConfig_p->SetValue("TERMTHRESH", config_p->GetValue("TERMTHRESH", termThresh));

  //Do we terminate each 1D distribution separately or the whole 2D simultaneously - 2D approach is the correct one, but we leave 1D for debug purposes + possible systematic
  const bool term1D = (bool)config_p->GetValue("TERM1D", 0);
  outUnfoldConfig_p->SetValue("TERM1D", config_p->GetValue("TERM1D", term1D));
  const bool term2D = (bool)config_p->GetValue("TERM2D", 0);
  outUnfoldConfig_p->SetValue("TERM2D", config_p->GetValue("TERM2D", term2D));
  const bool termBoth = (bool)config_p->GetValue("TERMBoth", 0);
  outUnfoldConfig_p->SetValue("TERMBoth", config_p->GetValue("TERMBoth", termBoth));

  //Make sure we have only picked an appropriate singular version of above
  if(term1D){
    if(term2D || termBoth){
      std::cout << "TERM1D \'" << term1D << "\' and 1 of TERM2D (" << term2D << ") or TERMBOTH (" << termBoth << ") selected. please pick 1. return 1" << std::endl;
      return 1;
    }
  }
  else if(term2D){
    if(term1D || termBoth){
      std::cout << "TERM2D \'" << term2D << "\' and 1 of TERM1D (" << term1D << ") or TERMBOTH (" << termBoth << ") selected. please pick 1. return 1" << std::endl;
      return 1;
    }
  }
  else if(!termBoth){
    std::cout << "All termination for best iterations 'TERM1D', 'TERM2D', 'TERMBOTH' are false. Please set one to true. Return 1" << std::endl;
    return 1;
  }

  //Continue grabbing parameters from config file
  const std::string saveTag = config_p->GetValue("SAVETAG", "");
  const std::string saveExt = config_p->GetValue("SAVEEXT", "");
  const bool doBinWidthAndSelfNorm = config_p->GetValue("DOBINWIDTHANDSELFNORM", 0);
  const bool doRelativeTerm = config_p->GetValue("DORELATIVETERM",  0);
  std::string deltaStr = "Rel";
  if(!doRelativeTerm) deltaStr = "Abs";

  const bool doLogX = config_p->GetValue("DOLOGX", 0);
  const bool doLogY = config_p->GetValue("DOLOGY", 0);

  const bool doJtPtCut = config_p->GetValue("DOJTPTCUT", 0);

  std::vector<std::string> validExt = {"pdf", "png"};

  if(!vectContainsStr(saveExt, &validExt)){
    std::cout << "Given SAVEEXT \'" << saveExt << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  const Int_t nIterCap = config_p->GetValue("NITERPLOTCAP", 100000);

  const Float_t jetVarMinGlobal = config_p->GetValue("JETVARMIN.GLOBAL", -1000.);
  const Float_t jetVarMaxGlobal = config_p->GetValue("JETVARMAX.GLOBAL", 1000.);


  const Float_t jetVarRatMinGlobal = config_p->GetValue("JETVARRATMIN.GLOBAL", -1000.0);
  const Float_t jetVarRatMaxGlobal = config_p->GetValue("JETVARRATMAX.GLOBAL", 1000.0);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //We define possibly multiple legends in the config file that we have complete control over
  //Gotta grab them all here
  std::vector<Float_t> jetVarLegXGlobal, jetVarLegYGlobal;
  std::vector<Int_t> jetVarLegNGlobal;
  for(Int_t i = 0; i < 100; ++i){
    Float_t tempLegX = config_p->GetValue(("JETVARLEGX.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegY = config_p->GetValue(("JETVARLEGY.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegN = config_p->GetValue(("JETVARLEGN.GLOBAL." + std::to_string(i)).c_str(), -1000.);

    if(tempLegX >= 0){
      jetVarLegXGlobal.push_back(tempLegX);
      jetVarLegYGlobal.push_back(tempLegY);
      if(tempLegN >= 0) jetVarLegNGlobal.push_back(tempLegN);
    }
  }
  if(jetVarLegNGlobal.size() != jetVarLegXGlobal.size()){
    std::cout << "Number of histograms per legend not specified manually for all legends - reverting to equal split" << std::endl;
    jetVarLegNGlobal.clear();
  }

  //Define label positions
  const Float_t jetVarLabelXGlobal = config_p->GetValue("JETVARLABELX.GLOBAL", -1000.);
  const Float_t jetVarLabelYGlobal = config_p->GetValue("JETVARLABELY.GLOBAL", -1000.);
  const Bool_t jetVarLabelAlignRightGlobal = config_p->GetValue("JETVARLABELALIGNRIGHT.GLOBAL", 0);

  //diag == diagnostic - we have an unfolding diagnostic plot that has similarly defined label and legend positionsing
  std::vector<Float_t> diagLegXGlobal, diagLegYGlobal;
  std::vector<Int_t> diagLegNGlobal;
  for(Int_t i = 0; i < 100; ++i){
    Float_t tempLegX = config_p->GetValue(("DIAGLEGX.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Float_t tempLegY = config_p->GetValue(("DIAGLEGY.GLOBAL." + std::to_string(i)).c_str(), -1000.);
    Int_t tempLegN = config_p->GetValue(("DIAGLEGN.GLOBAL." + std::to_string(i)).c_str(), -1000);

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
  const Float_t diagMinimum = config_p->GetValue("DIAGMINIMUM", 0.05);
  const Float_t diagMaximum = config_p->GetValue("DIAGMAXIMUM", 0.05);

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
  //Grabbing parameters from input config file complete

  //check that our given input files exist
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;

  //Grab our input unfolded histogram file, and internal configs + check for req. params
  TFile* inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");
  TEnv* inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");
  TEnv* inUnfoldFileLabel_p = (TEnv*)inUnfoldFile_p->Get("label");

  std::vector<std::string> necessaryUnfoldFileParams = {"ISMC",
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
							"NSUBJTPTBINS",
							"SUBJTPTBINSLOW",
							"SUBJTPTBINSHIGH",
							"SUBJTPTBINSDOLOG",
							"SUBJTPTBINSDOCUSTOM",
							"SUBJTGAMMAPTMIN",
							"SUBJTGAMMAPTMAX",
							"ISPP",
							"NITER",
							"JETR",
							"VARNAME"};

  //Global cut referential labels (dphi0 is the gamma+jet dphi cut)
  const std::string gammaJtDPhiStr = "DPhi0";
  //  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

  const Int_t nMaxSyst = 100;
  Int_t nSyst = inUnfoldFileConfig_p->GetValue("NSYST", -1);
  std::vector<std::string> systStrVect = strToVect(inUnfoldFileConfig_p->GetValue("SYST", ""));
  std::vector<std::string> systTypeStrVect = strToVect(inUnfoldFileConfig_p->GetValue("SYSTTYPE", ""));

  /*
  for(unsigned int sI = 0; sI < systStrVect.size(); ++sI){
    std::cout << systStrVect[sI] << ", " << systTypeStrVect[sI] << std::endl;
  }
  return 1;
  */

  const std::string inUnfoldName = inUnfoldFileConfig_p->GetValue("UNFOLD", "");
  std::vector<std::string> inUnfoldNames = strToVect(inUnfoldName);
  bool unfoldAll = false;
  for(unsigned int uI = 0; uI < inUnfoldNames.size(); ++uI){
    if(isStrSame(inUnfoldNames[uI], "ALL")){
      unfoldAll = true;
      break;
    }
  }


  //Additional params are required if it is Pb+Pb
  std::vector<std::string> pbpbParams = {"CENTBINS"};
  if(!checkEnvForParams(inUnfoldFileConfig_p, necessaryUnfoldFileParams)) return 1;

  std::string varName = inUnfoldFileConfig_p->GetValue("VARNAME", "");
  outUnfoldConfig_p->SetValue("VARNAME", varName.c_str());

  std::string varNameUpper = strLowerToUpper(varName);
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string varNameLabel = varNameToLabel(varName);


  std::string outUnfoldErrorName = outUnfoldConfigName.substr(0, outUnfoldConfigName.find("."));

  outUnfoldErrorName = "output/" + dateStr + "/" + outUnfoldErrorName + "_" + varNameUpper + "_SysError.root";
  TFile* outSysFile_p = new TFile(outUnfoldErrorName.c_str(), "RECREATE");


  const Bool_t isMultijet = isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");

  const Float_t jtPtBinsLow = inUnfoldFileConfig_p->GetValue("JTPTBINSLOW", -999.0);
  const Float_t jtPtBinsLowReco = inUnfoldFileConfig_p->GetValue("JTPTBINSLOWRECO", -999.0);
  const Float_t jtPtBinsHigh = inUnfoldFileConfig_p->GetValue("JTPTBINSHIGH", -999.0);

  //Add to global labels
  int jetR = inUnfoldFileConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  else{
    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  }

  std::string gammaJtDPhiCutLabel = inUnfoldFileConfig_p->GetValue("GAMMAJTDPHI", "");
  if(gammaJtDPhiCutLabel.size() == 0){
    std::cout << "GammaJtDPhi not found, check line L" << __LINE__ << ". return 1" << std::endl;
    return 1;
  }
  if(gammaJtDPhiCutLabel.find("pi") != std::string::npos) gammaJtDPhiCutLabel.replace(gammaJtDPhiCutLabel.find("pi"), 2, "#pi");

  std::string gammaMultiJtDPhiCut = inUnfoldFileConfig_p->GetValue("GAMMAMULTIJTDPHI", "");
  if(gammaMultiJtDPhiCut.find("pi") != std::string::npos) gammaMultiJtDPhiCut.replace(gammaMultiJtDPhiCut.find("pi"), 2, "#pi");

  //Add to the global labels
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
    if(!checkEnvForParams(inUnfoldFileConfig_p, pbpbParams)) return 1;
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
  //Define a bunch of parameters for the TCanvas
  const Int_t nPad = 3;

  Double_t height = 1200.0;
  Double_t width = 942.831;
  const Double_t leftMargin = 0.16;
  const Double_t bottomAbsFracOfLeft = 0.8;

  const Double_t bottomMarginFracHeight = bottomAbsFracOfLeft*leftMargin*width/height;

  const Double_t topPanelFrac = 0.48;
  const Double_t bottomPanelFrac = 0.39;
  const Double_t padSplit1 = 1.0 - topPanelFrac;
  //  const Double_t padSplit2 = padSplit1 - (padSplit1 - bottomMarginFracHeight)/2.0;
  const Double_t padSplit2 = bottomPanelFrac;

  const Double_t bottomMargin = bottomMarginFracHeight/padSplit2;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;

  const Double_t leftMargin2D = 0.12;
  const Double_t bottomMargin2D = leftMargin2D;
  const Double_t topMargin2D = 0.02;
  const Double_t rightMargin2D = 0.12;
  Double_t height2D = 900.0;
  Double_t width2D = height2D*(1.0 - topMargin2D* - bottomMargin2D)/(1.0 - leftMargin2D - rightMargin2D);

  const Int_t nMaxPtBins = 200;

  //GammaPt Bin Handling
  const Int_t nGammaPtBins = inUnfoldFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Float_t gammaPtBinsLowReco = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOWRECO", 80.0);
  const Float_t gammaPtBinsHighReco = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGHRECO", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inUnfoldFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inUnfoldFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string gammaPtBinsStr = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBins != nGammaPtBinsTemp){
      std::cout << "nGammaPtBins dont match - by stored value \'" << nGammaPtBins << "\', via vector \'" << nGammaPtBinsTemp << "\'" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  //SubJtPt Bin Handling
  const Int_t nSubJtPtBins = inUnfoldFileConfig_p->GetValue("NSUBJTPTBINS", 0);
  const Float_t subJtPtBinsLow = inUnfoldFileConfig_p->GetValue("SUBJTPTBINSLOW", -1.0);
  const Float_t subJtPtBinsHigh = inUnfoldFileConfig_p->GetValue("SUBJTPTBINSHIGH", -1.0);
  const Bool_t subJtPtBinsDoLog = (bool)inUnfoldFileConfig_p->GetValue("SUBJTPTBINSDOLOG", 1);
  const Bool_t subJtPtBinsDoCustom = (bool)inUnfoldFileConfig_p->GetValue("SUBJTPTBINSDOCUSTOM", 1);
  Double_t subJtPtBins[nMaxPtBins+1];
  Double_t subJtGammaPtMin = inUnfoldFileConfig_p->GetValue("SUBJTGAMMAPTMIN", -1000.0);
  Double_t subJtGammaPtMax = inUnfoldFileConfig_p->GetValue("SUBJTGAMMAPTMAX", -1000.0);
  Double_t subJtGammaPtBins[nMaxPtBins+1];

  if(subJtPtBinsDoLog) getLogBins(subJtPtBinsLow, subJtPtBinsHigh, nSubJtPtBins, subJtPtBins);
  else if(subJtPtBinsDoCustom){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string subJtPtBinsStr = inUnfoldFileConfig_p->GetValue("SUBJTPTBINSCUSTOM", "");
    if(subJtPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> subJtPtBinsTemp = strToVectF(subJtPtBinsStr);
    Int_t nSubJtPtBinsTemp = ((Int_t)subJtPtBinsTemp.size()) - 1;
    if(nSubJtPtBins != nSubJtPtBinsTemp){
      std::cout << "nSubJtPtBins dont match - by stored value \'" << nSubJtPtBins << "\', via vector \'" << nSubJtPtBinsTemp << "\'" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < subJtPtBinsTemp.size(); ++gI){
      subJtPtBins[gI] = subJtPtBinsTemp[gI];
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }
  else getLinBins(subJtPtBinsLow, subJtPtBinsHigh, nSubJtPtBins, subJtPtBins);

  binFlattener subJtGammaPtBinFlattener;
  //Now create your flattened bins using the binFlattener class
  if(!subJtGammaPtBinFlattener.Init("subJtGammaPtBinFlattener", nGammaPtBins, gammaPtBins, nSubJtPtBins, subJtPtBins)){
    std::cout << "Error gdjNTupleToHist: Failure to init binFlattener, L" << __LINE__ << ". return 1" << std::endl;
    return 1;
  }
  std::vector<double> subJtGammaPtBinsV = subJtGammaPtBinFlattener.GetFlattenedBins(subJtGammaPtMin, subJtGammaPtMax);
  for(unsigned int sgI = 0; sgI < subJtGammaPtBinsV.size(); ++sgI){
    subJtGammaPtBins[sgI] = subJtGammaPtBinsV[sgI];
  }



  if(doJtPtCut){
    bool jtPtCutFound = false;
    double deltaVal = 0.001;
    for(Int_t sI = 0; sI < nSubJtPtBins; ++sI){
      if(TMath::Abs(subJtPtBins[sI] - jtPtBinsLowReco) < deltaVal){
	jtPtCutFound = true;
	break;
      }
    }

    if(!jtPtCutFound){
      std::cout << "doJtPtCut is on but jtPtCut, " << jtPtBinsLowReco << " is not found in bins, ";
      for(Int_t sI = 0; sI < nSubJtPtBins; ++sI){
	std::cout << subJtPtBins[sI] << ", ";
      }
      std::cout << subJtPtBins[nSubJtPtBins] << ". return 1" << std::endl;
      return 1;
    }
  }

  gStyle->SetOptStat(0);

  //Jet Var Bin Handling
  std::vector<Float_t> jetVarMax, jetVarMin;
  std::string jtVar = inUnfoldFileConfig_p->GetValue("VARNAME", "");
  if(jtVar.size() == 0){
    std::cout << "Varname: \'" << jtVar << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  std::string jtVarLower = returnAllLowercaseString(jtVar);
  std::string jtVarUpper = returnAllCapsString(jtVar);

  std::string varPrefix = "";
  if(isStrSame(jtVarUpper, "PT")) varPrefix = "JT";

  std::string binVarStr = varPrefix + jtVarUpper;
  if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";
  else if(isStrSame(binVarStr, "DRJJ")) binVarStr = "DR";

   Int_t nVarBins = inUnfoldFileConfig_p->GetValue(("N" + binVarStr + "BINS").c_str(), -1);
  Float_t varBinsLow = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSLOW").c_str(), -1.0);
  Float_t varBinsHigh = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSHIGH").c_str(), -1.0);
  Bool_t varBinsDoLog = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSDOLOG").c_str(), 0);
  Bool_t varBinsDoCustom = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSDOCUSTOM").c_str(), 0);
  Float_t varBinsLowReco = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSLOWRECO").c_str(), varBinsLow);
  Float_t varBinsHighReco = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSHIGHRECO").c_str(), varBinsHigh);
  Double_t varBins[nMaxPtBins+1];

  /*
  std::cout << "VarbinsLow: " << (binVarStr + "BINSLOW") << ", " << varBinsLow << std::endl;
  std::cout << "VarbinsHigh: " << (binVarStr + "BINSHIGH") << ", " << varBinsHigh << std::endl;

  std::cout << "VarbinsLowReco: " << (binVarStr + "BINSLOWRECO") << ", " << varBinsLowReco << std::endl;
  std::cout << "VarbinsHighReco: " << (binVarStr + "BINSHIGHRECO") << ", " << varBinsHighReco << std::endl;
  return 1;
  */

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(varBinsDoLog) getLogBins(varBinsLow, varBinsHigh, nVarBins, varBins);
  else if(varBinsDoCustom){
    std::string varBinsStr = inUnfoldFileConfig_p->GetValue((binVarStr + "BINSCUSTOM").c_str(), "");
    if(varBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> varBinsTemp = strToVectF(varBinsStr);
    Int_t nVarBinsTemp = ((Int_t)varBinsTemp.size()) - 1;
    if(nVarBins != nVarBinsTemp){
      return 1;
    }

    for(unsigned int gI = 0; gI < varBinsTemp.size(); ++gI){
      varBins[gI] = varBinsTemp[gI];
    }
  }
  else getLinBins(varBinsLow, varBinsHigh, nVarBins, varBins);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //We need to create a map of bins to exclude in judging the iterative convergence
  //Map xBin, yBin, bool true false use in deltaAbs = deltaIter quadSum deltaStat
  std::map<int, std::map<int, bool> > goodBinMap;
  for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
     if(isMultijet){//Yaxis is a flattened 2-D array, handle as follows
      for(Int_t bIY = 0; bIY < nGammaPtBins*nSubJtPtBins; ++bIY){
	Int_t gammaBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
	Int_t subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);

	Double_t gammaBinCenter = (gammaPtBins[gammaBinPos] + gammaPtBins[gammaBinPos+1])/2.0;
	Double_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;

	bool isBinGood = true;
	Double_t binCenterX = (varBins[bIX] + varBins[bIX+1])/2.0;

	if(binCenterX < varBinsLowReco) isBinGood = false;
	else if(binCenterX > varBinsHighReco) isBinGood = false;

	//	std::cout << "bIX, bIY, varBinCenter, gammaBinCenter, subJtPtBinCenter, goodBin: " << bIX << ", " << bIY  << ", " << binCenterX << ", " << gammaBinCenter << ", " << subJtPtBinCenter << ", " << isBinGood << std::endl;

	if(isBinGood){
	  if(gammaBinCenter < gammaPtBinsLowReco) isBinGood = false;
	  else if(gammaBinCenter > gammaPtBinsHighReco) isBinGood = false;
	  else if(subJtPtBinCenter < jtPtBinsLowReco) isBinGood = false;
	}
	goodBinMap[bIX][bIY] = isBinGood;
      }
    }
    else{//y-axis defined by gamma pt bins
      for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
	Double_t binCenterY = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;

	bool isBinGood = true;
	Double_t binCenterX = (varBins[bIX] + varBins[bIX+1])/2.0;
	if(binCenterX < varBinsLowReco) isBinGood = false;
	else if(binCenterX > varBinsHighReco) isBinGood = false;

	if(isBinGood){
	  if(binCenterY < gammaPtBinsLowReco) isBinGood = false;
	  else if(binCenterY > gammaPtBinsHighReco) isBinGood = false;
	}
	goodBinMap[bIX][gI] = isBinGood;
      }
    }
  }

  /*
  std::cout << "GAMMA PT LOW HIGH : " << gammaPtBinsLowReco << "-" << gammaPtBinsHighReco << std::endl;
  std::cout << "SUBJET PT LOW HIGH : " << jtPtBinsLowReco << std::endl;
  std::cout << "VAR BINS LOW HIGH: " << varBinsLowReco << ", " << varBinsHighReco << std::endl;

  for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
    std::cout << "bIX: " << bIX << std::endl;
    for(Int_t bIY = 0; bIY < nGammaPtBins*nSubJtPtBins; ++bIY){
      std::cout << " bIY, value: " << bIY << ", " << goodBinMap[bIX][bIY] << std::endl;

    }
  }
  return 1;
  */

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  //Loop over all jet systematics
  for(Int_t systI = 0; systI < nSyst; ++systI){
    if(!unfoldAll){
      if(!vectContainsStr(returnAllCapsString(systStrVect[systI]), &inUnfoldNames)) continue;
    }

    //Photon-pt only; double vector for centrality, iterations
    std::vector<std::string> statsDeltaPhoPtNames, iterDeltaPhoPtNames, totalDeltaPhoPtNames;
    std::vector<std::string> recoHistPhoPtNames, truthHistPhoPtNames;
    std::vector<std::vector<std::string> > unfoldNamesPhoPt, refoldNamesPhoPt;

    //Photon-pt + Jet-var unfold
    std::vector<std::string> statsDeltaPhoPtJetVarNames, iterDeltaPhoPtJetVarNames, totalDeltaPhoPtJetVarNames;
    std::vector<std::string> recoHistPhoPtJetVarNames, truthHistPhoPtJetVarNames;
    std::vector<std::vector<std::string> > unfoldNamesPhoPtJetVar, refoldNamesPhoPtJetVar;


    std::string phoSystStrAlt = "NOMINAL";
    std::string phoSystStr = "Nominal";
    bool isIsoPur = isStrSame(systTypeStrVect[systI], "PHOISOANDPUR");
    bool isGammaSyst = isStrSame(systTypeStrVect[systI], "GESGER") || isIsoPur;
    if(isGammaSyst) phoSystStr = systStrVect[systI];
    if(isIsoPur) phoSystStrAlt = systStrVect[systI];


    //Fill out the name vectors
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
      unfoldNamesPhoPt.push_back({});
      refoldNamesPhoPt.push_back({});

      unfoldNamesPhoPtJetVar.push_back({});
      refoldNamesPhoPtJetVar.push_back({});

      recoHistPhoPtNames.push_back("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_" + phoSystStrAlt + "_PURCORR_COMBINED_h");
      truthHistPhoPtNames.push_back("photonPtTruth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h");
      statsDeltaPhoPtNames.push_back("statsDelta_PhoPt_" + centBinsStr[cI] + "_" + phoSystStr + "_PURCORR_COMBINED_h");
      iterDeltaPhoPtNames.push_back("iterDelta_PhoPt_" + centBinsStr[cI] + "_" + phoSystStr + "_PURCORR_COMBINED_h");
      totalDeltaPhoPtNames.push_back("totalDelta_PhoPt_" + centBinsStr[cI] + "_" + phoSystStr + "_PURCORR_COMBINED_h");

      std::string recoTruthString = "NOMINAL";
      if(systStrVect[systI].find("JTPTCUT") != std::string::npos) recoTruthString = "JTPTCUT";
      else if(isIsoPur) recoTruthString = phoSystStrAlt;

      recoHistPhoPtJetVarNames.push_back("photonPtJet" + jtVar + "Reco_PreUnfold_" + centBinsStr[cI] + "_" + recoTruthString + "_PURCORR_COMBINED_h");
      truthHistPhoPtJetVarNames.push_back("photonPtJet" + jtVar + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h");
      statsDeltaPhoPtJetVarNames.push_back("statsDelta_PhoPtJet" + jtVar + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");
      iterDeltaPhoPtJetVarNames.push_back("iterDelta_PhoPtJet" + jtVar + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");
      totalDeltaPhoPtJetVarNames.push_back("totalDelta_PhoPtJet" + jtVar + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");

      for(Int_t uI = 0; uI < nIterForLoop+1; ++uI){
	unfoldNamesPhoPt[cI].push_back("photonPtReco_Iter" + std::to_string(uI) + "_" + centBinsStr[cI] + "_" + phoSystStr + "_PURCORR_COMBINED_h");
	refoldNamesPhoPt[cI].push_back("photonPtReco_Iter" + std::to_string(uI) + "_" + centBinsStr[cI] + "_" + phoSystStr + "_PURCORR_COMBINED_Refolded_h");

	unfoldNamesPhoPtJetVar[cI].push_back("photonPtJet" + jtVar + "Reco_Iter" + std::to_string(uI) + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");
	refoldNamesPhoPtJetVar[cI].push_back("photonPtJet" + jtVar + "Reco_Iter" + std::to_string(uI) + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_Refolded_h");
      }
    }

    //two iterations of processing - one is pure photon, other is photon+jet 2d
    const Int_t nProc = 2;
    for(Int_t pI = 0; pI < nProc; ++pI){
      //pI = 0; Photon-pt processing
      //pI = 1; Photon+Jet Var processing
      //Skip all jet systematics on pure photon unfold

      bool isGammaSyst = isStrSame(systTypeStrVect[systI], "GESGER") || systI == 0 || isStrSame(systStrVect[systI], "ISO85") || isStrSame(systStrVect[systI], "ISO95");
      if(pI == 0 && !isGammaSyst) continue;

      std::vector<std::string> statsDeltaNames = statsDeltaPhoPtNames;
      std::vector<std::string> iterDeltaNames = iterDeltaPhoPtNames;
      std::vector<std::string> totalDeltaNames = totalDeltaPhoPtNames;
      std::vector<std::string> recoHistNames = recoHistPhoPtNames;
      std::vector<std::string> truthHistNames = truthHistPhoPtNames;
      std::vector<std::vector<std::string> > unfoldNames = unfoldNamesPhoPt;
      std::vector<std::vector<std::string> > refoldNames = refoldNamesPhoPt;

      if(pI == 1){
	statsDeltaNames = statsDeltaPhoPtJetVarNames;
	iterDeltaNames = iterDeltaPhoPtJetVarNames;
	totalDeltaNames = totalDeltaPhoPtJetVarNames;
	recoHistNames = recoHistPhoPtJetVarNames;
	truthHistNames = truthHistPhoPtJetVarNames;
	unfoldNames = unfoldNamesPhoPtJetVar;
	refoldNames = refoldNamesPhoPtJetVar;
      }

      for(unsigned int cI = 0; cI < statsDeltaNames.size(); ++cI){
	//photon pt alone unfolding result eval
	//Declare both new histograms and define/initialize null histograms we will grab from file
	TH1D* statsDelta_p = new TH1D(statsDeltaNames[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIterForLoop, 0.5, ((double)nIterForLoop)+0.5);
	TH1D* iterDelta_p = new TH1D(iterDeltaNames[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIterForLoop, 0.5, ((double)nIterForLoop)+0.5);
	TH1D* totalDelta_p = new TH1D(totalDeltaNames[cI].c_str(), (";Number of Iterations;#sqrt{#Sigma #delta^{2}_{" + deltaStr + "}}").c_str(), nIterForLoop, 0.5, ((double)nIterForLoop)+0.5);
	TH1D* unfold1D_p[nIterMax];
	TH1D* refold1D_p[nIterMax];
	TH2D* unfold2D_p[nIterMax];
	TH2D* refold2D_p[nIterMax];

	for(int i = 0; i < nIterMax; ++i){
	  unfold1D_p[i] = nullptr;
	  refold1D_p[i] = nullptr;
	  unfold2D_p[i] = nullptr;
	  refold2D_p[i] = nullptr;
	}

	//Global labels -> tempLabels in case I need to make variable specific alterations
	std::vector<std::string> tempLabels = globalLabels;
	if(!isMC && truthHistNames[cI].find("PURCORR") != std::string::npos) tempLabels = globalLabelsDataOverrides;
	else if(!isMC && truthHistNames[cI].find("HalfTo") != std::string::npos) tempLabels = globalLabelsMCOverrides;
	tempLabels.push_back("#bf{" + prettyString(gammaPtBins[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false) + "}");

	if(isMultijet){
	  tempLabels.push_back("|#Delta#phi_{JJ#gamma}| > " + gammaMultiJtDPhiCut);
	}

	if(!isPP) tempLabels.push_back(centBinsLabel[cI%centBinsLabel.size()]);
	else tempLabels.push_back("p+p");
	//EndLabelMods

	//Define tlatex label object
	TLine* line_p = new TLine();
	line_p->SetLineStyle(2);
	line_p->SetLineColor(1);

	TLatex* label_p = new TLatex();
	initLabel(label_p, titleFont, titleSize, diagLabelAlignRightGlobal);

	//Define leg var
	std::vector<TLegend*> leg_p, legBest_p;
	std::vector<TH1D*> histsForLeg_p = {iterDelta_p, statsDelta_p, totalDelta_p};
	std::vector<std::string> stringsForLeg_p = {"Iter. Change", "Statistical Unc.", "Quad. Sum"};
	if(!setupLegends(titleFont, titleSize, &leg_p, histsForLeg_p, stringsForLeg_p, {}, &diagLegNGlobal, diagLegXGlobal, diagLegYGlobal)) return 1;

	std::vector<TH1D*> unfoldedHists1D_p, refoldedHists1D_p;
	std::vector<TH2D*> unfoldedHists2D_p, refoldedHists2D_p;
	unfoldedHists1D_p.reserve(nIterForLoop);
	refoldedHists1D_p.reserve(nIterForLoop);
	unfoldedHists2D_p.reserve(nIterForLoop);
	refoldedHists2D_p.reserve(nIterForLoop);

	for(Int_t i = 1; i < nIterForLoop+1; ++i){
	  if(pI == 0){
	    unfold1D_p[i] = (TH1D*)inUnfoldFile_p->Get(unfoldNames[cI][i].c_str());
	    unfoldedHists1D_p.push_back(unfold1D_p[i]);

	    refold1D_p[i] = (TH1D*)inUnfoldFile_p->Get(refoldNames[cI][i].c_str());
	    refoldedHists1D_p.push_back(refold1D_p[i]);
	  }
	  else if(pI == 1){
	    unfold2D_p[i] = (TH2D*)inUnfoldFile_p->Get(unfoldNames[cI][i].c_str());
	    unfoldedHists2D_p.push_back(unfold2D_p[i]);

	    refold2D_p[i] = (TH2D*)inUnfoldFile_p->Get(refoldNames[cI][i].c_str());
	    refoldedHists2D_p.push_back(refold2D_p[i]);
	  }
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	bool doIterPrint = centBinsStr[cI].find("Cent0to10") != std::string::npos && pI == 1 && systI == 0 && false;

	TH1D* reco1D_p = nullptr;
	TH2D* reco2D_p = nullptr;

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;
	if(pI == 0){
	  std::cout << recoHistNames[cI] << std::endl;
	  reco1D_p = (TH1D*)inUnfoldFile_p->Get(recoHistNames[cI].c_str());
	  getIterativeHists(reco1D_p, unfoldedHists1D_p, statsDelta_p, iterDelta_p, totalDelta_p, doRelativeTerm, gammaPtBinsLowReco, gammaPtBinsHighReco);
	}
	else if(pI == 1){
	  reco2D_p = (TH2D*)inUnfoldFile_p->Get(recoHistNames[cI].c_str());

	  if(isMultijet) getIterativeHists2D(reco2D_p, unfoldedHists2D_p, statsDelta_p, iterDelta_p, totalDelta_p, doRelativeTerm, goodBinMap, doIterPrint);
 	  else getIterativeHists2D(reco2D_p, unfoldedHists2D_p, statsDelta_p, iterDelta_p, totalDelta_p, doRelativeTerm, goodBinMap);
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	TCanvas* canv_p = new TCanvas("canv_p", "", width, height*3.0/4.0);
	setMargins(canv_p, leftMargin, topMargin, rightMargin, leftMargin/2.0); //intentionally settng bottom to left margin value
	canv_p->cd();

	HIJet::Style::EquipHistogram(iterDelta_p, 2);
	HIJet::Style::EquipHistogram(statsDelta_p, 1);
	HIJet::Style::EquipHistogram(totalDelta_p, 0);

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	totalDelta_p->SetMinimum(diagMinimum);
	totalDelta_p->SetMaximum(diagMaximum);
	totalDelta_p->DrawCopy("HIST E1 P");
	statsDelta_p->DrawCopy("HIST E1 P SAME");
	iterDelta_p->DrawCopy("HIST E1 P SAME");
	totalDelta_p->DrawCopy("HIST E1 P SAME");

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	gPad->SetLogy();
	//Draw Legends
	for(unsigned int lI = 0; lI < leg_p.size(); ++lI){leg_p[lI]->Draw("SAME");}

	//Draw labels
	int labelPos = 0;
	for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	  if(pI == 0 && tempLabels[tI].find("jet") != std::string::npos) continue;

	  label_p->DrawLatex(diagLabelXGlobal, diagLabelYGlobal - labelPos*0.05, tempLabels[tI].c_str());
	  ++labelPos;
	}
	axisFix(canv_p);
	std::string saveName = "pdfDir/" + dateStr + "/diag";
	if(pI == 0) saveName = saveName + "PhoPt";
	else if(pI == 1) saveName = saveName + "PhoPtJet" + varName;
	saveName = saveName + "_" + centBinsStr[cI%centBinsStr.size()] + "_" + systStrVect[systI];
	saveName = saveName + "_Delta" + deltaStr + "_" + saveTag +  "." + saveExt;

	quietSaveAs(canv_p, saveName);
	delete canv_p;
	delete label_p;

	for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
	  if(leg_p[lI] == nullptr) continue;
	  delete leg_p[lI];
	}
	leg_p.clear();

	//Write photon pt termination point to our config file
	std::string configStr = "GAMMAPT_" + centBinsStr[cI%centBinsStr.size()] + "_" + systStrVect[systI];
	if(pI == 1) configStr = varNameUpper + "_" + configStr;

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	int termPos = -1;
	for(Int_t i = 1; i < nIterForLoop+1; ++i){
	  Double_t totalDelta = totalDelta_p->GetBinContent(i);
	  Double_t statsDelta = statsDelta_p->GetBinContent(i);

	  if(statsDelta/totalDelta  >= termThresh){
	    termPos = i;
	    break;
	  }
	}
	//MCSTAT our termination unfolding criteria breaks because the error bars are no longer from the measurement but from the MC response matrix
	//This variation should just take nominal termination
	if(isStrSame(systStrVect[systI], "MCSTAT")){
	  termPos = outUnfoldConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI%centBinsStr.size()] + "_Nominal").c_str(), -1);
	}

	outUnfoldConfig_p->SetValue(configStr.c_str(), termPos);

	TH1D* truth1D_p = nullptr;
	TH2D* truth2D_p = nullptr;

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	if(pI == 0){
	  if(truthHistNames[cI].find("TRUTH_COMBINED") == std::string::npos || isMC) truth1D_p = (TH1D*)inUnfoldFile_p->Get(truthHistNames[cI].c_str());
	  else if(!isMC) truth1D_p = (TH1D*)inUnfoldFile_p->Get(truthHistNames[cI].c_str());
	}
	else if(pI == 1){
	  if(truthHistNames[cI].find("TRUTH_COMBINED") == std::string::npos || isMC) truth2D_p = (TH2D*)inUnfoldFile_p->Get(truthHistNames[cI].c_str());
	  else if(!isMC) truth2D_p = (TH2D*)inUnfoldFile_p->Get(truthHistNames[cI].c_str());
	}

	int nGammaPtBinsForUnfold = 1;
	if(pI == 1) nGammaPtBinsForUnfold = nGammaPtBins;

	//	std::cout << "CHECKING RECO1D: " << reco1D_p->GetName() << std::endl;
	std::vector<TH1D*> reco_p = {reco1D_p};
	std::vector<TH1D*> truth_p = {truth1D_p};
	std::vector<std::vector<TH1D*> > unfold_p;
	std::vector<std::vector<TH1D*> > refold_p;


	//1D unfold still has to convert the unfold1d and refold1d to vector
	if(pI == 0){
	  unfold_p.push_back({});
	  refold_p.push_back({});

	  unfold_p[0].push_back(nullptr);
	  refold_p[0].push_back(nullptr);

	  for(int i = 1; i < nIterForLoop+1; ++i){
	    unfold_p[0].push_back(unfold1D_p[i]);
	    refold_p[0].push_back(refold1D_p[i]);
	  }
	}
	else{
	  reco_p.clear();
	  truth_p.clear();

	  for(Int_t gI = 0; gI < nGammaPtBinsForUnfold; ++gI){
	    reco_p.push_back(new TH1D(("reco_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBins, varBins));
	    truth_p.push_back(new TH1D(("truth_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBins, varBins));

	    reco_p[gI]->GetXaxis()->SetTitle(reco2D_p->GetXaxis()->GetTitle());
	    for(Int_t bIX = 0; bIX < reco2D_p->GetXaxis()->GetNbins(); ++bIX){
	      if(isMultijet){
		//Y-axis no longer necessarily gammaPtBins if isMultijet
		for(Int_t bIY = 0; bIY < reco2D_p->GetYaxis()->GetNbins(); ++bIY){
		  int gammaPtBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		  if(gammaPtBinPos != gI) continue;

		  if(doJtPtCut){
		    int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
		    Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;
		    if(subJtPtBinCenter < jtPtBinsLowReco) continue;
		  }

		  Float_t value = reco_p[gI]->GetBinContent(bIX+1);
		  Float_t error = reco_p[gI]->GetBinError(bIX+1);

		  value += reco2D_p->GetBinContent(bIX+1, bIY+1);
		  error = TMath::Sqrt(error*error + reco2D_p->GetBinError(bIX+1, bIY+1)*reco2D_p->GetBinError(bIX+1, bIY+1));

		  reco_p[gI]->SetBinContent(bIX+1, value);
		  reco_p[gI]->SetBinError(bIX+1, error);
		}
	      }
	      else{
		reco_p[gI]->SetBinContent(bIX+1, reco2D_p->GetBinContent(bIX+1, gI+1));
		reco_p[gI]->SetBinError(bIX+1, reco2D_p->GetBinError(bIX+1, gI+1));
	      }
	    }

	    //Repeat for Truth
	    for(Int_t bIX = 0; bIX < truth2D_p->GetXaxis()->GetNbins(); ++bIX){
	      if(isMultijet){
		//Y-axis no longer necessarily gammaPtBins if isMultijet
                for(Int_t bIY = 0; bIY < truth2D_p->GetYaxis()->GetNbins(); ++bIY){
		  int gammaPtBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		  if(gammaPtBinPos != gI) continue;

		  if(doJtPtCut){
		    int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
		    Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;
		    if(subJtPtBinCenter < jtPtBinsLowReco) continue;
		  }

		  Float_t value = truth_p[gI]->GetBinContent(bIX+1);
		  Float_t error = truth_p[gI]->GetBinError(bIX+1);

		  value += truth2D_p->GetBinContent(bIX+1, bIY+1);
		  error = TMath::Sqrt(error*error + truth2D_p->GetBinError(bIX+1, bIY+1)*truth2D_p->GetBinError(bIX+1, bIY+1));

		  truth_p[gI]->SetBinContent(bIX+1, value);
		  truth_p[gI]->SetBinError(bIX+1, error);
		}
	      }
	      else{
		truth_p[gI]->SetBinContent(bIX+1, truth2D_p->GetBinContent(bIX+1, gI+1));
		truth_p[gI]->SetBinError(bIX+1, truth2D_p->GetBinError(bIX+1, gI+1));
	      }
	    }

	    unfold_p.push_back({});
	    refold_p.push_back({});

	    unfold_p[gI].push_back(nullptr);
	    refold_p[gI].push_back(nullptr);

	    for(int i = 1; i < nIterForLoop+1; ++i){
	      unfold_p[gI].push_back(new TH1D(("unfold_GammaPt" + std::to_string(gI) + "_Iter" + std::to_string(i) + "_h").c_str(), "", nVarBins, varBins));
	      refold_p[gI].push_back(new TH1D(("refold_GammaPt" + std::to_string(gI) + "_Iter" + std::to_string(i) + "_h").c_str(), "", nVarBins, varBins));

	      for(Int_t bIX = 0; bIX < unfold2D_p[i]->GetXaxis()->GetNbins(); ++bIX){
		if(isMultijet){
		  for(Int_t bIY = 0; bIY < unfold2D_p[i]->GetYaxis()->GetNbins(); ++bIY){
		    int gammaPtBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		    if(gammaPtBinPos != gI) continue;


		    if(doJtPtCut){
		      int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
		      Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;
		      if(subJtPtBinCenter < jtPtBinsLowReco) continue;
		    }

		    Float_t value = unfold_p[gI][i]->GetBinContent(bIX+1);
		    Float_t error = unfold_p[gI][i]->GetBinError(bIX+1);

		    value += unfold2D_p[i]->GetBinContent(bIX+1, bIY+1);
		    error = TMath::Sqrt(error*error + unfold2D_p[i]->GetBinError(bIX+1, bIY+1)*unfold2D_p[i]->GetBinError(bIX+1, bIY+1));

		    unfold_p[gI][i]->SetBinContent(bIX+1, value);
		    unfold_p[gI][i]->SetBinError(bIX+1, error);
		  }
		}
		else{
		  unfold_p[gI][i]->SetBinContent(bIX+1, unfold2D_p[i]->GetBinContent(bIX+1, gI+1));
		  unfold_p[gI][i]->SetBinError(bIX+1, unfold2D_p[i]->GetBinError(bIX+1, gI+1));
		}
	      }

	      for(Int_t bIX = 0; bIX < refold2D_p[i]->GetXaxis()->GetNbins(); ++bIX){
		if(isMultijet){
		  for(Int_t bIY = 0; bIY < refold2D_p[i]->GetYaxis()->GetNbins(); ++bIY){
		    int gammaPtBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
		    if(gammaPtBinPos != gI) continue;

		    if(doJtPtCut){
		      int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
		      Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;
		      if(subJtPtBinCenter < jtPtBinsLowReco) continue;
		    }

		    Float_t value = refold_p[gI][i]->GetBinContent(bIX+1);
		    Float_t error = refold_p[gI][i]->GetBinError(bIX+1);

		    value += refold2D_p[i]->GetBinContent(bIX+1, bIY+1);
		    error = TMath::Sqrt(error*error + refold2D_p[i]->GetBinError(bIX+1, bIY+1)*refold2D_p[i]->GetBinError(bIX+1, bIY+1));

		    refold_p[gI][i]->SetBinContent(bIX+1, value);
		    refold_p[gI][i]->SetBinError(bIX+1, error);
		  }
		}
		else{
		  refold_p[gI][i]->SetBinContent(bIX+1, refold2D_p[i]->GetBinContent(bIX+1, gI+1));
		  refold_p[gI][i]->SetBinError(bIX+1, refold2D_p[i]->GetBinError(bIX+1, gI+1));
		}
	      }
	    }
	  }
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pI << std::endl;

	for(Int_t gI = 0; gI < nGammaPtBinsForUnfold; ++gI){
	  //First bin is left intentionally empty on the reco side so just skip it
	  if(pI == 1 && gI == 0) continue;
	  if(pI == 1 && gI == nGammaPtBinsForUnfold-1) continue;

	  //New TCanvs, set of pads, and set margins
	  canv_p = new TCanvas("canv_p", "", width, height);
	  TCanvas* canvBest_p = new TCanvas("canvBest_p", "", width, height);
	  TPad* pads_p[nPad];
	  TPad* padsBest_p[nPad];
	  setMargins(canv_p, noMargin, noMargin, noMargin, noMargin);
	  setMargins(canvBest_p, noMargin, noMargin, noMargin, noMargin);

	  pads_p[0] = new TPad("pad0", "", 0.0, padSplit1, 1.0, 1.0);
	  padsBest_p[0] = new TPad("padBest0", "", 0.0, padSplit1, 1.0, 1.0);
	  setMargins(pads_p[0], leftMargin, topMargin, rightMargin, noMargin);
	  setMargins(padsBest_p[0], leftMargin, topMargin, rightMargin, noMargin);

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  pads_p[1] = new TPad("pad1", "", 0.0, padSplit2, 1.0, padSplit1);
	  padsBest_p[1] = new TPad("padBest1", "", 0.0, padSplit2, 1.0, padSplit1);
	  setMargins(pads_p[1], leftMargin, noMargin, rightMargin, noMargin);
	  setMargins(padsBest_p[1], leftMargin, noMargin, rightMargin, noMargin);

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  pads_p[2] = new TPad("pad2", "", 0.0, 0.0, 1.0, padSplit2);
	  padsBest_p[2] = new TPad("padBest2", "", 0.0, 0.0, 1.0, padSplit2);
	  setMargins(pads_p[2], leftMargin, noMargin, rightMargin, bottomMargin);
	  setMargins(padsBest_p[2], leftMargin, noMargin, rightMargin, bottomMargin);

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  std::vector<TPad*> padsVect = {pads_p[0], pads_p[1], pads_p[2]};
	  std::vector<TPad*> padsBestVect = {padsBest_p[0], padsBest_p[1], padsBest_p[2]};

	  //Draw pads to canvas
	  canv_p->cd();
	  for(Int_t pI = 0; pI < nPad; ++pI){pads_p[pI]->Draw("SAME");}

	  canvBest_p->cd();
	  for(Int_t pI = 0; pI < nPad; ++pI){padsBest_p[pI]->Draw("SAME");}

	  canv_p->cd();
	  pads_p[0]->cd();

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  if(truth_p[gI] != nullptr){
	    HIJet::Style::EquipHistogram(truth_p[gI], 2); //TRUTH IS kAzure-3 lines (atlas style picked up here)
	    truth_p[gI]->SetMarkerSize(0.00001);
	  }
	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", pI=" << pI << std::endl;

	  reco_p[gI]->SetMarkerSize(0.00001);

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  label_p = new TLatex();
	  initLabel(label_p, titleFont, titleSize/(1.0 - padSplit1), jetVarLabelAlignRightGlobal);

	  std::vector<TH1D*> histsForLegBest_p;
	  std::vector<std::string> legFillStr, legFillStrBest, stringsForLegBest_p;
	  if(truth_p[gI] != nullptr){
	    if(isMC) stringsForLeg_p = {"Truth", "Reco."};
	    else stringsForLeg_p = {"Truth PYTHIA", "Reco. Data"};
	    histsForLeg_p = {truth_p[gI], reco_p[gI]};
	    legFillStr = {"L", "L"};
	  }
	  else{
	    stringsForLeg_p = {"Reco."};
	    histsForLeg_p = {reco_p[gI]};
	    legFillStr = {"L"};
	  }
	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  histsForLegBest_p = histsForLeg_p;
	  stringsForLegBest_p = stringsForLeg_p;
	  legFillStrBest = legFillStr;

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    stringsForLeg_p.push_back("Unfolded, " + std::to_string(i));
	    histsForLeg_p.push_back(unfold_p[gI][i]);
	    legFillStr.push_back("P L");

	    if(i == termPos){
	      stringsForLegBest_p.push_back("Unfolded, " + std::to_string(i));
	      histsForLegBest_p.push_back(unfold_p[gI][i]);
	      legFillStrBest.push_back("P L");
	    }
	  }
	  pads_p[0]->cd();

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  if(!setupLegends(titleFont, titleSize/(1.0 - padSplit1), &leg_p, histsForLeg_p, stringsForLeg_p, legFillStr, &jetVarLegNGlobal, jetVarLegXGlobal, jetVarLegYGlobal)) return 1;

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  padsBest_p[0]->cd();
	  if(!setupLegends(titleFont, titleSize/(1.0 - padSplit1), &legBest_p, histsForLegBest_p, stringsForLegBest_p, legFillStrBest, &jetVarLegNGlobal, jetVarLegXGlobal, jetVarLegYGlobal)) return 1;

	  pads_p[0]->cd();

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  if(doBinWidthAndSelfNorm) binWidthAndSelfNorm(reco_p[gI]);
	  prepTH1(reco_p[gI], titleFont, titleSize/(1.0-padSplit1), labelSize/(1.0 - padSplit1), 1, reco_p[gI]->GetMarkerStyle(), reco_p[gI]->GetMarkerSize(), 2, reco_p[gI]->GetXaxis()->GetTitleOffset(), yOffset);

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  reco_p[gI]->GetXaxis()->SetTitleSize(0.00001);
	  reco_p[gI]->GetXaxis()->SetLabelSize(0.00001);
	  reco_p[gI]->GetYaxis()->SetNdivisions(505);
	  reco_p[gI]->GetXaxis()->SetNdivisions(505);

	  Float_t recoMin = getMinGTZero(reco_p[gI]);
	  if(truth_p[gI] != nullptr){
	    if(doBinWidthAndSelfNorm) binWidthAndSelfNorm(truth_p[gI]);

	    Float_t tempMin = getMinGTZero(truth_p[gI]);
	    if(tempMin < recoMin) recoMin = tempMin;
	  }

	  if(pI == 1){

	    if(doBinWidthAndSelfNorm || true){
	      reco_p[gI]->SetMinimum(jetVarMinGlobal);
	      reco_p[gI]->SetMaximum(jetVarMaxGlobal);
	    }
	    else{//Quickly calculate a proper maximum
	      double maxVal = 0.0;
	      if(reco_p[gI]->GetMaximum() > maxVal) maxVal = reco_p[gI]->GetMaximum();
	      if(truth_p[gI] != nullptr){
		if(truth_p[gI]->GetMaximum() > maxVal) maxVal = truth_p[gI]->GetMaximum();
		for(Int_t i = 1; i < nIterForLoop+1; ++i){
		  if(unfold_p[gI][i]->GetMaximum() > maxVal) maxVal = unfold_p[gI][i]->GetMaximum();
		}
	      }

	      reco_p[gI]->SetMaximum(maxVal*1.1);
	    }
	  }
	  else reco_p[gI]->SetMinimum(recoMin/2.0);

	  if(doBinWidthAndSelfNorm){
	    if(pI == 0) reco_p[gI]->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN_{#gamma}}{dp_{T}}");
	    else if(varNameLower.find("jj") != std::string::npos) reco_p[gI]->GetYaxis()->SetTitle(("#frac{1}{N_{JJ#gamma}} #frac{dN_{JJ#gamma}}{d" + varNameLabel + "}").c_str());
	    else reco_p[gI]->GetYaxis()->SetTitle(("#frac{1}{N_{J#gamma}} #frac{dN_{J#gamma}}{d" + varNameLabel + "}").c_str());
	  }
	  else{
	    reco_p[gI]->GetYaxis()->SetTitle("Weighted Counts");
	  }

	  reco_p[gI]->DrawCopy("HIST E1");

	  if(doLogX || pI == 0) gPad->SetLogx();
	  if(doLogY || pI == 0) gPad->SetLogy();

	  if(truth_p[gI] != nullptr) truth_p[gI]->DrawCopy("HIST E1 SAME");

	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    HIJet::Style::EquipHistogram(unfold_p[gI][i], i-1);
	    if(doBinWidthAndSelfNorm) binWidthAndSelfNorm(unfold_p[gI][i]);

	    HIJet::Style::EquipHistogram(refold_p[gI][i], i-1);
	    if(doBinWidthAndSelfNorm) binWidthAndSelfNorm(refold_p[gI][i]);

	    refold_p[gI][i]->SetMarkerSize(0.00001);
	    refold_p[gI][i]->SetLineStyle(2);
	  }
	  for(Int_t i = 1; i < nIterForLoop+1; ++i){unfold_p[gI][i]->DrawCopy("HIST E1 P SAME");}

	  canvBest_p->cd();
	  padsBest_p[0]->cd();

	  reco_p[gI]->DrawCopy("HIST E1");
	  if(truth_p[gI] != nullptr) truth_p[gI]->DrawCopy("HIST E1 SAME");

	  if(doLogX || pI == 0) gPad->SetLogx();
	  if(doLogY || pI == 0) gPad->SetLogy();

	  std::vector<double> deltaVals;
	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    bool isGood = false;
	    if(i == termPos) isGood = true;
	    if(!isGood) continue;

	    HIJet::Style::EquipHistogram(unfold_p[gI][i], i-1);
	    HIJet::Style::EquipHistogram(refold_p[gI][i], i-1);
	    refold_p[gI][i]->SetMarkerSize(0.00001);
	    refold_p[gI][i]->SetLineStyle(2);

	    unfold_p[gI][i]->DrawCopy("HIST E1 P SAME");
	    refold_p[gI][i]->DrawCopy("HIST E1 SAME");

	    if(termPos == i){
	      Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);

	      Int_t minVal = TMath::Max(i-1, 1);
	      Int_t maxVal = TMath::Min(i+1, nIterForLoop-1);

	      //Just give it the same thing if there is no lower iteration/higher iterations via minVal/maxVal
	      drawBoxAdjacentIter(padsBest_p[0], {unfold_p[gI][minVal], unfold_p[gI][i], unfold_p[gI][maxVal]}, tempOpacity, &deltaVals);
	    }
	  }

	  canv_p->cd();
	  pads_p[1]->cd();

	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    TH1D* unfoldClone_p = (TH1D*)unfold_p[gI][i]->Clone("tempClone");

	    unfoldClone_p->Divide(reco_p[gI]);
	    unfoldClone_p->SetTitle("");

	    if(i == 1){
	      unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Reco.}");
	      prepTH1(unfoldClone_p, titleFont, titleSize/(padSplit1 - padSplit2), labelSize/(padSplit1 - padSplit2), unfoldClone_p->GetMarkerColor(), unfoldClone_p->GetMarkerStyle(), unfoldClone_p->GetMarkerSize(), unfoldClone_p->GetLineWidth(), 1.2, yOffset*(padSplit1 - padSplit2)/(1.0-padSplit1));

	      unfoldClone_p->GetYaxis()->SetNdivisions(404);
	      unfoldClone_p->SetMaximum(1.25);
	      unfoldClone_p->SetMinimum(0.75);
	      unfoldClone_p->DrawCopy("HIST E1 P");
	    }
	    else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	    delete unfoldClone_p;
	  }
	  if(doLogX || pI == 0) gPad->SetLogx();

	  canvBest_p->cd();
	  padsBest_p[1]->cd();

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  bool hasDrawn = false;
	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    bool isGood = false;
	    if(i == termPos) isGood = true;
	    if(!isGood) continue;

	    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << std::endl;

	    TH1D* unfoldClone_p = (TH1D*)unfold_p[gI][i]->Clone("tempClone");
	    unfoldClone_p->Divide(reco_p[gI]);

	    TH1D* refoldClone_p = (TH1D*)refold_p[gI][i]->Clone("tempClone2");
	    refoldClone_p->Divide(reco_p[gI]);

	    unfoldClone_p->SetTitle("");
	    refoldClone_p->SetTitle("");

	    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << std::endl;

	    if(!hasDrawn){
	      refoldClone_p->GetYaxis()->SetTitle("#frac{Refold}{Reco.}");

	      prepTH1(refoldClone_p, titleFont, titleSize/(padSplit1 - padSplit2), labelSize/(padSplit1 - padSplit2), refoldClone_p->GetMarkerColor(), refoldClone_p->GetMarkerStyle(), refoldClone_p->GetMarkerSize(), refoldClone_p->GetLineWidth(), 1.2, yOffset*(padSplit1 - padSplit2)/(1.0-padSplit1));

	      refoldClone_p->GetYaxis()->SetNdivisions(404);

	      refoldClone_p->SetMaximum(1.25);
	      refoldClone_p->SetMinimum(0.75);

	      refoldClone_p->DrawCopy("HIST E1 P");
	      hasDrawn = true;
	    }
	    else refoldClone_p->DrawCopy("HIST E1 P SAME");

	    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << std::endl;

	    if(i == termPos){
	      Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);

	      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << ", " << termPos << ", " << deltaVals.size() << std::endl;

	      drawBoxAdjacentIterRel(padsBest_p[1], refoldClone_p, tempOpacity, &deltaVals);
	      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << std::endl;
	    }

	    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << ", " << i << "/" << nIterForLoop+1 << std::endl;

	    delete refoldClone_p;
	    delete unfoldClone_p;
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE, gI/nGammaPtBinsForUnfold: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBinsForUnfold << std::endl;

	  if(doLogX || pI == 0) gPad->SetLogx();

	  //Third panel, divide by truth
	  canv_p->cd();
	  pads_p[2]->cd();

	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    TH1D* unfoldClone_p = (TH1D*)unfold_p[gI][i]->Clone("tempClone");

	    unfoldClone_p->Divide(truth_p[gI]);
	    unfoldClone_p->SetTitle("");

	    if(i == 1){
	      if(isMC) unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Truth}");
	      else unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded Data}{Truth PYTHIA}");

	      std::string xTitle = reco_p[gI]->GetXaxis()->GetTitle();
	      if(xTitle.find("Reco. ") != std::string::npos) xTitle.replace(xTitle.find("Reco. "), 6, "");
	      unfoldClone_p->GetXaxis()->SetTitle(xTitle.c_str());

	      prepTH1(unfoldClone_p, titleFont, titleSize/(padSplit2), labelSize/(padSplit2), unfoldClone_p->GetMarkerColor(), unfoldClone_p->GetMarkerStyle(), unfoldClone_p->GetMarkerSize(), unfoldClone_p->GetLineWidth(), 1.2, yOffset*(padSplit2)/(1.0-padSplit1));

	      unfoldClone_p->GetYaxis()->SetNdivisions(404);
	      unfoldClone_p->SetMaximum(jetVarRatMaxGlobal);
	      unfoldClone_p->SetMinimum(jetVarRatMinGlobal);
	      unfoldClone_p->DrawCopy("HIST E1 P");
	    }
	    else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	    delete unfoldClone_p;
	  }
	  if(doLogX || pI == 0) gPad->SetLogx();

	  if(pI == 1) line_p->DrawLine(varBins[0], 1.0, varBins[nVarBins], 1.0);
	  else line_p->DrawLine(gammaPtBins[0], 1.0, gammaPtBins[nGammaPtBins], 1.0);

	  canvBest_p->cd();
	  padsBest_p[2]->cd();

	  hasDrawn = false;
	  for(Int_t i = 1; i < nIterForLoop+1; ++i){
	    bool isGood = false;
	    if(i == termPos) isGood = true;
	    if(!isGood) continue;

	    TH1D* unfoldClone_p = (TH1D*)unfold_p[gI][i]->Clone(("unfoldOverTruthMC_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_MCNonClosure").c_str());
	    unfoldClone_p->Divide(truth_p[gI]);

	    //Write out the ratio of these two to file for systematic
	    if(systI == 0 && isMC){
	      outSysFile_p->cd();
	      unfoldClone_p->Write("", TObject::kOverwrite);
	      canvBest_p->cd();
	      padsBest_p[2]->cd();
	    }

	    if(!hasDrawn){
	      if(isMC) unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded}{Truth}");
	      else unfoldClone_p->GetYaxis()->SetTitle("#frac{Unfolded Data}{Truth PYTHIA}");

	      prepTH1(unfoldClone_p, titleFont, titleSize/(padSplit2), labelSize/(padSplit2), unfoldClone_p->GetMarkerColor(), unfoldClone_p->GetMarkerStyle(), unfoldClone_p->GetMarkerSize(), unfoldClone_p->GetLineWidth(), 1.2, yOffset*(padSplit2)/(1.0-padSplit1));

	      unfoldClone_p->GetYaxis()->SetNdivisions(404);

	      unfoldClone_p->SetMaximum(jetVarRatMaxGlobal);
	      unfoldClone_p->SetMinimum(jetVarRatMinGlobal);

	      std::string xTitle = reco_p[gI]->GetXaxis()->GetTitle();
	      if(xTitle.find("Reco. ") != std::string::npos) xTitle.replace(xTitle.find("Reco. "), 6, "");
	      unfoldClone_p->GetXaxis()->SetTitle(xTitle.c_str());
	      unfoldClone_p->DrawCopy("HIST E1 P");
	      hasDrawn = true;
	    }
	    else unfoldClone_p->DrawCopy("HIST E1 P SAME");

	    if(i == termPos){
	      Float_t tempOpacity = HIJet::Style::GetOpacity(i-1);
	      drawBoxAdjacentIterRel(padsBest_p[2], unfoldClone_p, tempOpacity, &deltaVals);
	    }

	    delete unfoldClone_p;
	  }
	  if(doLogX || pI == 0) gPad->SetLogx();

	  if(pI == 1) line_p->DrawLine(varBins[0], 1.0, varBins[nVarBins], 1.0);
	  else line_p->DrawLine(gammaPtBins[0], 1.0, gammaPtBins[nGammaPtBins], 1.0);

	  canv_p->cd();
	  pads_p[0]->cd();
	  for(unsigned int i = 0; i < leg_p.size(); ++i){
	    if(leg_p[i] != nullptr){
	      leg_p[i]->Draw("SAME");
	    }
	  }

	  int pos = 0;
	  for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	    if(pI == 0 && tempLabels[tI].find("jet") != std::string::npos) continue;

	    if(tempLabels[tI].find("Observable") != std::string::npos) continue;

	    if(tempLabels[tI].find("< p_{T,#gamma} <") != std::string::npos && pI == 1){
	      tempLabels[tI].replace(0, tempLabels[tI].find("<")+1, "");
	      tempLabels[tI].replace(tempLabels[tI].find("<"), tempLabels[tI].size(), "");

	      tempLabels[tI] = prettyString(gammaPtBins[gI], 1, false) + " <" + tempLabels[tI] + "< " + prettyString(gammaPtBins[gI+1], 1, false);
	    }


	    if(doJtPtCut){
	      if(tempLabels[tI].find("p_{T,jet}") != std::string::npos){
		tempLabels[tI].replace(0, tempLabels[tI].find("<"), "");
		tempLabels[tI] = prettyString(jtPtBinsLowReco, 1, false) + tempLabels[tI];
		//		std::cout << "RELEVANT LABEL: " << tempLabels[tI] << std::endl;
	      }
	    }

	    label_p->DrawLatex(jetVarLabelXGlobal, jetVarLabelYGlobal - pos*0.077*0.55/topPanelFrac, tempLabels[tI].c_str());
	    ++pos;
	  }

	  canvBest_p->cd();
	  padsBest_p[0]->cd();
	  for(unsigned int i = 0; i < legBest_p.size(); ++i){
	    if(legBest_p[i] != nullptr) legBest_p[i]->Draw("SAME");
	  }

	  pos = 0;
	  for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	    if(pI == 0 && tempLabels[tI].find("jet") != std::string::npos) continue;

	    label_p->DrawLatex(jetVarLabelXGlobal, jetVarLabelYGlobal - pos*0.077*0.55/topPanelFrac, tempLabels[tI].c_str());

	    ++pos;
	  }

	  tpadsAxisFix(padsVect);
	  std::string saveName = "pdfDir/" + dateStr + "/unfold";
	  if(pI == 0) saveName = saveName + "PhoPt";
	  else if(pI == 1) saveName = saveName + "PhoPtJet" + varName + "_GammaPt" + std::to_string(gI);
	  saveName = saveName + "_" + centBinsStr[cI%centBinsStr.size()];
	  if(pI == 1) saveName = saveName + "_" + systStrVect[systI];
	  saveName = saveName + "_Delta" + deltaStr + "_" + saveTag +  "." + saveExt;

	  quietSaveAs(canv_p, saveName);
	  for(Int_t i = 0; i < nPad; ++i){delete pads_p[i];}
	  delete canv_p;

	  tpadsAxisFix(padsBestVect);

	  saveName = "pdfDir/" + dateStr + "/unfold";
	  if(pI == 0) saveName = saveName + "PhoPt";
	  else if(pI == 1) saveName = saveName + "PhoPtJet" + varName + "_GammaPt" + std::to_string(gI);
	  saveName = saveName + "_BESTIterations_" + centBinsStr[cI%centBinsStr.size()];
	  if(pI == 1) saveName = saveName + "_" + systStrVect[systI];
	  saveName = saveName + "_Delta" + deltaStr + "_" + saveTag +  "." + saveExt;

	  quietSaveAs(canvBest_p, saveName);
	  for(Int_t i = 0; i < nPad; ++i){delete padsBest_p[i];}
	  delete canvBest_p;
	}

	delete statsDelta_p;
	delete iterDelta_p;
	delete totalDelta_p;

	if(pI == 1){
	  for(unsigned int i = 0; i < reco_p.size(); ++i){
	    delete reco_p[i];
	    delete truth_p[i];
	  }

	  for(unsigned int i = 0; i < unfold_p.size(); ++i){
	    for(unsigned int j = 0; j < unfold_p[i].size(); ++j){
	      delete unfold_p[i][j];
	      delete refold_p[i][j];
	    }
	  }

	  reco_p.clear();
	  truth_p.clear();
	}//End if(pI == 1)
      }//End for(unsigned int cI = 0; cI < statsDeltaNames.size();
    }//End for(Int_t pI = 0; pI < nProc;...
  }//End for(Int_t systI = 0; systI < nSyst;....

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

  outSysFile_p->Close();
  delete outSysFile_p;

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
