//Author: Chris McGinn (2022.11.21)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/varUtil.h"

int gdjHEPMCPlot(std::string inConfigFileName)
{
  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Global debug handler - do 'export DOGLOBALDEBUGROOT=1' at command line to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string pdfStr = "pdfDir/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INPPFILENAME",
					      "EXTENSION",
					      "MINY",
					      "MAXY",
					      "MINX",
					      "MAXX",
					      "RATMIN",
					      "RATMAX",
					      "LABELX",
					      "LABELY",
					      "LABELALIGNRIGHT",
					      "LEGX",
					      "LEGY"};
  
  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  const std::string inPPFileName = config_p->GetValue("INPPFILENAME", "");
  const std::string extStr = config_p->GetValue("EXTENSION", "");

  //Check extension given is valid
  std::vector<std::string> validExtensions = {"pdf", "png"};
  if(!vectContainsStr(extStr, &validExtensions)){
    std::cout << "Given extension \'" << extStr << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  //Grab the Pb+Pb files - we did not check for these initially so just check vector size at end
  const Int_t nMaxPbPbFiles = 10;
  std::vector<std::string> inPbPbFileNames;
  for(Int_t i = 0; i < nMaxPbPbFiles; ++i){
    std::string inPbPbFileName = config_p->GetValue(("INPBPBFILENAME." + std::to_string(i)).c_str(), "");
    if(check.checkFileExt(inPbPbFileName, ".root")) inPbPbFileNames.push_back(inPbPbFileName);
  }

  if(inPbPbFileNames.size() == 0){
    std::cout << "Config file \'" << inConfigFileName << "\' contains no valid INPBPBFILENAME.[0-9]. return 1" << std::endl;
    return 1;
  }

  //We will probably need to do some bin manipulations so grab these params
  float xMin = config_p->GetValue("MINX", -1000.0);
  float xMax = config_p->GetValue("MAXX", -10000.0);
  bool doXTrunc = xMin < xMax;
   
  //Check pp file given is valid
  if(!check.checkFileExt(inPPFileName, ".root")) return 1;

  TFile* inPPFile_p = new TFile(inPPFileName.c_str(), "READ");
  TEnv* ppConfig_p = (TEnv*)inPPFile_p->Get("config");
  TFile* inPbPbFile_p[nMaxPbPbFiles];
  TEnv* pbpbConfig_p[nMaxPbPbFiles];

  std::vector<std::string> matchParams = {"NGAMMAPTBINS",
					  "GAMMAPTBINSLOW",
					  "GAMMAPTBINSHIGH",
					  "GAMMAPTBINSDOLOG",
					  "GAMMAPTBINSDOCUSTOM",
					  "JTPTBINSLOWRECO",
					  "JTPTBINSHIGHRECO",
					  "JTETABINSLOW",
					  "JTETABINSHIGH",
					  "GAMMAEXCLUSIONDR",
					  "MIXJETEXCLUSIONDR",
					  "GAMMAJTDPHI",
					  "GAMMAMULTIJTDPHI",
					  "VARNAME"};
  
  for(unsigned int i = 0; i < inPbPbFileNames.size(); ++i){
    inPbPbFile_p[i] = new TFile(inPbPbFileNames[i].c_str(), "READ");
    pbpbConfig_p[i] = (TEnv*)inPbPbFile_p[i]->Get("config");

    if(!compEnvParams(ppConfig_p, pbpbConfig_p[i], matchParams)) return 1;
  }
  
  //Get gammaptbins, after defining a bin max
  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = ppConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = ppConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = ppConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Float_t gammaPtBinsLowReco = ppConfig_p->GetValue("GAMMAPTBINSLOWRECO", 80.0);
  const Float_t gammaPtBinsHighReco = ppConfig_p->GetValue("GAMMAPTBINSHIGHRECO", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)ppConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)ppConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];
 
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string gammaPtBinsStr = ppConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBins != nGammaPtBinsTemp){
      std::cout << "nGammaPtBins dont match - by stored value \'" << nGammaPtBins << "\', via vector \'\
" << nGammaPtBinsTemp << "\'" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  
  //Get Reco PT bounds
  const Float_t jtPtBinsLowReco = ppConfig_p->GetValue("JTPTBINSLOWRECO", -999.0);
  const Float_t jtPtBinsHighReco = ppConfig_p->GetValue("JTPTBINSHIGHRECO", -999.0);
  const Float_t jtEtaBinsLow = ppConfig_p->GetValue("JTETABINSLOW", -999.0);
  const Float_t jtEtaBinsHigh = ppConfig_p->GetValue("JTETABINSHIGH", -999.0);

  //Get Var name and derivatives
  std::string varName = ppConfig_p->GetValue("VARNAME", "");
  std::string varNameUpper = strLowerToUpper(varName);
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string varNameLabel = varNameToLabel(varName);

  std::string varPrefix = "";
  if(isStrSame(varNameUpper, "PT")) varPrefix = "JT";
  std::string binVarStr = varPrefix + varNameUpper;
  if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";
  else if(isStrSame(binVarStr, "DRJJ")) binVarStr = "DR";
  
  //Get Var Bins
  Int_t nVarBins = ppConfig_p->GetValue(("N" + binVarStr + "BINS").c_str(), -1);
  Float_t varBinsLow = ppConfig_p->GetValue((binVarStr + "BINSLOW").c_str(), -1.0);
  Float_t varBinsHigh = ppConfig_p->GetValue((binVarStr + "BINSHIGH").c_str(), -1.0);
  Bool_t varBinsDoLog = ppConfig_p->GetValue((binVarStr + "BINSDOLOG").c_str(), 0);
  Bool_t varBinsDoCustom = ppConfig_p->GetValue((binVarStr + "BINSDOCUSTOM").c_str(), 0);
  Float_t varBinsLowReco = ppConfig_p->GetValue((binVarStr + "BINSLOWRECO").c_str(), varBinsLow);
  Float_t varBinsHighReco = ppConfig_p->GetValue((binVarStr + "BINSHIGHRECO").c_str(), varBinsHigh);
  Double_t varBins[nMaxPtBins+1];
  Int_t nVarBinsTrunc = -1;
  Double_t varBinsTrunc[nMaxPtBins+1];
  std::string varBinsStr = "";
  
  if(varBinsDoLog) getLogBins(varBinsLow, varBinsHigh, nVarBins, varBins);
  else if(varBinsDoCustom){
    varBinsStr = ppConfig_p->GetValue((binVarStr + "BINSCUSTOM").c_str(), "");
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

  //Now that we have the varbins, quickly check that the requested truncation makes sense if doXTrunc
  if(doXTrunc){
    if(!checkHistContainsBins({xMin, xMax}, nVarBins, varBins, 0.0001)){
      std::cout << "gdjPlotResults ERROR: Requested plotting truncation in var \'" << varNameUpper << "\' to range " << xMin << "-" << xMax << " is not found in bins \'" << varBinsStr << "\'. return 1" << std::endl;
      return 1;
    }
  
    //If it passes construct the subbin array
    nVarBinsTrunc = truncateBinRange(nVarBins, varBins, xMin, xMax, 0.0001, varBinsTrunc);
  }

  //get specific gamma-jet cut values
  const Double_t gammaJtDRExclusionCut = ppConfig_p->GetValue("GAMMAEXCLUSIONDR", 100.0);
  const Double_t mixJtDRExclusionCut = ppConfig_p->GetValue("MIXJETEXCLUSIONDR", 100.0);
  const std::string gammaJtDPhiCutStr = ppConfig_p->GetValue("GAMMAJTDPHI", "");
  const Double_t gammaJtDPhiCut = mathStringToNum(gammaJtDPhiCutStr, doGlobalDebug);
  const std::string gammaMultiJtDPhiCutStr = ppConfig_p->GetValue("GAMMAMULTIJTDPHI", "");
  const Double_t gammaMultiJtDPhiCut = mathStringToNum(gammaMultiJtDPhiCutStr, doGlobalDebug);  
  
  //multijet defined by var name
  const Bool_t isMultijet = isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");
  
  const Int_t nMaxJetRs = 7;
  std::string jetRStr = ppConfig_p->GetValue("JETRS", "");
  std::vector<int> jetRVect = strToVectI(jetRStr);
  if(jetRVect.size() <= 0 || jetRVect.size() > (unsigned int)nMaxJetRs){
    std::cout << "Jet R query returned \'" << jetRStr << "\', either 0 valid jet Rs or in excess of cap, " << nMaxJetRs << ". return 1" << std::endl;
    return 1;
  }

  //Define a bunch of plotting parameters
  const Int_t nPad = 2;
  Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;
  Double_t height = 900.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);
  height = 1200.0;

  //Also define some fonts and sizes
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;

  //Define some global labels
  std::vector<std::string> labels = {"JEWEL", "Pb+Pb, #it{pp} at #sqrt{s_{NN}}=5.02 TeV"};    
  
  //Grab remaining config params for plotting
  float minVal = config_p->GetValue("MINY", 0.0);
  float maxVal = config_p->GetValue("MAXY", 0.0);
  float minRatVal = config_p->GetValue("RATMIN", 0.0);
  float maxRatVal = config_p->GetValue("RATMAX", 0.0);
  float labelX = config_p->GetValue("LABELX", 0.0);
  float labelY = config_p->GetValue("LABELY", 0.0);
  bool labelAlignRight = config_p->GetValue("LABELALIGNRIGHT", 0);
  float legX = config_p->GetValue("LEGX", 0.0);
  float legY = config_p->GetValue("LEGY", 0.0);
  
  
  //Loop to create basic plots, pbpb one centrality bin w/ pp, and ratio in second subpanel
  //loop over jet R values
  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    //loop over gamma pt bins, skipping those outside of acceptance    
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
      Float_t gammaBinCenter = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;
      if(gammaBinCenter < gammaPtBinsLowReco) continue;
      if(gammaBinCenter > gammaPtBinsHighReco) continue;

      //Grab pp hist + define the style
      
      TH1D* inPPHist_p = (TH1D*)inPPFile_p->Get((varNameLower + "_R" + std::to_string(jetRVect[rI]) + "_GammaPt" + std::to_string(gI) + "_h").c_str());
      std::string titleStr = ";" + ((std::string)inPPHist_p->GetXaxis()->GetTitle()) + ";" + ((std::string)inPPHist_p->GetYaxis()->GetTitle());
      TH1D* ppHist_p = nullptr;
      if(!doXTrunc) ppHist_p = new TH1D("ppHistTemp_p", titleStr.c_str(), nVarBins, varBins);
      else ppHist_p = new TH1D("ppHistTemp_p", titleStr.c_str(), nVarBinsTrunc, varBinsTrunc);
      fineHistToCoarseHist(inPPHist_p, ppHist_p);

      
      HIJet::Style::EquipHistogram(ppHist_p, 2);
      ppHist_p->GetXaxis()->SetTitleFont(titleFont);
      ppHist_p->GetYaxis()->SetTitleFont(titleFont);
      ppHist_p->GetXaxis()->SetTitleSize(0.00001);
      ppHist_p->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));
      
      ppHist_p->GetXaxis()->SetLabelFont(titleFont);
      ppHist_p->GetYaxis()->SetLabelFont(titleFont);
      ppHist_p->GetXaxis()->SetLabelSize(0.00001);
      ppHist_p->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));
      
      ppHist_p->GetYaxis()->SetTitleOffset(yOffset);
      
      ppHist_p->GetYaxis()->SetNdivisions(505);

      ppHist_p->SetMinimum(minVal);
      ppHist_p->SetMaximum(maxVal);
      

      //iterate over input files (which will define the centrality
      for(unsigned int i = 0; i < inPbPbFileNames.size(); ++i){
	std::vector<std::string> tempLabels = labels;
	tempLabels.push_back("anti-k_{t} R=0." + std::to_string(jetRVect[rI]) + " jets");
	tempLabels.push_back(prettyString(gammaPtBins[gI],1,false) + "<p_{T,#gamma}<" + prettyString(gammaPtBins[gI+1], 1, false));
	tempLabels.push_back(prettyString(jtPtBinsLowReco,1,false) + "<p_{T,Jet}<" + prettyString(jtPtBinsHighReco, 1, false));
	tempLabels.push_back("|#eta_{Jet}|<" + prettyString(jtEtaBinsHigh, 1, false));
	tempLabels.push_back("|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiCutStr);
	tempLabels.push_back("#DeltaR_{JJ} > " + prettyString(mixJtDRExclusionCut, 1, false));
	tempLabels.push_back("|#Delta#phi_{#gamma,JJ}| > " + gammaMultiJtDPhiCutStr);

	for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
	  if(tempLabels[tI].find("pi") != std::string::npos) tempLabels[tI].replace(tempLabels[tI].find("pi"), 2, "#pi");
	}
	
	//Grab the centrality + create a cent label
	std::string centStr = pbpbConfig_p[i]->GetValue("CENT", "");
	if(centStr.size() == 0){
	  std::cout << "CENT in file \'" << inPbPbFileNames[i] << "\' is not found! Please fix. return 1" << std::endl;
	  return 1;
	}	
	std::string centLabel = centStr.substr(4, centStr.find("to") - 4);
	centLabel = centLabel + "-" + centStr.substr(centStr.find("to")+2, centStr.size()) + "%";

	//grab pbpb histogram

	TH1D* inPbPbHist_p = (TH1D*)inPbPbFile_p[i]->Get((varNameLower + "_R" + std::to_string(jetRVect[rI]) + "_GammaPt" + std::to_string(gI) + "_h").c_str());
	TH1D* pbpbHist_p = nullptr;
	if(!doXTrunc) pbpbHist_p = new TH1D("pbpbHistTemp_p", titleStr.c_str(), nVarBins, varBins);
	else pbpbHist_p = new TH1D("pbpbHistTemp_p", titleStr.c_str(), nVarBinsTrunc, varBinsTrunc);
	fineHistToCoarseHist(inPbPbHist_p, pbpbHist_p);
	
	TCanvas* canv_p = new TCanvas("canv", "", width, height);
	TPad* pads_p[nPad];
	setMargins(canv_p, noMargin, noMargin, noMargin, noMargin);
	canv_p->cd();
	pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
	setMargins(pads_p[0], leftMargin, topMargin, rightMargin, noMargin);

	canv_p->cd();
	pads_p[0]->Draw("SAME");

	canv_p->cd();
	pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
	setMargins(pads_p[1], leftMargin, noMargin, rightMargin, bottomMargin);

	canv_p->cd();
	pads_p[1]->Draw("SAME");

	pads_p[0]->cd();
	HIJet::Style::EquipHistogram(pbpbHist_p, 1);

	ppHist_p->DrawCopy("HIST E1 P");
	pbpbHist_p->DrawCopy("HIST E1 P SAME");      
	
	gPad->SetTicks();
	gStyle->SetOptStat(0);

	//Dump labels and legends
	TLatex* label_p = new TLatex();
	label_p->SetNDC();
	label_p->SetTextFont(titleFont);
	label_p->SetTextSize(titleSize/(1.0-padSplit));
	if(labelAlignRight) label_p->SetTextAlign(31);

	for(unsigned int lI = 0; lI < tempLabels.size(); ++lI){
	  label_p->DrawLatex(labelX, labelY - 0.083*lI, tempLabels[lI].c_str());
	}
	delete label_p;
	
	TLegend* leg_p = new TLegend(legX, legY - 0.065*2, legX+0.25, legY);
	leg_p->SetTextFont(42);
	leg_p->SetTextSize(0.035/(1.0 - padSplit));
	leg_p->SetBorderSize(0);
	leg_p->SetFillStyle(0);
	
	leg_p->AddEntry(pbpbHist_p, ("Pb+Pb #bf{" + centLabel + "}").c_str(), "P L");
	leg_p->AddEntry(ppHist_p, "p+p", "P L");
	leg_p->Draw("SAME");

	
	//Now construct and plot ratio
	pads_p[1]->cd();

	pbpbHist_p->Divide(ppHist_p);

	pbpbHist_p->GetXaxis()->SetTitleFont(titleFont);
	pbpbHist_p->GetYaxis()->SetTitleFont(titleFont);
	pbpbHist_p->GetXaxis()->SetTitleSize(titleSize/padSplit);
	pbpbHist_p->GetYaxis()->SetTitleSize(titleSize/padSplit);
	
	pbpbHist_p->GetXaxis()->SetLabelFont(titleFont);
	pbpbHist_p->GetYaxis()->SetLabelFont(titleFont);
	pbpbHist_p->GetXaxis()->SetLabelSize(labelSize/padSplit);
	pbpbHist_p->GetYaxis()->SetLabelSize(labelSize/padSplit);
	
	pbpbHist_p->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0-padSplit));
	pbpbHist_p->GetXaxis()->SetTitleOffset(1.2);
	
	pbpbHist_p->GetYaxis()->SetNdivisions(505);
	
	pbpbHist_p->SetMinimum(minRatVal);
	pbpbHist_p->SetMaximum(maxRatVal);

	pbpbHist_p->GetYaxis()->SetTitle("Pb+Pb / pp");
	
	pbpbHist_p->DrawCopy("HIST E1 P");
	
	std::string saveName = pdfStr + "/" + varNameLower + "_R" + std::to_string(jetRVect[rI]) + "_GammaPt" + std::to_string(gI) + "_" + centStr + "." + extStr;
	quietSaveAs(canv_p, saveName);
	for(Int_t pI = 0; pI < nPad; ++pI){
	  delete pads_p[pI];
	}
	delete canv_p;
	delete leg_p;
	delete pbpbHist_p;
      }

      delete ppHist_p;
    }
  }
  
  

  for(unsigned int i = 0; i < inPbPbFileNames.size(); ++i){
    inPbPbFile_p[i]->Close();
    delete inPbPbFile_p[i];
  }
  
  inPPFile_p->Close();
  delete inPPFile_p;
  
  std::cout << "GDJHEPMCPLOT complete. return 0" << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHEPMCPlot.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHEPMCPlot(argv[1]);
  return retVal;
}
