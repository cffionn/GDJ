//Author: Chris McGinn (2021.02.18)
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

int gdjPlotResults(std::string inConfigFileName)
{
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;
  
  cppWatch globalTimer;
  globalTimer.start();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INPBPBFILENAME",
					      "INPPFILENAME",
					      "INPBPBTERMFILENAME",
					      "INPPTERMFILENAME",
					      "SAVETAG",
					      "SAVEEXT",
					      "PTMIN",
					      "PTMAX",
					      "PTRATMIN",
					      "PTRATMAX",
					      "PTDOLOGX",
					      "PTDOLOGY",
  					      "XJMIN",
					      "XJMAX",
  					      "XJRATMIN",
					      "XJRATMAX",
					      "XJDOLOGX",
					      "XJDOLOGY",
					      "XJJMIN",
					      "XJJMAX",
					      "XJJRATMIN",
					      "XJJRATMAX",
  					      "DOXJJMIN",
					      "DOXJJMAX"};

  
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  const std::string inPbPbFileName = config_p->GetValue("INPBPBFILENAME", "");
  const std::string inPPFileName = config_p->GetValue("INPPFILENAME", "");

  const std::string inPbPbTermFileName = config_p->GetValue("INPBPBTERMFILENAME", "");
  const std::string inPPTermFileName = config_p->GetValue("INPPTERMFILENAME", "");

  TEnv* inPbPbTermFileConfig_p = new TEnv(inPbPbTermFileName.c_str());
  TEnv* inPPTermFileConfig_p = new TEnv(inPPTermFileName.c_str());
  
  std::vector<std::string> termNecessaryParams = {"UNFOLDFILENAME",
						  "VARNAME"};

  if(!checkEnvForParams(inPbPbTermFileConfig_p, termNecessaryParams)) return 1;
  if(!checkEnvForParams(inPPTermFileConfig_p, termNecessaryParams)) return 1;  
  if(!compEnvParams(inPbPbTermFileConfig_p, inPPTermFileConfig_p, {"VARNAME"})) return 1;

  const std::string inPbPbTempFileName = inPbPbTermFileConfig_p->GetValue("UNFOLDFILENAME", "");
  const std::string inPPTempFileName = inPPTermFileConfig_p->GetValue("UNFOLDFILENAME", "");
     
  if(!isStrSame(inPbPbTempFileName, inPbPbFileName)){
    std::cout << "Given unfolding termination file for Pb+Pb \'" << inPbPbTermFileName << "\' is based on \'" << inPbPbTempFileName << "\', not matching given input \'" << inPbPbFileName << "\'. return 1" << std::endl;
    return 1;
  }
  if(!isStrSame(inPPTempFileName, inPPFileName)){
    std::cout << "Given unfolding termination file for p+p \'" << inPPTermFileName << "\' is based on \'" << inPPTempFileName << "\', not matching given input \'" << inPPFileName << "\'. return 1" << std::endl;
    return 1;
  }
  
  const std::string termVarName = returnAllLowercaseString(inPbPbTermFileConfig_p->GetValue("VARNAME", ""));
  
  const std::string saveTag = config_p->GetValue("SAVETAG", "");
  const std::string saveExt = config_p->GetValue("SAVEEXT", "");
  
  std::vector<std::string> validExt = {"pdf", "png"};

  if(!vectContainsStr(saveExt, &validExt)){
    std::cout << "Given SAVEEXT \'" << saveExt << "\' is not valid. return 1" << std::endl;
    return 1;
  }  
  
  //Get global labels
  std::vector<std::string> globalLabels;
  for(unsigned int i = 0; i < 100; ++i){
    std::string tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
    if(tempStr.size() != 0) globalLabels.push_back(tempStr);        
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //check that our given input files exist
  if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;
  if(!check.checkFileExt(inPPFileName, ".root")) return 1;
  
  TFile* inPbPbFile_p = new TFile(inPbPbFileName.c_str(), "READ");
  TEnv* inPbPbFileConfig_p = (TEnv*)inPbPbFile_p->Get("config");    
  TEnv* inPbPbFileLabel_p = (TEnv*)inPbPbFile_p->Get("label");    

  TFile* inPPFile_p = new TFile(inPPFileName.c_str(), "READ");
  TEnv* inPPFileConfig_p = (TEnv*)inPPFile_p->Get("config");    
  //  TEnv* inPPFileLabel_p = (TEnv*)inPPFile_p->Get("label");    
  
  std::vector<std::string> paramsToCompare = {"ISMC",
					      "NGAMMAPTBINS",
					      "GAMMAPTBINSLOW",
					      "GAMMAPTBINSHIGH",
					      "GAMMAPTBINSDOLOG",
					      "GAMMAPTBINSDOCUSTOM",
					      "NITER",
					      "JETR",
					      "VARNAME",
					      "ASSOCGENMINPT"};

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> necessaryFileParams = paramsToCompare;
  necessaryFileParams.push_back("ISPP");

  const std::string varName = returnAllLowercaseString(inPbPbFileConfig_p->GetValue("VARNAME", ""));
  if(!isStrSame(termVarName, varName)){
    std::cout << "Varname of input files \'" << varName << "\' does not match the varName of termination files \'" << termVarName << "\'. return 1" << std::endl;
    return 1;
  }


  const std::string varNameUpper = strLowerToUpper(varName);

  float minVal = config_p->GetValue((varNameUpper + "MIN").c_str(), 0.0);
  float maxVal = config_p->GetValue((varNameUpper + "MAX").c_str(), 0.0);
  float minRatVal = config_p->GetValue((varNameUpper + "RATMIN").c_str(), 0.0);
  float maxRatVal = config_p->GetValue((varNameUpper + "RATMAX").c_str(), 0.0);
  bool doLogX = config_p->GetValue((varNameUpper + "DOLOGX").c_str(), 0);
  bool doLogY = config_p->GetValue((varNameUpper + "DOLOGY").c_str(), 0);
  
  const std::string gammaJtDPhiStr = "DPhi0";
  //  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

  std::vector<std::string> pbpbParams = {"CENTBINS"};
  if(!checkEnvForParams(inPbPbFileConfig_p, necessaryFileParams)) return 1;
  if(!checkEnvForParams(inPPFileConfig_p, necessaryFileParams)) return 1;
  if(!compEnvParams(inPbPbFileConfig_p, inPPFileConfig_p, paramsToCompare)) return 1;
 
  //Add to global labels
  int jetR = inPbPbFileConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  else{
    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  }
  globalLabels.push_back(jetRStr);
  globalLabels.push_back(jtPtBinsGlobalStr);
  globalLabels.push_back(gammaJtDPhiStr);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  configParser labels(inPbPbFileLabel_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();

  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    if(labelMap.count(globalLabels[gI]) != 0){
      globalLabels[gI] = labelMap[globalLabels[gI]];
    }
  }
  
  //Grab first set of parameters defined in necessaryInFileParams
  const bool isPP = inPPFileConfig_p->GetValue("ISPP", 0);
  const bool isPbPb = !(inPbPbFileConfig_p->GetValue("ISPP", 0));

  if(!isPP){
    std::cout << "Given pp file is not pp return 1" << std::endl;
    return 1;
  }
  if(!isPbPb){
    std::cout << "Given pbpb file is not pbpb return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> tempCentBinsStr = strToVect(inPbPbFileConfig_p->GetValue("CENTBINS", ""));
  std::vector<std::string> centBinsStr ;
  std::vector<std::string> centBinsLabel;

  for(unsigned int cI = 0; cI < tempCentBinsStr.size()-1; ++cI){
    centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
    centBinsLabel.push_back("#bf{Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%}");
  }  
  
  //Define a bunch of parameters for the TCanvas
  const Int_t nPad = 2;
  Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;
  Double_t height = 900.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);

  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = inPbPbFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inPbPbFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inPbPbFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inPbPbFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inPbPbFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<Float_t> xJMax, xJMin;
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    std::string gammaPtBinsStr = inPbPbFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBins != nGammaPtBinsTemp){
      std::cout << "GAMMAPTBINS DONT MATCH: " << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  gStyle->SetOptStat(0);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t iterPPPhoPt = inPPTermFileConfig_p->GetValue("GAMMAPT_PP", -1); 
 
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    Int_t iterPP = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_PP").c_str(), -1);
    if(iterPP < 0){
      std::cout << "Terminating iteration is \'" << iterPP << "\'. return 1" << std::endl;
      return 1;
    }

    TH1F* ppHistPhoPt_p = (TH1F*)inPPFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPPPhoPt) + "_PP_PURCORR_COMBINED_h").c_str());
    Double_t scaledTotal = 0.0;
    for(Int_t bIX = 0; bIX < ppHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = ppHistPhoPt_p->GetBinCenter(bIX+1);
      if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

      scaledTotal += ppHistPhoPt_p->GetBinContent(bIX+1);
    }
    

    TH1F* ppHist_p = (TH1F*)inPPFile_p->Get((varName + "Reco_Iter" + std::to_string(iterPP) + "_GammaPt" + std::to_string(gI) + "_PP_PURCORR_COMBINED_h").c_str());

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    binWidthAndScaleNorm(ppHist_p, scaledTotal);
    //    ppHist_p->Scale(1./ppHist_p->Integral());
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
    
    ppHist_p->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}*}#frac{dN_{#vec{JJ}#gamma}}{d#vec{x}_{JJ#gamma}}");

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;
  
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){   
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
 
      const Int_t iterPbPbPhoPt = inPbPbTermFileConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI]).c_str(), -1);


      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
   
      Int_t iterPbPb = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_" + centBinsStr[cI]).c_str(), -1);
      if(iterPbPb < 0){
	std::cout << "Terminating iteration is \'" << iterPbPb << "\'. return 1" << std::endl;
	return 1;
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
       
      std::string centBins2Str = centBinsStr[cI];
      if(centBins2Str.find("to") != std::string::npos){
	centBins2Str.replace(0, 4, "");
	centBins2Str.replace(centBins2Str.find("to"), 2, "-");
	centBins2Str = centBins2Str + "%";
      }
      
      TH1F* pbpbHist_p = (TH1F*)inPbPbFile_p->Get((varName + "Reco_Iter" + std::to_string(iterPbPb) + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      TH1F* pbpbHistPhoPt_p = nullptr;
      if(iterPbPbPhoPt > 0) pbpbHistPhoPt_p = (TH1F*)inPbPbFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPbPbPhoPt) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      else pbpbHistPhoPt_p = (TH1F*)inPbPbFile_p->Get(("photonPtReco_Iter3_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << ", " <<  std::endl;
       
      //      std::cout << "NAME: photonPtReco_Iter" << std::to_string(iterPbPbPhoPt) << "_" << centBinsStr[cI] << "_PURCORR_COMBINED_h" << std::endl;
      
      scaledTotal = 0.0;
      for(Int_t bIX = 0; bIX < pbpbHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
	Float_t binCenter = pbpbHistPhoPt_p->GetBinCenter(bIX+1);
	if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;
	
	scaledTotal += pbpbHistPhoPt_p->GetBinContent(bIX+1);
      }
      
      std::cout << scaledTotal << std::endl;
      binWidthAndScaleNorm(pbpbHist_p, scaledTotal);      
      //      pbpbHist_p->Scale(1./pbpbHist_p->Integral());

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
	
      TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
      TPad* pads_p[nPad];
      canv_p->SetTopMargin(noMargin);
      canv_p->SetRightMargin(noMargin);
      canv_p->SetLeftMargin(noMargin);
      canv_p->SetBottomMargin(noMargin);
      canv_p->cd();
      
      pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
      pads_p[0]->SetTopMargin(topMargin);
      pads_p[0]->SetRightMargin(rightMargin);
      pads_p[0]->SetLeftMargin(leftMargin);
      pads_p[0]->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      pads_p[0]->cd();
      canv_p->cd();
      
      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
      pads_p[1]->SetTopMargin(noMargin);
      pads_p[1]->SetRightMargin(rightMargin);
      pads_p[1]->SetLeftMargin(leftMargin);
      pads_p[1]->SetBottomMargin(bottomMargin);
      
      canv_p->cd();
      pads_p[1]->Draw("SAME");
      pads_p[1]->cd();
      canv_p->cd();
      
      pads_p[0]->cd();
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      HIJet::Style::EquipHistogram(pbpbHist_p, 1); 

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      
      ppHist_p->SetMinimum(minVal);
      ppHist_p->SetMaximum(maxVal);
      
      ppHist_p->DrawCopy("HIST E1 P");
      pbpbHist_p->DrawCopy("HIST E1 P SAME");

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      
      std::vector<std::string> labelsAlt;
      std::vector<std::string> labelsTemp = getLabels(inPbPbFileConfig_p, pbpbHist_p, &labelMap, &labelsAlt);
      std::vector<std::string> labels = {"#bf{#it{ATLAS Internal}}", "Pb+Pb, #it{pp} at #sqrt{s_{NN}}=5.02 TeV"};

      for(unsigned int lI = 0; lI < labelsTemp.size(); ++lI){
	labels.push_back(labelsTemp[lI]);
      }
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      
      unsigned int pos = 0;
      while(pos < labels.size()){
	if(labels[pos].find("PURCORR") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("COMBINED") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("Iter") != std::string::npos) labels.erase(pos + labels.begin());
	else ++pos;
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      
      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(0.035/(1.0 - padSplit));
      label_p->SetTextAlign(31);
      
      for(unsigned int i = 0; i < labels.size(); ++i){
	label_p->DrawLatex(0.95, 0.81 - 0.075*i, labels[i].c_str());
      }
      delete label_p;
      
      TLegend* leg_p = new TLegend(0.25, 0.8, 0.5, 0.95);
      leg_p->SetTextFont(42);
      leg_p->SetTextSize(0.035/(1.0 - padSplit));
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
   
      leg_p->AddEntry(pbpbHist_p, ("Pb+Pb " + centBins2Str).c_str(), "P L");
      leg_p->AddEntry(ppHist_p, "p+p", "P L");


      leg_p->Draw("SAME");
      
      gPad->SetTicks();
      if(doLogX) pads_p[0]->SetLogx();
      if(doLogY) pads_p[0]->SetLogy();
      
      canv_p->cd();
      pads_p[1]->cd();
      
      pbpbHist_p->Divide(ppHist_p);
      
      pbpbHist_p->GetYaxis()->SetNdivisions(505);
      
      pbpbHist_p->GetXaxis()->SetTitleOffset(1.2);
      
      pbpbHist_p->GetXaxis()->SetTitleFont(titleFont);
      pbpbHist_p->GetYaxis()->SetTitleFont(titleFont);
      pbpbHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit));
      pbpbHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));
      
      pbpbHist_p->GetXaxis()->SetLabelFont(titleFont);
      pbpbHist_p->GetYaxis()->SetLabelFont(titleFont);
      pbpbHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit));
      pbpbHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));     
      
      pbpbHist_p->GetYaxis()->SetTitle("Pb+Pb/p+p");
      pbpbHist_p->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0 - padSplit));

      pbpbHist_p->SetMinimum(minRatVal);
      pbpbHist_p->SetMaximum(maxRatVal);

      gPad->SetTicks();
      
      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength()*(1.0 - padSplit)/padSplit);
      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength("Y")*(1.0 - padSplit)/padSplit, "Y");
      pbpbHist_p->DrawCopy("HIST E1 P");

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      
      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(pbpbHist_p->GetBinLowEdge(1), 1.0, pbpbHist_p->GetBinLowEdge(pbpbHist_p->GetXaxis()->GetNbins()+1), 1.0);
      delete line_p;

      canv_p->cd();
      
      TBox* box_p = new TBox();
      box_p->SetFillColor(0);
      box_p->DrawBox(0.11, 0.435, 0.1575, 0.46);
      delete box_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
         
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + dateStr + "." + saveExt);
      delete canv_p;
      delete leg_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
         
    }
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;
   }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  inPbPbFile_p->Close();
  delete inPbPbFile_p;

  inPPFile_p->Close();
  delete inPPFile_p;
  
  
  std::cout << "GDJPLOTRESULTS COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotResults.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotResults(argv[1]);
  return retVal;
}
