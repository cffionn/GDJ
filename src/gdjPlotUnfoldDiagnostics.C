//Author: Chris McGinn (2021.01.13)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
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
#include "include/treeUtil.h"

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

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INUNFOLDFILENAME",
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
  for(unsigned int i = 0; i < 100; ++i){
    std::string tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
    if(tempStr.size() != 0) globalLabels.push_back(tempStr);  
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
							"ISPP",
							"NITER",
							"JETR",
							"ASSOCGENMINPT"};

  const std::string gammaJtDPhiStr = "DPhi0";
  //  const std::string multiJtCutGlobalStr = "MultiJt0";
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> pbpbParams = {"CENTBINS"};
  if(!checkEnvForParams(inUnfoldFileConfig_p, necessaryUnfoldFileParams)) return 1;

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
  
  configParser labels(inUnfoldFileLabel_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();

  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    if(labelMap.count(globalLabels[gI]) != 0){
      globalLabels[gI] = labelMap[globalLabels[gI]];
    }
  }
  
  //Grab first set of parameters defined in necessaryInFileParams
  //  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  //  const bool isPP = inUnfoldFileConfig_p->GetValue("ISPP", 0);

  const Int_t nIterMax = 100;
  const Int_t nIter = inUnfoldFileConfig_p->GetValue("NITER", 0);

  if(nIter > nIterMax){
    std::cout << "nIter \'" << nIter << "\' exceeds maximum allowed nIterMax \'" << nIterMax << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nIterForLoop = TMath::Min(nIter, nIterCap);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Define a bunch of parameters for the TCanvas
  Double_t padSplit = 0.375;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit;
  const Double_t topMargin = 0.01;
  const Double_t rightMargin = 0.01;
  Double_t height = 450.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);

  const Double_t leftMargin2D = 0.12;
  const Double_t bottomMargin2D = leftMargin2D;
  const Double_t topMargin2D = 0.01;
  const Double_t rightMargin2D = 0.12;
  Double_t height2D = 450.0;
  Double_t width2D = height2D*(1.0 - topMargin2D*(1.0 - padSplit) - bottomMargin2D*padSplit)/(1.0 - leftMargin2D - rightMargin2D);
  
  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = inUnfoldFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inUnfoldFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  std::vector<Float_t> xJMax, xJMin;
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  gStyle->SetOptStat(0);

  TH2F* recoMatrix_p = (TH2F*)inUnfoldFile_p->Get("gammaPtJetXJReco_HalfResponse_h");
  TH2F* truthMatrix_p = (TH2F*)inUnfoldFile_p->Get("gammaPtJetXJTruth_HalfResponse_h");

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
    saveName = "pdfDir/" + saveName + ".pdf";
    quietSaveAs(canv_p, saveName);    
    delete canv_p;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;    
  
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){    
    TH1F* truth_p = (TH1F*)inUnfoldFile_p->Get(("xjTruth_HalfToUnfold_GammaPt" + std::to_string(gI) + "_h").c_str());
    TH1F* reco_p = (TH1F*)inUnfoldFile_p->Get(("xjReco_HalfToUnfold_GammaPt" + std::to_string(gI) + "_h").c_str());
    TH1F* unfold_p[nIterMax];

    TH1F* statsDelta_p = new TH1F(("statsDelta_GammaPt" + std::to_string(gI) + "_h").c_str(), ";Number of Iterations;#sqrt{#Sigma #delta^{2}}", nIter-1, 1.5, nIter+0.5);
    TH1F* iterDelta_p = new TH1F(("iterDelta_GammaPt" + std::to_string(gI) + "_h").c_str(), ";Number of Iterations;#sqrt{#Sigma #delta^{2}}", nIter-1, 1.5, ((double)nIter)+0.5);
    TH1F* totalDelta_p = new TH1F(("totalelta_GammaPt" + std::to_string(gI) + "_h").c_str(), ";Number of Iterations;#sqrt{#Sigma #delta^{2}}", nIter-1, 1.5, ((double)nIter)+0.5);
    
    for(Int_t i = 1; i < nIter; ++i){
      unfold_p[i] = (TH1F*)inUnfoldFile_p->Get(("xjReco_HalfToUnfold_Iter" + std::to_string(i) + "_GammaPt" + std::to_string(gI) + "_h").c_str());

      if(i > 1){
	Float_t deltaIter = 0.0;
	Float_t deltaStat = 0.0;
	for(Int_t bIX = 0; bIX < unfold_p[i]->GetXaxis()->GetNbins(); ++bIX){
	  if(unfold_p[i]->GetBinContent(bIX+1) <= TMath::Power(10,-50)) continue;

	  Float_t tempDeltaIter = (unfold_p[i]->GetBinContent(bIX+1) -  unfold_p[i-1]->GetBinContent(bIX+1))/unfold_p[i]->GetBinContent(bIX+1);
	  deltaIter += tempDeltaIter*tempDeltaIter;

	  Float_t tempDeltaStat = unfold_p[i]->GetBinError(bIX+1)/unfold_p[i]->GetBinContent(bIX+1);
	  deltaStat += tempDeltaStat*tempDeltaStat;
	}
	Float_t deltaTot = deltaIter + deltaStat;
	
	statsDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaStat));
	statsDelta_p->SetBinError(i-1, 0.0);

	iterDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaIter));
	iterDelta_p->SetBinError(i-1, 0.0);

	totalDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaTot));
	totalDelta_p->SetBinError(i-1, 0.0);
      }
    }

    std::vector<std::string> tempLabels = globalLabels;
    tempLabels.push_back(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false));

    TLatex* label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(titleFont);
    label_p->SetTextSize(titleSize);
    label_p->SetTextColor(1);
    if(diagLabelAlignRightGlobal) label_p->SetTextAlign(31);

    Int_t maxLegX = -1;
    std::vector<TLegend*> leg_p;
    std::vector<TH1F*> histsForLeg_p = {iterDelta_p, statsDelta_p, totalDelta_p};
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
      std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
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

    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
      leg_p[lI]->Draw("SAME");
    }

    for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
      label_p->DrawLatex(diagLabelXGlobal, diagLabelYGlobal - tI*0.06, tempLabels[tI].c_str());
    }    
    
    quietSaveAs(canv_p, "pdfDir/diag_GammaPt" + std::to_string(gI) + ".png");
    delete canv_p;

    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
      delete leg_p[lI];
    }
    leg_p.clear();
    
    Float_t xjMinTemp = config_p->GetValue(("XJMIN." + std::to_string(gI)).c_str(), -1000);
    Float_t xjMaxTemp = config_p->GetValue(("XJMAX." + std::to_string(gI)).c_str(), -1000);
    if(xjMaxTemp < 0){
      xjMinTemp = xjMinGlobal;
      xjMaxTemp = xjMaxGlobal;	
    }
    
    truth_p->SetMinimum(xjMinTemp);
    truth_p->SetMaximum(xjMaxTemp);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    canv_p = new TCanvas("canv_p", "", width, height);
    TPad* pads_p[2];
    canv_p->SetTopMargin(0.001);
    canv_p->SetRightMargin(0.001);
    canv_p->SetLeftMargin(0.001);
    canv_p->SetBottomMargin(0.001);
    canv_p->cd();
    
    pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
    pads_p[0]->SetTopMargin(topMargin);
    pads_p[0]->SetRightMargin(rightMargin);
    pads_p[0]->SetLeftMargin(leftMargin);
    pads_p[0]->SetBottomMargin(0.001);
    
    canv_p->cd();
    pads_p[0]->Draw("SAME");
    pads_p[0]->cd();
    canv_p->cd();

    pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
    pads_p[1]->SetTopMargin(0.001);
    pads_p[1]->SetRightMargin(rightMargin);
    pads_p[1]->SetLeftMargin(leftMargin);
    pads_p[1]->SetBottomMargin(bottomMargin);
    
    canv_p->cd();
    pads_p[1]->Draw("SAME");
    pads_p[1]->cd();
    canv_p->cd();
    
    
    pads_p[0]->cd();
   
    HIJet::Style::EquipHistogram(truth_p, 2); //TRUTH IS kAzure-3 lines (atlas style picked up here)
    HIJet::Style::EquipHistogram(reco_p, 1); //RECO IS kRed-4 lines (atlas style picked up here)
    reco_p->SetMarkerSize(0.00001);
    truth_p->SetMarkerSize(0.00001);    

    maxLegX = -1;
    //    Int_t nLeg = nIterForLoop + 2 - 1;

    label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(titleFont);
    label_p->SetTextSize(titleSize/(1.0 - padSplit));
    label_p->SetTextColor(1);
    if(xjLabelAlignRightGlobal) label_p->SetTextAlign(31);

    stringsForLeg_p = {"Truth", "Reco."};
    histsForLeg_p = {truth_p, reco_p};    
    std::vector<std::string> legFillStr = {"L", "L"};
    
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
      std::cout << "Number of histograms per legend requested custom \'" << totalCustom << "\' totals to less than required, \'" << numberHists << "\'. Reverting to equal splits." << std::endl;
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
	leg_p[lI] = new TLegend(xjLegXGlobal[lI], xjLegYGlobal[lI] - xjLegNGlobal[lI]*0.05, xjLegXGlobal[lI] + maxLegX*0.02, xjLegYGlobal[lI]);	
      }
      
      leg_p[lI]->SetTextFont(titleFont);
      leg_p[lI]->SetTextSize(titleSize/(1.0 - padSplit));
      leg_p[lI]->SetBorderSize(0);
      leg_p[lI]->SetFillColor(0);
      leg_p[lI]->SetFillStyle(0);
    }

    truth_p->GetXaxis()->SetTitleFont(titleFont);
    truth_p->GetYaxis()->SetTitleFont(titleFont);
    truth_p->GetXaxis()->SetTitleSize(0.00001);
    truth_p->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));

    truth_p->GetXaxis()->SetLabelFont(titleFont);
    truth_p->GetYaxis()->SetLabelFont(titleFont);
    truth_p->GetXaxis()->SetLabelSize(0.00001);
    truth_p->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));

    truth_p->GetYaxis()->SetTitleOffset(yOffset);

    truth_p->GetYaxis()->SetNdivisions(505);
    truth_p->GetXaxis()->SetNdivisions(505);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    truth_p->DrawCopy("HIST E1");
    reco_p->DrawCopy("HIST E1 SAME");

    for(Int_t i = 1; i < nIterForLoop; ++i){
      HIJet::Style::EquipHistogram(unfold_p[i], i-1);
      unfold_p[i]->DrawCopy("HIST E1 P SAME");
    }

    if(xjLegNGlobal.size() == 0){
      for(unsigned int lI = 0; lI < histsForLeg_p.size(); ++lI){
	leg_p[lI/histsPerLeg]->AddEntry(histsForLeg_p[lI], stringsForLeg_p[lI].c_str(), legFillStr[lI].c_str());
      }
    }
    else{
      Int_t legCounter = 0;
      for(unsigned int lI = 0; lI < xjLegNGlobal.size(); ++lI){
	for(unsigned int dI = 0; dI < xjLegNGlobal[lI]; ++dI){	
	  leg_p[lI]->AddEntry(histsForLeg_p[legCounter], stringsForLeg_p[legCounter].c_str(), legFillStr[legCounter].c_str());
	  ++legCounter;
	}
      }
    }

    
    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
      leg_p[lI]->Draw("SAME");
    }

    for(unsigned int tI = 0; tI < tempLabels.size(); ++tI){
      label_p->DrawLatex(xjLabelXGlobal, xjLabelYGlobal - tI*0.083, tempLabels[tI].c_str());
    }    
    
    canv_p->cd();
    pads_p[1]->cd();

    Float_t xjRatMinTemp = config_p->GetValue(("XJRATMIN." + std::to_string(gI)).c_str(), -1000);
    Float_t xjRatMaxTemp = config_p->GetValue(("XJRATMAX." + std::to_string(gI)).c_str(), -1000);
    if(xjRatMaxTemp < 0){
      xjRatMinTemp = xjRatMinGlobal;
      xjRatMaxTemp = xjRatMaxGlobal;	
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t i = 1; i < nIterForLoop; ++i){
      unfold_p[i]->Divide(truth_p);

      if(i == 1){
	unfold_p[i]->GetYaxis()->SetTitle("Unfolded/Truth");
	unfold_p[i]->SetMinimum(xjRatMinTemp);
	unfold_p[i]->SetMaximum(xjRatMaxTemp);

	unfold_p[i]->GetXaxis()->SetTitleOffset(1.2);
	unfold_p[i]->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0-padSplit));
	
	unfold_p[i]->GetXaxis()->SetTitleFont(titleFont);
	unfold_p[i]->GetYaxis()->SetTitleFont(titleFont);
	unfold_p[i]->GetXaxis()->SetTitleSize(titleSize/(padSplit));
	unfold_p[i]->GetYaxis()->SetTitleSize(titleSize/(padSplit));
	
	unfold_p[i]->GetXaxis()->SetLabelFont(titleFont);
	unfold_p[i]->GetYaxis()->SetLabelFont(titleFont);
	unfold_p[i]->GetXaxis()->SetLabelSize(labelSize/(padSplit));
	unfold_p[i]->GetYaxis()->SetLabelSize(labelSize/(padSplit));
	
	unfold_p[i]->GetYaxis()->SetNdivisions(505);	

	unfold_p[i]->DrawCopy("HIST E1 P");
      }
      else unfold_p[i]->DrawCopy("HIST E1 P SAME");
    }

    TLine* line_p = new TLine();
    line_p->SetLineStyle(2);
    line_p->DrawLine(reco_p->GetXaxis()->GetBinLowEdge(1), 1.0, reco_p->GetXaxis()->GetBinLowEdge(reco_p->GetXaxis()->GetNbins()+2), 1.0);
    delete line_p;
    delete label_p;
    
    quietSaveAs(canv_p, "pdfDir/test_GammaPt" + std::to_string(gI) + ".png");
    
    delete pads_p[0];
    delete canv_p;

    for(unsigned int lI = 0; lI < leg_p.size(); ++lI){
      delete leg_p[lI];      
    }
    leg_p.clear();
  }    

  inUnfoldFile_p->Close();
  delete inUnfoldFile_p;

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
