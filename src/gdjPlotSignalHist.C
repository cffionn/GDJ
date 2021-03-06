//Author: Chris McGinn (2020.08.07)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjPlotSignalHist(std::string inConfigFileName)
{
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"SIGNALJTPTYMIN",
					"SIGNALJTPTYMAX",
					"SIGNALJTPTLABELX",
					"SIGNALJTPTLABELY",
					"SIGNALJTPTLEGX",
					"SIGNALJTPTLEGY",
					"SIGNALVECTJTPTYMIN",
					"SIGNALVECTJTPTYMAX",
					"SIGNALVECTJTPTLABELX",
					"SIGNALVECTJTPTLABELY",
					"SIGNALVECTJTPTLEGX",
					"SIGNALVECTJTPTLEGY"};

  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir/" + dateStr);

  const Double_t signalJtPtYMin = inConfig_p->GetValue("SIGNALJTPTYMIN", -1.0);
  const Double_t signalJtPtYMax = inConfig_p->GetValue("SIGNALJTPTYMAX", 1.0);
  const Double_t signalJtPtLabelX = inConfig_p->GetValue("SIGNALJTPTLABELX", -1.0);
  const Double_t signalJtPtLabelY = inConfig_p->GetValue("SIGNALJTPTLABELY", -1.0);
  const Double_t signalJtPtLegX = inConfig_p->GetValue("SIGNALJTPTLEGX", -1.0);
  const Double_t signalJtPtLegY = inConfig_p->GetValue("SIGNALJTPTLEGY", -1.0);

  const Double_t signalVectJtPtYMin = inConfig_p->GetValue("SIGNALVECTJTPTYMIN", -1.0);
  const Double_t signalVectJtPtYMax = inConfig_p->GetValue("SIGNALVECTJTPTYMAX", 1.0);
  const Double_t signalVectJtPtLabelX = inConfig_p->GetValue("SIGNALVECTJTPTLABELX", -1.0);
  const Double_t signalVectJtPtLabelY = inConfig_p->GetValue("SIGNALVECTJTPTLABELY", -1.0);
  const Double_t signalVectJtPtLegX = inConfig_p->GetValue("SIGNALVECTJTPTLEGX", -1.0);
  const Double_t signalVectJtPtLegY = inConfig_p->GetValue("SIGNALVECTJTPTLEGY", -1.0);
  
  std::vector<std::string> globalLabels;
  for(Int_t gI = 0; gI < 20; ++gI){
    std::string globalLabel = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(globalLabel.size() == 0) continue;

    globalLabels.push_back(globalLabel);
  }

  
  const Double_t titleSizeX = 0.04;
  const Double_t titleSizeY = 0.04;
  const Double_t labelSizeX = 0.035;
  const Double_t labelSizeY = 0.035;
  const Int_t titleFont = 42;

  
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX);
  label_p->SetNDC();
 
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> reqFileParams = {"NGAMMAPTBINS",
					    "GAMMAPTBINSLOW",
					    "GAMMAPTBINSHIGH",
					    "NDPHICUTS",
					    "NTWOJETPTBINS",
					    "TWOJETPTBINSLOW",
					    "TWOJETPTBINSHIGH",
					    "JETR"};

  if(!checkEnvForParams(fileConfig_p, reqFileParams)) return 1;

  const int jetR = fileConfig_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\\
'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);

  const Int_t nMaxBins = 500;
  const Int_t nGammaPtBins = fileConfig_p->GetValue("NGAMMAPTBINS", -1);
  const Float_t gammaPtBinsLow = fileConfig_p->GetValue("GAMMAPTBINSLOW", -1);
  const Float_t gammaPtBinsHigh = fileConfig_p->GetValue("GAMMAPTBINSHIGH", -1);
  Double_t gammaPtBins[nMaxBins+1];
  getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  const Int_t nTwoJetPtBins = fileConfig_p->GetValue("NTWOJETPTBINS", -1);
  const Float_t twoJetPtBinsLow = fileConfig_p->GetValue("TWOJETPTBINSLOW", -1);
  const Float_t twoJetPtBinsHigh = fileConfig_p->GetValue("TWOJETPTBINSHIGH", -1);
  Double_t twoJetPtBins[nMaxBins+1];
  getLinBins(twoJetPtBinsLow, twoJetPtBinsHigh, nTwoJetPtBins, twoJetPtBins);
  
  const Int_t nDPhiCuts = fileConfig_p->GetValue("NDPHICUTS", 0);

  std::vector<std::string> dphiLabels;

  for(unsigned int dI = 0; dI < 100; ++dI){
    std::string labelStr = fileConfig_p->GetValue(("DPHICUTLABEL." + std::to_string(dI)).c_str(), "");
    if(labelStr.size() == 0) continue;

    dphiLabels.push_back(labelStr);
  }

  if(nDPhiCuts != (int)dphiLabels.size()){
    std::cout << "Mismatch between provided number of nDPhiCuts \'" << nDPhiCuts << "\' and number of provided labels \'" << dphiLabels.size() << "\'. return 1" << std::endl;
    return 1;
  }

  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.02);
    canv_p->SetRightMargin(0.02);
    canv_p->SetBottomMargin(0.15);
    canv_p->SetLeftMargin(0.15);

    TLegend* leg_p = new TLegend(signalJtPtLegX, signalJtPtLegY - 0.063*(double)dphiLabels.size(), signalJtPtLegX + 0.25, signalJtPtLegY);
    leg_p->SetTextFont(titleFont);
    leg_p->SetTextSize(titleSizeX);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);

    
    std::vector<TH1F*> hists_p;
    hists_p.reserve(nDPhiCuts);
    for(Int_t dI = 0; dI < nDPhiCuts; ++dI){
      hists_p.push_back((TH1F*)inFile_p->Get(("jetPtPerGammaPtDPhi_GammaPt" + std::to_string(gI) + "_DPhi" + std::to_string(dI) + "_h").c_str()));

      HIJet::Style::EquipHistogram(hists_p[dI], dI);

      leg_p->AddEntry(hists_p[dI], dphiLabels[dI].c_str(), "P L");

      hists_p[dI]->SetMinimum(signalJtPtYMin);
      hists_p[dI]->SetMaximum(signalJtPtYMax);

      hists_p[dI]->GetXaxis()->SetTitleSize(titleSizeX);
      hists_p[dI]->GetYaxis()->SetTitleSize(titleSizeY);
      hists_p[dI]->GetXaxis()->SetLabelSize(labelSizeX);
      hists_p[dI]->GetYaxis()->SetLabelSize(labelSizeY);

      hists_p[dI]->GetXaxis()->SetNdivisions(505);
      hists_p[dI]->GetYaxis()->SetNdivisions(505);
      
      if(dI == 0) hists_p[dI]->DrawCopy("HIST E1 P");
      else hists_p[dI]->DrawCopy("HIST E1 P SAME");
    }

    std::vector<std::string> tempLabels = globalLabels;
    tempLabels.push_back("p_{T}^{#gamma} > " + prettyString(gammaPtBins[gI], 1, false));
    for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
      label_p->DrawLatex(signalJtPtLabelX, signalJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
    }

    
    
    gStyle->SetOptStat(0);
    leg_p->Draw("SAME");

    std::string gIStr = std::to_string(gI);
    if(gI < 10) gIStr = "0" + gIStr;
    
    std::string saveName = "pdfDir/" + dateStr + "/jetPtPerGammaPtDPhi_GammaPt" + gIStr + "_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);

    delete leg_p;
    delete canv_p;
  }


  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    for(Int_t jI = 0; jI < nTwoJetPtBins; ++jI){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.02);
      canv_p->SetRightMargin(0.02);
      canv_p->SetBottomMargin(0.15);
      canv_p->SetLeftMargin(0.15);
      
      TLegend* leg_p = new TLegend(signalJtPtLegX, signalJtPtLegY - 0.063*(double)dphiLabels.size(), signalJtPtLegX + 0.25, signalJtPtLegY);
      leg_p->SetTextFont(titleFont);
      leg_p->SetTextSize(titleSizeX);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetFillStyle(0);
      
      std::vector<TH1F*> hists_p;
      hists_p.reserve(nDPhiCuts);
      for(Int_t dI = 0; dI < nDPhiCuts; ++dI){
	hists_p.push_back((TH1F*)inFile_p->Get(("jetPtPerGammaPtDPhi_GammaPt" + std::to_string(gI) + "_DPhi" + std::to_string(dI) + "_JetPt" + std::to_string(jI) + "_TwoJets_h").c_str()));

	HIJet::Style::EquipHistogram(hists_p[dI], dI);
	
	leg_p->AddEntry(hists_p[dI], dphiLabels[dI].c_str(), "P L");
	
	hists_p[dI]->SetMinimum(signalJtPtYMin);
	hists_p[dI]->SetMaximum(signalJtPtYMax);

	hists_p[dI]->GetXaxis()->SetTitleSize(titleSizeX);
	hists_p[dI]->GetYaxis()->SetTitleSize(titleSizeY);
	hists_p[dI]->GetXaxis()->SetLabelSize(labelSizeX);
	hists_p[dI]->GetYaxis()->SetLabelSize(labelSizeY);
	
	hists_p[dI]->GetXaxis()->SetNdivisions(505);
	hists_p[dI]->GetYaxis()->SetNdivisions(505);
	
	if(dI == 0) hists_p[dI]->DrawCopy("HIST E1 P");
	else hists_p[dI]->DrawCopy("HIST E1 P SAME");
      }

      std::vector<std::string> tempLabels = globalLabels;
      tempLabels.push_back("p_{T}^{#gamma} > " + prettyString(gammaPtBins[gI], 1, false));
      tempLabels.push_back("2 jets p_{T} > " + prettyString(twoJetPtBins[jI], 1, false));
      for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
	label_p->DrawLatex(signalJtPtLabelX, signalJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
      }
               
      gStyle->SetOptStat(0);
      leg_p->Draw("SAME");
      
      std::string gIStr = std::to_string(gI);
      if(gI < 10) gIStr = "0" + gIStr;

      std::string jIStr = std::to_string(jI);
      if(jI < 10) jIStr = "0" + jIStr;
      
      std::string saveName = "pdfDir/" + dateStr + "/jetPtPerGammaPtDPhi_GammaPt" + gIStr + "_JetPt" + jIStr + "_TwoJets_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
      quietSaveAs(canv_p, saveName);
      
      delete leg_p;
      delete canv_p;
    }
  }

  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    for(Int_t jI = 0; jI < nTwoJetPtBins; ++jI){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.02);
      canv_p->SetRightMargin(0.02);
      canv_p->SetBottomMargin(0.15);
      canv_p->SetLeftMargin(0.15);
      
      TLegend* leg_p = new TLegend(signalVectJtPtLegX, signalVectJtPtLegY - 0.063*(double)dphiLabels.size(), signalVectJtPtLegX + 0.25, signalVectJtPtLegY);
      leg_p->SetTextFont(titleFont);
      leg_p->SetTextSize(titleSizeX);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetFillStyle(0);
      
      std::vector<TH1F*> hists_p;
      hists_p.reserve(nDPhiCuts);
      for(Int_t dI = 0; dI < nDPhiCuts; ++dI){
	hists_p.push_back((TH1F*)inFile_p->Get(("jetVectPtPerGammaPtDPhi_GammaPt" + std::to_string(gI) + "_DPhi" + std::to_string(dI) + "_JetPt" + std::to_string(jI) + "_TwoJets_h").c_str()));

	HIJet::Style::EquipHistogram(hists_p[dI], dI);

	leg_p->AddEntry(hists_p[dI], dphiLabels[dI].c_str(), "P L");
	
	hists_p[dI]->SetMinimum(signalVectJtPtYMin);
	hists_p[dI]->SetMaximum(signalVectJtPtYMax);

	hists_p[dI]->GetXaxis()->SetTitleSize(titleSizeX);
	hists_p[dI]->GetYaxis()->SetTitleSize(titleSizeY);
	hists_p[dI]->GetXaxis()->SetLabelSize(labelSizeX);
	hists_p[dI]->GetYaxis()->SetLabelSize(labelSizeY);
	
	hists_p[dI]->GetXaxis()->SetNdivisions(505);
	hists_p[dI]->GetYaxis()->SetNdivisions(505);
	
	if(dI == 0) hists_p[dI]->DrawCopy("HIST E1 P");
	else hists_p[dI]->DrawCopy("HIST E1 P SAME");
      }

      std::vector<std::string> tempLabels = globalLabels;
      tempLabels.push_back("p_{T}^{#gamma} > " + prettyString(gammaPtBins[gI], 1, false));
      tempLabels.push_back("2 jets p_{T} > " + prettyString(twoJetPtBins[jI], 1, false));
      for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
	label_p->DrawLatex(signalVectJtPtLabelX, signalVectJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
      }
               
      gStyle->SetOptStat(0);
      leg_p->Draw("SAME");
      
      std::string gIStr = std::to_string(gI);
      if(gI < 10) gIStr = "0" + gIStr;

      std::string jIStr = std::to_string(jI);
      if(jI < 10) jIStr = "0" + jIStr;
      
      std::string saveName = "pdfDir/" + dateStr + "/jetVectPtPerGammaPtDPhi_GammaPt" + gIStr + "_JetPt" + jIStr + "_TwoJets_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
      quietSaveAs(canv_p, saveName);
      
      delete leg_p;
      delete canv_p;
    }
  }

  inFile_p->Close();
  delete inFile_p;
  delete inConfig_p;
  
  std::cout << "GDJPLOTSIGNALHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotSignalHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotSignalHist(argv[1]);
  return retVal;
}
