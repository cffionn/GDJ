//Author: Chris McGinn (2020.08.04)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/globalDebugHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjPlotMBHist(std::string inConfigFileName)
{
  const Int_t nStyles = 15;
  const Int_t styles[nStyles] = {20,21,33,34,29,24,25,27,28,30,23,20,21,33,34};

  const int nColors = 15;
  const int colors[nColors]={1,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kCyan+3,28,41,kGray};

  const int nSizes = 15;
  const double sizes[nSizes] = {1,1,1.6,1.2,1.6,1,1,1,1.6,1,1,1,1,1.6,1.2};

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"MIXEDJTPTYMIN",
					"MIXEDJTPTYMAX",
					"MIXEDJTPTLABELX",
					"MIXEDJTPTLABELY",
					"MIXEDJTPTLEGX",
					"MIXEDJTPTLEGY",
  					"PERCENTYMIN",
					"PERCENTYMAX",
					"PERCENTLABELX",
					"PERCENTLABELY",
					"PERCENTLEGX",
					"PERCENTLEGY"};
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const Double_t mixedJtPtYMin = inConfig_p->GetValue("MIXEDJTPTYMIN", -1.0);
  const Double_t mixedJtPtYMax = inConfig_p->GetValue("MIXEDJTPTYMAX", 1.0);

  const Double_t mixedJtPtLabelX = inConfig_p->GetValue("MIXEDJTPTLABELX", -1.0); 
  const Double_t mixedJtPtLabelY = inConfig_p->GetValue("MIXEDJTPTLABELY", -1.0);

  const Double_t mixedJtPtLegX = inConfig_p->GetValue("MIXEDJTPTLEGX", -1.0); 
  const Double_t mixedJtPtLegY = inConfig_p->GetValue("MIXEDJTPTLEGY", -1.0); 

  const Double_t percentYMin = inConfig_p->GetValue("PERCENTYMIN", -1.0);
  const Double_t percentYMax = inConfig_p->GetValue("PERCENTYMAX", 1.0);

  const Double_t percentLabelX = inConfig_p->GetValue("PERCENTLABELX", -1.0); 
  const Double_t percentLabelY = inConfig_p->GetValue("PERCENTLABELY", -1.0);

  const Double_t percentLegX = inConfig_p->GetValue("PERCENTLEGX", -1.0); 
  const Double_t percentLegY = inConfig_p->GetValue("PERCENTLEGY", -1.0); 
  
  std::vector<std::string> globalLabels;
  for(Int_t gI = 0; gI < 20; ++gI){
    std::string globalLabel = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(globalLabel.size() == 0) continue;

    globalLabels.push_back(globalLabel);
  }
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");  
  TEnv* prevConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> prevReqParams = {"NCENTBINS",
					    "CENTBINSLOW",
					    "CENTBINSHIGH",
					    "JETETABINSLOW",
					    "JETETABINSHIGH",
					    "JETETABINSDOABS",	
					    "NJETPTBINS",
					    "JETPTBINSLOW",
					    "JETPTBINSHIGH",
					    "JETR"};
  if(!checkEnvForParams(prevConfig_p, prevReqParams)) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TH1F* cent_p = (TH1F*)inFile_p->Get("cent_h");
  TH1F* centData_p = (TH1F*)inFile_p->Get("centData_h");
  
  const Int_t nMaxCentBins = 500;
  const Int_t nCentBins = prevConfig_p->GetValue("NCENTBINS", 0);
  if(nCentBins > nMaxCentBins) return 1;

  const Float_t centBinsLow = prevConfig_p->GetValue("CENTBINSLOW", 0);
  const Float_t centBinsHigh = prevConfig_p->GetValue("CENTBINSHIGH", 0);

  Double_t centBins[nMaxCentBins+1];
  getLinBins(centBinsLow, centBinsHigh, nCentBins, centBins);

  std::vector<std::string> centBinsStr, centBinsStr2;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centLow = std::to_string((Int_t)centBins[cI]);
    std::string centHigh = std::to_string((Int_t)centBins[cI+1]);

    centBinsStr.push_back("Cent" + centLow + "to" + centHigh);

    if(centBins[cI] < 10.0) centLow = "0" + centLow;
    if(centBins[cI+1] < 10.0) centHigh = "0" + centHigh;

    centBinsStr2.push_back("Cent" + centLow + "to" + centHigh);
  }

  const int jetR = prevConfig_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);

  const Int_t nJetPtBins = prevConfig_p->GetValue("NJETPTBINS", 0);
  if(nJetPtBins > nMaxCentBins) return 1;

  const Float_t jetPtBinsLow = prevConfig_p->GetValue("JETPTBINSLOW", 0);
  const Float_t jetPtBinsHigh = prevConfig_p->GetValue("JETPTBINSHIGH", 0);

  Double_t jetPtBins[nMaxCentBins+1];
  getLinBins(jetPtBinsLow, jetPtBinsHigh, nJetPtBins, jetPtBins);

  const Float_t jetEtaBinsLow = prevConfig_p->GetValue("JETETABINSLOW", 0.0);
  const Float_t jetEtaBinsHigh = prevConfig_p->GetValue("JETETABINSHIGH", 0.0);
  const Bool_t jetEtaBinsDoAbs = prevConfig_p->GetValue("JETETABINSDOABS", 0);
  std::string jetEtaLabelStr = "#eta_{jet}";

  if(jetEtaBinsDoAbs) jetEtaLabelStr = "|" + jetEtaLabelStr + "|";
  if(jetEtaBinsDoAbs && jetEtaBinsLow < 0.000001) jetEtaLabelStr = jetEtaLabelStr + " < " + prettyString(jetEtaBinsHigh, 1, false);  
  else jetEtaLabelStr = prettyString(jetEtaBinsLow, 1, false) + " < " + jetEtaLabelStr + " < " + prettyString(jetEtaBinsHigh, 1, false);

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  const Double_t titleSizeX = 0.04;
  const Double_t titleSizeY = 0.04;
  const Double_t labelSizeX = 0.035;
  const Double_t labelSizeY = 0.035;
  const Int_t titleFont = 42;
  
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX);
  label_p->SetNDC();
  
  std::vector<double> phiFractions = {2*TMath::Pi(), TMath::Pi(), TMath::Pi()/2., TMath::Pi()/3., TMath::Pi()/4.};
  std::vector<std::string> phiLabels = {"No #Delta#phi cut", "#Delta#phi_{#gamma,\\mathrm{Jet}}$>#pi/2", "#Delta#phi_{#gamma,\\mathrm{Jet}}$>3#pi/4", "#Delta#phi_{#gamma,\\mathrm{Jet}}$>5#pi/6", "#Delta#phi_{#gamma,\\mathrm{Jet}}$>7#pi/8"};
  std::vector<std::string> phiLabels2 = {"No $\\Delta\\phi$ cut", "$\\Delta\\phi_{\\gamma,\\mathrm{Jet}}$>$\\pi$/2", "$\\Delta\\phi_{\\gamma,\\mathrm{Jet}}$>3$\\pi$/4", "$\\Delta\\phi_{\\gamma,\\mathrm{Jet}}$>5$\\pi$/6", "$\\Delta\\phi_{\\gamma,\\mathrm{Jet}}$>7$\\pi$/8"};
  std::vector<double> phiCenters = {0.0, TMath::Pi()/2.0, 3.0*TMath::Pi()/4.0, 5.0*TMath::Pi()/6.0, 7.0*TMath::Pi()/8.0};

  std::vector<double> percentiles = {0.01, 0.05, 0.10};
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.02);
  canv_p->SetRightMargin(0.02);
  canv_p->SetLeftMargin(0.15);
  canv_p->SetBottomMargin(0.15);
  
  cent_p->GetXaxis()->SetNdivisions(505);
  cent_p->GetYaxis()->SetNdivisions(505);    
    
  cent_p->GetXaxis()->SetTitleFont(titleFont);
  cent_p->GetYaxis()->SetTitleFont(titleFont);
  cent_p->GetXaxis()->SetTitleSize(titleSizeX);
  cent_p->GetYaxis()->SetTitleSize(titleSizeY);
  
  cent_p->GetXaxis()->SetLabelFont(titleFont);
  cent_p->GetYaxis()->SetLabelFont(titleFont);
  cent_p->GetXaxis()->SetLabelSize(labelSizeX);
  cent_p->GetYaxis()->SetLabelSize(labelSizeY);

  cent_p->SetMarkerSize(sizes[0]);
  cent_p->SetMarkerStyle(styles[0]);
  cent_p->SetMarkerColor(colors[0]);
  cent_p->SetLineColor(colors[0]);      

  cent_p->Scale(1./cent_p->Integral());

  cent_p->SetMaximum(0.045);
  cent_p->SetMinimum(0.00);

  cent_p->DrawCopy("HIST E1 P");

  centData_p->SetMarkerSize(sizes[1]);
  centData_p->SetMarkerStyle(styles[1]);
  centData_p->SetMarkerColor(colors[1]);
  centData_p->SetLineColor(colors[1]);      
  centData_p->Scale(1./centData_p->Integral());
  centData_p->DrawCopy("HIST E1 P SAME");

  gStyle->SetOptStat(0);


  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  int lPos = 0;
  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    if(globalLabels[gI].find("Mixed Event") != std::string::npos) continue;
    label_p->DrawLatex(0.3, 0.93 - lPos*0.055, globalLabels[gI].c_str());
    ++lPos;
  }
    
  TLegend* leg_p = new TLegend(0.7, 0.85 - 0.063*2., 0.7 + 0.25, 0.85);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSizeX);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  leg_p->AddEntry(cent_p, "CC + PC", "P L");
  leg_p->AddEntry(centData_p, "#gamma Triggered", "P L");

  leg_p->Draw("SAME");
  
  std::string saveName = "pdfDir/" + dateStr + "/cent_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete leg_p;
  delete canv_p;    

  globalLabels.push_back("anti-k_{t} R=" + jetRStr);
  globalLabels.push_back(jetEtaLabelStr);  
  
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.02);
    canv_p->SetRightMargin(0.02);
    canv_p->SetLeftMargin(0.15);
    canv_p->SetBottomMargin(0.15);

    leg_p = new TLegend(mixedJtPtLegX, mixedJtPtLegY - 0.063*(double)phiLabels.size(), mixedJtPtLegX + 0.25, mixedJtPtLegY);
    leg_p->SetTextFont(titleFont);
    leg_p->SetTextSize(titleSizeX);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);


    std::string histName = "jtMultPerCent_" + centBinsStr[cI] + "_h";
    TH1F* inHist_p = (TH1F*)inFile_p->Get(histName.c_str());

    inHist_p->GetXaxis()->SetNdivisions(505);
    inHist_p->GetYaxis()->SetNdivisions(505);    

    inHist_p->SetMinimum(mixedJtPtYMin);
    inHist_p->SetMaximum(mixedJtPtYMax);

    inHist_p->GetXaxis()->SetTitleFont(titleFont);
    inHist_p->GetYaxis()->SetTitleFont(titleFont);
    inHist_p->GetXaxis()->SetTitleSize(titleSizeX);
    inHist_p->GetYaxis()->SetTitleSize(titleSizeY);

    inHist_p->GetXaxis()->SetLabelFont(titleFont);
    inHist_p->GetYaxis()->SetLabelFont(titleFont);
    inHist_p->GetXaxis()->SetLabelSize(labelSizeX);
    inHist_p->GetYaxis()->SetLabelSize(labelSizeY);

    std::vector<TH1F*> dummyHists_p;    
    for(unsigned int pI = 0; pI < phiFractions.size(); ++pI){
      inHist_p->SetMarkerSize(sizes[pI]);
      inHist_p->SetMarkerStyle(styles[pI]);
      inHist_p->SetMarkerColor(colors[pI]);
      inHist_p->SetLineColor(colors[pI]);    

      dummyHists_p.push_back(new TH1F(("dummy_" + std::to_string(pI)).c_str(), "", 1, 0, 1));

      dummyHists_p[pI]->SetMarkerSize(sizes[pI]);
      dummyHists_p[pI]->SetMarkerStyle(styles[pI]);
      dummyHists_p[pI]->SetMarkerColor(colors[pI]);
      dummyHists_p[pI]->SetLineColor(colors[pI]);    

      leg_p->AddEntry(dummyHists_p[pI], phiLabels[pI].c_str(), "P L");
      
      if(pI == 0) inHist_p->DrawCopy("HIST E1 P");
      else{
	inHist_p->Scale(phiFractions[pI]/phiFractions[pI-1]);
	inHist_p->DrawCopy("HIST E1 P SAME");
      }
      gStyle->SetOptStat(0);
    }
  
    inHist_p->Scale(phiFractions[0]/phiFractions[phiFractions.size()-1]);
    
    leg_p->Draw("SAME");
    
    std::vector<std::string> tempLabels = globalLabels;
    tempLabels.push_back("#bf{" + std::to_string((int)centBins[cI]) + "-" + std::to_string((int)centBins[cI+1]) + "%}");
    for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
      label_p->DrawLatex(mixedJtPtLabelX, mixedJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
    }
    
    histName = "jtMultPerCent_" + centBinsStr2[cI] + "_h";
    saveName = "pdfDir/" + dateStr + "/" + histName + "_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);

    for(unsigned int lI = 0; lI < dummyHists_p.size(); ++lI){
      delete dummyHists_p[lI];
    }

    delete leg_p;
    delete canv_p;    
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Lets look at the requested centrality combinations
  for(unsigned int cI = 0; cI < 20; ++cI){
    std::string centLowStr = inConfig_p->GetValue(("CENTCOMBOLOW." + std::to_string(cI)).c_str(), "");
    std::string centHighStr = inConfig_p->GetValue(("CENTCOMBOHIGH." + std::to_string(cI)).c_str(), "");

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    if(centLowStr.size() == 0 && centHighStr.size() == 0) continue;
    else if(centLowStr.size() == 0){
      std::cout << "CENTCOMBOHIGH." << cI << " IS GIVEN BUT CENTCOMBOLOW." + cI << " IS NOT. Skipping" << std::endl;
      continue;
    }
    else if(centHighStr.size() == 0){
      std::cout << "CENTCOMBOLOW." << cI << " IS GIVEN BUT CENTCOMBOHIGH." + cI << " IS NOT. Skipping" << std::endl;
      continue;
    }

    const Int_t centLow = std::stoi(centLowStr);
    const Int_t centHigh = std::stoi(centHighStr);

    TH1F* jetPt_CentCombo_p = new TH1F("jetPt_CentCombo_h", (";Reco. Jet p_{T} [GeV];N_{Jet,R=" + jetRStr + "}/N_{Event," + centLowStr + "-" + centHighStr + "%}").c_str(), nJetPtBins, jetPtBins);

    std::vector<TGraph*> graphs_p;

    
    for(unsigned int pI = 0; pI < percentiles.size(); ++pI){
      graphs_p.push_back(new TGraph(phiCenters.size()));
    }
    
    if(centLow < 10.0) centLowStr = "0" + centLowStr;
    if(centHigh < 10.0) centHighStr = "0" + centHighStr;

    std::string centBinsStrCombo = "Cent" + centLowStr + "to" + centHighStr;

    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.02);
    canv_p->SetRightMargin(0.02);
    canv_p->SetLeftMargin(0.15);
    canv_p->SetBottomMargin(0.15);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    leg_p = new TLegend(mixedJtPtLegX, mixedJtPtLegY - 0.063*(double)phiLabels.size(), mixedJtPtLegX + 0.25, mixedJtPtLegY);
    leg_p->SetTextFont(titleFont);
    leg_p->SetTextSize(titleSizeX);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    
    double centSum = 0.0;
    std::vector<double> centVals;
    std::vector<std::vector<double> > jtPtVals, jtPtErrs;

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(Int_t cI = centLow; cI < centHigh; ++cI){
      std::string histName = "jtMultPerCent_" + centBinsStr[cI] + "_h";
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << histName << ", " << inFile_p << std::endl;

      TH1F* inHist_p = (TH1F*)inFile_p->Get(histName.c_str());      
      jtPtVals.push_back({});
      jtPtErrs.push_back({});

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      centSum += centData_p->GetBinContent(cI+1);
      centVals.push_back(centData_p->GetBinContent(cI+1));
      
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t jI = 0; jI < inHist_p->GetXaxis()->GetNbins(); ++jI){
	jtPtVals[cI-centLow].push_back(inHist_p->GetBinContent(jI+1));
	jtPtErrs[cI-centLow].push_back(inHist_p->GetBinError(jI+1));
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t jI = 0; jI < nJetPtBins; ++jI){
      double binVal = 0.0;
      double errVal = 0.0;
      
      for(unsigned int cI = 0; cI < centVals.size(); ++cI){
	double weight = centVals[cI]/centSum;
	binVal += jtPtVals[cI][jI]*weight;
	errVal += jtPtErrs[cI][jI]*weight*jtPtErrs[cI][jI]*weight;
      }

      errVal = TMath::Sqrt(errVal);

      jetPt_CentCombo_p->SetBinContent(jI+1, binVal);
      jetPt_CentCombo_p->SetBinError(jI+1, errVal);
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    jetPt_CentCombo_p->GetXaxis()->SetNdivisions(505);
    jetPt_CentCombo_p->GetYaxis()->SetNdivisions(505);    

    jetPt_CentCombo_p->SetMinimum(mixedJtPtYMin);
    jetPt_CentCombo_p->SetMaximum(mixedJtPtYMax);

    jetPt_CentCombo_p->GetXaxis()->SetTitleFont(titleFont);
    jetPt_CentCombo_p->GetYaxis()->SetTitleFont(titleFont);
    jetPt_CentCombo_p->GetXaxis()->SetTitleSize(titleSizeX);
    jetPt_CentCombo_p->GetYaxis()->SetTitleSize(titleSizeY);

    jetPt_CentCombo_p->GetXaxis()->SetLabelFont(titleFont);
    jetPt_CentCombo_p->GetYaxis()->SetLabelFont(titleFont);
    jetPt_CentCombo_p->GetXaxis()->SetLabelSize(labelSizeX);
    jetPt_CentCombo_p->GetYaxis()->SetLabelSize(labelSizeY);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    std::vector<TH1F*> dummyHists_p;    
    for(unsigned int pI = 0; pI < phiFractions.size(); ++pI){
      jetPt_CentCombo_p->SetMarkerSize(sizes[pI]);
      jetPt_CentCombo_p->SetMarkerStyle(styles[pI]);
      jetPt_CentCombo_p->SetMarkerColor(colors[pI]);
      jetPt_CentCombo_p->SetLineColor(colors[pI]);    

          if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      dummyHists_p.push_back(new TH1F(("dummy_" + std::to_string(pI)).c_str(), "", 1, 0, 1));

          if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      dummyHists_p[pI]->SetMarkerSize(sizes[pI]);
      dummyHists_p[pI]->SetMarkerStyle(styles[pI]);
      dummyHists_p[pI]->SetMarkerColor(colors[pI]);
      dummyHists_p[pI]->SetLineColor(colors[pI]);    

          if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      leg_p->AddEntry(dummyHists_p[pI], phiLabels[pI].c_str(), "P L");

          if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(pI == 0) jetPt_CentCombo_p->DrawCopy("HIST E1 P");
      else{
	jetPt_CentCombo_p->Scale(phiFractions[pI]/phiFractions[pI-1]);
	jetPt_CentCombo_p->DrawCopy("HIST E1 P SAME");
      }
      unsigned int percentilesPos = 0;    
      for(Int_t bIX = jetPt_CentCombo_p->GetXaxis()->GetNbins()+1; bIX >= 0; --bIX){
	if(jetPt_CentCombo_p->GetBinContent(bIX+1) > percentiles.at(percentilesPos)){
	  graphs_p[percentilesPos]->SetPoint(pI, phiCenters[pI], jetPt_CentCombo_p->GetXaxis()->GetBinLowEdge(bIX+1));
	  ++percentilesPos;
	  if(percentilesPos >= percentiles.size()) break;
	}
      }

      while(percentilesPos < percentiles.size()){
	graphs_p[percentilesPos]->SetPoint(pI, phiCenters[pI], jetPt_CentCombo_p->GetXaxis()->GetBinLowEdge(1));	
	++percentilesPos;
      }
      
      gStyle->SetOptStat(0);
    }
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    leg_p->Draw("SAME");

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::vector<std::string> tempLabels = globalLabels;
    tempLabels.push_back("#bf{" + std::to_string(centLow) + "-" + std::to_string(centHigh) + "%}");
    for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
      label_p->DrawLatex(mixedJtPtLabelX, mixedJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
    }

    gStyle->SetOptStat(0);
    saveName = "pdfDir/" + dateStr + "/jtMultPerCent_" + centBinsStrCombo + "_R" + std::to_string(jetR) + "_" + dateStr +  ".pdf";
    quietSaveAs(canv_p, saveName);
    delete jetPt_CentCombo_p;
    for(unsigned int lI = 0; lI < dummyHists_p.size(); ++lI){
      delete dummyHists_p[lI];
    }
    delete leg_p;
    delete canv_p;  

    canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.02);
    canv_p->SetRightMargin(0.02);
    canv_p->SetLeftMargin(0.15);
    canv_p->SetBottomMargin(0.15);
    
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    leg_p = new TLegend(percentLegX, percentLegY - 0.063*(double)percentiles.size(), percentLegX + 0.25, percentLegY);
    leg_p->SetTextFont(titleFont);
    leg_p->SetTextSize(titleSizeX);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);

    TH1F* deltaPhi_p = new TH1F("deltaPhi_h", ";#Delta#phi{#gamma,Jet};Reco. p_{T} of Percentile", 10, -0.1, TMath::Pi()+0.1);

    deltaPhi_p->SetMaximum(percentYMax);
    deltaPhi_p->SetMinimum(percentYMin);

    deltaPhi_p->DrawCopy("HIST E1 P");

    for(unsigned int gI = 0; gI < graphs_p.size(); ++gI){
      graphs_p[gI]->SetLineColor(colors[gI]);
      graphs_p[gI]->SetMarkerColor(colors[gI]);
      graphs_p[gI]->SetMarkerSize(sizes[gI]);
      graphs_p[gI]->SetMarkerStyle(styles[gI]);

      graphs_p[gI]->Draw("P");

      leg_p->AddEntry(graphs_p[gI], (prettyString(percentiles[gI], 2, false) + " fake").c_str(), "P L");
    }

    for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
      label_p->DrawLatex(percentLabelX, percentLabelY - gI*0.055, tempLabels[gI].c_str());
    }
    
    leg_p->Draw("SAME");
    saveName = "pdfDir/" + dateStr + "/percentiles_" + centBinsStrCombo + "_R" + std::to_string(jetR) + "_" + dateStr +  ".pdf";
    quietSaveAs(canv_p, saveName);

    //OUTPUT TABULAR FOR CENTRALMOST
    if(centLow == 0){
      std::string lStr = "l | ";
      for(unsigned int lI = 0; lI < phiCenters.size()-1; ++lI){
	lStr = lStr + "l ";
      }
      lStr = lStr + "l";
	  
      std::cout << "\\begin{table}" << std::endl;
      std::cout << "\\fontsize{7}{7}\\selectfont" << std::endl;
      std::cout << "\\begin{center}" << std::endl;
      std::cout << "\\hspace*{-0.1cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
      std::string multiColumn = std::to_string(phiCenters.size()+1);

      std::cout << "\\multicolumn{" << multiColumn << "}{c}{\\pt [GeV] where given $\\Delta\\phi$ cut acheives X Fake Fraction (R=" << jetRStr << ")} \\\\" << std::endl;
      std::string columnStr = " X Fake Fraction &";
      for(unsigned int pI = 0; pI < phiCenters.size()-1; ++pI){
	columnStr = columnStr + " " + phiLabels2[pI] + " &";
      }
      columnStr = columnStr + " " + phiLabels2[phiLabels2.size()-1] + " \\\\ \\hline" ;
      
      std::cout << columnStr << std::endl;

      
      for(unsigned int cI = 0; cI < percentiles.size(); ++cI){
	std::string outStr = prettyString(percentiles[cI], 2, false);

	for(unsigned int pI = 0; pI < phiCenters.size(); ++pI){
	  outStr = outStr + " & " + graphs_p[cI]->GetPointY(pI);
	}	
	if(cI < percentiles.size()-1) std::cout << outStr << " \\\\" << std::endl;
	else std::cout << outStr << std::endl;
      }
      std::cout << "\\end{tabular}" << std::endl;
      std::cout << "\\end{center}" << std::endl;
      std::cout << "\\end{table}" << std::endl;
    }
    
    for(unsigned int gI = 0; gI < graphs_p.size(); ++gI){
      delete graphs_p[gI];
    }
    
    delete deltaPhi_p;
    delete leg_p;
    delete canv_p;  
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  delete label_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  delete inConfig_p;
  
  std::cout << "GDJPLOTMBHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotMBHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotMBHist(argv[1]);
  return retVal;
}
