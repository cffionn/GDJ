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
#include "TH1F.h"
#include "TLatex.h"
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
					"MIXEDJTPTLABELY"};
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const Double_t mixedJtPtYMin = inConfig_p->GetValue("MIXEDJTPTYMIN", -1.0);
  const Double_t mixedJtPtYMax = inConfig_p->GetValue("MIXEDJTPTYMAX", 1.0);

  const Double_t mixedJtPtLabelX = inConfig_p->GetValue("MIXEDJTPTLABELX", -1.0); 
  const Double_t mixedJtPtLabelY = inConfig_p->GetValue("MIXEDJTPTLABELY", -1.0);
 
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
					    "JETeTABINSDOABS",
					    "JETR"};
  if(!checkEnvForParams(prevConfig_p, prevReqParams)) return 1;
  
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
  globalLabels.push_back("anti-k_{t} R=" + jetRStr);

  const Float_t jetEtaBinsLow = prevConfig_p->GetValue("JETETABINSLOW", 0);
  const Float_t jetEtaBinsHigh = prevConfig_p->GetValue("JETETABINSHIGH", 0);
  const Bool_t jetEtaBinsDoAbs = prevConfig_p->GetValue("JETETABINSDOABS", 0);
  std::string jetEtaLabelStr = "#eta_{jet}";
  if(jetEtaBinsDoAbs) jetEtaLabelStr = "|" + jetEtaLabelStr + "|";
  if(jetEtaBinsDoAbs && jetEtaBinsLow < 0.000001) jetEtaLabelStr = jetEtaLabelStr + "<" + prettyString(jetEtaBinsHigh, 1, false);  
  else jetEtaLabelStr = prettyString(jetEtaBinsLow, 1, false) + "<" + jetEtaLabelStr + "<" + prettyString(jetEtaBinsHigh, 1, false);

  globalLabels.push_back(jetEtaLabelStr);
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  Double_t titleSizeX = 0.04;
  Double_t titleSizeY = 0.04;
  Double_t labelSizeX = 0.035;
  Double_t labelSizeY = 0.035;
  Int_t titleFont = 42;
  
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX);
  label_p->SetNDC();

  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.02);
    canv_p->SetRightMargin(0.02);
    canv_p->SetLeftMargin(0.15);
    canv_p->SetBottomMargin(0.15);

    std::string histName = "jtMultPerCent_" + centBinsStr[cI] + "_h";
    TH1F* inHist_p = (TH1F*)inFile_p->Get(histName.c_str());

    inHist_p->SetMarkerSize(sizes[0]);
    inHist_p->SetMarkerStyle(styles[0]);
    inHist_p->SetMarkerColor(colors[0]);
    inHist_p->SetLineColor(colors[0]);

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

    inHist_p->DrawCopy("HIST E1 P");
    gStyle->SetOptStat(0);

    std::vector<std::string> tempLabels = globalLabels;
    tempLabels.push_back("#bf{" + std::to_string((int)centBins[cI]) + "-" + std::to_string((int)centBins[cI+1]) + "%}");
    for(unsigned int gI = 0; gI < tempLabels.size(); ++gI){
      label_p->DrawLatex(mixedJtPtLabelX, mixedJtPtLabelY - gI*0.055, tempLabels[gI].c_str());
    }
    
    histName = "jtMultPerCent_" + centBinsStr2[cI] + "_h";
    std::string saveName = "pdfDir/" + dateStr + "/" + histName + "_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;    
  }

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
