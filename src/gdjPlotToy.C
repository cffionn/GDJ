//Author: Chris McGinn (2020.08.12)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjPlotToy(std::string inConfigFileName)
{
  const Int_t nStyles = 15;
  const Int_t styles[nStyles] = {20,21,33,34,29,24,25,27,28,30,23,20,21,33,34};

  const int nColors = 15;
  const int colors[nColors]={1,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kCyan+3,28,41,kGray};

  const int nSizes = 15;
  const double sizes[nSizes] = {1,1,1.6,1.2,1.6,1,1,1,1.6,1,1,1,1,1.6,1.2};

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"GLOBALTAG",
					"MIXEDYMIN",
					"MIXEDYMAX",
					"MIXEDYRATMIN",
					"MIXEDYRATMAX",
					"MIXEDLABELX",
					"MIXEDLABELY",
					"MIXEDLEGX",
					"MIXEDLEGY"};  

  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string globalTag = inConfig_p->GetValue("GLOBALTAG", "");
  
  const double mixedYMin = inConfig_p->GetValue("MIXEDYMIN", 100.0);
  const double mixedYMax = inConfig_p->GetValue("MIXEDYMAX", -1.0);
  const double mixedYRatMin = inConfig_p->GetValue("MIXEDYRATMIN", 100.0);
  const double mixedYRatMax = inConfig_p->GetValue("MIXEDYRATMAX", -1.0);

  const double mixedLabelX = inConfig_p->GetValue("MIXEDLABELX", 100.0);
  const double mixedLabelY = inConfig_p->GetValue("MIXEDLABELY", 100.0);
  const double mixedLegX = inConfig_p->GetValue("MIXEDLEGX", 100.0);
  const double mixedLegY = inConfig_p->GetValue("MIXEDLEGY", 100.0);

  std::vector<std::string> globalLabels;
  for(Int_t gI = 0; gI < 20; ++gI){
    std::string globalLabel = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(globalLabel.size() == 0) continue;
    globalLabels.push_back(globalLabel);
  }

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir/" + dateStr);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> fileReqParams = {"NEVT"};
  if(!checkEnvForParams(fileConfig_p, fileReqParams)) return 1;

  const ULong64_t nEvt = fileConfig_p->GetValue("NEVT", 0);

  globalLabels.push_back(prettyStringE(nEvt,1,false) + " Toys");
  
  std::vector<std::string> histNames = {"gammaHist_h",
					"jet1Hist_h", 
					"jet2Hist_h",
					"jetBkgdHist_h",
					"signalHist_h",
					"bkgdHist_h", 
					"mixedBkgdHist_h", 
					"pureBkgdHist_h", 
					"mixedHist_h",	
					"mixedHistCorrection_h",	
					"mixedHistTrue_h",	
					"signalAndBkgdHist_h"};

  const double topMargin = 0.04;
  const double rightMargin = topMargin;
  const double leftMargin = 0.16;
  const double bottomMargin = leftMargin;
  const double tinyMargin = 0.000001;

  const Double_t titleSizeX = 0.04;
  const Double_t titleSizeY = 0.04;
  const Double_t labelSizeX = 0.035;
  const Double_t labelSizeY = 0.035;

  const Int_t titleFont = 42;

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX);
  label_p->SetNDC();
  
  for(unsigned int hI = 0; hI < histNames.size(); ++hI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(topMargin);
    canv_p->SetRightMargin(rightMargin);
    canv_p->SetLeftMargin(leftMargin);
    canv_p->SetBottomMargin(bottomMargin);

    TH1F* hist_p = (TH1F*)inFile_p->Get(histNames[hI].c_str());

    hist_p->SetMarkerStyle(styles[0]);
    hist_p->SetMarkerSize(sizes[0]);
    hist_p->SetMarkerColor(colors[0]);
    hist_p->SetLineColor(colors[0]);

    hist_p->GetXaxis()->SetTitleFont(titleFont);
    hist_p->GetYaxis()->SetTitleFont(titleFont);
    hist_p->GetXaxis()->SetLabelFont(titleFont);
    hist_p->GetYaxis()->SetLabelFont(titleFont);

    hist_p->GetXaxis()->SetTitleSize(titleSizeX);
    hist_p->GetYaxis()->SetTitleSize(titleSizeY);
    hist_p->GetXaxis()->SetLabelSize(labelSizeX);
    hist_p->GetYaxis()->SetLabelSize(labelSizeY);

    hist_p->GetXaxis()->SetTitleOffset(1.2);
    hist_p->GetYaxis()->SetTitleOffset(1.2);
    
    hist_p->DrawCopy("HIST E1 P");
    
    gStyle->SetOptStat(0);
    
    std::string saveName = "pdfDir/" + dateStr + "/" + histNames[hI] + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  }
  
  TCanvas* canv_p =  new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(tinyMargin);
  canv_p->SetRightMargin(tinyMargin);
  canv_p->SetLeftMargin(tinyMargin);
  canv_p->SetBottomMargin(tinyMargin);

  const Int_t nPads = 2;
  const Double_t padSplit = 0.35;

  TLegend* leg_p = new TLegend(mixedLegX, mixedLegY - 0.063*3, mixedLegX + 0.25, mixedLegY);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSizeX/(1.0-padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  TPad* pads_p[nPads];
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetTopMargin(topMargin/(1.0-padSplit));
  pads_p[0]->SetRightMargin(rightMargin);
  pads_p[0]->SetBottomMargin(tinyMargin);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  
  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetTopMargin(tinyMargin);
  pads_p[1]->SetRightMargin(rightMargin);
  pads_p[1]->SetBottomMargin(leftMargin/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");

  canv_p->cd();
  pads_p[0]->cd();
  
  TH1F* bkgdHist_p = (TH1F*)inFile_p->Get("bkgdHist_h");
  TH1F* mixedHist_p = (TH1F*)inFile_p->Get("mixedHist_h");
  TH1F* mixedHistCorrection_p = (TH1F*)inFile_p->Get("mixedHistCorrection_h");
  TH1F* signalAndBkgdHist_p = (TH1F*)inFile_p->Get("signalAndBkgdHist_h");
  
  setSumW2({bkgdHist_p, mixedHist_p, mixedHistCorrection_p, signalAndBkgdHist_p});
  
  mixedHist_p->SetMarkerStyle(styles[1]);
  mixedHist_p->SetMarkerSize(sizes[1]);
  mixedHist_p->SetMarkerColor(colors[1]);
  mixedHist_p->SetLineColor(colors[1]);

  mixedHist_p->GetXaxis()->SetTitleFont(titleFont);
  mixedHist_p->GetYaxis()->SetTitleFont(titleFont);
  mixedHist_p->GetXaxis()->SetLabelFont(titleFont);
  mixedHist_p->GetYaxis()->SetLabelFont(titleFont);
  
  mixedHist_p->GetXaxis()->SetTitleSize(titleSizeX/(1.0 - padSplit));
  mixedHist_p->GetYaxis()->SetTitleSize(titleSizeY/(1.0 - padSplit));
  mixedHist_p->GetXaxis()->SetLabelSize(labelSizeX/(1.0 - padSplit));
  mixedHist_p->GetYaxis()->SetLabelSize(labelSizeY/(1.0 - padSplit));    

  mixedHist_p->Add(mixedHistCorrection_p, -1);
  
  bkgdHist_p->SetMarkerStyle(styles[0]);
  bkgdHist_p->SetMarkerSize(sizes[0]);
  bkgdHist_p->SetMarkerColor(colors[0]);
  bkgdHist_p->SetLineColor(colors[0]);

  signalAndBkgdHist_p->SetMarkerStyle(styles[2]);
  signalAndBkgdHist_p->SetMarkerSize(sizes[2]);
  signalAndBkgdHist_p->SetMarkerColor(colors[2]);
  signalAndBkgdHist_p->SetLineColor(colors[2]);

  mixedHist_p->GetYaxis()->SetTitleOffset(1.05);
  
  mixedHist_p->SetMinimum(mixedYMin);
  mixedHist_p->SetMaximum(mixedYMax);

  mixedHist_p->DrawCopy("HIST E1 P");
  signalAndBkgdHist_p->DrawCopy("HIST E1 P SAME");
  bkgdHist_p->DrawCopy("HIST E1 P SAME");

  label_p->SetTextSize(titleSizeX/(1.0 - padSplit));
  for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
    label_p->DrawLatex(mixedLabelX, mixedLabelY - gI*0.08, globalLabels[gI].c_str());
  }

  leg_p->AddEntry(mixedHist_p, "Mixed Event", "P L");
  leg_p->AddEntry(bkgdHist_p, "Pure Bkgd", "P L");
  leg_p->AddEntry(signalAndBkgdHist_p, "Signal+Bkgd", "P L");

  leg_p->Draw("SAME");
  gStyle->SetOptStat(0);
  
  canv_p->cd();
  pads_p[1]->cd();

  bkgdHist_p->GetXaxis()->SetTitleFont(titleFont);
  bkgdHist_p->GetYaxis()->SetTitleFont(titleFont);
  bkgdHist_p->GetXaxis()->SetLabelFont(titleFont);
  bkgdHist_p->GetYaxis()->SetLabelFont(titleFont);
  
  bkgdHist_p->GetXaxis()->SetTitleSize(titleSizeX/padSplit);
  bkgdHist_p->GetYaxis()->SetTitleSize(titleSizeY/padSplit);
  bkgdHist_p->GetXaxis()->SetLabelSize(labelSizeX/padSplit);
  bkgdHist_p->GetYaxis()->SetLabelSize(labelSizeY/padSplit);    
  
  bkgdHist_p->Divide(mixedHist_p);
  bkgdHist_p->GetYaxis()->SetNdivisions(505);
  bkgdHist_p->GetYaxis()->SetTitle("Ratio");
  bkgdHist_p->GetYaxis()->SetTitleOffset(1.05*padSplit/(1.0 - padSplit));
  bkgdHist_p->GetXaxis()->SetTitleOffset(1.05);

  bkgdHist_p->SetMinimum(mixedYRatMin);
  bkgdHist_p->SetMaximum(mixedYRatMax);

  bkgdHist_p->DrawCopy("HIST E1 P");
  
  std::string saveName = "pdfDir/" + dateStr + "/mixedEventToy_" + globalTag + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete pads_p[0];
  delete pads_p[1];
  delete canv_p;
  
  delete leg_p;
  delete label_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  delete inConfig_p;
  
  std::cout << "GDJPLOTTOY COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotToy.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotToy(argv[1]);
  return retVal;
}
