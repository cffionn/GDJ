//Author: Chris McGinn (2023.01.18)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TPad.h"
#include "TStyle.h"

//LOCAL
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjRunStabilityPlotter(std::string inConfigFileName)
{
  const int titleFont = 42;
  const double titleSize = 0.025;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  const std::string pdfStr = "pdfDir/" + dateStr;
  check.doCheckMakeDir(pdfStr);

  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  //Define necessary params for the input config file  
  std::vector<std::string> necessaryParams = {"INFILENAME",
					      "LABELALIGNRIGHT1",
					      "LABELALIGNRIGHT2",
					      "LABELX1",
					      "LABELX2",
					      "LABELY1",
					      "LABELY2",
                                              "SAVEEXT"};

  //Check all params exist
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Grab input params
  const std::string inFileName = config_p->GetValue("INFILENAME", "");

  const Double_t labelX1 = config_p->GetValue("LABELX1", -1.0);
  const Double_t labelX2 = config_p->GetValue("LABELX2", -1.0);
  const Double_t labelY1 = config_p->GetValue("LABELY1", -1.0);
  const Double_t labelY2 = config_p->GetValue("LABELY2", -1.0);

  const Bool_t labelAlignRight1 = (Bool_t)config_p->GetValue("LABELALIGNRIGHT1", 0);
  const Bool_t labelAlignRight2 = (Bool_t)config_p->GetValue("LABELALIGNRIGHT2", 0);
  
  const std::string saveExt = config_p->GetValue("SAVEEXT", "");

  //Check file exists
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  //Grab file, its internal config
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> inNecessaryParams = {"ISMC",
						"ISPP",
						"RUNSLIST"};

  //Check internal config has all req params
  if(!checkEnvForParams(inConfig_p, inNecessaryParams)) return 1;

  //Grab the req params
  Bool_t isMC = (Bool_t)inConfig_p->GetValue("ISMC", 0);
  Bool_t isPP = (Bool_t)inConfig_p->GetValue("ISPP", 0);
  std::string ppPbPbStr = "PP";
  std::string ppPbPbLabel = "pp 2017 Data";
  if(!isPP){
    ppPbPbStr = "PbPb";
    ppPbPbLabel = "Pb+Pb 2018 Data";
  }
  
  std::vector<std::string> runNumList = strToVect(inConfig_p->GetValue("RUNSLIST", ""));

  if(isMC){
    std::cout << "GDJRUNSTABILITYPLOTTER ERROR: Given file \'" << inFileName << "\' is MC, run list stability plot not defined. return 1" << std::endl;
    return 1;
  }
  
  TH1D* hist_p = (TH1D*)inFile_p->Get("eventsByRun_h");
  TH1D* hist_LumiScaled_p = (TH1D*)inFile_p->Get("eventsByRun_LumiScaled_h");

  std::vector<TH1D*> hists_p = {hist_p, hist_LumiScaled_p};
  std::vector<std::string> labels = {"NotLumiScaled", "LumiScaled"};
  std::vector<std::string> labels2 = {"Not Lumi. Scaled", "Lumi. Scaled"};
  std::vector<Double_t> labelXs = {labelX1, labelX2};
  std::vector<Double_t> labelYs = {labelY1, labelY2};
  std::vector<Bool_t> labelAligns = {labelAlignRight1, labelAlignRight2};
  
  //Define a bunch of parameters for the TCanvas
  const Int_t nPad = 2;
  Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit - 0.01;
  const Double_t topMargin = 0.05;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;
  Double_t height = 900.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);
  height = 1200.0;

  
  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
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

    HIJet::Style::EquipHistogram(hists_p[hI], 0);
    hists_p[hI]->GetXaxis()->SetTitleFont(titleFont);
    hists_p[hI]->GetXaxis()->SetLabelFont(titleFont);
    hists_p[hI]->GetYaxis()->SetTitleFont(titleFont);
    hists_p[hI]->GetYaxis()->SetLabelFont(titleFont);
 
    hists_p[hI]->GetXaxis()->SetTitleSize(0.00001);
    hists_p[hI]->GetXaxis()->SetLabelSize(0.00001);
    hists_p[hI]->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));
    hists_p[hI]->GetYaxis()->SetLabelSize(titleSize/(1.0-padSplit));

    hists_p[hI]->DrawCopy("HIST E1 P");

    
    TLatex* label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(titleFont);
    label_p->SetTextSize(titleSize/(1.0 - padSplit));
    if(labelAligns[hI]) label_p->SetTextAlign(31);

    std::vector<std::string> tempLabels = {"#bf{#it{ATLAS Internal}}", ppPbPbLabel, "After #gamma selection", labels2[hI]};
    for(Int_t tI = 0; tI < tempLabels.size(); ++tI){
      label_p->DrawLatex(labelXs[hI], labelYs[hI] - 0.06*tI, tempLabels[tI].c_str());
    }

    delete label_p;
    
    gPad->SetTicks();

    std::vector<Double_t> vals;
    for(Int_t bIX = 0; bIX < hists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
      vals.push_back(hists_p[hI]->GetBinContent(bIX+1));
    }

    if(vals.size() != runNumList.size()){
      std::cout << "GDJRUNSTABILITYPLOTTER ERROR: Values size \'" << vals.size() << "\' doesn't match runsList size \'" << runNumList.size() << "\'" << std::endl;
    }

    Double_t median = vals[vals.size()/2];
    if(vals.size()%2 == 0) median = (median + vals[vals.size()/2 -1])/2.0;

    hists_p[hI]->Scale(1.0/median);
    pads_p[1]->cd();
    std::string yTitle = hists_p[hI]->GetYaxis()->GetTitle();
    yTitle = "#frac{" + yTitle + "}{Median}";
    hists_p[hI]->GetYaxis()->SetTitle(yTitle.c_str());
    
    hists_p[hI]->GetXaxis()->SetTitleSize(titleSize/(padSplit));
    hists_p[hI]->GetXaxis()->SetLabelSize(titleSize/(padSplit));
    hists_p[hI]->GetYaxis()->SetTitleSize(titleSize/(padSplit));
    hists_p[hI]->GetYaxis()->SetLabelSize(titleSize/(padSplit));

    hists_p[hI]->DrawCopy("HIST E1 P");    

    gPad->SetTicks();
    
    gStyle->SetOptStat(0);
    
    std::string saveName = pdfStr + "/runStability_" + ppPbPbStr + "_" + labels[hI] + "_" + dateStr + "." + saveExt;
    quietSaveAs(canv_p, saveName);
    for(Int_t pI = 0; pI < nPad; ++pI){
      delete pads_p[pI];
    }   
    delete canv_p;
  }
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

  
int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjRunStabilityPlotter.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjRunStabilityPlotter(argv[1]);
  return retVal;
}
