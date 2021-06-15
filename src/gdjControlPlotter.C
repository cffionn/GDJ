//Author: Chris McGinn (2021.04.26)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <map>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjControlPlotter(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1;

  const std::string dateStr = getDateStr();

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> reqParams = {"INFILENAME", "SAVEEXT"};  
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, "root")) return 1;

  const std::string saveExtStr = inConfig_p->GetValue("SAVEEXT", "");
  std::vector<std::string> validExt = {"pdf", "png"};
  if(!vectContainsStr(saveExtStr, &validExt)){
    std::cout << "Given extension \'" << saveExtStr << "\' is not valid. Please pick from:" << std::endl;
    for(unsigned int eI = 0; eI < validExt.size(); ++eI){
      std::cout << " " << validExt[eI] << std::endl;
    }
    std::cout << "return 1. " << std::endl;
    return 1;
  }

  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> passingHistNames = {"vzPassing",
					       "centPassing",
					       "truPhoPt",
					       "leadingPhoPt",
					       "leadingPhoEta",
					       "leadingPhoIso",
					       "leadingPhoTight",
					       "jtEtaPassing",
					       "jtPtPassing",
					       "jtDPhiPassing",
					       "jtGammaDRPassing"};

  const Double_t leftMargin = 0.15;
  const Double_t rightMargin = 0.01;
  const Double_t topMargin = 0.05;
  const Double_t bottomMargin = leftMargin - (topMargin - rightMargin);
  const Double_t height = 600.0;
  const Double_t width = height*(1.0 - topMargin - bottomMargin)/(1.0 - leftMargin - rightMargin);
    
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int hI = 0; hI < passingHistNames.size(); ++hI){
    TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
    canv_p->SetTopMargin(topMargin);
    canv_p->SetLeftMargin(leftMargin);
    canv_p->SetRightMargin(rightMargin);
    canv_p->SetBottomMargin(bottomMargin);
    
    TH1F* hist_p = (TH1F*)inFile_p->Get((passingHistNames[hI] + "_h").c_str());

    HIJet::Style::EquipHistogram(hist_p, 0);
    hist_p->DrawCopy("HIST E1 P");

    gPad->SetTicks();
    gStyle->SetOptStat(0);

    std::string saveName = "pdfDir/" + dateStr + "/" + passingHistNames[hI] + "_" + dateStr + "." + saveExtStr;
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  } 

  inFile_p->Close();
  delete inFile_p;

  std::cout << "GDJCONTROLPLOTTER COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjControlPlotter.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjControlPlotter(argv[1]);
  return retVal;
}
