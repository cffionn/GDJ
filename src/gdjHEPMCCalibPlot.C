//Author: Chris McGinn (2022.12.02)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"

int gdjHEPMCCalibPlot(std::string inConfigFileName)
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
  //pdf of current job will all go to date tagged subdir
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INPBPBFILENAME",
					      "INPPFILENAME",
                                              "EXTENSION"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  const std::string inPbPbFileName = config_p->GetValue("INPBPBFILENAME", "");
  const std::string inPPFileName = config_p->GetValue("INPPFILENAME", "");
  const std::string extStr = config_p->GetValue("EXTENSION", "");

  //Check extension given is valid
  std::vector<std::string> validExtensions = {"pdf", "png"};
  if(!vectContainsStr(extStr, &validExtensions)){
    std::cout << "Given extension \'" << extStr << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  //Check pp file given is valid
  if(!check.checkFileExt(inPPFileName, ".root")) return 1;

  //Check Pb+Pb file given is valid
  if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;
  
  TFile* inPPFile_p = new TFile(inPPFileName.c_str(), "READ");
  TH1D* ppHist_p = (TH1D*)inPPFile_p->Get("jetPt_R4_h");

  TFile* inPbPbFile_p = new TFile(inPbPbFileName.c_str(), "READ");
  TH1D* pbpbHist_p = (TH1D*)inPbPbFile_p->Get("jetPt_R4_h");
  TH1D* raaHist_p = (TH1D*)inPbPbFile_p->Get("Hist1D_y1");
  TH1D* raaHist_SysUp_p = (TH1D*)inPbPbFile_p->Get("Hist1D_y1_e1plus");
  TH1D* raaHist_SysDown_p = (TH1D*)inPbPbFile_p->Get("Hist1D_y1_e1minus");
  TH1D* raaHist_Stat_p = (TH1D*)inPbPbFile_p->Get("Hist1D_y1_e2");

  for(Int_t bIX = 0; bIX < raaHist_p->GetNbinsX(); ++bIX){
    raaHist_p->SetBinError(bIX+1, raaHist_Stat_p->GetBinContent(bIX+1));
  }
  
  //Define some margins
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.16;
  const Double_t topMargin = 0.03;
  const Double_t rightMargin = 0.03;

    //Also define some fonts and sizes                                                                    
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;

  TCanvas* canv_p = new TCanvas("canv", "", 900, 900);
  setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin);

  HIJet::Style::EquipHistogram(raaHist_p, 2);
  raaHist_p->SetMaximum(1.1);
  raaHist_p->SetMinimum(0.0);
  raaHist_p->SetTitle("");
  raaHist_p->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  raaHist_p->GetXaxis()->SetTitleOffset(1.7);

  raaHist_p->GetXaxis()->SetTitleFont(titleFont);
  raaHist_p->GetYaxis()->SetTitleFont(titleFont);
  raaHist_p->GetXaxis()->SetTitleSize(titleSize);
  raaHist_p->GetYaxis()->SetTitleSize(titleSize);
 
  raaHist_p->GetXaxis()->SetLabelFont(titleFont);
  raaHist_p->GetYaxis()->SetLabelFont(titleFont);
  raaHist_p->GetXaxis()->SetLabelSize(labelSize);
  raaHist_p->GetYaxis()->SetLabelSize(labelSize);
 
  raaHist_p->DrawCopy("HIST E1 P");

  std::cout << "RAA OFFSET: " << raaHist_p->GetXaxis()->GetTitleOffset() << std::endl;
  
  pbpbHist_p->Divide(ppHist_p);
  HIJet::Style::EquipHistogram(pbpbHist_p, 1);
  pbpbHist_p->DrawCopy("HIST E1 P SAME");

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(pbpbHist_p->GetBinLowEdge(1), 1.0, pbpbHist_p->GetBinLowEdge(pbpbHist_p->GetNbinsX()+1), 1.0);
  delete line_p;

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSize);
  label_p->DrawLatex(185, -0.05, "200");
  label_p->DrawLatex(370, -0.05, "400");
  label_p->DrawLatex(740, -0.05, "800");

  delete label_p;
  
  gStyle->SetOptStat(0);
  gPad->SetTicks();
  gPad->SetLogx();
  
  std::string saveName = pdfStr + "/jetCalibration_" + dateStr + "." + extStr;
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  
  inPbPbFile_p->Close();
  delete inPbPbFile_p;

  inPPFile_p->Close();
  delete inPPFile_p;
  
  std::cout << "GDJHEPMCCALIBPLOT complete. return 0" << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHEPMCCalibPlot.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHEPMCCalibPlot(argv[1]);
  return retVal;
}
