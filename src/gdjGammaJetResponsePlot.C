//Author: Chris McGinn (2020.03.31)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <algorithm>
#include <iostream>
#include <map>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjGammaJetResponsePlot(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  const std::string dateStr = getDateStr();
  const std::string outDirStr = "pdfDir/" + dateStr;

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(42);
  label_p->SetTextSize(0.035);
  label_p->SetNDC();
  
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(outDirStr);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* configEnv_p = (TEnv*)inFile_p->Get("config");
  TEnv* labelEnv_p = (TEnv*)inFile_p->Get("label");
  configParser config(configEnv_p);
  configParser label(labelEnv_p);

  std::map<std::string, std::string> labelMap = label.GetConfigMap();
  std::vector<int> recoGammaPtVals, genGammaPtVals;
  std::vector<std::string> centBinsStr, recoGammaPtStr, genGammaPtStr;
  
  for(auto const & label : labelMap){
    if(label.first.find("RecoGammaPt") != std::string::npos || label.first.find("GenGammaPt") != std::string::npos){
      std::string testVal = label.first.substr(label.first.find("Pt")+2, label.first.size());

      if(label.first.find("RecoGammaPt") != std::string::npos) recoGammaPtVals.push_back(std::stoi(testVal));
      if(label.first.find("GenGammaPt") != std::string::npos) genGammaPtVals.push_back(std::stoi(testVal));
    }
    else if(label.first.find("Cent") != std::string::npos) centBinsStr.push_back(label.first);
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::sort(std::begin(recoGammaPtVals), std::end(recoGammaPtVals));
  std::sort(std::begin(genGammaPtVals), std::end(genGammaPtVals));

  recoGammaPtVals.erase(recoGammaPtVals.begin()+recoGammaPtVals.size()-1);
  genGammaPtVals.erase(genGammaPtVals.begin()+genGammaPtVals.size()-1);
  
  for(unsigned int rI = 0; rI < recoGammaPtVals.size(); ++rI){
    recoGammaPtStr.push_back("RecoGammaPt" + std::to_string(recoGammaPtVals[rI]));
  }
  for(unsigned int gI = 0; gI < genGammaPtVals.size(); ++gI){
    genGammaPtStr.push_back("GenGammaPt" + std::to_string(genGammaPtVals[gI]));
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  Int_t nMaxPtBins = 200;
  const Int_t nGenGammaPtBins = genGammaPtStr.size();
  const Int_t nRecoGammaPtBins = recoGammaPtStr.size();
  Double_t genGammaPtBins[nMaxPtBins+1];
  Double_t recoGammaPtBins[nMaxPtBins+1];

  
  for(int gI = 0; gI < nGenGammaPtBins; ++gI){
    std::string valStr = labelMap[genGammaPtStr[gI]];
    valStr.replace(0, valStr.rfind("<")+1, "");
    genGammaPtBins[gI+1] = std::stof(valStr);

    if(gI == 0){
      std::string valStr = labelMap[genGammaPtStr[gI]];
      valStr.replace(valStr.find("<"), valStr.size(), "");
      genGammaPtBins[0] = std::stof(valStr);
    }
  }

  for(int rI = 0; rI < nRecoGammaPtBins; ++rI){
    std::string valStr = labelMap[recoGammaPtStr[rI]];
    valStr.replace(0, valStr.rfind("<")+1, "");
    recoGammaPtBins[rI+1] = std::stof(valStr);

    if(rI == 0){
      std::string valStr = labelMap[recoGammaPtStr[rI]];
      valStr.replace(valStr.find("<"), valStr.size(), "");
      recoGammaPtBins[0] = std::stof(valStr);
    }
  }

  std::cout << "Gen: ";
  for(Int_t gI = 0; gI < nGenGammaPtBins; ++gI){
    std::cout << genGammaPtBins[gI] << ", ";
  }
  std::cout << genGammaPtBins[nGenGammaPtBins] << "." << std::endl;

  std::cout << "Reco: ";
  for(Int_t gI = 0; gI < nRecoGammaPtBins; ++gI){
    std::cout << recoGammaPtBins[gI] << ", ";
  }
  std::cout << recoGammaPtBins[nRecoGammaPtBins] << "." << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  double subMarginRight = 0.01;
  double subMarginTop = 0.045;
  double margin = 0.12;
  for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450*nRecoGammaPtBins, 450*nGenGammaPtBins);
    canv_p->SetTopMargin(subMarginTop);
    canv_p->SetRightMargin(subMarginRight);
    canv_p->SetLeftMargin(margin);
    canv_p->SetBottomMargin(margin);
    TH2F* hist_p = new TH2F("hist_h", ";Reco. p_{T,#gamma};Gen. p_{T,#gamma}", nRecoGammaPtBins, recoGammaPtBins, nGenGammaPtBins, genGammaPtBins);

    hist_p->GetXaxis()->SetTitleOffset(1.3);
    centerTitles(hist_p);   
    std::cout << "OFF: " << hist_p->GetXaxis()->GetTitleOffset() << std::endl;
    hist_p->DrawCopy("COL");    

    std::cout << "Y FONT, SIZE: " << hist_p->GetYaxis()->GetTitleFont() << ", " << hist_p->GetYaxis()->GetTitleSize() << std::endl;
    
    gPad->SetLogy();
    gPad->SetLogx();

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TPad* pad_p[nMaxPtBins][nMaxPtBins];
    
    for(int rI = 0; rI < nRecoGammaPtBins; ++rI){
      double xValLow = margin + rI*(1.0 - margin - subMarginRight)/(double)nRecoGammaPtBins;
      double xValHigh = margin + (rI+1)*(1.0 - margin - subMarginRight)/(double)nRecoGammaPtBins;

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(int gI = 0; gI < nGenGammaPtBins; ++gI){
	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGenGammaPtBins << std::endl;

	double yValLow = margin + gI*(1.0 - margin - subMarginTop)/(double)nRecoGammaPtBins;
	double yValHigh = margin + (gI+1)*(1.0 - margin - subMarginTop)/(double)nRecoGammaPtBins;

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGenGammaPtBins << std::endl;

	canv_p->cd();
	pad_p[rI][gI] = new TPad(("pad_" + std::to_string(rI) + "_" + std::to_string(gI)).c_str(), "", xValLow, yValLow, xValHigh, yValHigh);
	pad_p[rI][gI]->SetLeftMargin(0.06);
	pad_p[rI][gI]->SetBottomMargin(0.06);
	pad_p[rI][gI]->SetTopMargin(0.01);
	pad_p[rI][gI]->SetRightMargin(0.14);
	pad_p[rI][gI]->Draw("SAME");
	pad_p[rI][gI]->cd();
	
	TH2F* tempHist_p = (TH2F*)inFile_p->Get((centBinsStr[cI] + "/photonJtGenResVCentGenPtRecoPt_" + centBinsStr[cI] + "_" + genGammaPtStr[gI] + "_" + recoGammaPtStr[rI] + "_h").c_str());

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGenGammaPtBins << std::endl;
	if(gI == 4 && rI == 0) std::cout << tempHist_p->GetName() << std::endl;
	tempHist_p->DrawCopy("COLZ");
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetLogz();

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGenGammaPtBins << std::endl;
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    }
    
    gStyle->SetOptStat(0);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    canv_p->cd();
    gPad->SetTicks();
    gPad->RedrawAxis();

    std::string centStr = centBinsStr[cI].substr(4, centBinsStr[cI].size());
    centStr.replace(centStr.find("to"), 2, "-");
    centStr = centStr + "%";
    
    label_p->DrawLatex(margin+0.05, 1.0 - subMarginTop*4./5., ("ATLAS PYTHIA+Overlay, 5.02 TeV Pb+Pb, " + centStr).c_str());
    
    std::string saveName = outDirStr + "/gammaJetResponse_" + centBinsStr[cI] + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    for(int rI = 0; rI < nRecoGammaPtBins; ++rI){
      for(int gI = 0; gI < nGenGammaPtBins; ++gI){
	delete pad_p[rI][gI];
      }
    }
    delete canv_p;
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    delete hist_p;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->Close();
  delete inFile_p;

  std::cout << "GDJGAMMAJETRESPONSEPLOT COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjGammaJetResponsePlot.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjGammaJetResponsePlot(argv[1]);
  return retVal;
}
