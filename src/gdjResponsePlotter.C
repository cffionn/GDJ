//Author: Chris McGinn (2020.04.20)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
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
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

void plotMeanAndSigma(const bool doGlobalDebug, std::map<std::string, std::string>* labelMap, std::string dateStr, std::vector<std::string> centBins, TH1F* histMean_p[], TH1F* histWidth_p[])
{
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::string overMeanStr = histWidth_p[0]->GetName();
  std::string labelStr = histWidth_p[0]->GetName();

  std::cout << "LABEL STR: " << labelStr << std::endl;
  
  if(overMeanStr.find("OverMean") != std::string::npos) overMeanStr = "OverMean";
  else overMeanStr = "";

  labelStr.replace(0, labelStr.find("_")+1, "");
  labelStr.replace(labelStr.rfind("_h"), 2, "");
  std::vector<std::string> preLabels;
  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    if(centStr.find("Cent") == std::string::npos) preLabels.push_back(centStr);
    labelStr.replace(0, labelStr.find("_")+1, "");
  }
  if(labelStr.size() != 0 && labelStr.find("Cent") == std::string::npos) preLabels.push_back(labelStr);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  
  Double_t padSplit = 0.35;
  Double_t leftMargin = 0.11;

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  TPad* pads_p[2];
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);
  canv_p->cd();
  
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(0.01);
  pads_p[0]->SetRightMargin(0.01);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(0.001);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetTopMargin(0.001);
  pads_p[1]->SetRightMargin(0.01);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetBottomMargin(leftMargin/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  bool isPP = centBins[0].find("PP") != std::string::npos;
  kirchnerPalette kPal;
  const int nMarkers = 4;
  int markers[nMarkers] = {24,25,28,46};
  const int nColors = 4;
  int colors[nColors] = {0,1,3,4};

  TLegend* leg_p = new TLegend(0.35, 0.6, 0.9, 0.9);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(12);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  int nCentBins = 1;
  if(!isPP) nCentBins = centBins.size()-1;

  TF1* csnFit_p = new TF1("csnFit_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
  csnFit_p->SetLineStyle(2);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string legStr = "PP";
    if(!isPP) legStr = "PbPb, " + centBins[cI] + "-" + centBins[cI+1] + "%";
    
    canv_p->cd();
    pads_p[0]->cd();

    histWidth_p[cI]->SetMarkerSize(1);
    histWidth_p[cI]->SetMarkerStyle(markers[cI%nMarkers]);
    histWidth_p[cI]->SetMarkerColor(kPal.getColor(colors[cI%nColors]));
    histWidth_p[cI]->SetLineColor(kPal.getColor(colors[cI%nColors]));
    
    histWidth_p[cI]->SetMaximum(0.65);
    histWidth_p[cI]->SetMinimum(0.0);

    histWidth_p[cI]->GetXaxis()->SetTitleFont(43);
    histWidth_p[cI]->GetYaxis()->SetTitleFont(43);
    histWidth_p[cI]->GetXaxis()->SetTitleSize(14);
    histWidth_p[cI]->GetYaxis()->SetTitleSize(14);

    histWidth_p[cI]->GetXaxis()->SetLabelFont(43);
    histWidth_p[cI]->GetYaxis()->SetLabelFont(43);
    histWidth_p[cI]->GetXaxis()->SetLabelSize(12);
    histWidth_p[cI]->GetYaxis()->SetLabelSize(12);

    histWidth_p[cI]->GetYaxis()->SetNdivisions(505);
    histWidth_p[cI]->GetXaxis()->SetNdivisions(505);

    centerTitles({histWidth_p[cI], histMean_p[cI]});
    
    if(cI == 0) histWidth_p[cI]->DrawCopy("HIST E1 P");
    else histWidth_p[cI]->DrawCopy("HIST E1 P SAME");

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    csnFit_p->SetParameter(0, 0.05);
    csnFit_p->SetParameter(1, 1.0);
    csnFit_p->SetParameter(2, 0.00);

    csnFit_p->SetParLimits(0, 0.0, 1.0);
    csnFit_p->SetParLimits(1, 0.0, 10.0);
    csnFit_p->SetParLimits(2, 0.0, 100.0);
    
    histWidth_p[cI]->Fit("csnFit_p", "Q M N E", "", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    csnFit_p->SetMarkerColor(kPal.getColor(colors[cI%nColors]));
    csnFit_p->SetLineColor(kPal.getColor(colors[cI%nColors]));

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    csnFit_p->DrawCopy("SAME");
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    legStr = legStr + "; N=" + prettyString(csnFit_p->GetParameter(2),1,false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    legStr = legStr + ",S=" + prettyString(csnFit_p->GetParameter(1), 2, false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    std::cout << "CSN: " << csnFit_p->GetParameter(0) << std::endl;
    legStr = legStr + ",C=" + prettyString(csnFit_p->GetParameter(0),3,false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    leg_p->AddEntry(histMean_p[cI], legStr.c_str(), "P L");

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    canv_p->cd();
    pads_p[1]->cd();

    histMean_p[cI]->SetMarkerSize(1);
    histMean_p[cI]->SetMarkerStyle(markers[cI%nMarkers]);
    histMean_p[cI]->SetMarkerColor(kPal.getColor(colors[cI%nColors]));
    histMean_p[cI]->SetLineColor(kPal.getColor(colors[cI%nColors]));

    histMean_p[cI]->SetMaximum(1.15);
    histMean_p[cI]->SetMinimum(0.85);

    histMean_p[cI]->GetXaxis()->SetTitleFont(43);
    histMean_p[cI]->GetYaxis()->SetTitleFont(43);
    histMean_p[cI]->GetXaxis()->SetTitleSize(14);
    histMean_p[cI]->GetYaxis()->SetTitleSize(14);

    histMean_p[cI]->GetXaxis()->SetLabelFont(43);
    histMean_p[cI]->GetYaxis()->SetLabelFont(43);
    histMean_p[cI]->GetXaxis()->SetLabelSize(12);
    histMean_p[cI]->GetYaxis()->SetLabelSize(12);

    histMean_p[cI]->GetYaxis()->SetNdivisions(505);
    histMean_p[cI]->GetXaxis()->SetNdivisions(505);

    histMean_p[cI]->GetXaxis()->SetTitleOffset(2.8);
    
    if(cI == 0) histMean_p[cI]->DrawCopy("HIST E1 P");
    else histMean_p[cI]->DrawCopy("HIST E1 P SAME");    
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  canv_p->cd();
  pads_p[0]->cd();
  leg_p->Draw("SAME");

  labelStr = "";

  const double xPos = 0.18;
  double yPos = 0.9;
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    if(labelMap->count(preLabels[pI]) != 0) preLabels[pI] = (*labelMap)[preLabels[pI]];
    
    if(labelStr.size() + preLabels[pI].size() > 40){
      if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
      label_p->DrawLatex(xPos, yPos, labelStr.c_str());
      yPos -= 0.085;
      labelStr = "";
    }
    labelStr = labelStr + preLabels[pI] + "; ";
  }
  if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
  
  label_p->DrawLatex(xPos, yPos, labelStr.c_str());
  
  
  canv_p->cd();
  pads_p[1]->cd();

  line_p->DrawLine(histMean_p[0]->GetBinLowEdge(1), 1.0, histMean_p[0]->GetBinLowEdge(histMean_p[0]->GetXaxis()->GetNbins()+1), 1.0);
  
  canv_p->cd();
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::string saveName = "pdfDir/" + dateStr + "/meanAndSigma" + overMeanStr + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete csnFit_p;
  
  delete leg_p;
  delete pads_p[0];
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete pads_p[1];
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete canv_p;
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  delete line_p;
  delete label_p;
  
  return;
}

int gdjResponsePlotter(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, "root")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* config_p = (TEnv*)inFile_p->Get("config");
  TEnv* label_p = (TEnv*)inFile_p->Get("label");
  configParser labels(label_p);
  configParser configs(config_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();
  std::map<std::string, std::string> configMap = configs.GetConfigMap();
  
  std::vector<std::string> necessaryParams = {"ISMC",
					      "NJTPTBINS",
					      "JTPTBINSDOLOG",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "NGAMMAPTBINSSUB",
					      "ISPP",
					      "RECOJTPTMIN"};
  std::vector<std::string> pbpbParams = {"CENTBINS"};
  
  if(!configs.ContainsParamSet(necessaryParams)) return 1;

  const bool isMC = std::stoi(configs.GetConfigVal("ISMC"));
  if(!isMC){
    std::cout << "gdjResponsePlotter ERROR - This macro only takes MC input, but config gives \'isMC\'=" << isMC << ". return 1" << std::endl;
    return 1;
  }

  const int nMaxJtPtBins = 200;
  const int nJtPtBins = std::stoi(configs.GetConfigVal("NJTPTBINS"));
  if(nJtPtBins > nMaxJtPtBins){
    std::cout << "Given nJtPtBins \'" << nJtPtBins << "\' is greater than max value \'" << nMaxJtPtBins << "\'. return 1" << std::endl;
    return 1;
  }
  const bool jtPtBinsDoLog = std::stoi(configs.GetConfigVal("JTPTBINSDOLOG"));
  const Float_t jtPtBinsLow = std::stof(configs.GetConfigVal("JTPTBINSLOW"));
  const Float_t jtPtBinsHigh = std::stof(configs.GetConfigVal("JTPTBINSHIGH"));
  Double_t jtPtBins[nMaxJtPtBins+1];
  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  
  const int nGammaPtBinsSub = std::stoi(configs.GetConfigVal("NGAMMAPTBINSSUB"));
  const bool isPP = std::stoi(configs.GetConfigVal("ISPP"));

  const Float_t recoJtPtMin = std::stof(configs.GetConfigVal("RECOJTPTMIN"));
  
  const int nMaxCentBins = 20;
  TH1F* recoOverGenVCent_FitMean_p[nMaxCentBins];
  TH1F* recoOverGenVCent_FitWidth_p[nMaxCentBins];
  TH1F* recoOverGenVCent_FitWidthOverMean_p[nMaxCentBins];
  
  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!configs.ContainsParamSet(pbpbParams)) return 1;    
    
    centBins = strToVect(configs.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;
    if(nCentBins > nMaxCentBins){
      std::cout << "Given nCentBins \'" << nCentBins << "\' is greater than max value \'" << nMaxCentBins << "\'. return 1" << std::endl;
      return 1;
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];

    recoOverGenVCent_FitMean_p[cI] = new TH1F(("recoOverGenFitMeanVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#LT Reco./Gen #GT;", nJtPtBins, jtPtBins);
    recoOverGenVCent_FitWidth_p[cI] = new TH1F(("recoOverGenFitWidthVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#sigma(Reco./Gen);", nJtPtBins, jtPtBins);
    recoOverGenVCent_FitWidthOverMean_p[cI] = new TH1F(("recoOverGenFitWidthOverMeanVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#sigma(Reco./Gen/)/#LT Reco./Gen #GT;", nJtPtBins, jtPtBins);
  }

  TF1* fit_p = new TF1("fit_p", "gaus", 0, 3.0);
  
  int nXVal = 0;
  int nYVal = 0;
  getNXNYPanels(nJtPtBins, &nXVal, &nYVal);
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TCanvas* canv_p = new TCanvas("canv_p", "", nXVal*450, nYVal*450);
    canv_p->SetTopMargin(0.001);
    canv_p->SetLeftMargin(0.001);
    canv_p->SetRightMargin(0.001);
    canv_p->SetBottomMargin(0.001);

    canv_p->Divide(nXVal, nYVal);

    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];

    for(Int_t jI = 0; jI < nJtPtBins; ++jI){
      canv_p->cd();
      canv_p->cd(jI+1);

      gPad->SetTopMargin(0.001);
      gPad->SetRightMargin(0.001);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12);
      
      std::string jtPtStr = "JtPt";
      if(jI < 10) jtPtStr = jtPtStr + "0";
      jtPtStr = jtPtStr + std::to_string(jI);
	 
      std::string histStr = centStr + "/photonJtRecoOverGenVCentJtPt_" + centStr + "_" + jtPtStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h";
      TH1F* hist_p = (TH1F*)inFile_p->Get(histStr.c_str());

      hist_p->SetMarkerSize(1);
      hist_p->SetMarkerColor(1);
      hist_p->SetLineColor(1);
      hist_p->SetMarkerStyle(24);
      
      hist_p->DrawCopy("HIST E1 P");
      gStyle->SetOptStat(0);

      Double_t minVal = TMath::Max((Double_t)hist_p->GetBinLowEdge(1), (Double_t)(recoJtPtMin/jtPtBinsLow));
      hist_p->Fit("fit_p", "M E Q N", "", minVal, hist_p->GetBinLowEdge(hist_p->GetXaxis()->GetNbins()+1));

      Double_t mean1 = fit_p->GetParameter(1);
      Double_t sigma1 = fit_p->GetParameter(2);
      
      minVal = TMath::Max(mean1 - sigma1*1.5, (Double_t)(recoJtPtMin/jtPtBinsLow));
      Double_t maxVal = TMath::Min(mean1 + sigma1*1.5, (Double_t)hist_p->GetBinLowEdge(hist_p->GetXaxis()->GetNbins()+1));
      
      hist_p->Fit("fit_p", "M E Q N", "", minVal, maxVal);      

      //      fit_p->SetMarkerSize(0.5);
      //      fit_p->SetLineWidth(0.5);
      fit_p->SetLineStyle(2);
      fit_p->DrawCopy("SAME");
      
      Double_t mean2 = fit_p->GetParameter(1);
      Double_t sigma2 = fit_p->GetParameter(2);
      Double_t meanErr2 = fit_p->GetParError(1);
      Double_t sigmaErr2 = fit_p->GetParError(2);
      
      recoOverGenVCent_FitMean_p[cI]->SetBinContent(jI+1, mean2);
      recoOverGenVCent_FitMean_p[cI]->SetBinError(jI+1, meanErr2);

      recoOverGenVCent_FitWidth_p[cI]->SetBinContent(jI+1, sigma2);
      recoOverGenVCent_FitWidth_p[cI]->SetBinError(jI+1, sigmaErr2);

      Double_t val = sigma2/mean2;
      Double_t err = TMath::Sqrt(meanErr2*meanErr2/(mean2*mean2) + sigmaErr2*sigmaErr2/(sigma2*sigma2));

      recoOverGenVCent_FitWidthOverMean_p[cI]->SetBinContent(jI+1, val);
      recoOverGenVCent_FitWidthOverMean_p[cI]->SetBinError(jI+1, val*err);
    }    
    
    std::string saveName = "pdfDir/" + dateStr + "/photonJtRecoOverGen_Cent" + std::to_string(cI) + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  }
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  plotMeanAndSigma(doGlobalDebug, &labelMap, dateStr, centBins, recoOverGenVCent_FitMean_p, recoOverGenVCent_FitWidth_p);
  plotMeanAndSigma(doGlobalDebug, &labelMap, dateStr, centBins, recoOverGenVCent_FitMean_p, recoOverGenVCent_FitWidthOverMean_p);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "GDJRESPONSEPLOTTER COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjResponsePlotter.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjResponsePlotter(argv[1]);
  return retVal;
}
