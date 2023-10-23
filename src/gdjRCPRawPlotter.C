//Author: Chris McGinn (2020.05.06)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

//Order is data, MC-Reco, MC-Truth
void plotRCP(std::string saveName, std::vector<TH1F*> hists_p, std::vector<std::string> legLabels, std::vector<std::string> texLabels, std::vector<bool> permaTex, TEnv* config_p, std::string envStr, bool doRenorm, bool doGlobalDebug)
{
  bool doAllLeg = config_p->GetValue("DOALLLEG", 1);
  bool goodLeg = false;
  for(Int_t lI = 0; lI < 10; ++lI){
    std::string label = config_p->GetValue(("PERMALEG." + std::to_string(lI)).c_str(), "");
    if(label.size() == 0) continue;
    if(saveName.find(label) != std::string::npos){
      goodLeg = true;
      break;
    }
  } 
  
  const Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.18;

  const Double_t max = config_p->GetValue((envStr + "MAX").c_str(), 0.5);
  const Double_t min = config_p->GetValue((envStr + "MIN").c_str(), 0.5);
  const Double_t ratMax = config_p->GetValue((envStr + "RATMAX").c_str(), 0.5);
  const Double_t ratMin = config_p->GetValue((envStr + "RATMIN").c_str(), 0.5);

  const Double_t xLegPos = config_p->GetValue((envStr + "LEGX").c_str(), 0.7);
  Double_t yLegPos = config_p->GetValue((envStr + "LEGY").c_str(), 0.7);

  const Double_t xPos = config_p->GetValue((envStr + "LABELX").c_str(), 0.7);
  Double_t yPos = config_p->GetValue((envStr + "LABELY").c_str(), 0.7);

  const Bool_t alignRight = config_p->GetValue((envStr + "ALIGNRIGHT").c_str(), 0);

  const Bool_t doLogX = config_p->GetValue((envStr + "LOGX").c_str(), 0);
  const Bool_t doLogY = config_p->GetValue((envStr + "LOGY").c_str(), 0);
   
  std::vector<std::vector<float> > generalBoxes;
  for(int gI = 0; gI < 10; ++gI){
    std::string tempStr = config_p->GetValue((envStr + "GENBOX." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    generalBoxes.push_back(strToVectF(tempStr));
  }

  TCanvas* canv_p = new TCanvas("canv_p", "", 450*2, 450*2);
  canv_p->SetTopMargin(0.001);
  canv_p->SetBottomMargin(0.001);
  canv_p->SetLeftMargin(0.001);
  canv_p->SetRightMargin(0.001);

  TPad* pads_p[2];
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(0.01);
  pads_p[0]->SetRightMargin(0.02);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(0.001);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetTopMargin(0.001);
  pads_p[1]->SetRightMargin(0.02);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetBottomMargin(leftMargin*0.8/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  pads_p[0]->cd();
  
  Double_t titleSizeX = 0.04;//histData_p->GetXaxis()->GetTitleSize();
  Double_t titleSizeY = 0.04;//histData_p->GetYaxis()->GetTitleSize();
  Double_t labelSizeX = 0.035;//histData_p->GetXaxis()->GetLabelSize();
  Double_t labelSizeY = 0.035;//histData_p->GetYaxis()->GetLabelSize();  

  Int_t titleFont = hists_p[0]->GetXaxis()->GetTitleFont();
  
  hists_p[0]->GetXaxis()->SetLabelSize(0.0001);
  hists_p[0]->GetYaxis()->SetLabelSize(labelSizeY/(1.0 - padSplit));
  hists_p[0]->GetXaxis()->SetTitleSize(0.0001);
  hists_p[0]->GetYaxis()->SetTitleSize(titleSizeY/(1.0 - padSplit));
  
  gStyle->SetOptStat(0);

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  int nHistsInLeg = hists_p.size();
  
  TLegend* leg_p = new TLegend(xLegPos, yLegPos - 0.07*(double)nHistsInLeg, xLegPos + 0.25, yLegPos);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSizeX/(1.0 - padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX/(1.0 - padSplit));
  label_p->SetNDC();

  if(doRenorm){
    for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
      hists_p[hI]->Scale(1./hists_p[hI]->Integral());
      hists_p[hI]->GetYaxis()->SetTitle("Unity Norm.");
    }
  }
  
  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    HIJet::Style::EquipHistogram(hists_p[hI], hI);
  
    hists_p[hI]->GetYaxis()->SetNdivisions(505);
    hists_p[hI]->GetXaxis()->SetNdivisions(505);
    
    hists_p[hI]->SetMaximum(max);
    hists_p[hI]->SetMinimum(min);
    
    hists_p[hI]->GetYaxis()->SetTitleOffset(1.0);
    
    if(hI == 0) hists_p[hI]->DrawCopy("HIST E1 P");
    else hists_p[hI]->DrawCopy("HIST E1 P SAME");
    gPad->SetTicks();
    leg_p->AddEntry(hists_p[hI], legLabels[hI].c_str(), "P L");
  }
  
  if(doLogX) gPad->SetLogx();
  if(doLogY) gPad->SetLogy();
  
  if(doAllLeg || goodLeg) leg_p->Draw("SAME");

  if(alignRight){
    /*
    double maxSize = 0.0;
    for(unsigned int pI = 0; pI < texLabels.size(); ++pI){    
      label_p->SetText(xPos, yPos, texLabels[pI].c_str());
      if(label_p->GetXsize() > maxSize) maxSize = label_p->GetXsize();
    }

    for(unsigned int pI = 0; pI < texLabels.size(); ++pI){    
      label_p->SetText(xPos, yPos, texLabels[pI].c_str());
      while(label_p->GetXsize() < maxSize){
	texLabels[pI] = " " + texLabels[pI];
	label_p->SetText(xPos, yPos, texLabels[pI].c_str());
      }
    }
    */
    label_p->SetTextAlign(31);
  }

  if(!doAllLeg && !goodLeg){
    for(unsigned int gI = 0; gI < texLabels.size(); ++gI){
      if(permaTex[gI]) continue;
      texLabels[gI] = "";
    }
  }

  //  const double xPos = 0.2;
  //  double yPos = 0.94;
  std::string labelName = "";
  for(unsigned int pI = 0; pI < texLabels.size(); ++pI){    
    /*
    if(labelName.size() + texLabels[pI].size() > 1){//1 to force each one to its own line
      if(labelName.find(";") != std::string::npos) labelName.replace(labelName.rfind(";"), labelName.size(), "");
      label_p->DrawLatex(xPos, yPos, labelName.c_str());
      yPos -= 0.065;
      labelName = "";
    }
    labelName = labelName + texLabels[pI] + "; ";
    */

    label_p->DrawLatex(xPos, yPos, texLabels[pI].c_str());
    yPos -= 0.09;
  }
  
  /*
  if(labelName.find(";") != std::string::npos) labelName.replace(labelName.rfind(";"), labelName.size(), "");
  label_p->DrawLatex(xPos, yPos, labelName.c_str());
  */
  
  if(min < 0.0) line_p->DrawLine(hists_p[0]->GetXaxis()->GetBinLowEdge(1), 0.0, hists_p[0]->GetXaxis()->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);

  canv_p->cd();
  pads_p[1]->cd();

  TH1D* histDataClone_p = (TH1D*)hists_p[0]->Clone("histDataClone_h");
  
  histDataClone_p->Sumw2();
  histDataClone_p->Divide(hists_p[1]);

  histDataClone_p->SetMaximum(ratMax);
  histDataClone_p->SetMinimum(ratMin);

  histDataClone_p->GetYaxis()->SetNdivisions(505);
  histDataClone_p->GetXaxis()->SetNdivisions(505);
  
  histDataClone_p->GetXaxis()->SetTitleSize(titleSizeX/padSplit);
  histDataClone_p->GetYaxis()->SetTitleSize(titleSizeY/padSplit);
  histDataClone_p->GetXaxis()->SetLabelSize(labelSizeX/padSplit);
  histDataClone_p->GetYaxis()->SetLabelSize(labelSizeY/padSplit);
  
  histDataClone_p->GetYaxis()->SetTitleOffset(histDataClone_p->GetYaxis()->GetTitleOffset()*padSplit/(1.0 - padSplit));

  histDataClone_p->GetYaxis()->SetTitle("R_{CP}");
  histDataClone_p->DrawCopy("HIST E1 P");

  gPad->SetTicks();  
  if(doLogX) gPad->SetLogx();

  if(ratMax > 1.0 && ratMin < 1.0){
    line_p->DrawLine(histDataClone_p->GetXaxis()->GetBinLowEdge(1), 1.0, histDataClone_p->GetXaxis()->GetBinLowEdge(histDataClone_p->GetXaxis()->GetNbins()+1), 1.0);
  }

  for(unsigned int gI = 0; gI < generalBoxes.size(); ++gI){
    canv_p->cd();
    drawWhiteBoxNDC(canv_p, generalBoxes[gI][0], generalBoxes[gI][1], generalBoxes[gI][2], generalBoxes[gI][3]);
  }

  quietSaveAs(canv_p, saveName.c_str());
  delete canv_p;

  delete histDataClone_p;

  delete line_p;
  delete leg_p;
  delete label_p;
  
  return;
}

int gdjRCPRawPlotter(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

 
  TEnv* plotConfig_p = new TEnv(inConfigFileName.c_str());
  
  configParser config(plotConfig_p);
  std::vector<std::string> necessaryParams = {"INPBPBFILENAME",
					      "PBPBLABEL",
					      "SAVETAG",
					      "SAVEEXT",
					      "DOGAMMAEFFCORR",					      
					      "GAMMAPTMIN",
					      "GAMMAPTMAX",
					      "GAMMAPTRATMIN",
					      "GAMMAPTRATMAX",
					      "JTDPHIMIN",
					      "JTDPHIMAX",
					      "JTDPHIRATMIN",
					      "JTDPHIRATMAX",
					      "JTETAMAX",
					      "JTETAMIN",
					      "JTETARATMAX",
					      "JTETARATMIN",
					      "JTPTMIN",
					      "JTPTMAX",
					      "JTPTRATMIN",
					      "JTPTRATMAX",
					      "JTXJMIN",
					      "JTXJMAX",
					      "JTXJRATMIN",
					      "JTXJRATMAX",
					      "JTMULTMIN",
					      "JTMULTMAX",
					      "JTMULTRATMIN",
					      "JTMULTRATMAX",
  					      "MULTIJTPTMIN",
					      "MULTIJTPTMAX",
					      "MULTIJTPTRATMIN",
					      "MULTIJTPTRATMAX",
					      "MULTIJTXJMIN",
					      "MULTIJTXJMAX",
					      "MULTIJTXJRATMIN",
					      "MULTIJTXJRATMAX",
  					      "MULTIJTXJJMIN",
					      "MULTIJTXJJMAX",
					      "MULTIJTXJJRATMIN",
					      "MULTIJTXJJRATMAX",
  					      "MULTIJTDPHIJJMIN",
					      "MULTIJTDPHIJJMAX",
					      "MULTIJTDPHIJJRATMIN",
					      "MULTIJTDPHIJJRATMAX",
					      "DOALLLEG"};

  if(!config.ContainsParamSet(necessaryParams)) return 1;

  const std::string inPbPbFileName = config.GetConfigVal("INPBPBFILENAME");

  if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;

  std::vector<std::string> globalLabels;
  for(Int_t gI = 0; gI < 10; ++gI){
    std::string plotStr = plotConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(plotStr.size() == 0) break;
    globalLabels.push_back(plotStr);
  }

  bool hasJetLabel = false;
  for(auto const & label : globalLabels){
    if(label.find("R=") != std::string::npos || label.find("R =") != std::string::npos){
      hasJetLabel = true;
      break;
    }
    else if(label.find("anti-k") != std::string::npos){
      hasJetLabel = true;
      break;
    }
  }

  std::string pbpbLabel = plotConfig_p->GetValue("PBPBLABEL", "");
  std::string ppLabel = plotConfig_p->GetValue("PPLABEL", "");

  const std::string saveTag = plotConfig_p->GetValue("SAVETAG", "");
  const std::string saveExt = plotConfig_p->GetValue("SAVEEXT", "pdf");
  std::vector<std::string> validExts = {"pdf", "png"};
  if(!vectContainsStr(saveExt, &validExts)){
    std::cout << "Given SAVEEXT \'" << saveExt << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  const bool doGammaEffCorr = (bool)plotConfig_p->GetValue("DOGAMMAEFFCORR", 0);
  std::vector<std::string> reqGammaEffParams = {"PBPBGAMMAEFFCORRFILE"};
  
  std::string pbpbGammaEffCorrFileName = "";
  std::string ppGammaEffCorrFileName = "";
  TFile* pbpbGammaEffCorrFile_p = nullptr;
  std::vector<TH1F*> pbpbEffCorr_p;
  std::vector<Int_t> centBinsGammaEffVect;
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(doGammaEffCorr){
    if(!config.ContainsParamSet(reqGammaEffParams)) return 1;

    pbpbGammaEffCorrFileName = plotConfig_p->GetValue("PBPBGAMMAEFFCORRFILE", "");

    if(!check.checkFileExt(pbpbGammaEffCorrFileName, ".root")) return 1;

    pbpbGammaEffCorrFile_p = new TFile(pbpbGammaEffCorrFileName.c_str(), "READ");
    
    TEnv* tempEffConfig_p = (TEnv*)pbpbGammaEffCorrFile_p->Get("config");
    std::string centBinsStr = tempEffConfig_p->GetValue("CENTBINS", "");
    centBinsGammaEffVect = strToVectI(centBinsStr);

    for(unsigned int cI = 0; cI < centBinsGammaEffVect.size()-1; ++cI){
      pbpbEffCorr_p.push_back(nullptr);
      std::string centStr = "Cent" + std::to_string(centBinsGammaEffVect[cI]) + "to" + std::to_string(centBinsGammaEffVect[cI+1]);
      pbpbEffCorr_p[cI] = (TH1F*)pbpbGammaEffCorrFile_p->Get(("photonEff_TOT_CentDep_Eff_" + centStr + "_Eta0p00to2p37_h").c_str());
    }
  }  
  
  TFile* inPbPbFile_p = new TFile(inPbPbFileName.c_str(), "READ");
  TEnv* pbpbEnv_p = (TEnv*)inPbPbFile_p->Get("config");
  TEnv* pbpbLabels_p = (TEnv*)inPbPbFile_p->Get("label");
  configParser configPbPb(pbpbEnv_p);
  configParser labelPbPb(pbpbLabels_p);

  std::vector<std::string> paramsToCheck = {"ISMC",
					    "NGAMMAPTBINS",
					    "GAMMAPTBINSLOW",
					    "GAMMAPTBINSHIGH",
					    "GAMMAPTBINSDOLOG",
					    "GAMMAPTBINSDOCUSTOM",
					    "NJTPTBINS",
					    "JTPTBINSLOW",
					    "JTPTBINSHIGH",
					    "JTPTBINSDOLOG",
					    "JTPTBINSDOCUSTOM",
					    "JETR"};

  if(!configPbPb.ContainsParamSet(paramsToCheck)) return 1;

  const bool isPP = pbpbEnv_p->GetValue("ISPP", 0);
  plotConfig_p->SetValue("ISPP", isPP);
  if(isPP){
    std::cout << "RCP isPP: " << isPP << " is not valid. return 1" << std::endl;
    return 1;
  }
  
  const int jetR = std::stoi(configPbPb.GetConfigVal("JETR"));
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);

  if(!hasJetLabel) globalLabels.push_back("anti-k_{t} #it{R}=" + jetRStr + " jets");

  const std::string gammaExclusionDRStr = pbpbEnv_p->GetValue("GAMMAEXCLUSIONDR", "");
  const std::string gammaJtDPhiCutStr = pbpbEnv_p->GetValue("GAMMAJTDPHI", "");
  
  std::string gammaDRStr = "";
  if(gammaExclusionDRStr.size() != 0) gammaDRStr = "#DeltaR_{#gammaJ} > " + gammaExclusionDRStr;

  std::string gammaJtDPhiStr = "";
  if(gammaJtDPhiCutStr.size() != 0){
    gammaJtDPhiStr = "#Delta#phi_{#gammaJ} > " + gammaJtDPhiCutStr;
    if(gammaJtDPhiStr.find("pi") != std::string::npos){
      gammaJtDPhiStr.replace(gammaJtDPhiStr.find("pi"), 2, "#pi");
    }
  }

  
  const Bool_t doMix = std::stoi(configPbPb.GetConfigVal("DOMIX"));
  if(!doMix && !isPP){
    std::cout << "GDJRCPRAWPLOTTER - DOMIX must be true if Pb+Pb. return 1" << std::endl;
    return 1;
  }

  const Int_t nMaxBins = 100;
  const std::string centBinsStr = configPbPb.GetConfigVal("CENTBINS");
  std::vector<int> centBins = strToVectI(configPbPb.GetConfigVal("CENTBINS"));
  const Int_t nCentBins = centBins.size()-1;
  const Int_t nGammaEtaBinsSub = std::stoi(configPbPb.GetConfigVal("NGAMMAETABINSSUB"));
  const Int_t nGammaPtBinsSub = std::stoi(configPbPb.GetConfigVal("NGAMMAPTBINSSUB"))+1;

  //CFM Edit 2021.11.29 - for an RCP we have to identify the peripheral most bin
  int peripheralPos = nCentBins-1;
  int peripheralLower = centBins[peripheralPos];
  int peripheralUpper = centBins[peripheralPos+1];

  std::string centPeripheralStr = "Cent" + std::to_string(centBins[peripheralPos]) + "to" + std::to_string(centBins[peripheralPos+1]);
  std::cout << "PERIPHERALPOSITIONS: " << peripheralPos << ", " << peripheralLower << "-" << peripheralUpper << std::endl;  
  
  //CFM Edit 2021.11.12 - need to do a massive update on this macro to run on latest ntuples, particular the mixMachine MIXMODE0, MIXMODE1, MIXMODE2 output histograms

  std::string barrelECStr = "BarrelAndEC"; //So we can later put it in an nBarrelAndECLoop
  std::vector<std::string> observables = {"JtPt"};
  std::vector<std::string> observablesXTitle = {"Jet p_{T} [GeV]"};  
  std::vector<std::string> observablesYTitle = {"p_{T}^{Jet}"};

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(cI == peripheralPos) break;
    std::string centStr = "Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]);
    
    for(unsigned int oI = 0; oI < observables.size(); ++oI){      
      bool isMultijet = isStrSame(observables[oI], "JtXJJ") || isStrSame(observables[oI], "JtAJJ") || isStrSame(observables[oI], "JtDPhiJJG") || isStrSame(observables[oI], "JtDPhiJJ") || isStrSame(observables[oI], "JtDRJJ");
      
      std::string mixModeStr = "MIXMODE1";
      if(isMultijet) mixModeStr = "MIXMODE2";      

      std::string titleStr = ";" + observablesXTitle[oI] + ";#frac{1}{N_{#gamma}}#frac{dN_{Jet}}{d" + observablesYTitle[oI] + "}";
      std::string upperObsStr = strLowerToUpper(observables[oI]);

      TH1F* histPhoPt_PbPbPeripheral_p = (TH1F*)inPbPbFile_p->Get((centPeripheralStr + "/photonPtVCent_" + centPeripheralStr + "_" + barrelECStr + "_PURCORR_h").c_str());            
      TH2F* hist_PbPbPeripheral_p = (TH2F*)inPbPbFile_p->Get((centPeripheralStr + "/photonPt" + observables[oI] + "VCent_" + centPeripheralStr + "_" + barrelECStr + "_DPhi0_PURCORR_h" ).c_str());

      TH1F* histPhoPt_PbPb_p = (TH1F*)inPbPbFile_p->Get((centStr + "/photonPtVCent_" + centStr + "_" + barrelECStr + "_PURCORR_h").c_str());            
      TH2F* hist_PbPb_p = (TH2F*)inPbPbFile_p->Get((centStr + "/photonPt" + observables[oI] + "VCent_" + centStr + "_" + barrelECStr + "_DPhi0_PURCORR_h" ).c_str());
      
      //Declare an array and a maxbins value for rebinning purposes
      const Int_t nMaxBins = 1000;
      Int_t nBinsX=-1;
      Double_t newBinsX[nMaxBins+1];
      for(Int_t bIX = 0; bIX < hist_PbPb_p->GetXaxis()->GetNbins()+1; ++bIX){
	newBinsX[bIX] = hist_PbPb_p->GetXaxis()->GetBinLowEdge(bIX+1);
	++nBinsX;
      }
      
      //We need to process these along the Y-Axis
      for(Int_t bIY = 0; bIY < hist_PbPb_p->GetYaxis()->GetNbins()+1; ++bIY){	
	TH1F* hist1D_PbPbPeripheral_p = new TH1F("hist1D_PbPbPeripheral_p", titleStr.c_str(), nBinsX, newBinsX);
	TH1F* hist1D_PbPb_p = new TH1F("hist1D_PbPb_p", titleStr.c_str(), nBinsX, newBinsX);
	/*
 	std::vector<TH1F*> histsToConstruct = {hist1D_PP_p, hist1D_PPTRUTH_p, hist1D_PbPb_p};
 	std::vector<TH2F*> histsIn = {hist_PP_p, hist_PPTRUTH_p, hist_PbPb_p};
	std::vector<TH1F*> gammaHists = {histPhoPt_PP_p, histPhoPt_PPTRUTH_p, histPhoPt_PbPb_p};
	std::vector<Float_t> gammaNormVal = {0.0, 0.0, 0.0};
	*/
	std::vector<TH1F*> histsToConstruct = {hist1D_PbPbPeripheral_p, hist1D_PbPb_p};
 	std::vector<TH2F*> histsIn = {hist_PbPbPeripheral_p, hist_PbPb_p};
	std::vector<TH1F*> gammaHists = {histPhoPt_PbPbPeripheral_p, histPhoPt_PbPb_p};
	std::vector<Float_t> gammaNormVal = {0.0, 0.0};
      	
	if(bIY == hist_PbPbPeripheral_p->GetYaxis()->GetNbins()){
	  //We need to add an efficiency correction!
	  //Match your centrality bins
	  std::vector<TH1F*> effCorr;
	
	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, bIY: " << __FILE__ << ", " << __LINE__ << ", " << bIY << std::endl;
	  if(doGammaEffCorr){
	    Int_t centEffPos = -1;
	    Int_t centPeripheralEffPos = -1;
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    for(unsigned int cI2 = 0; cI2 < centBinsGammaEffVect.size()-1; ++cI2){
	      std::string centEffStr = "Cent" + std::to_string(centBinsGammaEffVect[cI2]) + "to" +  std::to_string(centBinsGammaEffVect[cI2+1]);
	      
	      if(isStrSame(centEffStr, centStr)){
		centEffPos = cI2;
		//		break;
	      }

	      if(isStrSame(centEffStr, centPeripheralStr)){
		centPeripheralEffPos = cI2;
		//		break;
	      }
	    }
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    if(centEffPos < 0){
	      std::cout << "Skipping all gamma-pt for centrality: " << centStr << ", no matching efficiency bin" << std::endl;
	    
	      delete hist1D_PbPbPeripheral_p;
	      delete hist1D_PbPb_p;
	      
	      continue;
	    }

	    if(centPeripheralEffPos < 0){
	      std::cout << "Skipping all gamma-pt for peripheral centrality in RCP: " << centPeripheralStr << ", no matching efficiency bin" << std::endl;
	    
	      delete hist1D_PbPbPeripheral_p;
	      delete hist1D_PbPb_p;
	      
	      continue;
	    }
	    
	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    effCorr = {pbpbEffCorr_p[centPeripheralEffPos], pbpbEffCorr_p[centEffPos]};
	  }
	  
	  Float_t lowPhoPtVal = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(1);
	  Float_t highPhoPtVal = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(hist_PbPbPeripheral_p->GetYaxis()->GetNbins()+1);	  

	    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  for(unsigned int hI = 0; hI < histsIn.size(); ++hI){	  
	    //Get our normalizing values of nGamma
	    for(Int_t bIX = 0; bIX < gammaHists[hI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = gammaHists[hI]->GetXaxis()->GetBinCenter(bIX+1);

	      if(center < lowPhoPtVal) continue;
	      if(center >= highPhoPtVal) continue;	      

	      gammaNormVal[hI] += gammaHists[hI]->GetBinContent(bIX+1);
	    }
	    
	    //Construct our histograms
	    for(Int_t bIX = 0; bIX < hist_PbPbPeripheral_p->GetXaxis()->GetNbins(); ++bIX){
	      for(Int_t bIY2 = 0; bIY2 < hist_PbPbPeripheral_p->GetYaxis()->GetNbins(); ++bIY2){	
		Float_t corrFactor = 1.0;
		if(doGammaEffCorr){
		  //we need to find the bin for the correction factor
		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		  Int_t binPos = -1;

		  Float_t binEdgeLow1 = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(bIY2+1);
		  Float_t binEdgeHigh1 = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(bIY2+2);

		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, bIY2: " << __FILE__ << ", " << __LINE__ << ", " << bIY2 << std::endl;
		  
		  for(Int_t bIY3 = 0; bIY3 < effCorr[hI]->GetXaxis()->GetNbins(); ++bIY3){
		    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, bIY3: " << __FILE__ << ", " << __LINE__ << ", " << bIY3 << std::endl;

		    Float_t binEdgeLow2 = effCorr[hI]->GetXaxis()->GetBinLowEdge(bIY3+1);
		    Float_t binEdgeHigh2 = effCorr[hI]->GetXaxis()->GetBinLowEdge(bIY3+2);

		    if(TMath::Abs(binEdgeLow1 - binEdgeLow2) > 0.0001) continue;
		    if(TMath::Abs(binEdgeHigh1 - binEdgeHigh2) > 0.0001) continue;

		    binPos = bIY3;
		    break;
		  }

		  if(binPos < 0){
		    std::cout << "WARNING WE HAVE NO FOUND BINPOS FOR BIN " << binEdgeLow1 << "-" << binEdgeHigh1 << std::endl;
		    std::cout << "Possible matches are: " << std::endl;
		    for(Int_t bIY3 = 0; bIY3 < effCorr[hI]->GetXaxis()->GetNbins(); ++bIY3){
		      Float_t binEdgeLow2 = effCorr[hI]->GetXaxis()->GetBinLowEdge(bIY3+1);
		      Float_t binEdgeHigh2 = effCorr[hI]->GetXaxis()->GetBinLowEdge(bIY3+2);

		      std::cout << " " << binEdgeLow2 << "-" << binEdgeHigh2 << std::endl;
		    }
		    effCorr[hI]->Print("ALL");
		  }
		  
		  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		  
		  corrFactor = 1.0/effCorr[hI]->GetBinContent(binPos+1);
		  //		  std::cout << "CORR FACTOR: " << corrFactor << std::endl;
		}
		
		Float_t newVal = histsIn[hI]->GetBinContent(bIX+1, bIY2+1)*corrFactor;
		Float_t newErr = histsIn[hI]->GetBinError(bIX+1, bIY2+1)*corrFactor;

		newVal = newVal + histsToConstruct[hI]->GetBinContent(bIX+1);
		newErr = TMath::Sqrt(newErr*newErr + histsToConstruct[hI]->GetBinError(bIX+1)*histsToConstruct[hI]->GetBinError(bIX+1));
		
		histsToConstruct[hI]->SetBinContent(bIX+1, newVal);
		histsToConstruct[hI]->SetBinError(bIX+1, newErr);
	      }
	    }

	    //And Normalize
	    histsToConstruct[hI]->Sumw2();
	    histsToConstruct[hI]->Scale(1.0/gammaNormVal[hI]);	    
	  }
	}
	else{
	  for(unsigned int hI = 0; hI < histsIn.size(); ++hI){
	    Float_t lowPhoPtVal = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(bIY+1);
	    Float_t highPhoPtVal = hist_PbPbPeripheral_p->GetYaxis()->GetBinLowEdge(bIY+2);	  
	  
	    //Get our normalizing values of nGamma
	    for(Int_t bIX = 0; bIX < gammaHists[hI]->GetXaxis()->GetNbins(); ++bIX){
	      Float_t center = gammaHists[hI]->GetXaxis()->GetBinCenter(bIX+1);

	      if(center < lowPhoPtVal) continue;
	      if(center >= highPhoPtVal) continue;	      

	      gammaNormVal[hI] += gammaHists[hI]->GetBinContent(bIX+1);
	    }
	    
	    //Construct our histograms
	    for(Int_t bIX = 0; bIX < hist_PbPbPeripheral_p->GetXaxis()->GetNbins(); ++bIX){
	      histsToConstruct[hI]->SetBinContent(bIX+1, histsIn[hI]->GetBinContent(bIX+1, bIY+1));
	      histsToConstruct[hI]->SetBinError(bIX+1, histsIn[hI]->GetBinError(bIX+1, bIY+1));
	    }

	    //Normalize our histograms
	    histsToConstruct[hI]->Sumw2();
	    histsToConstruct[hI]->Scale(1.0/gammaNormVal[hI]);	    
	  }
	}

	//Quickly Rescale for bin width
	for(unsigned int hI = 0; hI < histsToConstruct.size(); ++hI){
	  for(Int_t bIX = 0; bIX < histsToConstruct[hI]->GetXaxis()->GetNbins(); ++bIX){
	    Float_t newValue = histsToConstruct[hI]->GetBinContent(bIX+1);
	    Float_t newError = histsToConstruct[hI]->GetBinError(bIX+1);
	    Float_t binWidth = histsToConstruct[hI]->GetBinWidth(bIX+1);

	    histsToConstruct[hI]->SetBinContent(bIX+1, newValue/binWidth);
	    histsToConstruct[hI]->SetBinError(bIX+1, newError/binWidth);
	  }
	}
	
	std::vector<TH1F*> hists_p = {hist1D_PbPb_p, hist1D_PbPbPeripheral_p};
	std::vector<std::string> legLabels = {pbpbLabel + ", " + labelPbPb.GetConfigVal(centStr), pbpbLabel + ", " + labelPbPb.GetConfigVal(centPeripheralStr)};
	std::vector<std::string> tempGlobalLabels = globalLabels;
	std::vector<bool> permaTex;
	for(unsigned int gI = 0; gI < tempGlobalLabels.size(); ++gI){
	  permaTex.push_back(false);
	}
      
	tempGlobalLabels.push_back(labelPbPb.GetConfigVal("GammaPt" + std::to_string(bIY)));
	permaTex.push_back(true);

	tempGlobalLabels.push_back(labelPbPb.GetConfigVal("GlobalJtPt0"));
	permaTex.push_back(false);        

	if(gammaDRStr.size() != 0){
	  tempGlobalLabels.push_back(gammaDRStr);
	  permaTex.push_back(false);
	}          	

	if(gammaJtDPhiStr.size() != 0){
          tempGlobalLabels.push_back(gammaJtDPhiStr);
          permaTex.push_back(false);
        }
	
	std::string canvStr = "pdfDir/" + dateStr + "/photonPt" + observables[oI] + "VCent_" + centStr + "_GammaPt" + std::to_string(bIY) + "_R" + std::to_string(jetR) + "_RCP";
	if(saveTag.size() != 0) canvStr = canvStr + "_" + saveTag;
	canvStr = canvStr + "_" + dateStr + "." + saveExt;
	//Note that doRenorm is listed false
	plotRCP(canvStr, hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, upperObsStr, false, doGlobalDebug);       	
	
	delete hist1D_PbPbPeripheral_p;
	delete hist1D_PbPb_p;
      }
    }
  }

  inPbPbFile_p->Close();
  delete inPbPbFile_p;

  delete plotConfig_p;

  if(doGammaEffCorr){
    pbpbGammaEffCorrFile_p->Close();

    delete pbpbGammaEffCorrFile_p;
  }
  
  std::cout << "GDJRCPRAWPLOTTER COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjRCPRawPlotter.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjRCPRawPlotter(argv[1]);
  return retVal;
}
