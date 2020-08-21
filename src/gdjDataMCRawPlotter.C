//Author: Chris McGinn (2020.05.06)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
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
void plotDataMC(std::string saveName, std::vector<TH1F*> hists_p, std::vector<std::string> legLabels, std::vector<std::string> texLabels, std::vector<bool> permaTex, TEnv* config_p, std::string envStr, bool doRenorm, bool doGlobalDebug)
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

  /*
  const std::string isPPStr = config_p->GetValue("ISPP", "");
  if(isPPStr.size() != 1){
    std::cout << "ISPP VALUE \'" << isPPStr << "\' is invalid. return" << std::endl;
    return;
  }
  const bool isPP = std::stoi(isPPStr);
  */
  
  const Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.125;

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


  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
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
  pads_p[1]->SetBottomMargin(leftMargin/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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

  TLegend* leg_p = new TLegend(xLegPos, yLegPos - 0.07*(double)hists_p.size(), xLegPos + 0.25, yLegPos);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSizeX/(1.0 - padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
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
    
    hists_p[hI]->GetYaxis()->SetTitleOffset(0.7);
    
    if(hI == 0) hists_p[hI]->DrawCopy("HIST E1 P");
    else hists_p[hI]->DrawCopy("HIST E1 P SAME");

    leg_p->AddEntry(hists_p[hI], legLabels[hI].c_str(), "P L");
  }
  
  if(doLogX) gPad->SetLogx();
  if(doLogY) gPad->SetLogy();
    
  
  if(doAllLeg || goodLeg) leg_p->Draw("SAME");

  if(alignRight){
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
    yPos -= 0.0975;
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

  histDataClone_p->GetYaxis()->SetTitle("Data / MC (Reco.)");
  histDataClone_p->DrawCopy("HIST E1 P");

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

int gdjDataMCRawPlotter(std::string inConfigFileName)
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
  std::vector<std::string> necessaryParams = {"INDATAFILENAME",
					      "INMCFILENAME",
					      "DATALABEL",
					      "MCLABEL",
					      "MCTRUTHLABEL",
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!config.ContainsParamSet(necessaryParams)) return 1;

  const std::string inDataFileName = config.GetConfigVal("INDATAFILENAME");
  const std::string inMCFileName = config.GetConfigVal("INMCFILENAME");

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


  std::string dataLabel = plotConfig_p->GetValue("DATALABEL", "");
  std::string mcLabel = plotConfig_p->GetValue("MCLABEL", "");
  std::string mcTruthLabel = plotConfig_p->GetValue("MCTRUTHLABEL", "");

  TFile* inDataFile_p = new TFile(inDataFileName.c_str(), "READ");
  TEnv* dataEnv_p = (TEnv*)inDataFile_p->Get("config");
  TEnv* dataLabels_p = (TEnv*)inDataFile_p->Get("label");
  configParser configData(dataEnv_p);
  configParser labelData(dataLabels_p);

  TFile* inMCFile_p = new TFile(inMCFileName.c_str(), "READ");
  TEnv* mcEnv_p = (TEnv*)inMCFile_p->Get("config");
  configParser configMC(mcEnv_p);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> paramsToCheck = {"CENTBINS",
					    "JETR",
					    "DOMIX",
					    "GAMMAJTDPHI",
					    "GAMMAPTBINSDOLOG",
					    "GAMMAPTBINSHIGH",
					    "GAMMAPTBINSLOW",
					    "GAMMAPTBINSSUBDOLOG",
					    "GAMMAPTBINSSUBHIGH",
					    "GAMMAPTBINSSUBLOW",
					    "ISPP",
					    "JTPTBINSDOLOG",
					    "JTPTBINSHIGH",
					    "JTPTBINSLOW",
					    "NDPHIBINS",
					    "NGAMMAETABINSSUB",
					    "NGAMMAPTBINS",
					    "NGAMMAPTBINSSUB",
					    "NJTPTBINS",
					    "NPHIBINS",
					    "NXJBINS",
					    "RECOJTPTMIN",
					    "XJBINSHIGH",
					    "XJBINSLOW"};

  std::vector<std::string> mcParamsToCheck = {"GAMMAEXCLUSIONDR"};
  
  if(!configData.ContainsParamSet(paramsToCheck)) return 1;
  if(!configMC.ContainsParamSet(paramsToCheck)) return 1;
  if(!configMC.ContainsParamSet(mcParamsToCheck)) return 1;

  if(!configData.CheckConfigParams(&configMC, paramsToCheck)) return 1;

  plotConfig_p->SetValue(dataEnv_p->GetValue("ISPP", ""));
  
  const int jetR = std::stoi(configMC.GetConfigVal("JETR"));
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);

  if(!hasJetLabel) globalLabels.push_back("anti-k_{t} #it{R}=" + jetRStr + " jets");

  const std::string gammaExclusionDRStr = mcEnv_p->GetValue("GAMMAEXCLUSIONDR", "");

  std::string gammaDRStr = "";
  if(gammaExclusionDRStr.size() != 0) gammaDRStr = "#DeltaR_{#gammaJ} > " + gammaExclusionDRStr;
  
  const Bool_t doMix = std::stoi(configData.GetConfigVal("DOMIX"));
  const Bool_t isPP = std::stoi(configData.GetConfigVal("ISPP"));
  if(!doMix && !isPP){
    std::cout << "GDJDATAMCRAWPLOTTER - DOMIX must be true if Pb+Pb. return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxBins = 100;
  const std::string centBinsStr = configData.GetConfigVal("CENTBINS");
  std::vector<int> centBins = strToVectI(configData.GetConfigVal("CENTBINS"));
  const Int_t nCentBins = centBins.size()-1;
  const Int_t nGammaEtaBinsSub = std::stoi(configData.GetConfigVal("NGAMMAETABINSSUB"));

  const Int_t nGammaPtBinsSub = std::stoi(configData.GetConfigVal("NGAMMAPTBINSSUB"))+1;

  TH1F* photonPtVCentEta_MC_p[nMaxBins][nMaxBins+1];

  TH1F* photonSubJtDPhiVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtPtVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtPtVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtEtaVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtXJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtXJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtXJJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtDPhiJJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtMultModVCentPt_MC_p[nMaxBins][nMaxBins+1];

  TH1F* photonGenJtXJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonGenMultiJtXJJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonGenMultiJtDPhiJJVCentPt_MC_p[nMaxBins][nMaxBins+1];
  TH1F* photonGenJtDPhiVCentPt_MC_p[nMaxBins][nMaxBins+1];

  TH1F* photonPtVCentEta_Data_p[nMaxBins][nMaxBins+1];

  TH1F* photonSubJtDPhiVCentPt_Data_p[nMaxBins][nMaxBins+1];  
  TH1F* photonSubJtPtVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtPtVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtEtaVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtXJVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtXJVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtXJJVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubMultiJtDPhiJJVCentPt_Data_p[nMaxBins][nMaxBins+1];
  TH1F* photonSubJtMultModVCentPt_Data_p[nMaxBins][nMaxBins+1];

  /*
  const double xPos1 = 0.2;
  const double yPos1 = 0.94;
  const double xPos2 = 0.7;
  const double yPos2 = 0.82;

  const double xPos3 = 0.2;
  const double yPos3 = 0.54;
  */

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]);

    std::cout << (centStr + "/photonPtVCentEta_" + centStr + "_AbsEta" + std::to_string(nGammaEtaBinsSub) + "_h") << std::endl;

    photonPtVCentEta_MC_p[cI][0] = (TH1F*)inMCFile_p->Get((centStr + "/photonPtVCentEta_" + centStr + "_AbsEta" + std::to_string(nGammaEtaBinsSub) + "_h").c_str());

    photonPtVCentEta_MC_p[cI][0]->Sumw2();
    photonPtVCentEta_MC_p[cI][0]->Scale(1./photonPtVCentEta_MC_p[cI][0]->Integral());
    photonPtVCentEta_MC_p[cI][0]->GetYaxis()->SetTitle("Unity normalization");
    
    photonPtVCentEta_Data_p[cI][0] = (TH1F*)inDataFile_p->Get((centStr + "/photonPtVCentEta_" + centStr + "_AbsEta" + std::to_string(nGammaEtaBinsSub) + "_h").c_str());

    photonPtVCentEta_Data_p[cI][0]->Sumw2();
    photonPtVCentEta_Data_p[cI][0]->Scale(1./photonPtVCentEta_Data_p[cI][0]->Integral());
    photonPtVCentEta_Data_p[cI][0]->GetYaxis()->SetTitle("Unity normalization");
    
    for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      photonSubJtDPhiVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubJtDPhiVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_h").c_str());
      if(isPP) photonGenJtDPhiVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonGenJtDPhiVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_h").c_str());

      photonSubJtDPhiVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubJtDPhiVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_h").c_str());

      photonSubJtEtaVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubJtEtaVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());
      photonSubJtEtaVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubJtEtaVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());

      photonSubJtXJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());
      if(isPP) photonGenJtXJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonGenJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());
      photonSubJtXJVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());

      photonSubMultiJtXJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubMultiJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());
      photonSubMultiJtXJVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubMultiJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());

      photonSubJtMultModVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubJtMultModVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());
      photonSubJtMultModVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubJtMultModVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_h").c_str());
      
      photonSubJtPtVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_DPhi0_h").c_str());
      photonSubJtPtVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_DPhi0_h").c_str());

      photonSubMultiJtPtVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubMultiJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_DPhi0_MultiJt0_h").c_str());
      photonSubMultiJtPtVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubMultiJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_DPhi0_MultiJt0_h").c_str());

      photonSubMultiJtXJJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubMultiJtXJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());
      if(isPP) photonGenMultiJtXJJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonGenMultiJtXJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());

      photonSubMultiJtXJJVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubMultiJtXJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());

      photonSubMultiJtDPhiJJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonSubMultiJtDPhiJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());
      if(isPP) photonGenMultiJtDPhiJJVCentPt_MC_p[cI][pI] = (TH1F*)inMCFile_p->Get((centStr + "/photonGenMultiJtDPhiJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());

      photonSubMultiJtDPhiJJVCentPt_Data_p[cI][pI] = (TH1F*)inDataFile_p->Get((centStr + "/photonSubMultiJtDPhiJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_GlobalJtPt0_DPhi0_MultiJt0_h").c_str());

      std::vector<std::string> legLabels = {dataLabel, mcLabel, mcTruthLabel};
      std::vector<TH1F*> hists_p;
      
      std::vector<std::string> tempGlobalLabels = globalLabels;
      std::vector<bool> permaTex;
      for(unsigned int gI = 0; gI < tempGlobalLabels.size(); ++gI){
	permaTex.push_back(false);
      }

      if(!isPP){
	tempGlobalLabels.push_back(labelData.GetConfigVal(centStr));
	permaTex.push_back(true);
      }

      if(pI == 0){
	hists_p = {photonPtVCentEta_Data_p[cI][0], photonPtVCentEta_MC_p[cI][0]};
	plotDataMC("pdfDir/" + dateStr + "/photonPtVCentEta_" + centStr + "_AbsEta" + std::to_string(nGammaEtaBinsSub) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "GAMMAPT", false, doGlobalDebug);
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      tempGlobalLabels.push_back(labelData.GetConfigVal("GammaPt" + std::to_string(pI)));
      permaTex.push_back(true);

      tempGlobalLabels.push_back(labelData.GetConfigVal("GlobalJtPt0"));
      permaTex.push_back(false);

      if(gammaDRStr.size() != 0){
	tempGlobalLabels.push_back(gammaDRStr);
	permaTex.push_back(false);
      }
      
      if(isPP) hists_p = {photonSubJtDPhiVCentPt_Data_p[cI][pI], photonSubJtDPhiVCentPt_MC_p[cI][pI], photonGenJtDPhiVCentPt_MC_p[cI][pI]}; 
      else hists_p = {photonSubJtDPhiVCentPt_Data_p[cI][pI], photonSubJtDPhiVCentPt_MC_p[cI][pI]}; 
      plotDataMC("pdfDir/" + dateStr + "/photonSubJtDPhiVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTDPHI", false, doGlobalDebug);

      //      plotDataMC("pdfDir/" + dateStr + "/photonSubJtDPhiVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_UNITYNORM_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTDPHI", true, doGlobalDebug);


      if(gammaDRStr.size() != 0){
	tempGlobalLabels[tempGlobalLabels.size()-1] = labelData.GetConfigVal("DPhi0");
      }
      else{
	tempGlobalLabels.erase(tempGlobalLabels.begin() + tempGlobalLabels.size()-1);

	tempGlobalLabels.push_back(labelData.GetConfigVal("DPhi0"));
	//	permaTex.push_back(false);
      }
      
      hists_p = {photonSubJtEtaVCentPt_Data_p[cI][pI], photonSubJtEtaVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubJtEtaVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTETA", false, doGlobalDebug);

      photonSubJtMultModVCentPt_Data_p[cI][pI]->Sumw2();
      photonSubJtMultModVCentPt_MC_p[cI][pI]->Sumw2();
      
      photonSubJtMultModVCentPt_Data_p[cI][pI]->Scale(1./photonSubJtMultModVCentPt_Data_p[cI][pI]->Integral());
      photonSubJtMultModVCentPt_MC_p[cI][pI]->Scale(1./photonSubJtMultModVCentPt_MC_p[cI][pI]->Integral());

      photonSubJtMultModVCentPt_Data_p[cI][pI]->GetYaxis()->SetTitle("Unity Normalization");
      photonSubJtMultModVCentPt_MC_p[cI][pI]->GetYaxis()->SetTitle("Unity Normalization");
      
      hists_p = {photonSubJtMultModVCentPt_Data_p[cI][pI], photonSubJtMultModVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubJtMultModVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTMULT", false, doGlobalDebug);

      if(isPP) hists_p = {photonSubJtXJVCentPt_Data_p[cI][pI], photonSubJtXJVCentPt_MC_p[cI][pI], photonGenJtXJVCentPt_MC_p[cI][pI]};
      else hists_p = {photonSubJtXJVCentPt_Data_p[cI][pI], photonSubJtXJVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTXJ", false, doGlobalDebug);

      tempGlobalLabels.push_back(labelData.GetConfigVal("MultiJt0"));
      permaTex.push_back(false);

      hists_p = {photonSubMultiJtXJVCentPt_Data_p[cI][pI], photonSubMultiJtXJVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubMultiJtXJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "MULTIJTXJ", false, doGlobalDebug);
      
      tempGlobalLabels = globalLabels;
      for(unsigned int gI = 0; gI < tempGlobalLabels.size(); ++gI){
	permaTex.push_back(false);
      }
      if(!isPP){
	tempGlobalLabels.push_back(labelData.GetConfigVal(centStr));
	permaTex.push_back(true);
      }
      tempGlobalLabels.push_back(labelData.GetConfigVal("GammaPt" + std::to_string(pI)));
      permaTex.push_back(true);
      tempGlobalLabels.push_back(labelData.GetConfigVal("GlobalJtPt0"));
      permaTex.push_back(false);
      tempGlobalLabels.push_back(labelData.GetConfigVal("DPhi0"));
      permaTex.push_back(false);

      hists_p = {photonSubJtPtVCentPt_Data_p[cI][pI], photonSubJtPtVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "JTPT", false, doGlobalDebug);

      tempGlobalLabels.push_back(labelData.GetConfigVal("MultiJt0"));
      permaTex.push_back(false);

      hists_p = {photonSubMultiJtPtVCentPt_Data_p[cI][pI], photonSubMultiJtPtVCentPt_MC_p[cI][pI]};
      plotDataMC("pdfDir/" + dateStr + "/photonSubMultiJtPtVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "MULTIJTPT", false, doGlobalDebug);

      hists_p = {photonSubMultiJtXJJVCentPt_Data_p[cI][pI], photonSubMultiJtXJJVCentPt_MC_p[cI][pI]};
      if(isPP) hists_p.push_back(photonGenMultiJtXJJVCentPt_MC_p[cI][pI]);
      plotDataMC("pdfDir/" + dateStr + "/photonSubMultiJtXJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "MULTIJTXJJ", false, doGlobalDebug);

      
      
      hists_p = {photonSubMultiJtDPhiJJVCentPt_Data_p[cI][pI], photonSubMultiJtDPhiJJVCentPt_MC_p[cI][pI]};
      if(isPP) hists_p.push_back(photonGenMultiJtDPhiJJVCentPt_MC_p[cI][pI]);
      plotDataMC("pdfDir/" + dateStr + "/photonSubMultiJtDPhiJJVCentPt_" + centStr + "_GammaPt" + std::to_string(pI) + "_R" + std::to_string(jetR) + "_DataMC_" + dateStr + ".pdf", hists_p, legLabels, tempGlobalLabels, permaTex, plotConfig_p, "MULTIJTDPHIJJ", false, doGlobalDebug);
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  inMCFile_p->Close();
  delete inMCFile_p;

  inDataFile_p->Close();
  delete inDataFile_p;

  delete plotConfig_p;

  std::cout << "GDJDATAMCRAWPLOTTER COMPLETE. return 0." << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjDataMCRawPlotter.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjDataMCRawPlotter(argv[1]);
  return retVal;
}
