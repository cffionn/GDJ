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
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

void plotMixClosure(const bool doGlobalDebug, std::map<std::string, std::string>* labelMap, TEnv* plotConfig_p, std::string envStr, std::string dateStr,  std::vector<TH1F*> hists_p)
{
  const double max = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "MAX").c_str(), 0.5);
  const double min = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "MIN").c_str(), -0.05);
  const double zoomMax = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "ZOOMMAX").c_str(), 0.05);
  const double zoomMin = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "ZOOMMIN").c_str(), -0.05);

  const bool isMC = plotConfig_p->GetValue("ISMC", 0);
  const bool doLogX = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "DOLOGX").c_str(), 0);
  const bool doLogY = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "DOLOGY").c_str(), 0);

  const bool alignRight = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LABELALIGNRIGHT").c_str(), 0);

  std::vector<std::vector<float> > generalBoxes;

  for(int gI = 0; gI < 10; ++gI){
    std::string tempStr = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "GENBOX." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    generalBoxes.push_back(strToVectF(tempStr));
  }


  std::vector<std::string> legStrs = {"Raw", "Mixed", "Raw - Mixed", "Gen. Matched"};
  if(legStrs.size() < hists_p.size()) return;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string preStr = hists_p[0]->GetName();
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  while(preStr.find("/") != std::string::npos){preStr.replace(0, preStr.find("/")+1, "");}
  preStr.replace(preStr.find("_"), preStr.size(), "");

  int titleFont = 42;
  double titleSize = 0.035;
  double labelSize = titleSize*0.9;

  std::string labelStr = hists_p[0]->GetName();  
  labelStr.replace(0, labelStr.find("_")+1, "");
  labelStr.replace(labelStr.rfind("_h"), 2, "");
  std::vector<std::string> preLabels;
  int nGlobalLabels = 0;
  for(Int_t gI = 0; gI < 10; ++gI){
    std::string tempStr = plotConfig_p->GetValue(("MIXEDEVTPLOT.GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    while(tempStr.find("\"") != std::string::npos){tempStr.replace(tempStr.find("\""), 1, "");}
    preLabels.push_back(tempStr);
    ++nGlobalLabels;
  }

  int jetR = plotConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  else{
    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  }
  preLabels.push_back(jetRStr);
  
  /*
  if(isMC) preLabels.push_back("#bf{ATLAS Internal} Monte Carlo");
  else preLabels.push_back("#bf{ATLAS Internal} Data");
  */

  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    preLabels.push_back(centStr);
    labelStr.replace(0, labelStr.find("_")+1, "");
  }
  if(labelStr.size() != 0) preLabels.push_back(labelStr);
    
  Double_t padSplit = 0.375;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit;
  const Double_t topMargin = 0.01;
  const Double_t rightMargin = 0.01;
  Double_t height = 450.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);


  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSize/(1.0 - padSplit));

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
  TPad* pads_p[2];
  canv_p->SetTopMargin(0.001);
  canv_p->SetRightMargin(0.001);
  canv_p->SetLeftMargin(0.001);
  canv_p->SetBottomMargin(0.001);
  canv_p->cd();
  
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(topMargin);
  pads_p[0]->SetRightMargin(rightMargin);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(0.001);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetTopMargin(0.001);
  pads_p[1]->SetRightMargin(rightMargin);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetBottomMargin(bottomMargin);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  const double legX = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LEGX").c_str(), 0.7);
  const double legY = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LEGY").c_str(), 0.9);

  Int_t maxLegX = 0;
  Int_t nLeg = hists_p.size();
  for(Int_t lI = 0; lI < nLeg; ++lI){
    if(legStrs[lI].size() > (unsigned int)maxLegX) maxLegX = legStrs[lI].size();
  }

  TLegend* leg_p = new TLegend(legX, legY - nLeg*0.04/(1.0 - padSplit), legX + maxLegX*0.02, legY);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSize/(1.0-padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  double yOffset = 1.2;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int cI = 0; cI < hists_p.size(); ++cI){
    leg_p->AddEntry(hists_p[cI], legStrs[cI].c_str(), "P L");
    
    canv_p->cd();
    pads_p[0]->cd();

    HIJet::Style::EquipHistogram(hists_p[cI], cI);

    hists_p[cI]->SetMaximum(max);
    hists_p[cI]->SetMinimum(min);
    
    hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
    hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
    hists_p[cI]->GetXaxis()->SetTitleSize(0.00001);
    hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));

    hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
    hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
    hists_p[cI]->GetXaxis()->SetLabelSize(0.00001);
    hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));
    
    hists_p[cI]->GetYaxis()->SetNdivisions(505);
    hists_p[cI]->GetXaxis()->SetNdivisions(505);

    hists_p[cI]->GetYaxis()->SetTitle("N_{#gamma,jet}/N_{#gamma}");

    hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset);
    
    if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
    else hists_p[cI]->DrawCopy("HIST E1 P SAME");
    gStyle->SetOptStat(0);
    if(doLogX) gPad->SetLogx();
    if(doLogY) gPad->SetLogy();
    
    canv_p->cd();
    pads_p[1]->cd();

    hists_p[cI]->SetMaximum(zoomMax);
    hists_p[cI]->SetMinimum(zoomMin);
    if(!isMC || hists_p.size() == 3){
      std::string yTitle = hists_p[cI]->GetYaxis()->GetTitle();
      yTitle = yTitle + " (ZOOM)";
      hists_p[cI]->GetYaxis()->SetTitle(yTitle.c_str());
      hists_p[cI]->GetXaxis()->SetTitleOffset(1.2);
      hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0-padSplit));
      
      hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetXaxis()->SetTitleSize(titleSize/(padSplit));
      hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(padSplit));
      
      hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetXaxis()->SetLabelSize(labelSize/(padSplit));
      hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(padSplit));
      
      hists_p[cI]->GetYaxis()->SetNdivisions(505);
      
      if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
      else hists_p[cI]->DrawCopy("HIST E1 P SAME");  
    }
    else if(hists_p.size() == 4){
      if(cI < 2) continue;
      if(cI == 3) continue;

      hists_p[cI]->Divide(hists_p[3]);
      
      hists_p[cI]->GetYaxis()->SetTitle("(Raw - Mixed)/Gen.");
      hists_p[cI]->GetXaxis()->SetTitleOffset(1.2);
      hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0-padSplit));
      
      hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetXaxis()->SetTitleSize(titleSize/(padSplit));
      hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(padSplit));
      
      hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetXaxis()->SetLabelSize(labelSize/(padSplit));
      hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(padSplit));
      
      hists_p[cI]->GetYaxis()->SetNdivisions(505);
      
      hists_p[cI]->DrawCopy("HIST E1 P");
    }

      
    if(doLogX) gPad->SetLogx();
    if(doLogY) gPad->SetLogy();  
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  canv_p->cd();
  pads_p[0]->cd();
  leg_p->Draw("SAME");

  labelStr = "";

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const double xPos = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LABELX").c_str(), 0.18);
  double yPos = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LABELY").c_str(), 0.9);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string preLabelSaveStr = "";
  for(unsigned int pI = nGlobalLabels; pI < preLabels.size(); ++pI){
    if(preLabels[pI].find("anti-k") != std::string::npos) continue;
    if(preLabels[pI].find("ATLAS") != std::string::npos) continue;
    else preLabelSaveStr = preLabelSaveStr + preLabels[pI] + "_";
  }

  if(isMC) preLabelSaveStr = preLabelSaveStr + "MC_";
  else preLabelSaveStr = preLabelSaveStr + "DATA_";
  preLabelSaveStr = preLabelSaveStr + "R" + std::to_string(jetR) + "_";

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  double maxX = 0.0;
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    canv_p->cd();
    pads_p[0]->cd();

    std::string preStr = "";
    if(preLabels[pI].find("Cent") != std::string::npos) preStr = "";//"Pb+Pb";
    if(labelMap->count(preLabels[pI]) != 0){
      if(preStr.size() != 0) preLabels[pI] = preStr + " " + (*labelMap)[preLabels[pI]];
      else preLabels[pI] = (*labelMap)[preLabels[pI]];	
    }

    if(alignRight){
      canv_p->cd();
      pads_p[0]->cd();

      label_p->SetText(0.1, 0.1, preLabels[pI].c_str());
      if(label_p->GetXsize() > maxX) maxX = label_p->GetXsize();
    }
  }

  std::cout << "FINAL MAXX: " << maxX << std::endl;
  
  if(isMC){
    if(hists_p.size() >= 4){
      std::string assocGenMinPtStr = plotConfig_p->GetValue("ASSOCGENMINPT", "");
      
      if(assocGenMinPtStr.size() != 0) preLabels.push_back("p_{T}^{Gen. Match} > " + assocGenMinPtStr);
    }
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << maxX << ", " << preLabels.size() << std::endl;

  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){

    if(alignRight){
      label_p->SetText(xPos, yPos, preLabels[pI].c_str());
      while(label_p->GetXsize() < maxX){
	preLabels[pI] = " " + preLabels[pI];
	label_p->SetText(xPos, yPos, preLabels[pI].c_str());
      }
    }

    label_p->DrawLatex(xPos, yPos, preLabels[pI].c_str());
    yPos -= 0.083;
    /*
    if(labelStr.size() + preLabels[pI].size() > 0){
      if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
      label_p->DrawLatex(xPos, yPos, labelStr.c_str());
      yPos -= 0.085;
      labelStr = "";
    }
    labelStr = labelStr + preLabels[pI] + "; ";
    */
  }
  //  if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
  
  //  label_p->DrawLatex(xPos, yPos, labelStr.c_str());
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  canv_p->cd();
  pads_p[0]->cd();

  if(!doLogX){
     line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 0.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);
  }

  canv_p->cd();
  pads_p[1]->cd();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  if(!doLogX){
    if(!isMC || hists_p.size() == 3) line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 0.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);
    else if(hists_p.size() == 4) line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 1.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 1.0);
  }

  canv_p->cd();

  for(unsigned int gI = 0; gI < generalBoxes.size(); ++gI){
    canv_p->cd();
    drawWhiteBoxNDC(canv_p, generalBoxes[gI][0], generalBoxes[gI][1], generalBoxes[gI][2], generalBoxes[gI][3]);
  }

  std::string saveName = "pdfDir/" + dateStr + "/" + preStr + "_" + preLabelSaveStr + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  delete leg_p;
  delete pads_p[0];
  delete pads_p[1];
  delete canv_p;

  delete line_p;
  delete label_p;
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  return;
}

int gdjMixedEventPlotter(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  TEnv* plotConfig_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> reqConfigParams = {"MIXEDEVTPLOT.INFILENAME"};
  if(!checkEnvForParams(plotConfig_p, reqConfigParams)) return 1;

  std::string inFileName = plotConfig_p->GetValue("MIXEDEVTPLOT.INFILENAME", "");
  if(!check.checkFileExt(inFileName, "root")) return 1;

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
					      "NGAMMAPTBINSSUB",
					      "ISPP",
					      "JETR"};

  std::vector<std::string> pbpbParams = {"CENTBINS"};

  std::vector<std::string> mcParams = {"ASSOCGENMINPT"};
  
  if(!configs.ContainsParamSet(necessaryParams)) return 1;

  const bool isMC = config_p->GetValue("ISMC", 0);

  if(isMC){
    if(!configs.ContainsParamSet(mcParams)) return 1;

    plotConfig_p->SetValue("ASSOCGENMINPT", config_p->GetValue("ASSOCGENMINPT", ""));
  }
  plotConfig_p->SetValue("ISMC", isMC);

  const int nGammaPtBinsSub = std::stoi(configs.GetConfigVal("NGAMMAPTBINSSUB"));
  const bool isPP = std::stoi(configs.GetConfigVal("ISPP"));
  
  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!configs.ContainsParamSet(pbpbParams)) return 1;    
    
    centBins = strToVect(configs.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;
  }

  std::vector<std::string> observables1 = {"JtDPhi", "JtXJ", "JtPt", "JtMult", "MultiJtPt", "MultiJtXJ", "MultiJtXJJ"};//, "MultiJtDPhiJJ"};
  std::vector<std::string> backStr1 = {"_GlobalJtPt0", "_GlobalJtPt0_DPhi0", "_DPhi0", "_GlobalJtPt0_DPhi0", "_DPhi0_MultiJt0", "_GlobalJtPt0_DPhi0_MultiJt0", "_GlobalJtPt0_DPhi0_MultiJt0"};//, "_GlobalJtPt0_DPhi0_MultiJt0"};

  plotConfig_p->SetValue("JETR", config_p->GetValue("JETR", 100));
  int jetR = plotConfig_p->GetValue("JETR", 100);
  if(jetR != 2 && jetR != 4){
    std::cout << "JETR given \'" << jetR << "\' is not valid. Please pick \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }

  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
    
   for(Int_t gI = 0; gI < nGammaPtBinsSub+1; ++gI){
      for(unsigned int oI = 0; oI < observables1.size(); ++oI){     
	std::string rawName = centStr + "/photon" + observables1[oI] + "VCentPt_" + centStr + "_GammaPt" + std::to_string(gI) + backStr1[oI] + "_h";

	TH1F* raw_p = (TH1F*)inFile_p->Get(rawName.c_str());
	rawName.replace(rawName.find("photon"), std::string("photon").size(), "photonSub");
	TH1F* sub_p = (TH1F*)inFile_p->Get(rawName.c_str());

	rawName.replace(rawName.find("photonSub"), std::string("photonSub").size(), "photonMix");
	
	if(isStrSame(observables1[oI], "MultiJtXJJ") || isStrSame(observables1[oI], "MultiJtDPhiJJ")){
	  rawName.replace(rawName.find("photonMix"), std::string("photonMix").size(), "photonMixCorrected");
	}
	TH1F* mix_p = (TH1F*)inFile_p->Get(rawName.c_str());

	TH1F* mc_p = nullptr;
	if(isMC){
	  rawName = centStr + "/photonGenMatched" + observables1[oI] + "VCentPt_" + centStr + "_GammaPt" + std::to_string(gI) + backStr1[oI] + "_h";
	  mc_p = (TH1F*)inFile_p->Get(rawName.c_str());	  
	}
	
	if(doGlobalDebug) std::cout << "OBSERVABLE: " << observables1[oI] << std::endl;
	std::vector<TH1F*> hists_p = {raw_p, mix_p, sub_p};
	if(isMC && mc_p != nullptr) hists_p.push_back(mc_p);

	plotMixClosure(doGlobalDebug, &labelMap, plotConfig_p, strLowerToUpper(observables1[oI]), dateStr, hists_p);
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "GDJMIXEDEVENTPLOTTER COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjMixedEventPlotter.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjMixedEventPlotter(argv[1]);
  return retVal;
}
