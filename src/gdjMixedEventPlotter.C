//Author: Chris McGinn (2020.04.20)
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
#include "include/binUtils.h"
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

//Based on file provided by YJ 2021.04.27
//See input/yjDPhi for pdf and C file

void fillYJHist(TH1F* inHist_p)
{
  const Int_t nBins = 16;
  const Double_t binVals[nBins] = {0.04494337, 0.05285863, 0.06226765, 0.09007777, 0.08434847, 0.07955648, 0.08816509, 0.09671815, 0.1008858, 0.1199762, 0.1605424, 0.2304045, 0.3904447, 0.6729873, 1.310041, 3.538677};
  const Double_t errVals[nBins] = {0.002512662, 0.003018931, 0.002655415, 0.002963426, 0.003064175, 0.002952764, 0.003097857, 0.00307974, 0.003240865, 0.003212133, 0.003609592, 0.003765436, 0.004419775, 0.005151995, 0.007508904, 0.01069313};

  if(inHist_p->GetNbinsX() != nBins){
    std::cout << "FillYJHist given an input w/ differing numbers of bins than needed. return" << std::endl;
    return;
  }

  for(Int_t i = 0; i < nBins; ++i){
    inHist_p->SetBinContent(i+1, binVals[i]);
    inHist_p->SetBinError(i+1, errVals[i]);
  }
  
  return;
}


/*
std::vector<std::string> getLabels(TEnv* plotConfig_p, TH1F* histForLabels_p, std::map<std::string, std::string>* labelMap, std::vector<std::string>* labelsForSaveStr = nullptr)
{
  std::vector<std::string> labelVect;
  Int_t nGlobalLabels = 0;
  for(Int_t gI = 0; gI < 10; ++gI){
    std::string tempStr = plotConfig_p->GetValue(("MIXEDEVTPLOT.GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    while(tempStr.find("\"") != std::string::npos){tempStr.replace(tempStr.find("\""), 1, "");}

    labelVect.push_back(tempStr);
    ++nGlobalLabels;
  }

  int jetR = plotConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  else{
    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  }
  labelVect.push_back(jetRStr);

  std::string histName = histForLabels_p->GetName();
  std::string labelStr = histName;
  if(labelStr.find("_") != std::string::npos) labelStr.replace(0, labelStr.find("_")+1, "");
  if(labelStr.rfind("_h") != std::string::npos) labelStr.replace(labelStr.rfind("_h"), 2, "");

  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    labelStr.replace(0, labelStr.find("_")+1, "");

    if(histName.find("DPhiVCent") != std::string::npos || histName.find("DPhiJJGVCent") != std::string::npos){
      if(labelStr.find("DPhi") == std::string::npos) continue;
    }
    
    if(centStr.find("Alt") == std::string::npos) labelVect.push_back(centStr);
    if(isStrSame("RAW_h", labelStr)) break;
  }
  if(labelStr.size() != 0) labelVect.push_back(labelStr);
  
  
  if(labelsForSaveStr != nullptr){
    for(unsigned int pI = nGlobalLabels; pI < labelVect.size(); ++pI){
      if(labelVect[pI].find("anti-k") != std::string::npos) continue;
      if(labelVect[pI].find("ATLAS") != std::string::npos) continue;      

      labelsForSaveStr->push_back(labelVect[pI]);
    }       
  }
  
  for(unsigned int pI = 0; pI < labelVect.size(); ++pI){
    if(labelMap->count(labelVect[pI]) != 0){
      labelVect[pI] = (*labelMap)[labelVect[pI]];
    }
  }


  std::string jtPtLowStr = plotConfig_p->GetValue("JTPTBINSLOW", "");
  std::string jtPtHighStr = plotConfig_p->GetValue("JTPTBINSHIGH", "");
  std::string jtEtaLowStr = plotConfig_p->GetValue("JTETABINSLOW", "");
  std::string jtEtaHighStr = plotConfig_p->GetValue("JTETABINSHIGH", "");
  std::string mixJJDR = plotConfig_p->GetValue("MIXJETEXCLUSIONDR", "");
  //Add the three relevant jet cuts
  labelVect.push_back(jtPtLowStr + "<p_{T,Jet}<" + jtPtHighStr);
  labelVect.push_back(jtEtaLowStr + "<#eta_{Jet}<" + jtEtaHighStr);
  labelVect.push_back("#DeltaR_{JJ}>" + mixJJDR);

  for(unsigned int pI = 0; pI < labelVect.size(); ++pI){
    if(labelVect[pI].find("Barrel") != std::string::npos || labelVect[pI].find("EC") != std::string::npos) labelVect[pI] = labelVect[pI] + " #gamma";
  }

  if(histName.find("XJJVCent") != std::string::npos || histName.find("AJJVCent") != std::string::npos){
    std::string multiDPhiJJStr = "|#Delta#phi_{#vec{JJ},#gamma}| > " + std::string(plotConfig_p->GetValue("GAMMAMULTIJTDPHI", ""));
    //    std::cout << "MULTIJTDPHI 1: " << multiDPhiJJStr << std::endl;
    if(multiDPhiJJStr.find("pi") != std::string::npos) multiDPhiJJStr.replace(multiDPhiJJStr.find("pi"), 2, "#pi");
    //    std::cout << "MULTIJTDPHI 2: " << multiDPhiJJStr << std::endl;
    labelVect.push_back(multiDPhiJJStr);    
  }

  const bool isMC = plotConfig_p->GetValue("ISMC", 0);
  if(isMC){
    std::string assocGenMinPtStr = plotConfig_p->GetValue("ASSOCGENMINPT", "");
    if(assocGenMinPtStr.size() != 0) labelVect.push_back("p_{T,Jet}^{Truth Match} > " + assocGenMinPtStr);
  }

  return labelVect;
}
*/

std::string plotMixClosure(const bool doGlobalDebug, std::map<std::string, std::string>* labelMap, TEnv* plotConfig_p, std::string envStr, std::string dateStr, std::string saveTag,  std::vector<TH1F*> hists_p, std::vector<std::string> legStrs, std::vector<TH1F*> refHist_p, std::vector<std::string> refLegStr, std::vector<std::string> yLabels, std::vector<bool> yLogs, bool doReducedLabel)
{
  if(hists_p.size() == 0){
    std::cout << "No hists given to plotMixClosure - return" << std::endl;
    return "";
  }
  if(refHist_p.size() >= 3){
    std::cout << "Max 2 references allowed. return" << std::endl;
    return "";
  }
  if(refHist_p.size() == 0){
    std::cout << "Must have at least 1 reference allowed. return" << std::endl;
    return "";
  }
  for(unsigned int hI = 0; hI < refHist_p.size(); ++hI){
    hists_p.push_back(refHist_p[hI]);
    legStrs.push_back(refLegStr[hI]);
  }
  
  //Fix the y-labels
  std::string yLabelNominal = "N_{J,#gamma}/N_{#gamma}";
  std::string firstHistName = hists_p[0]->GetName();
  if(firstHistName.find("XJJVCent") != std::string::npos){
    yLabelNominal = "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{d#vec{x}_{JJ#gamma}}";        
  }
  else if(firstHistName.find("AJJVCent") != std::string::npos){
    yLabelNominal = "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{dA_{JJ#gamma}}";        
  }
  else if(firstHistName.find("XJVCent") != std::string::npos){
    yLabelNominal = "#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{dx_{J#gamma}}";        
  }
  else if(firstHistName.find("JJVCent") != std::string::npos){
    yLabelNominal = "N_{JJ,#gamma}/N_{#gamma}";        
  }
  else if(firstHistName.find("DPhiVCent") != std::string::npos){
     bool yjHistFound = false;
    for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
      std::string histName = hists_p[hI]->GetName();

      if(histName.find("yjHist") != std::string::npos){
	yjHistFound = true;
	break;
      }
    }

    if(yjHistFound){
      yLabelNominal = "#frac{1}{N_{#gamma}} #frac{dN_{J#gamma}}{d#Delta#phi_{J#gamma}}";
    }
  }
  
  const bool doVaryMax = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "DOVARYMAX").c_str(), 0);
  const double max = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "MAX").c_str(), 0.5);
  const double min = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "MIN").c_str(), -0.05);
  const double zoomMax = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "ZOOMMAX").c_str(), 0.05);
  const double zoomMin = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "ZOOMMIN").c_str(), -0.05);

  /*
  std::cout << "ENVSTR: " << envStr << std::endl;
  std::cout << " doVaryMax: " << doVaryMax << std::endl;
  std::cout << " max: " << max << std::endl;
  std::cout << " min: " << min << std::endl;
  */  

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
  if(legStrs.size() != hists_p.size()) return "";

  const int jetR = plotConfig_p->GetValue("JETR", 100);  

  std::string preStr = hists_p[0]->GetName();
  while(preStr.find("/") != std::string::npos){preStr.replace(0, preStr.find("/")+1, "");}
  preStr.replace(preStr.find("_"), preStr.size(), "");

  int titleFont = 42;
  double titleSize = 0.035;
  double labelSize = titleSize*0.9;

  std::vector<std::string> preLabelsAlt;
  std::vector<std::string> preLabels = getLabels(plotConfig_p, hists_p[0], labelMap, &preLabelsAlt);
  
  const Int_t nMaxBins = 1000;
  const Int_t nBinsTemp = hists_p[0]->GetXaxis()->GetNbins();
  if(nBinsTemp > nMaxBins){
    std::cout << "Number of bins needed \'" << nBinsTemp << "\' exceeds max \'" << nMaxBins << "\'. return 1" << std::endl;
    return "";
  }
  Double_t binsTemp[nMaxBins+1];
  for(Int_t bIX = 0; bIX < nBinsTemp+1; ++bIX){
    binsTemp[bIX] = hists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1);
  }

  Double_t padSplit = 0.32;
  const Double_t leftMargin2Pad = 0.16;
  const Double_t bottomMargin2Pad = 0.125/padSplit;
  const Double_t topMargin2Pad = 0.01;
  const Double_t rightMargin2Pad = 0.01;
  Double_t height2Pad = 900.0;
  Double_t width2Pad = height2Pad*(1.0 - topMargin2Pad*(1.0 - padSplit) - bottomMargin2Pad*padSplit)/(1.0 - leftMargin2Pad - rightMargin2Pad);

  //Re-write the 3 pad code for a longer plot
  Double_t height3Pad = 1200.0;
  Double_t width3Pad = 942.831;
  const Double_t leftMargin3Pad = 0.20;
  const Double_t bottomAbsFracOfLeft = 0.8;
  const Double_t bottomMarginFracHeight = bottomAbsFracOfLeft*leftMargin3Pad*width3Pad/height3Pad;
  const Double_t topPanelFrac = 0.48;
  const Double_t bottomPanelFrac = 0.39;
  const Double_t padSplit1 = 1.0 - topPanelFrac;
  const Double_t padSplit2 = bottomPanelFrac;

  const Double_t bottomMargin3Pad = bottomMarginFracHeight/padSplit2;
  const Double_t topMargin3Pad = 0.01;
  const Double_t rightMargin3Pad = 0.01;
  
  Int_t nPads = 3;
  Double_t leftMargin = leftMargin3Pad;
  Double_t rightMargin = rightMargin3Pad;
  Double_t bottomMargin = bottomMargin3Pad;
  Double_t topMargin = topMargin3Pad;
  Double_t height = height3Pad;
  Double_t width = width3Pad;
  
  if(refHist_p.size() == 1){
    nPads = 2;
    leftMargin = leftMargin2Pad;
    rightMargin = rightMargin2Pad;
    bottomMargin = bottomMargin2Pad;
    topMargin = topMargin2Pad;
    height = height2Pad;
    width = width2Pad;
  }

  if(nPads != (Int_t)yLabels.size()){
    std::cout << "Number of pads \'" << nPads << "\' does not match the yLables.size() \'" << yLabels.size() << "\'. return" << std::endl;
    return "";
  }
  
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  if(nPads == 3) label_p->SetTextSize(titleSize/(1.0-padSplit1));
  else label_p->SetTextSize(titleSize/(1.0-padSplit));

  TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
  TPad* pads_p[3];
  canv_p->SetTopMargin(0.001);
  canv_p->SetRightMargin(0.001);
  canv_p->SetLeftMargin(0.001);
  canv_p->SetBottomMargin(0.001);
  canv_p->cd();  
  
  if(nPads == 3) pads_p[0] = new TPad("pad0", "", 0.0, padSplit1, 1.0, 1.0);
  else pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(topMargin);
  pads_p[0]->SetRightMargin(rightMargin);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(0.001);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

  if(nPads == 3) pads_p[1] = new TPad("pad1", "", 0.0, padSplit2, 1.0, padSplit1);
  else pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);

  pads_p[1]->SetTopMargin(0.001);
  pads_p[1]->SetRightMargin(rightMargin);
  pads_p[1]->SetLeftMargin(leftMargin);
  if(nPads == 3) pads_p[1]->SetBottomMargin(0.001);
  else pads_p[1]->SetBottomMargin(bottomMargin);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  if(nPads == 3){
    pads_p[2] = new TPad("pad2", "", 0.0, 0.0, 1.0, padSplit2);
    pads_p[2]->SetTopMargin(0.001);
    pads_p[2]->SetRightMargin(rightMargin3Pad);
    pads_p[2]->SetLeftMargin(leftMargin3Pad);
    pads_p[2]->SetBottomMargin(bottomMargin3Pad);

    canv_p->cd();
    pads_p[2]->Draw("SAME");
    pads_p[2]->cd();
    canv_p->cd();
  }
  
  const double legX = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LEGX").c_str(), 0.7);
  const double legY = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LEGY").c_str(), 0.9);

  Int_t maxLegX = 0;
  Int_t nLeg = hists_p.size();
  
  for(Int_t lI = 0; lI < nLeg; ++lI){
    if(legStrs[lI].size() > (unsigned int)maxLegX) maxLegX = legStrs[lI].size();
  }
  
  TLegend* leg_p = nullptr;
  if(nPads == 3) leg_p = new TLegend(legX, legY - nLeg*0.04/(1.0 - padSplit1), legX + maxLegX*0.016, legY);
  else leg_p = new TLegend(legX, legY - nLeg*0.04/(1.0 - padSplit), legX + maxLegX*0.016, legY);

  leg_p->SetTextFont(titleFont);

  if(nPads == 3) leg_p->SetTextSize(titleSize/(1.0-padSplit1));
  else leg_p->SetTextSize(titleSize/(1.0-padSplit));

  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  double yOffset = 0.9;

  //Gotta do a quick pre-process to get nice ranges
  std::vector<Double_t> refMaxIfLog, refMinIfLog;
  for(unsigned int hI = 0; hI < refHist_p.size(); ++hI){
    refMaxIfLog.push_back(-1);
    refMinIfLog.push_back(-1);

    if(yLogs[hI+1]){
      for(unsigned int cI = 0; cI < hists_p.size(); ++cI){
	std::string histName = hists_p[cI]->GetName();
	if(histName.find("yjHist_p") != std::string::npos) continue;
	if(refHist_p[hI] == hists_p[cI]) continue;
       	
	TH1F* tempHist_p = new TH1F("tempHist_p", "", nBinsTemp, binsTemp);	
	tempHist_p->Divide(hists_p[cI], refHist_p[hI]);

	
	Double_t tempMin = getMinGTZero(tempHist_p);
	Double_t tempMax = getMax(tempHist_p);

	if(refMaxIfLog[hI] < tempMax) refMaxIfLog[hI] = tempMax;

	if(refMinIfLog[hI] < 0) refMinIfLog[hI] = tempMin;
	else if(refMinIfLog[hI] > tempMin) refMinIfLog[hI] = tempMin;

	delete tempHist_p;
      }
    }
  }

  double varyMax = -1.0;
  if(doVaryMax){
    for(unsigned int cI = 0; cI < hists_p.size(); ++cI){
      double tempMax = hists_p[cI]->GetMaximum();
      if(tempMax > varyMax) varyMax = tempMax;
    }    
  }
  
  for(unsigned int cI = 0; cI < hists_p.size(); ++cI){
    leg_p->AddEntry(hists_p[cI], legStrs[cI].c_str(), "P L");
    
    canv_p->cd();
    pads_p[0]->cd();

    HIJet::Style::EquipHistogram(hists_p[cI], cI);

    if(doVaryMax) hists_p[cI]->SetMaximum(varyMax*1.1);
    else hists_p[cI]->SetMaximum(max);
    hists_p[cI]->SetMinimum(min);
    
    hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
    hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
    hists_p[cI]->GetXaxis()->SetTitleSize(0.00001);

    if(nPads == 3) hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit1));
    else hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));

    hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
    hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
    hists_p[cI]->GetXaxis()->SetLabelSize(0.00001);
    if(nPads == 3) hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit1));
    else hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));
    
    hists_p[cI]->GetYaxis()->SetNdivisions(505);
    hists_p[cI]->GetXaxis()->SetNdivisions(505);
  
    //    if(yLabels[0].size() == 0) hists_p[cI]->GetYaxis()->SetTitle(yLabelNominal.c_str());
    if(yLabels[0].size() != 0) hists_p[cI]->GetYaxis()->SetTitle(yLabels[0].c_str());

    //    if(yLabelNominal.find("#frac") != std::string::npos) hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset - 0.2);
    
    if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
    else hists_p[cI]->DrawCopy("HIST E1 P SAME");

    gPad->SetTicks();
    gStyle->SetOptStat(0);
    if(doLogX) gPad->SetLogx();
    if(doLogY) gPad->SetLogy();
    
    //editing here

    for(unsigned int hI = 0; hI < refHist_p.size(); ++hI){
      if(refHist_p[hI] == hists_p[cI]) continue;
      std::string histName = hists_p[cI]->GetName();
      if(histName.find("yjHist_p") != std::string::npos) continue;

      canv_p->cd();
      pads_p[hI+1]->cd();

      TH1F* tempHist_p = new TH1F("tempHist_p", "", nBinsTemp, binsTemp);
      tempHist_p->GetXaxis()->SetTitle(hists_p[cI]->GetXaxis()->GetTitle());

      
      tempHist_p->Divide(hists_p[cI], refHist_p[hI]);
    
      if(yLogs[hI+1]){
	tempHist_p->SetMaximum(refMaxIfLog[hI]*10.);
	tempHist_p->SetMinimum(refMinIfLog[hI]/10.);
      }
      else{
	tempHist_p->SetMaximum(zoomMax);
	tempHist_p->SetMinimum(zoomMin);
      }
      
      std::string yTitle = tempHist_p->GetYaxis()->GetTitle();
      yTitle = yTitle + " (ZOOM)";
      tempHist_p->GetXaxis()->SetTitleOffset(yOffset);


      //      if(yLabels[hI+1].find("#frac") == std::string::npos){

      if(nPads == 3){
	if(hI == 0) tempHist_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit1 - padSplit2)/(1.0 - padSplit1));
	else tempHist_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit2)/(1.0 - padSplit1));
      }
      else tempHist_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit)/(1.0 - padSplit));
      /*
    }
      else{
	if(nPads == 3){
	  if(hI == 0) tempHist_p->GetYaxis()->SetTitleOffset((yOffset-0.2)*(padSplit1 - padSplit2)/(1.0 - padSplit1));
	  else tempHist_p->GetYaxis()->SetTitleOffset((yOffset-0.2)*(padSplit2)/(1.0 - padSplit1));
	}
	else tempHist_p->GetYaxis()->SetTitleOffset((yOffset-0.2)*(padSplit)/(1.0 - padSplit));
      }
      */
      tempHist_p->GetXaxis()->SetTitleFont(titleFont);
      tempHist_p->GetYaxis()->SetTitleFont(titleFont);
      
      tempHist_p->GetXaxis()->SetLabelFont(titleFont);
      tempHist_p->GetYaxis()->SetLabelFont(titleFont);

      if(nPads == 3){
	if(hI == 0){
	  tempHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  tempHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
	  tempHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	  tempHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
	}
	else{
	  tempHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit2));
	  tempHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit2));
	  tempHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit2));
	  tempHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit2));
	}
      }
      else{
	tempHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit));
	tempHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));
	tempHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit));
	tempHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));
      }
      
      tempHist_p->SetMarkerColor(hists_p[cI]->GetMarkerColor());
      tempHist_p->SetMarkerSize(hists_p[cI]->GetMarkerSize());
      tempHist_p->SetMarkerStyle(hists_p[cI]->GetMarkerStyle());
      tempHist_p->SetLineColor(hists_p[cI]->GetLineColor());
    
      if(yLabels[hI+1].size() == 0) tempHist_p->GetYaxis()->SetTitle(yTitle.c_str());
      else tempHist_p->GetYaxis()->SetTitle(yLabels[hI+1].c_str());
        
      if(hI == 0) tempHist_p->GetYaxis()->SetNdivisions(404);
      else tempHist_p->GetYaxis()->SetNdivisions(505);

      tempHist_p->GetXaxis()->SetNdivisions(505);
      
      if(cI == 0) tempHist_p->DrawCopy("HIST E1 P");
      else tempHist_p->DrawCopy("HIST E1 P SAME");  

      gPad->SetTicks();
      
      if(yLogs[1+hI]) gPad->SetLogy();
      else{
	if(zoomMax > 1.0 && zoomMin < 1.0) line_p->DrawLine(tempHist_p->GetBinLowEdge(1), 1.0, tempHist_p->GetBinLowEdge(tempHist_p->GetXaxis()->GetNbins()+1), 1.0);
;
      }
      
      delete tempHist_p;
    }
    /*
    //EDITING HERE 2 - ADVANCED
    else if(hists_p.size() == 4){
      if(cI < 2) continue;
      if(cI == 3) continue;

      TH1F* tempHist_p = new TH1F("tempHist_p", "", nBinsTemp, binsTemp);
      tempHist_p->GetXaxis()->SetTitle(hists_p[cI]->GetXaxis()->GetTitle());            
      
      
      tempHist_p->Divide(hists_p[2], hists_p[1]);

      tempHist_p->SetMaximum(100);
      tempHist_p->SetMinimum(0.01);
      
      HIJet::Style::EquipHistogram(tempHist_p, cI);
      
      tempHist_p->GetYaxis()->SetTitle("Sig/Bkgd");
      tempHist_p->GetXaxis()->SetTitleOffset(1.2);
      tempHist_p->GetYaxis()->SetTitleOffset(yOffset*(padSplit1 - padSplit2)/(1.0-padSplit1));
      
      tempHist_p->GetXaxis()->SetTitleFont(titleFont);
      tempHist_p->GetYaxis()->SetTitleFont(titleFont);
      tempHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
      tempHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit1 - padSplit2));
      
      tempHist_p->GetXaxis()->SetLabelFont(titleFont);
      tempHist_p->GetYaxis()->SetLabelFont(titleFont);
      tempHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
      tempHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit1 - padSplit2));
      
      tempHist_p->GetYaxis()->SetNdivisions(505);
      
      tempHist_p->DrawCopy("HIST E1 P");
      gPad->SetLogy();
      
      delete tempHist_p;
    }

    if(doLogX) gPad->SetLogx();
  
    canv_p->cd();
    pads_p[2]->cd();

    if(!isMC || hists_p.size() == 3){
      hists_p[cI]->SetMaximum(zoomMax);
      hists_p[cI]->SetMinimum(zoomMin);
      
    
      std::string yTitle = hists_p[cI]->GetYaxis()->GetTitle();
      yTitle = yTitle + " (ZOOM)";
      hists_p[cI]->GetYaxis()->SetTitle(yTitle.c_str());
      hists_p[cI]->GetXaxis()->SetTitleOffset(1.2);
      hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset*(padSplit2)/(1.0 - padSplit1));
      
      hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetXaxis()->SetTitleSize(titleSize/(padSplit2));
      hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(padSplit2));
      
      hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetXaxis()->SetLabelSize(labelSize/(padSplit2));
      hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(padSplit2));
      
      hists_p[cI]->GetYaxis()->SetNdivisions(505);
      
      if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
      else hists_p[cI]->DrawCopy("HIST E1 P SAME");  
    }
    else if(hists_p.size() == 4){
      if(cI < 2) continue;
      if(cI == 3) continue;

      hists_p[cI]->SetMaximum(2.0);
      hists_p[cI]->SetMinimum(0.0);
      
      hists_p[cI]->Divide(hists_p[3]);
      
      hists_p[cI]->GetYaxis()->SetTitle("Sub./Gen.");
      hists_p[cI]->GetXaxis()->SetTitleOffset(1.2);
      hists_p[cI]->GetYaxis()->SetTitleOffset(yOffset*(padSplit2)/(1.0-padSplit1));
      
      hists_p[cI]->GetXaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetYaxis()->SetTitleFont(titleFont);
      hists_p[cI]->GetXaxis()->SetTitleSize(titleSize/(padSplit2));
      hists_p[cI]->GetYaxis()->SetTitleSize(titleSize/(padSplit2));
      
      hists_p[cI]->GetXaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetYaxis()->SetLabelFont(titleFont);
      hists_p[cI]->GetXaxis()->SetLabelSize(labelSize/(padSplit2));
      hists_p[cI]->GetYaxis()->SetLabelSize(labelSize/(padSplit2));
      
      hists_p[cI]->GetYaxis()->SetNdivisions(505);
      
      hists_p[cI]->DrawCopy("HIST E1 P");
    }

    if(doLogX) gPad->SetLogx();
    */
  }

  canv_p->cd();
  pads_p[0]->cd();
  leg_p->Draw("SAME");

  const double xPos = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LABELX").c_str(), 0.18);
  double yPos = plotConfig_p->GetValue(("MIXEDEVTPLOT." + envStr + "LABELY").c_str(), 0.9);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string preLabelSaveStr = "";
  for(unsigned int pI = 0; pI < preLabelsAlt.size(); ++pI){
    if(preLabelsAlt[pI].find("anti-k") != std::string::npos) continue;
    if(preLabelsAlt[pI].find("ATLAS") != std::string::npos) continue;
    else preLabelSaveStr = preLabelSaveStr + preLabelsAlt[pI] + "_";
  }

  if(isMC) preLabelSaveStr = preLabelSaveStr + "MC_";
  else preLabelSaveStr = preLabelSaveStr + "DATA_";
  preLabelSaveStr = preLabelSaveStr + "R" + std::to_string(jetR) + "_";
  
  if(alignRight) label_p->SetTextAlign(31);
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){    
    if(preLabels[pI].find("Alt") != std::string::npos) continue;

    if(!doReducedLabel || preLabels[pI].find("ATLAS") != std::string::npos || preLabels[pI].find("PYT") != std::string::npos || preLabels[pI].find("Pb+Pb") != std::string::npos) label_p->DrawLatex(xPos, yPos, preLabels[pI].c_str());
    yPos -= 0.083;
  }

  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  canv_p->cd();
  pads_p[0]->cd();

  if(!doLogX){
    line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 0.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);
  }

  if(zoomMax > 1.0 && zoomMin < 1.0 && !doLogX){
    for(unsigned int hI = 0; hI < refHist_p.size(); ++hI){
      canv_p->cd();
      pads_p[hI+1]->cd();

      line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 1.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 1.0);      
    }
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  
  if(nPads == 3){
    canv_p->cd();
    pads_p[2]->cd();
    
    if(!doLogX){
      if(hists_p.size() == 4) line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 1.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 1.0);
    }
  }

  canv_p->cd();

  for(unsigned int gI = 0; gI < generalBoxes.size(); ++gI){
    canv_p->cd();
    drawWhiteBoxNDC(canv_p, generalBoxes[gI][0], generalBoxes[gI][1], generalBoxes[gI][2], generalBoxes[gI][3]);
  }

  std::string extStr = plotConfig_p->GetValue("SAVEEXT", "");
  if(extStr.substr(0,1).find(".") == std::string::npos) extStr = "." + extStr;
  std::string saveName = "pdfDir/" + dateStr + "/" + preStr + "_" + preLabelSaveStr;

  bool barrelECStr = false;
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    if(preLabels[pI].find("Barrel") != std::string::npos){
      if(preLabels[pI].find("EC") != std::string::npos){
	barrelECStr = true;
	break;
      }
    }
  }
  if(barrelECStr && false){
    saveName.replace(saveName.find("Barrel"), 6, "BarrelAndEC");
  }
  
  if(saveName.substr(saveName.size()-1, 1).find("_") != std::string::npos) saveName = saveName.substr(0, saveName.size()-1);
  if(saveTag.size() != 0){
    saveName = saveName + "_" + saveTag;
  }
  saveName = saveName + "_" + dateStr + extStr;
  quietSaveAs(canv_p, saveName);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  delete leg_p;
  delete pads_p[0];
  delete pads_p[1];
 
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(nPads == 3) delete pads_p[2];
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete canv_p;
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  delete line_p;
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete label_p;
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  return saveName;
}

int gdjMixedEventPlotter(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  TEnv* plotConfig_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> reqConfigParams = {"MIXEDEVTPLOT.INFILENAME", "SAVETAG", "SAVEEXT", "DOREDUCEDLABEL"};
  if(!checkEnvForParams(plotConfig_p, reqConfigParams)) return 1;

  const std::string globalSaveTag = plotConfig_p->GetValue("SAVETAG", "");
  const bool doReducedLabel = plotConfig_p->GetValue("DOREDUCEDLABEL", 0);
  
  std::string inFileName = plotConfig_p->GetValue("MIXEDEVTPLOT.INFILENAME", "");
  if(!check.checkFileExt(inFileName, "root")) return 1;

  std::string mcFileName = plotConfig_p->GetValue("MIXEDEVTPLOT.MCFILENAME", "");
  if(mcFileName.size() != 0){
    if(!check.checkFileExt(mcFileName, "root")) return 1;
  }

  const std::string extStr = plotConfig_p->GetValue("SAVEEXT", "");
  std::vector<std::string> validExt = {"pdf", "png"};
  if(!vectContainsStr(extStr, &validExt)){
    std::cout << "Given extension \'" << extStr << "\' is not valid. Please pick from:" << std::endl;
    for(unsigned int eI = 0; eI < validExt.size(); ++eI){
      std::cout << " " << validExt[eI] << std::endl;
    }
    std::cout << "return 1. " << std::endl;
    return 1;
  }
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  TFile* mcFile_p = nullptr;
  if(mcFileName.size() != 0) mcFile_p = new TFile(mcFileName.c_str(), "READ");
  
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
  
  if(!configs.ContainsParamSet(necessaryParams)) return 1;
  const bool isMC = config_p->GetValue("ISMC", 0);
  
  plotConfig_p->SetValue("ISMC", isMC);
 
  //  const int nGammaPtBinsSub = std::stoi(configs.GetConfigVal("NGAMMAPTBINSSUB"));
  const bool isPP = std::stoi(configs.GetConfigVal("ISPP"));
  
  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!configs.ContainsParamSet(pbpbParams)) return 1;    
    
    centBins = strToVect(configs.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;
  }

  //  std::vector<std::string> observables1 = {"JtXJ", "JtXJJ", "JtAJJ", "JtDPhiJJG", "JtDPhiJJ", "JtDRJJ", "JtPt", "JtDPhi"};
  //  std::vector<std::string> mixStrings = {"MIX", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED", "MIX", "MIX"};
  //  std::vector<std::string> barrelECStr = {"Barrel", "EC", "BarrelAndEC"};
  std::vector<std::string> barrelECStr = {"BarrelAndEC"};

  /*
  std::vector<std::string> normYAxisTitle = {"#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{dx_{J#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{dx_{JJ#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{dA_{JJ#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{d#Delta#phi_{JJ#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ}}{d#Delta#phi_{JJ}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ}}{d#DeltaR_{JJ}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{dp_{T,J#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{d#Delta#phi_{J#gamma}}"};
  */
  
  std::vector<std::string> observables1 = {"JtPt", "JtXJ", "JtDPhi", "JtXJJ", "JtDRJJ", "JtDPhiJJG", "JtDPhiJJ", "JtAJJ"};
  std::vector<std::string> mixStrings = {"MIX", "MIX", "MIX", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED", "MIXCORRECTED"};
  std::vector<std::string> normYAxisTitle = {"#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{dp_{T}^{Jet}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{dx_{J#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{J#gamma}}{d#Delta#phi_{J#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{dx_{JJ#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ}}{d#DeltaR_{JJ}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{d#Delta#phi_{JJ#gamma}}",
					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ}}{d#Delta#phi_{JJ}}",
  					     "#frac{1}{N_{#gamma}}#frac{dN_{JJ#gamma}}{dA_{JJ#gamma}}"};

  //  std::vector<std::string> observables1 = {"JtXJJ", "JtDRJJ"};
  //  std::vector<std::string> mixStrings = {"MIXCORRECTED", "MIXCORRECTED"};
  //  std::vector<std::string> barrelECStr = {"Barrel"};

  std::vector<std::vector<std::string> > globalLabels;
  
  plotConfig_p->SetValue("JETR", config_p->GetValue("JETR", 100));
  plotConfig_p->SetValue("GAMMAMULTIJTDPHI", config_p->GetValue("GAMMAMULTIJTDPHI", ""));
  int jetR = plotConfig_p->GetValue("JETR", 100);
  if(jetR != 2 && jetR != 4){
    std::cout << "JETR given \'" << jetR << "\' is not valid. Please pick \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }

  plotConfig_p->SetValue("JTPTBINSLOW", config_p->GetValue("JTPTBINSLOW", ""));
  plotConfig_p->SetValue("JTPTBINSHIGH", config_p->GetValue("JTPTBINSHIGH", ""));
  plotConfig_p->SetValue("JTPTBINSLOWRECO", config_p->GetValue("JTPTBINSLOWRECO", ""));
  plotConfig_p->SetValue("JTPTBINSHIGHRECO", config_p->GetValue("JTPTBINSHIGHRECO", ""));
  plotConfig_p->SetValue("JTETABINSLOW", config_p->GetValue("JTETABINSLOW", ""));
  plotConfig_p->SetValue("JTETABINSHIGH", config_p->GetValue("JTETABINSHIGH", ""));
  plotConfig_p->SetValue("MIXJETEXCLUSIONDR", config_p->GetValue("MIXJETEXCLUSIONDR", ""));

  //We need gamma pt bins to know which y axis positions to exclude
  const Float_t gammaPtBinsLowReco = config_p->GetValue("GAMMAPTBINSLOWRECO", 1000.0);
  const Float_t gammaPtBinsHighReco = config_p->GetValue("GAMMAPTBINSHIGHRECO", -1000.0);
  
  const Int_t nMaxBins = 100;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
    
    for(unsigned int bI = 0; bI < barrelECStr.size(); ++bI){
      //THIS MUST BE FIXED TO BE PURCORR (ABOVE)
      std::string mainBarrelECStr = barrelECStr[bI];
      bool isBarrelAndEC = barrelECStr[bI].find("BarrelAndEC") != std::string::npos;
      //      if(isBarrelAndEC) mainBarrelECStr = "Barrel";	

      std::string phoName = centStr + "/photonPtVCent_" + centStr + "_" + mainBarrelECStr + "_RAW_h";
      if(doGlobalDebug) std::cout << "FILE, LINE, PHONAME: " << __FILE__ << ", " << __LINE__ << ", " << phoName << std::endl;
      TH1F* photonHist_p = (TH1F*)inFile_p->Get(phoName.c_str());
      TH1F* photonHistMC_p = nullptr;
      if(mcFileName.size() != 0) photonHistMC_p = (TH1F*)mcFile_p->Get(phoName.c_str());

      //Commenting out all references to EC
      //      TH1F* photonHistEC_p = nullptr;
      //      TH1F* photonHistMCEC_p = nullptr;
      /*
      if(isBarrelAndEC){	
	if(!strReplace(&phoName, "Barrel", "EC")) return 1;
	//	phoName.replace(phoName.find("Barrel"), 6, "EC");
	//        photonHistEC_p = (TH1F*)inFile_p->Get(phoName.c_str());
	//	photonHist_p->Add(photonHistEC_p);

	if(mcFileName.size() != 0){
	  //	  photonHistMCEC_p = (TH1F*)mcFile_p->Get(phoName.c_str());
	  //	  photonHistMC_p->Add(photonHistMCEC_p);
	}	
      }      
      */
      
      //Start grabbing histograms for plots
      for(unsigned int oI = 0; oI < observables1.size(); ++oI){     
	std::string rebinStr = strLowerToUpper("MIXEDEVTPLOT." + observables1[oI]);
	
	bool doTempCorr = plotConfig_p->GetValue(strLowerToUpper("MIXEDEVTPLOT." + observables1[oI] + "TEMPCORR").c_str(), 0);
	
	bool doRebinX = (bool)plotConfig_p->GetValue((rebinStr + "DOREBINX").c_str(), 0);
	bool doRebinY = (bool)plotConfig_p->GetValue((rebinStr + "DOREBINY").c_str(), 0);
	bool doRebin = doRebinX || doRebinY;
	
	bool isMultijet = isStrSame(observables1[oI], "JtXJJ") || isStrSame(observables1[oI], "JtAJJ") || isStrSame(observables1[oI], "JtDPhiJJG") || isStrSame(observables1[oI], "JtDPhiJJ") || isStrSame(observables1[oI], "JtDRJJ");

	std::string mixMode = "MIXMODE1";
	if(isMultijet) mixMode = "MIXMODE2";

	std::string rawName = centStr + "/photonPt" + observables1[oI] + "VCent_" + centStr + "_" + mainBarrelECStr + "_DPhi0_Fine_" + mixMode + "_RAW_h";
	if(!doRebin) rawName.replace(rawName.find("_Fine"), 5, "");

	std::string mixName = rawName;
	if(!strReplace(&mixName, "RAW", mixStrings[oI])) return 1;
	std::string subName = rawName;
	if(!strReplace(&subName, "RAW", "SUB")) return 1;

	/*
	std::string rawNameEC = rawName;
	std::string mixNameEC = mixName;
	std::string subNameEC = subName;
	if(isBarrelAndEC){
	  if(!strReplace(&rawNameEC, "Barrel", "EC")) return 1;
	  if(!strReplace(&mixNameEC, "Barrel", "EC")) return 1;
	  if(!strReplace(&subNameEC, "Barrel", "EC")) return 1;
	}
	*/
	
	//Big block declaring all histograms
	TH2F* raw_p = nullptr;
	TH2F* mix_p = nullptr;
	TH2F* sub_p = nullptr;
      
	/*
	TH2F* raw_p = (TH2F*)inFile_p->Get(rawName.c_str());
	TH2F* mix_p = (TH2F*)inFile_p->Get(mixName.c_str());
	TH2F* sub_p = (TH2F*)inFile_p->Get(subName.c_str());
	*/
	TH2F* mc_p = nullptr;
       
	TH2F* mixUncorrected_p = nullptr;
	TH2F* mixCorrection_p = nullptr;
	
	TH2F* pureBkgdMC_p = nullptr;
	TH2F* pureBkgd_p = nullptr;

	TH2F* mixBkgdMC_p = nullptr;
	TH2F* mixBkgd_p = nullptr;
	TH2F* mixBkgdUncorr_p = nullptr;
	  
	TH2F* rawEC_p = nullptr;
	TH2F* mixEC_p = nullptr;
	TH2F* subEC_p = nullptr;
	TH2F* mcEC_p = nullptr;

	TH2F* subMC_p = nullptr;
	TH2F* subMCEC_p = nullptr;

	TH2F* singleTruthToMultiFake_p = nullptr;
	TH2F* singleTruthToMultiFakeMix_p = nullptr;
	
	std::string mixUncorrectedName = mixName;
	std::string mixCorrectionName = mixName;
	std::string pureBkgdName = mixName;
	std::string pureBkgdMCName = mixName;
	std::string mixBkgdName = mixName;
	std::string mixBkgdUncorrName = mixName;
	std::string mixBkgdMCName = mixName;

	std::string singleTruthToMultiFakeName = mixName;
	std::string singleTruthToMultiFakeMixName = mixName;
	
	std::string mcName = rawName;
	std::string mcNameEC;
	if(isMC){
	  if(!strReplace(&mcName, "RAW", "RAWWITHTRUTHMATCH")) return 1;
	  mcNameEC = mcName;

	  if(isBarrelAndEC){
	    if(!strReplace(&mcNameEC, "Barrel", "EC")) return 1;
	  }
	}

	if(isMultijet){
	  if(!strReplace(&mixUncorrectedName, "MIXCORRECTED", "MIX")) return 1;
	  if(!strReplace(&mixCorrectionName, "MIXCORRECTED", "MIXCORRECTION")) return 1;
	  if(!strReplace(&pureBkgdName, "MIXCORRECTED", "MIXPURE")) return 1;
	  if(!strReplace(&pureBkgdMCName, "MIXCORRECTED", "PUREBKGD")) return 1;
	  if(!strReplace(&mixBkgdName, "MIXCORRECTED", "MIXMIXCORRECTED")) return 1;
	  if(!strReplace(&mixBkgdUncorrName, "MIXCORRECTED", "MIXMIX")) return 1;
	  if(!strReplace(&mixBkgdMCName, "MIXCORRECTED", "MIXBKGD")) return 1;

	  if(!strReplace(&singleTruthToMultiFakeName, "MIXCORRECTED", "SINGLETRUTHTOMULTIFAKE")) return 1;
	  if(!strReplace(&singleTruthToMultiFakeMixName, "MIXCORRECTED", "SINGLETRUTHTOMULTIFAKEMIX")) return 1;
	}

	//Declare an array and a maxbins value for rebinning purposes
	const Int_t nMaxBins = 1000;
	Int_t nBinsX, nBinsY;
	Double_t newBinsX[nMaxBins+1], newBinsY[nMaxBins+1];
	
	//Check the rebin
	if(doRebin){
	  std::string rebinXStr = plotConfig_p->GetValue((rebinStr + "REBINX").c_str(), "");
	  std::string rebinYStr = plotConfig_p->GetValue((rebinStr + "REBINY").c_str(), "");

	  if(rebinXStr.size() == 0 && doRebinX){
	    std::cout << "ERROR: '" << rebinStr << "DOREBIN' is " << doRebin << ", but rebinXStr is empty. return 1" << std::endl;
	    return 1;
	  }
	  if(rebinYStr.size() == 0 && doRebinY){
	    std::cout << "ERROR: '" << rebinStr << "DOREBIN' is " << doRebin << ", but rebinYStr is empty. return 1" << std::endl;
	    return 1;
	  }

	  //Grab a histogram for bin checks 
	  TH2F* temp_p = (TH2F*)inFile_p->Get(rawName.c_str());
	  std::vector<float> rebinXVect = strToVectF(rebinXStr);
	  std::vector<float> rebinYVect = strToVectF(rebinYStr);
	  Float_t deltaValue = 0.001;

	  //Pre-check the bins requested are unique at the level of 2xDeltaValue
	  if(doRebinX && !checkMinBinWidth(rebinXVect, 2*deltaValue)) return 1;
	  if(doRebinY && !checkMinBinWidth(rebinYVect, 2*deltaValue)) return 1;

	  //Check that the macro histogram bins align with requested new binning
	  if(doRebinX && !checkHistContainsBins(rebinXVect, temp_p, deltaValue, false)) return 1;
	  if(doRebinY && !checkHistContainsBins(rebinYVect, temp_p, deltaValue, true)) return 1;

	  if(doRebinX){
	    nBinsX = rebinXVect.size() - 1;
	    for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	      newBinsX[bIX] = rebinXVect[bIX];
	    }
	  }
	  else{
	    nBinsX = temp_p->GetXaxis()->GetNbins();
	    for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	      newBinsX[bIX] = temp_p->GetXaxis()->GetBinLowEdge(bIX+1);
	    }
	  }
	    
	  if(doRebinY){
	    nBinsY = rebinYVect.size() - 1;
	    for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
	      newBinsY[bIY] = rebinYVect[bIY];
	    }
	  }
	  else{
	    nBinsY = temp_p->GetYaxis()->GetNbins();
	    for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
	      newBinsY[bIY] = temp_p->GetYaxis()->GetBinLowEdge(bIY+1);
	    }
	  }

	
	  //CFM EDITING HERE 2021.11.04
	  //Double pointer since we are newing in the loop
	  std::vector<TH2F**> histsToRebin_p = {&raw_p, &mix_p, &sub_p};	  
	  std::vector<std::string> rebinHistsName = {rawName, mixName, subName};

	  if(mcFileName.size() != 0){
	    histsToRebin_p.push_back(&subMC_p);
	    rebinHistsName.push_back(subName);
	  }

	  /*
	  if(isBarrelAndEC){
	    histsToRebin_p.push_back(&rawEC_p);
	    histsToRebin_p.push_back(&mixEC_p);
	    histsToRebin_p.push_back(&subEC_p);

	    rebinHistsName.push_back(rawNameEC);
	    rebinHistsName.push_back(mixNameEC);
	    rebinHistsName.push_back(subNameEC);

	    if(mcFileName.size() != 0){
	      histsToRebin_p.push_back(&subMCEC_p);
	      rebinHistsName.push_back(subNameEC);
	    }
	  }
	  */
	  if(isMC){
	    histsToRebin_p.push_back(&mc_p);
	    rebinHistsName.push_back(mcName);

	    /*
	    if(isBarrelAndEC){
	      histsToRebin_p.push_back(&mcEC_p);
	      rebinHistsName.push_back(mcNameEC);
	    }
	    */
	  }
	  
	  if(isMultijet){
	    
	    histsToRebin_p.push_back(&mixUncorrected_p);
	    histsToRebin_p.push_back(&mixCorrection_p);
	    //	    histsToRebin_p.push_back(&pureBkgd_p);
	    //	    histsToRebin_p.push_back(&mixBkgd_p);
	    //	    histsToRebin_p.push_back(&mixBkgdUncorr_p);

	    rebinHistsName.push_back(mixUncorrectedName);
	    rebinHistsName.push_back(mixCorrectionName);
	    //	    rebinHistsName.push_back(pureBkgdName);
	    //	    rebinHistsName.push_back(mixBkgdName);
	    //	    rebinHistsName.push_back(mixBkgdUncorrName);
	    
	    
	    if(isMC){
	      //	      histsToRebin_p.push_back(&pureBkgdMC_p);
	      //	      histsToRebin_p.push_back(&mixBkgdMC_p);

	      //	      rebinHistsName.push_back(pureBkgdMCName);
	      //	      rebinHistsName.push_back(mixBkgdMCName);
	   
	      if(isStrSame(observables1[oI], "JtDRJJ")){		
		histsToRebin_p.push_back(&singleTruthToMultiFake_p);
		histsToRebin_p.push_back(&singleTruthToMultiFakeMix_p);

		rebinHistsName.push_back(singleTruthToMultiFakeName);
		rebinHistsName.push_back(singleTruthToMultiFakeMixName);
	      }	      
	    }
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	      
	  for(unsigned int i = 0; i < histsToRebin_p.size(); ++i){
	    std::string tempName = rebinHistsName[i];
	    TH2F* temp2_p = (TH2F*)inFile_p->Get(tempName.c_str());
	    if(!strReplace(&tempName, "_h", "_Rebin_h")) return 1;
	    while(tempName.find("/") != std::string::npos){tempName.replace(0, tempName.find("/")+1, "");}
	    std::string xTitle = temp2_p->GetXaxis()->GetTitle();
	    std::string yTitle = temp2_p->GetYaxis()->GetTitle();
	    std::string titleStr = ";" + xTitle + ";" + yTitle;
	    
	    (*histsToRebin_p[i]) = new TH2F(tempName.c_str(), titleStr.c_str(), nBinsX, newBinsX, nBinsY, newBinsY);	    
	    fineHistToCoarseHist(temp2_p, *(histsToRebin_p[i]));  
	  }
	}
	else{
	  raw_p = (TH2F*)inFile_p->Get(rawName.c_str());
	  mix_p = (TH2F*)inFile_p->Get(mixName.c_str());
	  sub_p = (TH2F*)inFile_p->Get(subName.c_str());
	  if(mcFileName.size() != 0) subMC_p = (TH2F*)mcFile_p->Get(subName.c_str());	  	  

	  /*	    
	  if(isBarrelAndEC){
	    rawEC_p = (TH2F*)inFile_p->Get(rawNameEC.c_str());
	    mixEC_p = (TH2F*)inFile_p->Get(mixNameEC.c_str());
	    subEC_p = (TH2F*)inFile_p->Get(subNameEC.c_str());

	    if(mcFileName.size() != 0) subMCEC_p = (TH2F*)mcFile_p->Get(subNameEC.c_str());
	  }
	  */
	  
	  if(isMC){
	    if(doGlobalDebug) std::cout << "FILE, LINE, mcNAME: " << __FILE__ << ", " << __LINE__ << ", " << mcName << std::endl;
	    mc_p = (TH2F*)inFile_p->Get(mcName.c_str());
	    //	    if(isBarrelAndEC) mcEC_p = (TH2F*)inFile_p->Get(mcNameEC.c_str());	      
	  }

	  if(isMultijet){	  
	    mixUncorrected_p = (TH2F*)inFile_p->Get(mixUncorrectedName.c_str());
	    mixCorrection_p = (TH2F*)inFile_p->Get(mixCorrectionName.c_str());
	    
	    //	    pureBkgd_p = (TH2F*)inFile_p->Get(pureBkgdName.c_str());
	    //	    mixBkgd_p = (TH2F*)inFile_p->Get(mixBkgdName.c_str());
	    //      mixBkgdUncorr_p = (TH2F*)inFile_p->Get(mixBkgdUncorrName.c_str());
	    
	    if(isMC){
	      //	      pureBkgdMC_p = (TH2F*)inFile_p->Get(pureBkgdMCName.c_str());	  
	      //	      mixBkgdMC_p = (TH2F*)inFile_p->Get(mixBkgdMCName.c_str());

	      if(isStrSame(observables1[oI], "JtDRJJ")){
		singleTruthToMultiFake_p = (TH2F*)inFile_p->Get(singleTruthToMultiFakeName.c_str());
		singleTruthToMultiFakeMix_p = (TH2F*)inFile_p->Get(singleTruthToMultiFakeMixName.c_str());
	      }
	    }
	  }	  
	}

	/*
	if(isBarrelAndEC){
	  raw_p->Add(rawEC_p);
	  mix_p->Add(mixEC_p);
	  sub_p->Add(subEC_p);
	  
	  if(mcFileName.size() != 0) subMC_p->Add(subMCEC_p);
	  if(isMC) mc_p->Add(mcEC_p);
	}
	*/
 		        
	const Int_t nBinsTemp = raw_p->GetXaxis()->GetNbins();
	if(doGlobalDebug) std::cout << "FILE, LINE, nBinsTemp: " << __FILE__ << ", " << __LINE__ << ", " << nBinsTemp << std::endl;

	if(nBinsTemp > nMaxBins){
	  std::cout << "Number of bins needed \'" << nBinsTemp << "\' exceeds max \'" << nMaxBins << "\'. return 1" << std::endl;
	  return 1;
	}
	Double_t binsTemp[nMaxBins+1];
	for(Int_t bIX = 0; bIX < nBinsTemp+1; ++bIX){
	  binsTemp[bIX] = raw_p->GetXaxis()->GetBinLowEdge(bIX+1);
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	for(Int_t bIY = 0; bIY < raw_p->GetYaxis()->GetNbins(); ++bIY){	
	  Float_t binCenter = raw_p->GetYaxis()->GetBinCenter(bIY+1);

	  if(binCenter < gammaPtBinsLowReco) continue;
	  if(binCenter >= gammaPtBinsHighReco) continue;

	  //Calculate photon integral of this bin
	  Float_t lowEdge = raw_p->GetYaxis()->GetBinLowEdge(bIY+1);
	  Float_t highEdge = raw_p->GetYaxis()->GetBinLowEdge(bIY+2);

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  if(doGlobalDebug) std::cout << "FILE, LINE, photonHistName: " << __FILE__ << ", " << __LINE__ << ", " << photonHist_p->GetName() << std::endl;

	  Float_t photonIntegral = 0.0;
	  for(Int_t pI = 0; pI < photonHist_p->GetXaxis()->GetNbins(); ++pI){
	    Float_t binCenter = photonHist_p->GetXaxis()->GetBinCenter(pI+1);
	    if(binCenter < lowEdge) continue;
	    if(binCenter >= highEdge) continue;

	    photonIntegral += photonHist_p->GetBinContent(pI+1);
	  }
	  Float_t photonIntegralMC = 0.0;
	  if(mcFileName.size() != 0){
	    for(Int_t pI = 0; pI < photonHistMC_p->GetXaxis()->GetNbins(); ++pI){
	      Float_t binCenter = photonHistMC_p->GetXaxis()->GetBinCenter(pI+1);
	      if(binCenter < lowEdge) continue;
	      if(binCenter >= highEdge) continue;
	      
	      photonIntegralMC += photonHistMC_p->GetBinContent(pI+1);
	    }	    
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  std::string axisStr = ";" + std::string(raw_p->GetXaxis()->GetTitle()) + ";Counts";
	  std::string newName = rawName.substr(0, rawName.rfind("_RAW"));
	  newName = newName + "_GammaPt" + std::to_string(bIY) + "_h";

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  TH1F* rawTemp_p = new TH1F(newName.c_str(), axisStr.c_str(), nBinsTemp, binsTemp);
	  TH1F* mixTemp_p = new TH1F("mixTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
	  TH1F* subTemp_p = new TH1F("subTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
	  TH1F* mcTemp_p = nullptr;
	  TH1F* mixUncorrectedTemp_p = nullptr;
	  TH1F* mixCorrectionTemp_p = nullptr;

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  TH1F* pureBkgdTemp_p = nullptr;
	  TH1F* pureBkgdMCTemp_p = nullptr;
	  TH1F* mixBkgdTemp_p = nullptr;
	  TH1F* mixBkgdUncorrTemp_p = nullptr;
	  TH1F* mixBkgdMCTemp_p = nullptr;

	  TH1F* singleTruthToMultiFakeTemp_p = nullptr;
	  TH1F* singleTruthToMultiFakeMixTemp_p = nullptr;
	  
	  TH1F* subMCTemp_p = nullptr;
	  if(mcFileName.size() != 0) subMCTemp_p = new TH1F("subMCTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);

	  if(isMC) mcTemp_p = new TH1F("mcTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
	  
	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  //	  if(isStrSame(observables1[oI], "JtXJJ") || isStrSame(observables1[oI], "JtDPhiJJG")){
	  if(isMultijet){
	    newName = mixUncorrectedName.substr(0, mixUncorrectedName.rfind("_MIX"));
	    newName = newName + "_Alt_GammaPt" + std::to_string(bIY) + "_h";

	    mixUncorrectedTemp_p = new TH1F(newName.c_str(), axisStr.c_str(), nBinsTemp, binsTemp);
 	    mixCorrectionTemp_p = new TH1F("mixCorrectionTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);

	    if(!strReplace(&newName, "Alt", "Alt2")) return 1;
	    //	    pureBkgdTemp_p = new TH1F(newName.c_str(), axisStr.c_str(), nBinsTemp, binsTemp);

	    if(!strReplace(&newName, "Alt2", "Alt3")) return 1;
	    //	    mixBkgdTemp_p = new TH1F(newName.c_str(), axisStr.c_str(), nBinsTemp, binsTemp);
	    //	    mixBkgdUncorrTemp_p = new TH1F("mixBkgdUncorrTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);

	    if(isMC){
	      //	      mixBkgdMCTemp_p = new TH1F("mixBkgdMCTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
	      //	      pureBkgdMCTemp_p = new TH1F("pureBkgdMCTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);

	      if(isStrSame(observables1[oI], "JtDRJJ")){
		singleTruthToMultiFakeTemp_p = new TH1F("singleTruthToMultiFakeTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
		singleTruthToMultiFakeMixTemp_p = new TH1F("singleTruthToMultiFakeMixTemp_h", axisStr.c_str(), nBinsTemp, binsTemp);
	      }
	    }	    
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  for(Int_t bIX = 0; bIX < raw_p->GetXaxis()->GetNbins(); ++bIX){	
	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    rawTemp_p->SetBinContent(bIX+1, raw_p->GetBinContent(bIX+1, bIY+1));
	    mixTemp_p->SetBinContent(bIX+1, mix_p->GetBinContent(bIX+1, bIY+1));
	    subTemp_p->SetBinContent(bIX+1, sub_p->GetBinContent(bIX+1, bIY+1));

	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    rawTemp_p->SetBinError(bIX+1, raw_p->GetBinError(bIX+1, bIY+1));
	    mixTemp_p->SetBinError(bIX+1, mix_p->GetBinError(bIX+1, bIY+1));
	    subTemp_p->SetBinError(bIX+1, sub_p->GetBinError(bIX+1, bIY+1));

	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    if(isMC){
	      mcTemp_p->SetBinContent(bIX+1, mc_p->GetBinContent(bIX+1, bIY+1));
	      mcTemp_p->SetBinError(bIX+1, mc_p->GetBinError(bIX+1, bIY+1));
	    }

	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    if(mcFileName.size() != 0){	
	      subMCTemp_p->SetBinContent(bIX+1, subMC_p->GetBinContent(bIX+1, bIY+1));
	      subMCTemp_p->SetBinError(bIX+1, subMC_p->GetBinError(bIX+1, bIY+1));      
	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    }
	    
	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    //	    if(isStrSame(observables1[oI], "JtXJJ") || isStrSame(observables1[oI], "JtDPhiJJG")){
	    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    if(isMultijet){

	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	      mixUncorrectedTemp_p->SetBinContent(bIX+1, mixUncorrected_p->GetBinContent(bIX+1, bIY+1));
              mixUncorrectedTemp_p->SetBinError(bIX+1, mixUncorrected_p->GetBinError(bIX+1, bIY+1));

	      mixCorrectionTemp_p->SetBinContent(bIX+1, mixCorrection_p->GetBinContent(bIX+1, bIY+1));
              mixCorrectionTemp_p->SetBinError(bIX+1, mixCorrection_p->GetBinError(bIX+1, bIY+1));

	      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	      //	      mixBkgdTemp_p->SetBinContent(bIX+1, mixBkgd_p->GetBinContent(bIX+1, bIY+1));
	      //	      mixBkgdTemp_p->SetBinError(bIX+1, mixBkgd_p->GetBinError(bIX+1, bIY+1));
	      
	      //	      mixBkgdUncorrTemp_p->SetBinContent(bIX+1, mixBkgdUncorr_p->GetBinContent(bIX+1, bIY+1));
	      //	      mixBkgdUncorrTemp_p->SetBinError(bIX+1, mixBkgdUncorr_p->GetBinError(bIX+1, bIY+1));
	      
	      //	      pureBkgdTemp_p->SetBinContent(bIX+1, pureBkgd_p->GetBinContent(bIX+1, bIY+1));
	      //	      pureBkgdTemp_p->SetBinError(bIX+1, pureBkgd_p->GetBinError(bIX+1, bIY+1));
	      
	      if(isMC){
		//		mixBkgdMCTemp_p->SetBinContent(bIX+1, mixBkgdMC_p->GetBinContent(bIX+1, bIY+1));
		//		mixBkgdMCTemp_p->SetBinError(bIX+1, mixBkgdMC_p->GetBinError(bIX+1, bIY+1));

		//		pureBkgdMCTemp_p->SetBinContent(bIX+1, pureBkgdMC_p->GetBinContent(bIX+1, bIY+1));
		//		pureBkgdMCTemp_p->SetBinError(bIX+1, pureBkgdMC_p->GetBinError(bIX+1, bIY+1));

		if(isStrSame(observables1[oI], "JtDRJJ")){
		  singleTruthToMultiFakeTemp_p->SetBinContent(bIX+1, singleTruthToMultiFake_p->GetBinContent(bIX+1, bIY+1));
		  singleTruthToMultiFakeTemp_p->SetBinError(bIX+1, singleTruthToMultiFake_p->GetBinError(bIX+1, bIY+1));
		  
		  singleTruthToMultiFakeMixTemp_p->SetBinContent(bIX+1, singleTruthToMultiFakeMix_p->GetBinContent(bIX+1, bIY+1));
		  singleTruthToMultiFakeMixTemp_p->SetBinError(bIX+1, singleTruthToMultiFakeMix_p->GetBinError(bIX+1, bIY+1));
		}
	      }
	    }
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  //Temp check of the remaining non-closure
	  if(isStrSame(observables1[oI], "JtDRJJ")){
	    if(isMC){
	      if(doTempCorr){
		mixTemp_p->Add(singleTruthToMultiFakeTemp_p, 1.0);
		mixTemp_p->Add(singleTruthToMultiFakeMixTemp_p, -1.0);
		
		subTemp_p->Add(singleTruthToMultiFakeTemp_p, -1.0);
		subTemp_p->Add(singleTruthToMultiFakeMixTemp_p, 1.0);
	      }
	    }	    
	  }
		      

	  
	  std::vector<TH1F*> hists_p = {rawTemp_p, subTemp_p};
	  std::vector<std::string> legStrs = {"Raw", "Raw - Mixed"};

	  int equipPos = 2;
	  TH1F* yjHist_p = nullptr;

	  
	  
	  if(isStrSame(observables1[oI], "JtDPhi") && bIY == 3 && jetR == 4){
	    //YJ properly normalized the histogram by the bin widths so I must as well for the comparison
	    /*
       	    for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
	      for(int bIX = 0; bIX < hists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
		Double_t val = hists_p[hI]->GetBinContent(bIX+1);
		Double_t err = hists_p[hI]->GetBinError(bIX+1);
		Double_t width = hists_p[hI]->GetBinWidth(bIX+1);

		hists_p[hI]->SetBinContent(bIX+1, val/width);
		hists_p[hI]->SetBinError(bIX+1, err/width);
	      }
	    }

	    hists_p[0]->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}#frac{dN_{J#gamma}}{d#Delta#phi_{J#gamma}}");
	    */
	    
	    yjHist_p = new TH1F("yjHist_p", "", 16, hists_p[0]->GetXaxis()->GetBinLowEdge(1), hists_p[0]->GetXaxis()->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1));
	    HIJet::Style::EquipHistogram(yjHist_p, equipPos);
	    fillYJHist(yjHist_p);
	    ++equipPos;
	    
	    hists_p.push_back(yjHist_p);
	    legStrs.push_back("YJ Raw - Mixed");
	  }

	  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
	    for(int bIX = 0; bIX < hists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
	      Double_t val = hists_p[hI]->GetBinContent(bIX+1);
	      Double_t err = hists_p[hI]->GetBinError(bIX+1);
	      Double_t width = hists_p[hI]->GetBinWidth(bIX+1);
	      
	      hists_p[hI]->SetBinContent(bIX+1, val/width);
	      hists_p[hI]->SetBinError(bIX+1, err/width);
	    }

	    hists_p[hI]->GetYaxis()->SetTitle(normYAxisTitle[oI].c_str());
	  }	  
	  
	  if(bIY == 5 && false){
	    std::cout << "PRINTING rawTemp_p" << std::endl;
	    std::cout << "SCALING PHOTON INTEGRAL '" << photonIntegral << "' (" << lowEdge << "-" << highEdge << ")" << std::endl;
	    std::cout << "PHOTON HIST: " << photonHist_p->GetName() << std::endl;
	    
	    rawTemp_p->Print("ALL");
	  }

	  HIJet::Style::EquipHistogram(mixTemp_p, equipPos);    
	  ++equipPos;
	  
	  std::vector<TH1F*> refHists_p = {mixTemp_p};
	  std::vector<std::string> refLegStrs = {"Mixed"};
	  std::vector<std::string> yLabels = {"", ("X/#color[" + std::to_string(mixTemp_p->GetMarkerColor()) + "]{Bkgd}").c_str()};
	  std::vector<bool> yLogs = {false, true};
	  if(isMC && mcTemp_p != nullptr){
	    HIJet::Style::EquipHistogram(mcTemp_p, equipPos);
	    ++equipPos;
	    refHists_p.push_back(mcTemp_p);
	    refLegStrs.push_back("Truth Matched");

	    yLabels.push_back("#frac{Reco.}{#color[" + std::to_string(mcTemp_p->GetMarkerColor()) + "]{Truth-Matched}}");
	    yLogs.push_back(false);
	  }

	  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
	    std::string histName = hists_p[hI]->GetName();
	    if(histName.find("yjHist_p") != std::string::npos) continue;
	    
	    hists_p[hI]->Scale(1.0/photonIntegral);
	  }

	  for(unsigned int hI = 0; hI < refHists_p.size(); ++hI){
	    refHists_p[hI]->Scale(1.0/photonIntegral);
	  }

	  if(!isMC && mcFileName.size() != 0){ //EDITING HERE FOR PHOTON INTEGRAL, PUSHBACK, AND SCALING
 	    subMCTemp_p->Scale(1.0/photonIntegralMC);

	    HIJet::Style::EquipHistogram(subMCTemp_p, equipPos);

	    
	    refHists_p.push_back(subMCTemp_p);
	    refLegStrs.push_back("Raw - Mixed MC");
	    yLabels.push_back("Data/#color[" + std::to_string(subMCTemp_p->GetMarkerColor()) + "]{MC}");	    
	  }
	  
	  if(yjHist_p != nullptr){
	    //YJ properly normalized the histogram by the bin widths so I must as well for the comparison
	    /*
       	    for(unsigned int hI = 0; hI < refHists_p.size(); ++hI){
	      for(int bIX = 0; bIX < refHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
		Double_t val = refHists_p[hI]->GetBinContent(bIX+1);
		Double_t err = refHists_p[hI]->GetBinError(bIX+1);
		Double_t width = refHists_p[hI]->GetBinWidth(bIX+1);
		refHists_p[hI]->SetBinContent(bIX+1, val/width);
		refHists_p[hI]->SetBinError(bIX+1, err/width);
	      }
	    }
	    */
	  }
		  

	  for(unsigned int hI = 0; hI < refHists_p.size(); ++hI){
	    for(int bIX = 0; bIX < refHists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
	      Double_t val = refHists_p[hI]->GetBinContent(bIX+1);
	      Double_t err = refHists_p[hI]->GetBinError(bIX+1);
	      Double_t width = refHists_p[hI]->GetBinWidth(bIX+1);
	      
	      refHists_p[hI]->SetBinContent(bIX+1, val/width);
	      refHists_p[hI]->SetBinError(bIX+1, err/width);
	    }
	    
	    refHists_p[hI]->GetYaxis()->SetTitle(normYAxisTitle[oI].c_str());
	  }

	  std::map<std::string, std::string> labelMapTemp = labelMap;
	  if(isBarrelAndEC){
	    labelMapTemp["Barrel"] = "Barrel+EC";
	    labelMapTemp["BarrelAndEC"] = "Barrel+EC";
	  }

	  /*
	  std::cout << "STARTING PLOTMIXCLOSURE (GammaPt: " << lowEdge << "-" << highEdge << ")" << std::endl;

	  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
	    std::cout << "VECTOR HIST " << hI << std::endl;
	    hists_p[hI]->Print("ALL");
	  }
	  */	  

	  std::string firstPlotName = plotMixClosure(doGlobalDebug, &labelMapTemp, plotConfig_p, strLowerToUpper(observables1[oI]), dateStr, globalSaveTag, hists_p, legStrs, refHists_p, refLegStrs, yLabels, yLogs, doReducedLabel);
	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  //	  std::cout << "ENDING PLOTMIXCLOSURE" << std::endl;
	  
	  if(yjHist_p != nullptr) delete yjHist_p;
	  //	  if(bIY == 5) return 1;
	  
	  
	  if((isStrSame(observables1[oI], "JtXJJ") || isStrSame(observables1[oI], "JtDPhiJJG")) && !isBarrelAndEC && false){
	    hists_p.clear();
	    legStrs.clear();
	    hists_p = {mixUncorrectedTemp_p, mixCorrectionTemp_p};
 	    legStrs = {"Uncorrected Mix", "Mix Correction", "Final Mix"};

	    for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
	      hists_p[hI]->Scale(1.0/photonIntegral);
	    }
	    hists_p.push_back(mixTemp_p);
	    
	    //plotMixClosure(doGlobalDebug, &labelMap, plotConfig_p, strLowerToUpper(observables1[oI]), dateStr, "MultiMixCheck", hists_p, legStrs);

	    if(isMC){
	      hists_p.clear();
	      legStrs.clear();
	      hists_p = {mixBkgdTemp_p, mixBkgdUncorrTemp_p};
	      legStrs = {"Corr. Sig.+Fake Bkgd.", "Uncorr. Sig.+Fake Bkgd."};
	      for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
		hists_p[hI]->Scale(1.0/photonIntegral);
	      }
	      mixBkgdMCTemp_p->Scale(1.0/photonIntegral);
	      HIJet::Style::EquipHistogram(mixBkgdMCTemp_p, 2);
	      
	      
	      plotMixClosure(doGlobalDebug, &labelMap, plotConfig_p, strLowerToUpper(observables1[oI]), dateStr, "MixBkgdCheck_" + globalSaveTag, hists_p, legStrs, {mixBkgdMCTemp_p}, {"Truth Sig.+Fake Bkgd."}, {"", ("Reco./#color[" + std::to_string(mixBkgdMCTemp_p->GetMarkerColor()) + "]{Truth}").c_str()}, {false, false}, doReducedLabel);
	    
	      hists_p.clear();
	      legStrs.clear();
	      hists_p = {pureBkgdTemp_p};
	      refHists_p = {pureBkgdMCTemp_p};
	      legStrs = {"Fake+Fake Bkgd."};
	      refLegStrs = {"Fake+Fake Truth Bkgd."};
	      HIJet::Style::EquipHistogram(pureBkgdMCTemp_p, 1);
	      for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
		hists_p[hI]->Scale(1.0/photonIntegral);
		refHists_p[hI]->Scale(1.0/photonIntegral);
	      }


	    
	      plotMixClosure(doGlobalDebug, &labelMap, plotConfig_p, strLowerToUpper(observables1[oI]), dateStr, "PureBkgdCheck_" + globalSaveTag, hists_p, legStrs, refHists_p, refLegStrs, {"", ("Reco./#color[" + std::to_string(pureBkgdMCTemp_p->GetMarkerColor()) + "]{Truth}").c_str()}, {false, false}, doReducedLabel);

	    }
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  //EDITING HERE - LETS MAKE THIS SUPER NICE AND AUTOMATED - N Labels -> XxY, Lengths well placed, etc.
	  //Making one of these for every single hist is dumb come up with a better system
	  if(doReducedLabel){
	    std::vector<std::string> preLabels = getLabels(plotConfig_p, rawTemp_p, &labelMapTemp);
	 
	    bool isInGlobal = false;
	    for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
	      std::vector<std::string> globalLabelTemp = globalLabels[gI];

	      if(preLabels.size() != globalLabelTemp.size()) continue;

	      bool allLabelsMatch = true;
	      for(unsigned int gI2 = 0; gI2 < globalLabelTemp.size(); ++gI2){		  

		bool isLabelFound = false;
		for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
		  if(isStrSame(globalLabelTemp[gI2], preLabels[pI])){
		    isLabelFound = true;
		    break;
		  }       	  
		}

		if(!isLabelFound){
		  allLabelsMatch = false;
		  break;
		}
	      }
	      
	      if(allLabelsMatch){
		isInGlobal = true;
		break;
	      }
	    }
	    if(!isInGlobal){
	      globalLabels.push_back(preLabels);
	    
	      std::vector<std::string> fourLabels = {"All", "NoPtGamma", "NoCent", "NoDeltaR"};
	      int titleFont = 42;
	      double titleSize = 0.1875;
	      const Int_t nMaxX = 4;
	      const Int_t nMaxY = 3;
	      
	      for(unsigned int tI = 0; tI < fourLabels.size(); ++tI){	      	    
		std::vector<std::string> tempLabels = preLabels;

		unsigned int iter = 0;
		while(tempLabels.size() > iter){
		  if(tempLabels[iter].find("ATLAS") != std::string::npos) tempLabels.erase(tempLabels.begin()+iter);
		  else if(tempLabels[iter].find("PYT") != std::string::npos) tempLabels.erase(tempLabels.begin()+iter);
		  else if(tempLabels[iter].find("Pb+") != std::string::npos) tempLabels.erase(tempLabels.begin()+iter);
		  else ++iter;
		}
		
		if(fourLabels[tI].find("All") == std::string::npos){
		  for(unsigned int fI = 0; fI < tempLabels.size(); ++fI){
		    if(fourLabels[tI].find("NoPtGamma") != std::string::npos && tempLabels[fI].find("p_{T,#gamma}") != std::string::npos){
		      tempLabels.erase(tempLabels.begin()+fI);
		      break;
		    }
		    
		    if(fourLabels[tI].find("NoCent") != std::string::npos && tempLabels[fI].find("%") != std::string::npos){
		      tempLabels.erase(tempLabels.begin()+fI);
		      break;
		    }
		    
		    if(fourLabels[tI].find("NoDeltaR") != std::string::npos && tempLabels[fI].find("#DeltaR") != std::string::npos){
		      tempLabels.erase(tempLabels.begin()+fI);
		      break;
		    }
		    
		  }
		}
		
		Int_t nLabels = tempLabels.size();
		
		Double_t xWidth[nMaxX];
		Double_t yHeight = -1.0;
		
		for(Int_t xI = 0; xI < nMaxX; ++xI){
		  xWidth[xI] = -1.0;
		}
		
		Int_t nX = 3;
		Int_t nY = 2;
		if(nLabels >= 7) nY = 3;
		if(nLabels >= 10) nX = 4;
		
		if(nX > nMaxX){
		  std::cout << "nX \'" << nX << "\' for total labels \'" << tempLabels.size() << "\' exceeds maxX \'" << nMaxX << "\'. return 1" << std::endl;
		  return 1;
		}
		if(nY > nMaxY){
		  std::cout << "nY \'" << nY << "\' for total labels \'" << tempLabels.size() << "\' exceeds maxY \'" << nMaxY << "\'. return 1" << std::endl;
		  return 1;
		}
		
		const Int_t width = 700 + 475*(nX-1);
		const Int_t height = 120*nY;
		
		TCanvas* labelCanv_p = new TCanvas("labelCanv_p", "", width, height);
		labelCanv_p->SetTopMargin(0.001);
		labelCanv_p->SetLeftMargin(0.001);
		labelCanv_p->SetRightMargin(0.001);
		labelCanv_p->SetBottomMargin(0.001);
		
		labelCanv_p->cd();
		
		TLatex* label_p = new TLatex();
		label_p->SetNDC();
		label_p->SetTextFont(titleFont);
		label_p->SetTextSize(titleSize);	    
		
		for(unsigned int pI = 0; pI < tempLabels.size(); ++pI){
		  label_p->SetText(0.01, 0.01, tempLabels[pI].c_str());
		  
		  Double_t tempX = label_p->GetXsize();
		  Double_t tempY = label_p->GetYsize();
		  
		  if(tempY > yHeight) yHeight = tempY;
		  Int_t xPos = pI/nY;
		  
		  if(tempX > xWidth[xPos]) xWidth[xPos] = tempX;
		}
		
		for(unsigned int eI = 0; eI < tempLabels.size(); ++eI){
		  if(tempLabels[eI].find("Alt") != std::string::npos){
		    tempLabels.erase(tempLabels.begin() + eI);
		    break;
		  }
		}
		
		Double_t xStart = 0.01;
		
		bool doBreak = false;
		for(Int_t xI = 0; xI < nX; ++xI){
		  Double_t yStart = 0.8;
		  
		  for(Int_t yI = 0; yI < nY; ++yI){
		    unsigned int labelPos = xI*nY + yI;
		    
		    if(labelPos >= tempLabels.size()){
		      doBreak = true;
		      break;
		    }
		    
		    label_p->DrawLatex(xStart, yStart - yI*yHeight, tempLabels[labelPos].c_str());
		    
		  }
		  if(doBreak) break;
		  
		  xStart += 0.04 + xWidth[xI];
		}
				
		std::string saveName = firstPlotName;
		saveName.replace(saveName.find(dateStr + "."), dateStr.size()+1, fourLabels[tI] + "_" + dateStr + ".");
		//		std::cout << "SAVENAME: " << saveName << std::endl;
		quietSaveAs(labelCanv_p, saveName);
		delete labelCanv_p;
	      }
	    }
	  }
	  
	  delete rawTemp_p;
	  delete mixTemp_p;
	  delete subTemp_p;	  
	  if(isMC) delete mcTemp_p;	  
	  if(mcFileName.size() != 0) delete subMCTemp_p;

	  if(isMultijet){
	    delete mixUncorrectedTemp_p;
	    delete mixCorrectionTemp_p;
	    
	    //	    delete mixBkgdTemp_p;
	    //	    delete mixBkgdUncorrTemp_p;
	    //	    delete pureBkgdTemp_p;
    
	    if(isMC){
	      //	      delete mixBkgdMCTemp_p;
	      //	      delete pureBkgdMCTemp_p;
	      //Quickly we are gonna dump to a test canvas
	      if(isStrSame(observables1[oI], "JtDRJJ")){
		//Define some margins
		const Double_t padSplit = 0.5;
		const Double_t height = 1200.0;
		const Double_t width = 942.831;

		const Double_t noMargin = 0.001;
		const Double_t leftMargin = 0.20;
		const Double_t topMargin = 0.02;
		const Double_t rightMargin = topMargin;
		const Double_t bottomMargin = leftMargin;

		const int titleFont = 42;
		const double titleSize = 0.035;
		const double labelSize = titleSize*0.9;

		const Double_t yOffset = 0.9;
		
		TCanvas* canv_p = new TCanvas("testCanv_p", "", height, width);
		canv_p->SetTopMargin(noMargin);
		canv_p->SetLeftMargin(noMargin);
		canv_p->SetRightMargin(noMargin);
		canv_p->SetBottomMargin(noMargin);
		canv_p->cd();

		TPad* pads_p[2];
		pads_p[0] = new TPad("pads0", "", 0.0, padSplit, 1.0, 1.0);
		pads_p[0]->Draw("SAME");
		pads_p[0]->SetBottomMargin(noMargin);
		pads_p[0]->SetTopMargin(topMargin);
		pads_p[0]->SetRightMargin(rightMargin);
		pads_p[0]->SetLeftMargin(leftMargin);
		canv_p->cd();

		pads_p[1] = new TPad("pads1", "", 0.0, 0.0, 1.0, padSplit);
		pads_p[1]->SetTopMargin(noMargin);
		pads_p[1]->SetRightMargin(rightMargin);
		pads_p[1]->SetBottomMargin(bottomMargin);
		pads_p[1]->SetLeftMargin(leftMargin);
		pads_p[1]->Draw("SAME");

		canv_p->cd();
		pads_p[0]->cd();
		
		for(int bIX = 0; bIX < singleTruthToMultiFakeTemp_p->GetXaxis()->GetNbins(); ++bIX){
		  Double_t val = singleTruthToMultiFakeTemp_p->GetBinContent(bIX+1);
		  Double_t err = singleTruthToMultiFakeTemp_p->GetBinError(bIX+1);
		  Double_t width = singleTruthToMultiFakeTemp_p->GetBinWidth(bIX+1);
		  
		  singleTruthToMultiFakeTemp_p->SetBinContent(bIX+1, val/width);
		  singleTruthToMultiFakeTemp_p->SetBinError(bIX+1, err/width);

		  val = singleTruthToMultiFakeMixTemp_p->GetBinContent(bIX+1);
		  err = singleTruthToMultiFakeMixTemp_p->GetBinError(bIX+1);

		  singleTruthToMultiFakeMixTemp_p->SetBinContent(bIX+1, val/width);
		  singleTruthToMultiFakeMixTemp_p->SetBinError(bIX+1, err/width);
		}
		  
		singleTruthToMultiFakeTemp_p->Scale(1.0/photonIntegral);
		singleTruthToMultiFakeMixTemp_p->Scale(1.0/photonIntegral);
		
		HIJet::Style::EquipHistogram(singleTruthToMultiFakeTemp_p, 0);
		HIJet::Style::EquipHistogram(singleTruthToMultiFakeMixTemp_p, 2);

		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetTitleFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetLabelFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetLabelFont(titleFont);
		
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleSize(titleSize/(1.0 - padSplit));
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetTitleSize(0.000001);
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetLabelSize(labelSize/(1.0 - padSplit));
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetLabelSize(0.000001);
		
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitle(normYAxisTitle[oI].c_str());
		singleTruthToMultiFakeTemp_p->SetMaximum(singleTruthToMultiFakeTemp_p->GetMaximum()*1.5);

		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleOffset(yOffset);
		
		singleTruthToMultiFakeTemp_p->DrawCopy("HIST E1 P");
		singleTruthToMultiFakeMixTemp_p->DrawCopy("HIST E1 P SAME");

 		TLegend* leg_p = new TLegend(leftMargin + 0.02, 0.85, 0.6, 0.95);
		leg_p->SetTextFont(titleFont);	       
		leg_p->SetTextSize(titleSize/padSplit);
		leg_p->SetBorderSize(0);
		leg_p->SetFillColor(0);
		leg_p->SetFillStyle(0);

		leg_p->AddEntry(singleTruthToMultiFakeTemp_p, "Multijet w/ 1 Truth-match", "P L");
		leg_p->AddEntry(singleTruthToMultiFakeMixTemp_p, "Multijet w/ 1 Truth-match mix", "P L");
		leg_p->Draw("SAME");
		
		singleTruthToMultiFakeTemp_p->Add(singleTruthToMultiFakeMixTemp_p, -1.0);
		canv_p->cd();
		pads_p[1]->cd();
		HIJet::Style::EquipHistogram(singleTruthToMultiFakeTemp_p, 1);

		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitle("Raw - Mix");
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetTitleFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetLabelFont(titleFont);
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetLabelFont(titleFont);

		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetTitleSize(titleSize/padSplit);
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));
		singleTruthToMultiFakeTemp_p->GetXaxis()->SetLabelSize(labelSize/padSplit);
		
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetTitleOffset(yOffset*(1.0-padSplit)/padSplit);
		singleTruthToMultiFakeTemp_p->GetYaxis()->SetNdivisions(404);
		singleTruthToMultiFakeTemp_p->DrawCopy("HIST E1 P");
		
		HIJet::Style::EquipHistogram(singleTruthToMultiFakeTemp_p, 0);
				
		quietSaveAs(canv_p, "pdfDir/" + dateStr + "/quickTest_Cent" + std::to_string(cI) + "_" + barrelECStr[bI] + "_BIY" + std::to_string(bIY) + "_Obs" + observables1[oI] + ".png");
	      		
		delete pads_p[0];
 		delete pads_p[1];		
		delete canv_p;
		delete singleTruthToMultiFakeTemp_p;
		delete singleTruthToMultiFakeMixTemp_p;
	      }
	    }
	  }	
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	//If we did a rebin we new'd and therefore need to delete
	if(doRebin){
	  delete raw_p;
	  delete mix_p;
	  delete sub_p;

	  if(mcFileName.size() != 0) delete subMC_p;

	  /*
	  if(isBarrelAndEC){
	    delete rawEC_p;
	    delete mixEC_p;
	    delete subEC_p;

	    if(mcFileName.size() != 0) delete subMCEC_p;	   
	  }
	  */
	  
	  if(isMC){
	    delete mc_p;
	    //  if(isBarrelAndEC) delete mcEC_p;	    
	  }

	  if(isMultijet){
	    delete mixUncorrected_p;
	    delete mixCorrection_p;
	    //	    delete pureBkgd_p;
	    //	    delete mixBkgd_p;
	    //	    delete mixBkgdUncorr_p;

	    if(isMC){
	      //	      delete pureBkgdMC_p;
	      //	      delete mixBkgdMC_p;

	      delete singleTruthToMultiFake_p;
	      delete singleTruthToMultiFakeMix_p;
	    }
	  }

	  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	}
      
      }
    }
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(mcFileName.size() != 0){  
    mcFile_p->Close();
    delete mcFile_p;
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
