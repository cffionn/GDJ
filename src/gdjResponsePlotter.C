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
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

const Double_t titleSizeX = 0.04;
const Double_t titleSizeY = 0.04;
const Double_t labelSizeX = 0.035;
const Double_t labelSizeY = 0.035;
const Int_t titleFont = 42;

const Double_t padSplit = 0.35;
const Double_t leftMargin = 0.15;
const Double_t rightMargin = 0.12;
const Double_t bottomMargin = leftMargin;
const Double_t tinyMargin = 0.0001;
const Double_t otherMargin = 0.02;

void plotMeanAndSigma(const bool doGlobalDebug, TEnv* plotConfig_p, std::map<std::string, std::string>* labelMap, std::string dateStr, std::vector<std::string> centBins, TH1F* histMean_p[], TH1F* histWidth_p[])
{
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::string overMeanStr = histWidth_p[0]->GetName();
  std::string labelStr = histWidth_p[0]->GetName();

  
  if(overMeanStr.find("OverMean") != std::string::npos) overMeanStr = "OverMean";
  else overMeanStr = "";

  labelStr.replace(0, labelStr.find("_")+1, "");
  labelStr.replace(labelStr.rfind("_h"), 2, "");
  std::vector<std::string> preLabels;
  for(Int_t gI = 0; gI < 20; ++gI){
    std::string globalLabel = plotConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(globalLabel.size() == 0) continue;    
    preLabels.push_back(globalLabel);
  }

  const int jetR = plotConfig_p->GetValue("JETR", 10);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return" << std::endl;
    return;
  }
  preLabels.push_back("anti-k_{T} R=0." + std::to_string(jetR) + " jets");

  const Double_t widthMax = plotConfig_p->GetValue("ERESWIDTHMAX", 1.0);
  const Double_t widthMin = plotConfig_p->GetValue("ERESWIDTHMIN", 0.0);
  const Double_t meanMax = plotConfig_p->GetValue("ERESMEANMAX", 1.5);
  const Double_t meanMin = plotConfig_p->GetValue("ERESMEANMIN", 0.5);

  const Bool_t doLogX = plotConfig_p->GetValue("ERESDOLOGX", 0);
  
  const Double_t labelX = plotConfig_p->GetValue("ERESLABELX", 0.5);
  const Double_t labelY = plotConfig_p->GetValue("ERESLABELY", 0.5);
  const Bool_t labelAlignRight = plotConfig_p->GetValue("ERESLABELALIGNRIGHT", 0);

  const Double_t legX = plotConfig_p->GetValue("ERESLEGX", 0.5);
  const Double_t legY = plotConfig_p->GetValue("ERESLEGY", 0.5);

  std::vector<std::vector<float> > generalBoxes;
  for(int gI = 0; gI < 10; ++gI){
    std::string tempStr = plotConfig_p->GetValue(("ERESBOX." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    generalBoxes.push_back(strToVectF(tempStr));    
  }
     
  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    if(centStr.find("Cent") == std::string::npos) preLabels.push_back(centStr);
    labelStr.replace(0, labelStr.find("_")+1, "");
  }
  if(labelStr.size() != 0 && labelStr.find("Cent") == std::string::npos) preLabels.push_back(labelStr);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
  const Double_t yOffset = 1.0;
  const Double_t xOffset = 1.0;
  
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSizeX/(1.0-padSplit));

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  TPad* pads_p[2];
  canv_p->SetTopMargin(tinyMargin);
  canv_p->SetRightMargin(tinyMargin);
  canv_p->SetLeftMargin(tinyMargin);
  canv_p->SetBottomMargin(tinyMargin);
  canv_p->cd();
  
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(otherMargin/(1.0-padSplit));
  pads_p[0]->SetRightMargin(otherMargin);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(tinyMargin);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetTopMargin(tinyMargin);
  pads_p[1]->SetRightMargin(otherMargin);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetBottomMargin(bottomMargin/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  bool isPP = centBins[0].find("PP") != std::string::npos;
  std::string inFileNamePP = plotConfig_p->GetValue("INFILENAMEPP", "");
  const bool addPP = inFileNamePP.size() != 0;
  int nCentBins = 1;
  if(!isPP){
    nCentBins = centBins.size()-1;
    if(addPP) ++nCentBins;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    if(labelMap->count(preLabels[pI]) != 0) preLabels[pI] = (*labelMap)[preLabels[pI]];
  }    
  double maxLabelXSize = -0.1;
  if(labelAlignRight){
    for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
      label_p->SetText(0.1, 0.1, preLabels[pI].c_str());

      if(label_p->GetXsize() > maxLabelXSize) maxLabelXSize = label_p->GetXsize();
    }

    for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
      label_p->SetText(0.1, 0.1, preLabels[pI].c_str());

      while(label_p->GetXsize() < maxLabelXSize){
	preLabels[pI] = " " + preLabels[pI];
	label_p->SetText(0.1, 0.1, preLabels[pI].c_str());
      }    
    }
  }
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  double maxLegXSize = -0.1;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string legStr = "PP";
    if(!isPP){
      if(addPP && cI == nCentBins-1) legStr = "p+p";
      else legStr = "Pb+Pb, " + centBins[cI] + "-" + centBins[cI+1] + "%";     
    }

    label_p->SetText(0.1, 0.1, legStr.c_str());
    if(label_p->GetXsize() > maxLegXSize) maxLegXSize = label_p->GetXsize();    
  }      

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TLegend* leg_p = new TLegend(legX, legY - nCentBins*0.06, legX + maxLegXSize, legY);
  leg_p->SetTextFont(titleFont);
  leg_p->SetTextSize(titleSizeX/(1.0-padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TF1* csnPP_p = nullptr;
  TF1* csnFit_p = nullptr;
  
  if(addPP){
    csnPP_p = new TF1("csnPP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));

    csnPP_p->SetParameter(0, 0.05);
    csnPP_p->SetParameter(1, 1.0);
    
    histWidth_p[nCentBins-1]->Fit("csnPP_p", "Q M N E I", "", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));    
  }
 
  
  if(isPP) csnFit_p = new TF1("csnFit_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
  else{
    if(addPP) csnFit_p = new TF1("csnFit_p", ("TMath::Sqrt(" + std::to_string(csnPP_p->GetParameter(0)*csnPP_p->GetParameter(0)) + " + " + std::to_string(csnPP_p->GetParameter(1)*csnPP_p->GetParameter(1)) + "/x + [0]*[0]/(x*x))").c_str(), histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
    else csnFit_p = new TF1("csnFit_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
  }
  
  

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string legStr = "PP";
    if(!isPP){
      if(addPP && cI == nCentBins-1) legStr = "p+p";
      else legStr = "Pb+Pb, " + centBins[cI] + "-" + centBins[cI+1] + "%";     
    }
    
    canv_p->cd();
    pads_p[0]->cd();

    HIJet::Style::EquipHistogram(histWidth_p[cI], cI);
    
    histWidth_p[cI]->SetMaximum(widthMax);
    histWidth_p[cI]->SetMinimum(widthMin);

    histWidth_p[cI]->GetXaxis()->SetTitleFont(titleFont);
    histWidth_p[cI]->GetYaxis()->SetTitleFont(titleFont);
    histWidth_p[cI]->GetXaxis()->SetTitleSize(titleSizeX/(1.0 - padSplit));
    histWidth_p[cI]->GetYaxis()->SetTitleSize(titleSizeY/(1.0 - padSplit));

    histWidth_p[cI]->GetXaxis()->SetLabelFont(titleFont);
    histWidth_p[cI]->GetYaxis()->SetLabelFont(titleFont);
    histWidth_p[cI]->GetXaxis()->SetLabelSize(labelSizeX/(1.0 - padSplit));
    histWidth_p[cI]->GetYaxis()->SetLabelSize(labelSizeY/(1.0 - padSplit));

    histWidth_p[cI]->GetYaxis()->SetNdivisions(505);
    histWidth_p[cI]->GetXaxis()->SetNdivisions(505);

    histWidth_p[cI]->GetYaxis()->SetTitleOffset(yOffset);

    if(cI == 0) histWidth_p[cI]->DrawCopy("HIST E1 P");
    else histWidth_p[cI]->DrawCopy("HIST E1 P SAME");


    //Fit sequence
    if(isPP){		
      csnFit_p->SetParameter(0, 0.05);
      csnFit_p->SetParameter(1, 1.0);
	
      csnFit_p->SetParLimits(0, 0.0, 1.0);
      csnFit_p->SetParLimits(1, 0.0, 10.0);
	
      histWidth_p[cI]->Fit("csnFit_p", "Q M N E I", "", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
      
      HIJet::Style::EquipTF1(csnFit_p, cI);
      csnFit_p->SetLineStyle(2);    
      
      csnFit_p->DrawCopy("SAME");
    }
    else{
      if(addPP){
	if(cI == nCentBins-1){
	  HIJet::Style::EquipTF1(csnPP_p, cI);
	  csnPP_p->SetLineStyle(2);
	  
	  csnPP_p->DrawCopy("SAME");
	}
	else{
	  csnFit_p->SetParameter(0, 0.00);
	
	  csnFit_p->SetParLimits(0, 0.0, 100.0);
	
	  histWidth_p[cI]->Fit("csnFit_p", "M N E I", "", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
	  HIJet::Style::EquipTF1(csnFit_p, cI);
	  csnFit_p->SetLineStyle(2);    
	  
	  csnFit_p->DrawCopy("SAME"); 
	}	
      }
      else{
	csnFit_p->SetParameter(0, 0.05);
	csnFit_p->SetParameter(1, 1.0);
	csnFit_p->SetParameter(2, 0.00);
	
	csnFit_p->SetParLimits(0, 0.0, 1.0);
	csnFit_p->SetParLimits(1, 0.0, 10.0);
	csnFit_p->SetParLimits(2, 0.0, 100.0);
	
	histWidth_p[cI]->Fit("csnFit_p", "Q M N E I", "", histWidth_p[0]->GetBinLowEdge(1), histWidth_p[0]->GetBinLowEdge(histWidth_p[0]->GetXaxis()->GetNbins()+1));
	
	HIJet::Style::EquipTF1(csnFit_p, cI);
	csnFit_p->SetLineStyle(2);    
	
	csnFit_p->DrawCopy("SAME");
      }
    }
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    /*
    legStr = legStr + "; N=" + prettyString(csnFit_p->GetParameter(2),1,false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    legStr = legStr + ",S=" + prettyString(csnFit_p->GetParameter(1), 2, false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    std::cout << "CSN: " << csnFit_p->GetParameter(0) << std::endl;
    legStr = legStr + ",C=" + prettyString(csnFit_p->GetParameter(0),3,false);
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    */
  
    if(addPP){
      if(cI == nCentBins-1) legStr = legStr + "; C=" + prettyString(csnPP_p->GetParameter(0), 3, false) + "; S=" + prettyString(csnPP_p->GetParameter(1), 2, false);
      else legStr = legStr + "; N=" + prettyString(csnFit_p->GetParameter(0), 1, false);
    }

    leg_p->AddEntry(histMean_p[cI], legStr.c_str(), "P L");
      
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    canv_p->cd();
    pads_p[1]->cd();

    HIJet::Style::EquipHistogram(histMean_p[cI], cI);
    
    histMean_p[cI]->SetMaximum(meanMax);
    histMean_p[cI]->SetMinimum(meanMin);

    histMean_p[cI]->GetXaxis()->SetTitleFont(titleFont);
    histMean_p[cI]->GetYaxis()->SetTitleFont(titleFont);
    histMean_p[cI]->GetXaxis()->SetTitleSize(titleSizeX/(padSplit));
    histMean_p[cI]->GetYaxis()->SetTitleSize(titleSizeY/(padSplit));

    histMean_p[cI]->GetXaxis()->SetLabelFont(titleFont);
    histMean_p[cI]->GetYaxis()->SetLabelFont(titleFont);
    histMean_p[cI]->GetXaxis()->SetLabelSize(labelSizeX/(padSplit));
    histMean_p[cI]->GetYaxis()->SetLabelSize(labelSizeY/(padSplit));

    histMean_p[cI]->GetYaxis()->SetNdivisions(505);
    histMean_p[cI]->GetXaxis()->SetNdivisions(505);

    histMean_p[cI]->GetXaxis()->SetTitleOffset(xOffset);
    histMean_p[cI]->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0-padSplit));
    
    if(cI == 0) histMean_p[cI]->DrawCopy("HIST E1 P");
    else histMean_p[cI]->DrawCopy("HIST E1 P SAME");    
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  canv_p->cd();
  pads_p[0]->cd();
  leg_p->Draw("SAME");
  if(doLogX) pads_p[0]->SetLogx();
  
  labelStr = "";

  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    label_p->DrawLatex(labelX, labelY - pI*0.085, preLabels[pI].c_str());
  }
  
  
  canv_p->cd();
  pads_p[1]->cd();
  if(doLogX) pads_p[1]->SetLogx();
  

  line_p->DrawLine(histMean_p[0]->GetBinLowEdge(1), 1.0, histMean_p[0]->GetBinLowEdge(histMean_p[0]->GetXaxis()->GetNbins()+1), 1.0);
  
  canv_p->cd();
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int gI = 0; gI < generalBoxes.size(); ++gI){
    canv_p->cd();
    drawWhiteBoxNDC(canv_p, generalBoxes[gI][0], generalBoxes[gI][1], generalBoxes[gI][2], generalBoxes[gI][3]);
    generalBoxes[gI].clear();
  }
  generalBoxes.clear();
    
  std::string saveName = "pdfDir/" + dateStr + "/meanAndSigma" + overMeanStr + "_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  delete csnFit_p;
  if(addPP) delete csnPP_p;
  
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

int gdjResponsePlotter(const std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"ERESWIDTHMAX",
					"ERESWIDTHMIN",
					"ERESMEANMAX",
					"ERESMEANMIN",
					"ERESDOLOGX",
					"ERESLABELX",
					"ERESLABELY",
					"ERESLABELALIGNRIGHT",
					"ERESLEGX",
					"ERESLEGY"};
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string inFileNamePP = inConfig_p->GetValue("INFILENAMEPP", "");
  const bool addPP = inFileNamePP.size() != 0;
  if(addPP){
    if(!check.checkFileExt(inFileNamePP, ".root")) return 1;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> globalLabels;
  for(Int_t gI = 0; gI < 20; ++gI){
    std::string globalLabel = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(globalLabel.size() == 0) continue;    
    globalLabels.push_back(globalLabel);
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TFile* inFilePP_p = nullptr;
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");
  TEnv* fileLabel_p = (TEnv*)inFile_p->Get("label");
  configParser labels(fileLabel_p);
  configParser configs(fileConfig_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();
  std::map<std::string, std::string> configMap = configs.GetConfigMap();
  
  std::vector<std::string> fileReqParams = {"ISMC",
					    "NJTPTBINS",
					    "JTPTBINSDOLOG",
					    "JTPTBINSLOW",
					    "JTPTBINSHIGH",
					    "NJTETABINSSUB",
					    "JTETABINSSUBLOW",
					    "JTETABINSSUBHIGH",
					    "NGAMMAPTBINSSUB",
					    "ISPP",
					    "JETR",
					    "RECOJTPTMIN"};  

  std::vector<std::string> fileMatchParams = {"ISMC",
					      "NJTPTBINS",
					      "JTPTBINSDOLOG",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "NGAMMAPTBINSSUB",
					      "JETR",
					      "RECOJTPTMIN"};  
  
  std::vector<std::string> pbpbParams = {"CENTBINS"};

  if(!checkEnvForParams(fileConfig_p, fileReqParams)) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(addPP){
    inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
    TEnv* fileConfigPP_p = (TEnv*)inFilePP_p->Get("config");

    if(!checkEnvForParams(fileConfigPP_p, fileReqParams)) return 1;
    if(!compEnvParams(fileConfig_p, fileConfigPP_p, fileMatchParams)) return 1;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  inConfig_p->SetValue("JETR", fileConfig_p->GetValue("JETR", ""));
  
  const bool isMC = fileConfig_p->GetValue("ISMC", 0);
  if(!isMC){   
    std::cout << "gdjResponsePlotter ERROR - This macro only takes MC input, but config gives \'isMC\'=" << isMC << ". return 1" << std::endl;
    return 1;
  }

  const int nMaxJtPtBins = 200;
  const int nJtPtBins = fileConfig_p->GetValue("NJTPTBINS", 0);
  if(nJtPtBins > nMaxJtPtBins){
    std::cout << "Given nJtPtBins \'" << nJtPtBins << "\' is greater than max value \'" << nMaxJtPtBins << "\'. return 1" << std::endl;
    return 1;
  }
  const bool jtPtBinsDoLog = fileConfig_p->GetValue("JTPTBINSDOLOG", 0);
  const Float_t jtPtBinsLow = fileConfig_p->GetValue("JTPTBINSLOW", 10.0);
  const Float_t jtPtBinsHigh = fileConfig_p->GetValue("JTPTBINSHIGH", 100.0);
  Double_t jtPtBins[nMaxJtPtBins+1];
  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  
  const int nJtEtaBinsSub = fileConfig_p->GetValue("NJTETABINSSUB", 0);
  if(nJtEtaBinsSub > nMaxJtPtBins){
    std::cout << "Given nJtEtaBinsSub \'" << nJtEtaBinsSub << "\' is greater than max value \'" << nMaxJtPtBins << "\'. return 1" << std::endl;
    return 1;
  }
  const Float_t jtEtaBinsSubLow = fileConfig_p->GetValue("JTETABINSSUBLOW", 10.0);
  const Float_t jtEtaBinsSubHigh = fileConfig_p->GetValue("JTETABINSSUBHIGH", 100.0);
  
  const Bool_t jtEtaBinsSubDoAbs = fileConfig_p->GetValue("JTETABINSSUBDOABS", 0);
  const Bool_t jtEtaBinsSubDoLog = fileConfig_p->GetValue("JTETABINSSUBDOLOG", 0);
  const Bool_t jtEtaBinsSubDoLin = fileConfig_p->GetValue("JTETABINSSUBDOLIN", 0);
  const Bool_t jtEtaBinsSubDoCustom = fileConfig_p->GetValue("JTETABINSSUBDOCUSTOM", 0);
  Double_t jtEtaBinsSub[nMaxJtPtBins+1];
  if(jtEtaBinsSubDoLin) getLinBins(jtEtaBinsSubLow, jtEtaBinsSubHigh, nJtEtaBinsSub, jtEtaBinsSub);
  else if(jtEtaBinsSubDoLog) getLogBins(jtEtaBinsSubLow, jtEtaBinsSubHigh, nJtEtaBinsSub, jtEtaBinsSub);
  else if(jtEtaBinsSubDoCustom){
    std::vector<float> jtEtaBinsCustom = strToVectF(fileConfig_p->GetValue("JTETABINSSUBCUSTOM", ""));

    for(unsigned int jI = 0; jI < jtEtaBinsCustom.size(); ++jI){
      jtEtaBinsSub[jI] = jtEtaBinsCustom[jI];
    }   
  }
  

  const int nGammaPtBinsSub = fileConfig_p->GetValue("NGAMMAPTBINSSUB", 0);
  const bool isPP = fileConfig_p->GetValue("ISPP", 0);
  const Float_t recoJtPtMin = fileConfig_p->GetValue("RECOJTPTMIN", -100);
  
  const int nMaxCentBins = 20;
  TH1F* recoOverGenVCent_FitMean_p[nMaxCentBins];
  TH1F* recoOverGenVCent_FitWidth_p[nMaxCentBins];
  TH1F* recoOverGenVCent_FitWidthOverMean_p[nMaxCentBins];

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!configs.ContainsParamSet(pbpbParams)) return 1;    
    centBins = strToVect(fileConfig_p->GetValue("CENTBINS", ""));
    nCentBins = centBins.size()-1;

    if(addPP){
      ++nCentBins;
    }

    if(nCentBins > nMaxCentBins){
      std::cout << "Given nCentBins \'" << nCentBins << "\' is greater than max value \'" << nMaxCentBins << "\'. return 1" << std::endl;
      return 1;
    }
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    if(!isPP){
      if(addPP && cI == nCentBins-1) centStr = "PP";
      else centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
    }

    recoOverGenVCent_FitMean_p[cI] = new TH1F(("recoOverGenFitMeanVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#LT Reco./Gen #GT;", nJtPtBins, jtPtBins);
    recoOverGenVCent_FitWidth_p[cI] = new TH1F(("recoOverGenFitWidthVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#sigma(Reco./Gen);", nJtPtBins, jtPtBins);
    recoOverGenVCent_FitWidthOverMean_p[cI] = new TH1F(("recoOverGenFitWidthOverMeanVCent_" + centStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h").c_str(), ";Gen. Jet p_{T} [GeV];#sigma(Reco./Gen/)/#LT Reco./Gen #GT;", nJtPtBins, jtPtBins);
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TF1* fit_p = new TF1("fit_p", "gaus", 0, 3.0);

  const double tinyCanvMargin = 0.0001;
  const double tinyPadMargin = 0.01;
  const double bottomMargin = 0.14;
  const double leftMargin = 0.2;

  const int jetR = fileConfig_p->GetValue("JETR", 10);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  globalLabels.push_back("anti-k_{T} R=0." + std::to_string(jetR) + " jets");
  
  int nXVal = 0;
  int nYVal = 0;
  getNXNYPanels(nJtPtBins, &nXVal, &nYVal);
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;
  
    TCanvas* canv_p = new TCanvas("canv_p", "", nXVal*450, nYVal*450);
    canv_p->SetTopMargin(tinyCanvMargin);
    canv_p->SetLeftMargin(tinyCanvMargin);
    canv_p->SetRightMargin(tinyCanvMargin);
    canv_p->SetBottomMargin(tinyCanvMargin);

    canv_p->Divide(nXVal, nYVal, 0.0, 0.0);
    
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;
    
    std::string centStr = "PP";
    std::string centLabel = "p+p";
    if(!isPP){
      if(addPP && cI == nCentBins-1){
	centStr = "PP";
	centLabel = "p+p";
      }
      else{
	centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
	centLabel = "Pb+Pb, " + centBins[cI] + "-" + centBins[cI+1] + "%";
      }
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;
    
    for(Int_t jI = 0; jI < nJtPtBins; ++jI){
      canv_p->cd();
      canv_p->cd(jI+1);

      gPad->SetTopMargin(tinyPadMargin);
      gPad->SetRightMargin(tinyPadMargin);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetLeftMargin(leftMargin);
      
      std::string jtPtStr = "JtPt";
      if(jI < 10) jtPtStr = jtPtStr + "0";
      jtPtStr = jtPtStr + std::to_string(jI);

      std::string jtPtLabel = prettyString(jtPtBins[jI], 1, false) + " < p_{T}^{Gen.} < " + prettyString(jtPtBins[jI+1], 1, false);
      
      std::string histStr = centStr + "/photonJtRecoOverGenVCentJtPt_" + centStr + "_" + jtPtStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h";
      TH1F* hist_p = nullptr;
      if(addPP && cI == nCentBins-1) hist_p = (TH1F*)inFilePP_p->Get(histStr.c_str());
      else hist_p = (TH1F*)inFile_p->Get(histStr.c_str());
   
      hist_p->Scale(1.0/hist_p->Integral());
      hist_p->GetYaxis()->SetTitle("Unity Norm.");
      hist_p->GetXaxis()->SetTitleFont(titleFont);
      hist_p->GetYaxis()->SetTitleFont(titleFont);
      hist_p->GetXaxis()->SetTitleSize(titleSizeX*(double)TMath::Min(((double)nXVal)/2., ((double)nYVal)/2.));
      hist_p->GetYaxis()->SetTitleSize(titleSizeY*(double)TMath::Min(((double)nXVal)/2., ((double)nYVal)/2.));

      hist_p->GetXaxis()->SetLabelFont(titleFont);
      hist_p->GetYaxis()->SetLabelFont(titleFont);
      hist_p->GetXaxis()->SetLabelSize(labelSizeX*(double)TMath::Min(((double)nXVal)*2./5., ((double)nYVal)*2./5.));
      hist_p->GetYaxis()->SetLabelSize(labelSizeY*(double)TMath::Min(((double)nXVal)*2./5., ((double)nYVal)*2./5.)); 
     
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", cI=" << cI << "/" << nCentBins << ", jI=" << jI << "/" << nJtPtBins << std::endl;
      
      HIJet::Style::EquipHistogram(hist_p, 0);

      hist_p->GetXaxis()->SetNdivisions(505);
      hist_p->GetYaxis()->SetNdivisions(505);      
      hist_p->GetXaxis()->SetTitleOffset(0.7);
      hist_p->GetYaxis()->SetTitleOffset(1.1);
      hist_p->SetMinimum(0.0);
      hist_p->DrawCopy("HIST E1 P");
      gStyle->SetOptStat(0);

      label_p->SetTextSize(labelSizeX*(double)TMath::Min(((double)nXVal)*2./5., ((double)nYVal)*2./5.));

      if(jI == 0){
	for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
	  label_p->DrawLatex(0.25, 0.92 - gI*0.07, globalLabels[gI].c_str());
	}
	label_p->DrawLatex(0.25, 0.92 - ((double)(globalLabels.size()))*0.07, centLabel.c_str());
      }

      label_p->DrawLatex(0.25, 0.92 - ((double)(globalLabels.size()+1))*0.07, jtPtLabel.c_str());
      
      
      Double_t minVal = TMath::Max((Double_t)hist_p->GetBinLowEdge(1), (Double_t)(recoJtPtMin/jtPtBinsLow));
      hist_p->Fit("fit_p", "M E Q N", "", minVal, hist_p->GetBinLowEdge(hist_p->GetXaxis()->GetNbins()+1));

      Double_t mean1 = fit_p->GetParameter(1);
      Double_t sigma1 = fit_p->GetParameter(2);
      
      minVal = TMath::Max(mean1 - sigma1*1.5, (Double_t)(recoJtPtMin/jtPtBinsLow));
      Double_t maxVal = TMath::Min(mean1 + sigma1*1.5, (Double_t)hist_p->GetBinLowEdge(hist_p->GetXaxis()->GetNbins()+1));
      
      hist_p->Fit("fit_p", "M E Q N", "", minVal, maxVal);      

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << std::endl;
      
      //      fit_p->SetMarkerSize(0.5);
      //      fit_p->SetLineWidth(0.5);
      TF1* fitTemp_p = new TF1("fitTemp_p", "gaus", minVal, maxVal);
      fitTemp_p->SetLineStyle(2);
      fitTemp_p->SetParameter(0, fit_p->GetParameter(0));
      fitTemp_p->SetParameter(1, fit_p->GetParameter(1));
      fitTemp_p->SetParameter(2, fit_p->GetParameter(2));
      fitTemp_p->DrawCopy("SAME");
      delete fitTemp_p;
      
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
    
    std::string saveName = "pdfDir/" + dateStr + "/photonJtRecoOverGen_Cent" + std::to_string(cI) + "_R" + std::to_string(jetR) + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  }
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  delete fit_p;
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  plotMeanAndSigma(doGlobalDebug, inConfig_p, &labelMap, dateStr, centBins, recoOverGenVCent_FitMean_p, recoOverGenVCent_FitWidth_p);
  plotMeanAndSigma(doGlobalDebug, inConfig_p, &labelMap, dateStr, centBins, recoOverGenVCent_FitMean_p, recoOverGenVCent_FitWidthOverMean_p);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Lets look at the corrections

  int nPanelY = 450;
  int nPanelX = nPanelY*(1 - otherMargin - bottomMargin)/(1 - leftMargin - rightMargin);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    std::string centLabel = "p+p";
    if(!isPP){
      if(addPP && cI == nCentBins-1){
	centStr = "PP";
	centLabel = "p+p";
      }
      else{
	centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
	centLabel = "Pb+Pb, " + centBins[cI] + "-" + centBins[cI+1] + "%";
      }
    }

    for(int jI = 0; jI < nJtEtaBinsSub; ++jI){
      std::string jtEtaStr = std::to_string(jI);
      if(jI < 10) jtEtaStr = "0" + jtEtaStr;

      std::string jtLabelStr = "#eta_{Jet}";
      if(jtEtaBinsSubDoAbs) jtLabelStr = "|" + jtLabelStr + "|";
      jtLabelStr = jtLabelStr + " < " + prettyString(jtEtaBinsSub[jI+1],1, false);
      if(!jtEtaBinsSubDoAbs || TMath::Abs(jtEtaBinsSub[jI]) > 0.00001) jtLabelStr = prettyString(jtEtaBinsSub[jI],1, false) + " < " + jtLabelStr;			    
			      
      TCanvas* canv_p = new TCanvas("canv_p", "", nPanelX, nPanelY);
      canv_p->SetTopMargin(otherMargin);
      canv_p->SetRightMargin(rightMargin);
      canv_p->SetBottomMargin(bottomMargin);
      canv_p->SetLeftMargin(leftMargin);

      const std::string histName = centStr + "/photonJtEMScaleOverConstituentVCentJtEta_" + centStr + "_JtEta" + jtEtaStr + "_GammaPt" + std::to_string(nGammaPtBinsSub) + "_DPhi0_h";
      
      TH1F* tempHist_p = (TH1F*)inFile_p->Get(histName.c_str());

      tempHist_p->GetXaxis()->SetTitleFont(titleFont);
      tempHist_p->GetYaxis()->SetTitleFont(titleFont);
      tempHist_p->GetXaxis()->SetLabelFont(titleFont);
      tempHist_p->GetYaxis()->SetLabelFont(titleFont);

      tempHist_p->GetXaxis()->SetTitleSize(titleSizeX);
      tempHist_p->GetYaxis()->SetTitleSize(titleSizeY);
      tempHist_p->GetXaxis()->SetLabelSize(labelSizeX);
      tempHist_p->GetYaxis()->SetLabelSize(labelSizeY);
      
      tempHist_p->DrawCopy("COLZ");      

      for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
	label_p->DrawLatex(0.25, 0.56 - gI*0.07, globalLabels[gI].c_str());
      }
      label_p->DrawLatex(0.25, 0.56 - ((double)(globalLabels.size()))*0.07, centLabel.c_str());
          
      label_p->DrawLatex(0.25, 0.56 - ((double)(globalLabels.size()+1))*0.07, jtLabelStr.c_str());
      
      std::string saveName = "pdfDir/" + dateStr + "/emScalOverConstituentScale_JtEta" + jtEtaStr + "_" +  centStr + "_" + dateStr + ".pdf";
      quietSaveAs(canv_p, saveName);
      delete canv_p;
    }
  }
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    delete recoOverGenVCent_FitWidth_p[cI];
    delete recoOverGenVCent_FitMean_p[cI];
    delete recoOverGenVCent_FitWidthOverMean_p[cI];
  }
  
  inFile_p->Close();
  delete inFile_p;

  if(addPP){
    inFilePP_p->Close();
    delete inFilePP_p;
  }
  
  delete inConfig_p;
  
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
