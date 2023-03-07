//Author: Chris McGinn (2020.02.10)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

//Local
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

bool recursiveHistSearch(std::string dateStr, TFile* inFile_p, std::map<std::string, std::string>* labelMap, std::string filterStr = "", std::string topDir = "", TFile* matchedFile_p=nullptr, std::vector<std::string> legVect = {})
{
  globalDebugHandler gDebug;
  bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  bool retVal = true;
  TIter* next;
  
  if(topDir.size() == 0){
    next = new TIter(inFile_p->GetListOfKeys());
  }
  else{
    TDirectoryFile* dir_p = (TDirectoryFile*)inFile_p->Get(topDir.c_str());
    next = new TIter(dir_p->GetListOfKeys());
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t markerSize = 1;
  const Int_t markerStyle = 24;
  const Int_t markerColor = 1;
  const Int_t lineColor = 1;

  const Int_t markerSize2 = 1;
  const Int_t markerStyle2 = 25;
  const Int_t markerColor2 = 2;
  const Int_t lineColor2 = 2;
  
  const Float_t yOffset = 1.8;

  const Float_t topMargin = 0.06;
  const Float_t leftMargin = 0.16;
  const Float_t rightMargin = 0.12;
  const Float_t bottomMargin = 0.12;

  const Float_t height = 450;
  const Float_t width = height*(1.0 - topMargin - bottomMargin)/(1.0 - leftMargin - rightMargin);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TLine* line_p = new TLine();
  line_p->SetLineWidth(2);
  line_p->SetLineColor(1);
  line_p->SetLineStyle(2);
  
  TKey* key;
  while((key=(TKey*)((*next)()))){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();
    
    if(isStrSame(className, "TDirectory") || isStrSame(className, "TDirectoryFile")){
      if(topDir.size() == 0) retVal = retVal && recursiveHistSearch(dateStr, inFile_p, labelMap, filterStr, name, matchedFile_p, legVect);
      else retVal = retVal && recursiveHistSearch(dateStr, inFile_p, labelMap, filterStr, topDir + "/" + name, matchedFile_p, legVect);
    }

    if(!isStrSame(className, "TH1D") && !isStrSame(className, "TH2D") && !isStrSame(className, "TH1F") && !isStrSame(className, "TH2F")) continue;

    //TEMP TO GET THE TH2 only!
    if(!isStrSame(className, "TH2D") && !isStrSame(className, "TH2F")) continue;

    
    else if(name.find("GenRes") != std::string::npos && name.find("GenPt") != std::string::npos && name.find("RecoPt") != std::string::npos) continue;

    if(filterStr.size() != 0){
      if(name.find(filterStr) == std::string::npos) continue;
    }
    
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
    canv_p->SetLeftMargin(leftMargin);
    canv_p->SetTopMargin(topMargin);
    canv_p->SetRightMargin(rightMargin);
    canv_p->SetBottomMargin(bottomMargin);

    std::cout << "NAME: " << name << std::endl;
    std::cout << "NAME2: " << topDir << "/" << name << std::endl;

    TLegend* leg_p = nullptr;
    
    if(isStrSame(className, "TH1D")){
      TH1D* tempHist_p = (TH1D*)key->ReadObj();
      
      TH1D* matchedHist_p = nullptr;
      if(matchedFile_p != nullptr){
	if(topDir.size() != 0) matchedHist_p = (TH1D*)matchedFile_p->Get((topDir + "/" + name).c_str());
	else matchedHist_p = (TH1D*)matchedFile_p->Get(name.c_str());
      }
    
      tempHist_p->SetMarkerSize(markerSize);
      tempHist_p->SetMarkerStyle(markerStyle);
      tempHist_p->SetMarkerColor(markerColor);
      tempHist_p->SetLineColor(lineColor);

      //Temp!!!!
      binWidthAndSelfNorm(tempHist_p);
      if(matchedFile_p != nullptr){
	binWidthAndSelfNorm(matchedHist_p);
	matchedHist_p->SetMarkerSize(markerSize2);
	matchedHist_p->SetMarkerStyle(markerStyle2);
	matchedHist_p->SetMarkerColor(markerColor2);
	matchedHist_p->SetLineColor(lineColor2);
      }
      std::string xTitle = tempHist_p->GetXaxis()->GetTitle();
      xTitle.replace(0,xTitle.find("Jet")+4, "");
      if(xTitle.find(" [") != std::string::npos) xTitle.replace(xTitle.rfind(" ["), xTitle.size(), "");

      if(name.find("DPhiJJ") != std::string::npos) tempHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{JJ}} #frac{dN_{JJ}}{d" + xTitle + "}").c_str());
      else tempHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{Jet}} #frac{dN_{Jet}}{d" + xTitle + "}").c_str());
      
      tempHist_p->GetYaxis()->SetTitleOffset(yOffset);

      if(name.find("JtDPhi") != std::string::npos){
	tempHist_p->SetMaximum(0.5);
	tempHist_p->SetMinimum(-0.02);
      }

      Double_t tempMax = getMax(tempHist_p);
      if(matchedHist_p != nullptr) tempMax = TMath::Max(tempMax, getMax(matchedHist_p));
      tempHist_p->SetMaximum(1.15*tempMax);
      if(name.find("multijetPt") == std::string::npos) tempHist_p->SetMinimum(0.0);
      
      tempHist_p->DrawCopy("HIST E1 P");
      if(matchedFile_p != nullptr){
	matchedHist_p->DrawCopy("HIST E1 P SAME");

	Float_t legX = 0.7;
	Float_t legY = 0.8;
	if(name.find("DPhi") != std::string::npos){
	  legX -= 0.3;
	}
	else if(name.find("Phi") != std::string::npos){
	  legY -= 0.3;
	}
	
	leg_p = new TLegend(legX, legY, legX+0.25, legY+0.1);

	leg_p->SetTextFont(42);
        leg_p->SetTextSize(0.035);
        leg_p->SetBorderSize(0);
        leg_p->SetFillStyle(0);

	leg_p->AddEntry(tempHist_p, legVect[0].c_str(), "P L");
	leg_p->AddEntry(matchedHist_p, legVect[1].c_str(), "P L");

	leg_p->Draw("SAME");
      }
      if(name.find("multijetPt") != std::string::npos) gPad->SetLogy();
      
      if(name.find("JtDPhi") != std::string::npos) line_p->DrawLine(tempHist_p->GetBinLowEdge(1), 0.0, tempHist_p->GetBinLowEdge(tempHist_p->GetXaxis()->GetNbins()+1), 0.0);
    }
    else if(isStrSame(className, "TH2F")){
      TH1F* tempHist_p = (TH1F*)key->ReadObj();
      tempHist_p->SetMaximum(1.15*getMax(tempHist_p));
      tempHist_p->SetMinimum(0.0);

      tempHist_p->SetMarkerSize(markerSize);
      tempHist_p->SetMarkerStyle(markerStyle);
      tempHist_p->SetMarkerColor(markerColor);
      tempHist_p->SetLineColor(lineColor);

      tempHist_p->GetYaxis()->SetTitleOffset(yOffset);

      if(name.find("JtDPhi") != std::string::npos){
	tempHist_p->SetMaximum(0.5);
	tempHist_p->SetMinimum(-0.02);
      }

      tempHist_p->DrawCopy("HIST E1 P");

      if(name.find("JtDPhi") != std::string::npos) line_p->DrawLine(tempHist_p->GetBinLowEdge(1), 0.0, tempHist_p->GetBinLowEdge(tempHist_p->GetXaxis()->GetNbins()+1), 0.0);
    }
    else if(isStrSame(className, "TH2D")){
      TH2D* tempHist_p = (TH2D*)key->ReadObj();
      tempHist_p->DrawCopy("COLZ");

      if(name.find("rooResJtPtGammaM") != std::string::npos){
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetLogz();
      }
      else if(name.find("rooResJtXJJGammaM") != std::string::npos){
	gPad->SetLogz();
      }
    }
    else{
      TH2F* tempHist_p = (TH2F*)key->ReadObj();
      tempHist_p->DrawCopy("COLZ");
      
      if(name.find("rooResJtPtGammaM") != std::string::npos){
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetLogz();
      }
      else if(name.find("rooResJtXJJGammaM") != std::string::npos){
	gPad->SetLogz();
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string saveName = name;
    while(saveName.find("/") != std::string::npos){
      saveName.replace(0, saveName.find("/")+1, "");
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string labelName = saveName.substr(saveName.find("_")+1, saveName.size());//First part of the name should be clear from the figure itself so to simplify chop it off
    if(labelName.find("_h") != std::string::npos) labelName.replace(labelName.rfind("_h"), 2, "");
    std::vector<std::string> preLabels;
    while(labelName.find("_") != std::string::npos){
      preLabels.push_back(labelName.substr(0, labelName.find("_")));
      labelName.replace(0, labelName.find("_")+1, "");
    }
    if(labelName.size() != 0) preLabels.push_back(labelName);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
    labelName = "";

    const double xPos = 0.22;
    double yPos = 0.965;
    for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
      if(labelMap->count(preLabels[pI]) != 0) preLabels[pI] = (*labelMap)[preLabels[pI]];

      if(labelName.size() + preLabels[pI].size() > 60){
	if(labelName.find(";") != std::string::npos) labelName.replace(labelName.rfind(";"), labelName.size(), "");
	label_p->DrawLatex(xPos, yPos, labelName.c_str());	
	yPos -= 0.065;
	labelName = "";
      }
      labelName = labelName + preLabels[pI] + "; ";
    }
    if(labelName.find(";") != std::string::npos) labelName.replace(labelName.rfind(";"), labelName.size(), "");

    label_p->DrawLatex(xPos, yPos, labelName.c_str());
    gStyle->SetOptStat(0);

    saveName = "pdfDir/" + dateStr + "/" + saveName + "_" + dateStr + ".png";
    quietSaveAs(canv_p, saveName);

    if(leg_p != nullptr) delete leg_p;
    delete canv_p;
  }

  delete line_p;
  delete label_p;
  
  return retVal;
}

int gdjHistDumper(std::string inFileName, std::string filterStr = "", std::string matchFileName = "", std::string commaSepLegStr = "")
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, "root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* label_p = (TEnv*)inFile_p->Get("label");
  configParser config(label_p);
  std::map<std::string, std::string> labelMap = config.GetConfigMap();

  TFile* matchedFile_p = nullptr;
  std::vector<std::string> legVect;
  if(matchFileName.size() != 0){
    if(check.checkFileExt(matchFileName, "root")){

      matchedFile_p = new TFile(matchFileName.c_str(), "READ");
      
      legVect = strToVect(commaSepLegStr);
      if(legVect.size() != 2){
	std::cout << "Not enough leg: " << commaSepLegStr << ", " << legVect.size() << std::endl;
	return 1;
      }
    }
  }

  recursiveHistSearch(dateStr, inFile_p, &labelMap, filterStr, "", matchedFile_p, legVect);  
  
  inFile_p->Close();
  delete inFile_p;
  
  std::cout << "HISTDUMPING COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3 && argc != 5){
    std::cout << "Usage: ./bin/gdjHistDumper.exe <inFileName> <filterStr=''> <matchedFile=''> <commaSepLegStr>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  if(argc == 2) retVal += gdjHistDumper(argv[1]);
  else if(argc == 3) retVal += gdjHistDumper(argv[1], argv[2]);
  else if(argc == 5) retVal += gdjHistDumper(argv[1], argv[2], argv[3], argv[4]);
  return retVal;
}
