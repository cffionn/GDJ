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
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

bool recursiveHistSearch(std::string dateStr, TFile* inFile_p, std::string topDir = "")
{
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
  
  TKey* key;
  while((key=(TKey*)((*next)()))){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    if(isStrSame(className, "TDirectory") || isStrSame(className, "TDirectoryFile")){
      if(topDir.size() == 0) retVal = retVal && recursiveHistSearch(dateStr, inFile_p, name);
      else retVal = retVal && recursiveHistSearch(dateStr, inFile_p, topDir + "/" + name);
    }
    if(!isStrSame(className, "TH1D") && !isStrSame(className, "TH2D") && !isStrSame(className, "TH1F") && !isStrSame(className, "TH2F")) continue;
    
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetTopMargin(0.05);
    canv_p->SetRightMargin(0.12);
    canv_p->SetBottomMargin(0.12);

    if(isStrSame(className, "TH1D")){
      TH1D* tempHist_p = (TH1D*)key->ReadObj();
      tempHist_p->DrawCopy("HIST E1 P");
    }
    else if(isStrSame(className, "TH1F")){
      TH1F* tempHist_p = (TH1F*)key->ReadObj();
      tempHist_p->DrawCopy("HIST E1 P");
    }
    else if(isStrSame(className, "TH1D")){
      TH2D* tempHist_p = (TH2D*)key->ReadObj();
      tempHist_p->DrawCopy("COLZ");
    }
    else{
      TH2F* tempHist_p = (TH2F*)key->ReadObj();
      tempHist_p->DrawCopy("COLZ");
    }
   
    std::string saveName = name;
    while(saveName.find("/") != std::string::npos){
      saveName.replace(0, saveName.find("/")+1, "");
    }

    std::string labelName = saveName.substr(saveName.find("_")+1, saveName.size());//First part of the name should be clear from the figure itself so to simplify chop it off
    labelName.replace(labelName.rfind("_h"), 2, "");
    while(labelName.find("_") != std::string::npos){labelName.replace(labelName.find("_"), 1, "; ");}    
    label_p->DrawLatex(0.20, 0.965, labelName.c_str());
    gStyle->SetOptStat(0);

    saveName = "pdfDir/" + dateStr + "/" + saveName + "_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  }

  delete label_p;
  
  return retVal;
}

int HistDumping(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, "root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  recursiveHistSearch(dateStr, inFile_p);  
  
  inFile_p->Close();
  delete inFile_p;
  
  std::cout << "HISTDUMPING COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/HistDumping.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += HistDumping(argv[1]);
  return retVal;
}