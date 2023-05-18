//Author: Chris McGinn (2021.06.21)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/configParser.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/unfoldingUtil.h"
// #include "include/treeUtil.h"

std::vector<std::string> getDirContents(TDirectoryFile* dir_p, std::vector<std::string> classFilter, std::vector<std::string>* classVect=nullptr)
{
  std::vector<std::string> contentsToReturn;
  TIter next(dir_p->GetListOfKeys());
  TKey* key = nullptr;
  while((key=(TKey*)next())){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    bool isGoodClass = false;
    for(unsigned int cI = 0; cI < classFilter.size(); ++cI){
      if(isStrSame(classFilter[cI], className)){isGoodClass = true; break;}
    }

    if(!isGoodClass){
      if(className.find("TDirectory") != std::string::npos){
	TDirectoryFile* recursiveDir_p = (TDirectoryFile*)key->ReadObj();
	std::vector<std::string> recursiveClassVect;	
	std::vector<std::string> dirContents;
	if(classVect == nullptr) getDirContents(recursiveDir_p, classFilter);
	else getDirContents(recursiveDir_p, classFilter, &recursiveClassVect);

	contentsToReturn.insert(contentsToReturn.end(), dirContents.begin(), dirContents.end());
	if(classVect != nullptr) classVect->insert(classVect->end(), recursiveClassVect.begin(), recursiveClassVect.end());
      }
    }
    else{
      std::string finalName = dir_p->GetName();
      finalName = finalName + "/" + name;
      contentsToReturn.push_back(finalName);
      if(classVect != nullptr) classVect->push_back(className);
    }    
  }

  return contentsToReturn;
}

int gdjHistDQM(std::string inConfigFileName)
{
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INOLDFILENAME",
					      "INNEWFILENAME",
					      "PRECISION"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  const std::string inOldFileName = config_p->GetValue("INOLDFILENAME", "");
  const std::string inNewFileName = config_p->GetValue("INNEWFILENAME", "");
  const double precision = config_p->GetValue("PRECISION", -1.0);

  std::vector<std::string> stringsToSubOut, stringsToSubIn;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(int i = 0; i < 100; ++i){
    std::string stringToSubOut = config_p->GetValue(("STRINGTOSUBOUT." + std::to_string(i)).c_str(), "");
    std::string stringToSubIn = config_p->GetValue(("STRINGTOSUBIN." + std::to_string(i)).c_str(), "");

    if(stringToSubOut.size() == 0) continue;
    if(stringToSubIn.size() == 0) continue;

    stringsToSubOut.push_back(stringToSubOut);
    stringsToSubIn.push_back(stringToSubIn);
  }
    
  if(!check.checkFileExt(inOldFileName, ".root")) return 1;
  if(!check.checkFileExt(inNewFileName, ".root")) return 1;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<std::string> classFilter = {"TH1F", "TH1D", "TH2F", "TH2D"};  
  std::vector<std::string> goodObjectsOld, goodObjectsNew;
  std::vector<std::string> objectClassOld, objectClassNew;
     
  TFile* oldFile_p = new TFile(inOldFileName.c_str(), "READ");
  TIter nextOld(oldFile_p->GetListOfKeys());
  TKey* key = nullptr;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  while((key=(TKey*)nextOld())){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    bool isGoodClass = vectContainsStr(className, &classFilter);
    
    if(!isGoodClass){
      if(className.find("TDirectory") != std::string::npos){
	TDirectoryFile* dir_p = (TDirectoryFile*)key->ReadObj();
	std::vector<std::string> classContents;
	std::vector<std::string> dirContents = getDirContents(dir_p, classFilter, &classContents);
	goodObjectsOld.insert(goodObjectsOld.end(), dirContents.begin(), dirContents.end());
	objectClassOld.insert(objectClassOld.end(),  classContents.begin(),  classContents.end());
      }
      
      continue;
    }
    else{
      goodObjectsOld.push_back(name);
      objectClassOld.push_back(className);
    }
  }
	       
  oldFile_p->Close();
  delete oldFile_p;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TFile* newFile_p = new TFile(inNewFileName.c_str(), "READ");
  TIter nextNew(newFile_p->GetListOfKeys());
  key = nullptr;
  while((key=(TKey*)nextNew())){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    bool isGoodClass = vectContainsStr(className, &classFilter);
    
    if(!isGoodClass){
      if(className.find("TDirectory") != std::string::npos){
	TDirectoryFile* dir_p = (TDirectoryFile*)key->ReadObj();
	std::vector<std::string> classContents;
	std::vector<std::string> dirContents = getDirContents(dir_p, classFilter, &classContents);
	goodObjectsNew.insert(goodObjectsNew.end(), dirContents.begin(), dirContents.end());
	objectClassNew.insert(objectClassNew.end(),  classContents.begin(),  classContents.end());
      }
      
      continue;
    }
    else{      
      goodObjectsNew.push_back(name);
      objectClassNew.push_back(className);	
    }
  }
	       
  newFile_p->Close();
  delete newFile_p;

  newFile_p = new TFile(inNewFileName.c_str(), "READ");
  oldFile_p = new TFile(inOldFileName.c_str(), "READ");

  const Int_t nPad1D = 2;
  const Float_t noMargin = 0.0001;

  const Float_t padSplit1D = 0.36;
  
  const Float_t topMargin1D = 0.01;
  const Float_t leftMargin1D = 0.16;
  const Float_t rightMargin1D = 0.01;
  const Float_t bottomMargin1D = 0.125/padSplit1D;
  Double_t height1D = 900.0;
  Double_t width1D = height1D*(1.0 - topMargin1D*(1.0 - padSplit1D) - bottomMargin1D*padSplit1D)/(1.0 - leftMargin1D - rightMargin1D);

  const Int_t nPad2D = 3;
  const Float_t topMargin2D = 0.06;
  const Float_t bottomMargin2D = 0.10;
  const Float_t leftMargin2D = 0.16;
  const Float_t rightMargin2D = 0.16;
  const Float_t width2D = 1400;
  const Float_t height2D = 500;

  const Float_t padSplit2D1 = 0.333;
  const Float_t padSplit2D2 = 0.667;
  
  
  const Int_t nPadMax = TMath::Max(nPad1D, nPad2D);

  
  const int titleFont = 42;
  const double titleSize = 0.035;
  //  const double labelSize = titleSize*0.9;
  const double yOffset = 1.5;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<std::string> oldHistSkipped;
  std::vector<std::string> histFailingBins;
  std::vector<std::string> histFailingPrecision;
  std::map<std::string, int> newHistUsedCounter;  
  for(unsigned int sI = 0; sI < goodObjectsNew.size(); ++sI){
    newHistUsedCounter[goodObjectsNew[sI]] = 0;
  }
  
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int gI = 0; gI < goodObjectsOld.size(); ++gI){
    std::string newStrToFind = goodObjectsOld[gI];
    for(unsigned int sI = 0; sI < stringsToSubOut.size(); ++sI){
      if(newStrToFind.find(stringsToSubOut[sI]) != std::string::npos){
	std::string tempNewStr = newStrToFind;

	tempNewStr.replace(tempNewStr.find(stringsToSubOut[sI]), stringsToSubOut[sI].size(), stringsToSubIn[sI]);

	if(vectContainsStr(tempNewStr, &goodObjectsNew)){newStrToFind = tempNewStr; break;}
      }
    }


    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
    
    int pos = vectContainsStrPos(newStrToFind, &goodObjectsNew);
    if(pos < 0){
      oldHistSkipped.push_back(goodObjectsOld[gI]);
      std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' in old file \'" << inOldFileName << "\' had no corresponding partner in new file" << std::endl;
      continue;
    }

    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
    
    ++(newHistUsedCounter[goodObjectsNew[pos]]);
    
    TCanvas* canv_p = nullptr;
    TPad* pads_p[nPadMax];

    Int_t nPadForThisHist = -1;
    if(objectClassOld[gI].find("TH1") != std::string::npos){
      nPadForThisHist = 2;

      canv_p = new TCanvas("canv_p", "", width1D, height1D);
      canv_p->SetTopMargin(noMargin);
      canv_p->SetLeftMargin(noMargin);
      canv_p->SetRightMargin(noMargin);
      canv_p->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0] = new TPad("pad0", "", 0.0, padSplit1D, 1.0, 1.0);
      pads_p[0]->SetTopMargin(topMargin1D);
      pads_p[0]->SetLeftMargin(leftMargin1D);
      pads_p[0]->SetRightMargin(rightMargin1D);
      pads_p[0]->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      
      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit1D);
      pads_p[1]->SetTopMargin(noMargin);
      pads_p[1]->SetLeftMargin(leftMargin1D);
      pads_p[1]->SetRightMargin(rightMargin1D);
      pads_p[1]->SetBottomMargin(bottomMargin1D);
      
      pads_p[1]->Draw("SAME");
    }
    else{
      nPadForThisHist = 3;

      canv_p = new TCanvas("canv_p", "", width2D, height2D);
      canv_p->SetTopMargin(noMargin);
      canv_p->SetLeftMargin(noMargin);
      canv_p->SetRightMargin(noMargin);
      canv_p->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0] = new TPad("pad0", "", 0.0, 0.0, padSplit2D1, 1.0);
      pads_p[0]->SetTopMargin(topMargin2D);
      pads_p[0]->SetLeftMargin(leftMargin2D);
      pads_p[0]->SetRightMargin(rightMargin2D);
      pads_p[0]->SetBottomMargin(bottomMargin2D);
      
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      
      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", padSplit2D1, 0.0, padSplit2D2, 1.0);
      pads_p[1]->SetTopMargin(topMargin2D);
      pads_p[1]->SetLeftMargin(leftMargin2D);
      pads_p[1]->SetRightMargin(rightMargin2D);
      pads_p[1]->SetBottomMargin(bottomMargin2D);
      
      pads_p[1]->Draw("SAME");

      canv_p->cd();
      pads_p[2] = new TPad("pad2", "", padSplit2D2, 0.0, 1.0, 1.0);
      pads_p[2]->SetTopMargin(topMargin2D);
      pads_p[2]->SetLeftMargin(leftMargin2D);
      pads_p[2]->SetRightMargin(rightMargin2D);
      pads_p[2]->SetBottomMargin(bottomMargin2D);
      
      pads_p[2]->Draw("SAME");      
    }
    
    canv_p->cd();
    pads_p[0]->cd();

    if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
    
    TH1F* histOldTH1F_p = nullptr;
    TH1F* histNewTH1F_p = nullptr;
    TH1D* histOldTH1D_p = nullptr;
    TH1D* histNewTH1D_p = nullptr;

    TH2F* histOldTH2F_p = nullptr;
    TH2F* histNewTH2F_p = nullptr;
    TH2D* histOldTH2D_p = nullptr;
    TH2D* histNewTH2D_p = nullptr;

    TLatex* label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(42);
    label_p->SetTextSize(0.0225/(1.0 - padSplit1D));
    label_p->SetTextAlign(31);
      
    TLine* line_p = new TLine();
    line_p->SetLineStyle(2);
    line_p->SetLineColor(1);

    TLegend* leg_p = new TLegend(0.65, 0.6, 0.95, 0.75);
    leg_p->SetTextFont(42);
    leg_p->SetTextSize(0.035/(1.0 - padSplit1D));
    leg_p->SetBorderSize(0);
    leg_p->SetFillStyle(0);
        
    if(isStrSame(objectClassOld[gI], "TH1F")){
      histOldTH1F_p = (TH1F*)oldFile_p->Get(goodObjectsOld[gI].c_str());
      histNewTH1F_p = (TH1F*)newFile_p->Get(goodObjectsNew[pos].c_str());

      if(histOldTH1F_p->GetEntries() < TMath::Power(10,-100) || histNewTH1F_p->GetEntries() < TMath::Power(10,-100)){
	std::cout << "GDJHISTDQM: Skipping \'" << goodObjectsOld[gI] << "\' as empty..." <<  std::endl;

	for(Int_t pI = 0; pI < nPadForThisHist; ++pI){
	  delete pads_p[pI];
	}    
	delete line_p;
	delete canv_p;
	delete leg_p;
	delete label_p;
	
	continue;
      }

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
      
      HIJet::Style::EquipHistogram(histOldTH1F_p, 0);
      HIJet::Style::EquipHistogram(histNewTH1F_p, 1);
      
      canv_p->cd();
      pads_p[0]->cd();

      histNewTH1F_p->GetYaxis()->SetTitleOffset(yOffset);
      
      histNewTH1F_p->DrawCopy("HIST E1 P");
      histOldTH1F_p->DrawCopy("HIST E1 P SAME");

      leg_p->AddEntry(histNewTH1F_p, "New", "P L");
      leg_p->AddEntry(histOldTH1F_p, "Old", "P L");
      leg_p->Draw("SAME");

      label_p->DrawLatex(0.96, 0.93, ("New: " + goodObjectsNew[pos]).c_str()); 
      label_p->DrawLatex(0.96, 0.885, ("Old: " + goodObjectsOld[gI]).c_str());
     
      gPad->SetTicks();

      if(histNewTH1F_p->GetSumw2()->fN == 0) histNewTH1F_p->Sumw2();

      bool allGood = true;
      bool allGoodBins = true;
      Int_t newNBinsX = histNewTH1F_p->GetXaxis()->GetNbins();
      Int_t oldNBinsX = histOldTH1F_p->GetXaxis()->GetNbins();

      if(newNBinsX != oldNBinsX) allGoodBins = false;

      if(allGoodBins){
	for(Int_t bIX = 0; bIX < histOldTH1F_p->GetXaxis()->GetNbins()+2; ++bIX){
	  //First check bin edges are identical within precision
	  Float_t binLowEdgeOld = histOldTH1F_p->GetXaxis()->GetBinLowEdge(bIX+1);
	  Float_t binLowEdgeNew = histNewTH1F_p->GetXaxis()->GetBinLowEdge(bIX+1);

	  if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
	    std::cout << "Bin edges dont match. Precision test failure" << std::endl;
	    allGoodBins = false;
	    break;
	  }
	  
	  if(bIX < histOldTH1F_p->GetXaxis()->GetNbins()+1){
	    double delta = TMath::Abs(histOldTH1F_p->GetBinContent(bIX+1) - histOldTH1F_p->GetBinContent(bIX+1));
	    if(delta > precision){allGood = false; break;}
	  }

	}
      }

      if(!allGoodBins){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check of binning" << std::endl;
	histFailingBins.push_back(goodObjectsOld[gI]);
      }
      
      if(!allGood){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check at precision \'" << precision << "\'" << std::endl;
	histFailingPrecision.push_back(goodObjectsOld[gI]);
      }

      canv_p->cd();
      pads_p[1]->cd();

      if(allGoodBins){
	histNewTH1F_p->Divide(histOldTH1F_p);	
	histNewTH1F_p->GetYaxis()->SetTitle("New/Old");
	histNewTH1F_p->SetMaximum(1.25);
	histNewTH1F_p->SetMinimum(0.75);
	
	histNewTH1F_p->GetXaxis()->SetTitleFont(titleFont);
	histNewTH1F_p->GetYaxis()->SetTitleFont(titleFont);
	histNewTH1F_p->GetXaxis()->SetLabelFont(titleFont);
	histNewTH1F_p->GetYaxis()->SetLabelFont(titleFont);
	
	histNewTH1F_p->GetXaxis()->SetTitleSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1F_p->GetYaxis()->SetTitleSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1F_p->GetXaxis()->SetLabelSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1F_p->GetYaxis()->SetLabelSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	
	histNewTH1F_p->GetYaxis()->SetTitleOffset(yOffset*padSplit1D/(1.0 - padSplit1D));
	
	histNewTH1F_p->GetYaxis()->SetNdivisions(505);
	
	histNewTH1F_p->DrawCopy("HIST E1 P");

	line_p->DrawLine(histNewTH1F_p->GetBinLowEdge(1), 1.0, histNewTH1F_p->GetBinLowEdge(histNewTH1F_p->GetXaxis()->GetNbins()+1), 1.0);
      }
      else{
	label_p->SetTextAlign(21);
	label_p->DrawLatex(0.5, 0.5, "Bin mismatch; No division performed");
	label_p->SetTextAlign(11);
      }
    }
    else if(isStrSame(objectClassOld[gI], "TH1D")){
      histOldTH1D_p = (TH1D*)oldFile_p->Get(goodObjectsOld[gI].c_str());
      histNewTH1D_p = (TH1D*)newFile_p->Get(goodObjectsNew[pos].c_str());

      if(histOldTH1D_p->GetEntries() < TMath::Power(10,-100) || histNewTH1D_p->GetEntries() < TMath::Power(10,-100)){
	std::cout << "GDJHISTDQM: Skipping \'" << goodObjectsOld[gI] << "\' as empty..." <<  std::endl;

	for(Int_t pI = 0; pI < nPadForThisHist; ++pI){
	  delete pads_p[pI];
	}    
	delete line_p;
	delete canv_p;
	delete leg_p;
	delete label_p;

	continue;
      }
      
      HIJet::Style::EquipHistogram(histOldTH1D_p, 0);
      HIJet::Style::EquipHistogram(histNewTH1D_p, 1);
     
      canv_p->cd();
      pads_p[0]->cd();

      histNewTH1D_p->GetYaxis()->SetTitleOffset(yOffset);

      leg_p->AddEntry(histNewTH1D_p, "New", "P L");
      leg_p->AddEntry(histOldTH1D_p, "Old", "P L");
      leg_p->Draw("SAME");

      label_p->DrawLatex(0.96, 0.93, ("New: " + goodObjectsNew[pos]).c_str()); 
      label_p->DrawLatex(0.96, 0.885, ("Old: " + goodObjectsOld[gI]).c_str());
      
      histOldTH1D_p->DrawCopy("HIST E1 P");
      histNewTH1D_p->DrawCopy("HIST E! P SAME");
 
      gPad->SetTicks();

      if(histNewTH1D_p->GetSumw2()->fN == 0) histNewTH1D_p->Sumw2();

      bool allGood = true;
      bool allGoodBins = true;
      Int_t newNBinsX = histNewTH1D_p->GetXaxis()->GetNbins();
      Int_t oldNBinsX = histOldTH1D_p->GetXaxis()->GetNbins();

      if(newNBinsX != oldNBinsX) allGoodBins = false;

      if(allGoodBins){
	for(Int_t bIX = 0; bIX < histOldTH1D_p->GetXaxis()->GetNbins()+2; ++bIX){
	  Float_t binLowEdgeOld = histOldTH1D_p->GetXaxis()->GetBinLowEdge(bIX+1);
	  Float_t binLowEdgeNew = histNewTH1D_p->GetXaxis()->GetBinLowEdge(bIX+1);

	  if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
	    std::cout << "Bin edges dont match. Precision test failure" << std::endl;
	    allGoodBins = false;
	    break;
	  }

	  if(bIX < histOldTH1D_p->GetXaxis()->GetNbins()+1){
	    double delta = TMath::Abs(histOldTH1D_p->GetBinContent(bIX+1) - histOldTH1D_p->GetBinContent(bIX+1));
	    if(delta > precision){allGood = false; break;}
	  }
	}
      }
      
      if(!allGoodBins){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check of binning'" << std::endl;
	histFailingBins.push_back(goodObjectsOld[gI]);
      }

      if(!allGood){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check at precision \'" << precision << "\'" << std::endl;
	histFailingPrecision.push_back(goodObjectsOld[gI]);
      }
      
      canv_p->cd();
      pads_p[1]->cd();

      if(allGoodBins){
	histNewTH1D_p->Divide(histOldTH1D_p);
	
	histNewTH1D_p->GetYaxis()->SetTitle("New/Old");
	histNewTH1D_p->SetMaximum(1.25);
	histNewTH1D_p->SetMinimum(0.75);
	
	histNewTH1D_p->GetXaxis()->SetTitleFont(titleFont);
	histNewTH1D_p->GetYaxis()->SetTitleFont(titleFont);
	histNewTH1D_p->GetXaxis()->SetLabelFont(titleFont);
	histNewTH1D_p->GetYaxis()->SetLabelFont(titleFont);
	
	histNewTH1D_p->GetXaxis()->SetTitleSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1D_p->GetYaxis()->SetTitleSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1D_p->GetXaxis()->SetLabelSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	histNewTH1D_p->GetYaxis()->SetLabelSize(titleSize*(1.0 - padSplit1D)/padSplit1D);
	
	histNewTH1D_p->GetYaxis()->SetTitleOffset(yOffset*padSplit1D/(1.0 - padSplit1D));
	
	histNewTH1D_p->GetYaxis()->SetNdivisions(505);
	
	histNewTH1D_p->DrawCopy("HIST E1 P");    
	line_p->DrawLine(histNewTH1D_p->GetBinLowEdge(1), 1.0, histNewTH1D_p->GetBinLowEdge(histNewTH1D_p->GetXaxis()->GetNbins()+1), 1.0);
      }
      else{
	label_p->SetTextAlign(21);
	label_p->DrawLatex(0.5, 0.5, "Bin mismatch; No division performed");
	label_p->SetTextAlign(11);
      }
    }
    else if(isStrSame(objectClassOld[gI], "TH2F")){
      histOldTH2F_p = (TH2F*)oldFile_p->Get(goodObjectsOld[gI].c_str());
      histNewTH2F_p = (TH2F*)newFile_p->Get(goodObjectsNew[pos].c_str());

      if(histOldTH2F_p->GetEntries() < TMath::Power(10,-100) || histNewTH2F_p->GetEntries() < TMath::Power(10,-100)){
	std::cout << "GDJHISTDQM: Skipping \'" << goodObjectsOld[gI] << "\' as empty..." <<  std::endl;

	for(Int_t pI = 0; pI < nPadForThisHist; ++pI){
	  delete pads_p[pI];
	}    
	delete line_p;
	delete canv_p;
	delete leg_p;
	delete label_p;
	
	continue;
      }
      
      HIJet::Style::EquipHistogram(histOldTH2F_p, 0);
      HIJet::Style::EquipHistogram(histNewTH2F_p, 1);
     
      canv_p->cd();
      pads_p[0]->cd();

      histOldTH2F_p->GetYaxis()->SetTitleOffset(yOffset);

      std::string titleOld = goodObjectsOld[gI];
      while(titleOld.find("/") != std::string::npos){
	titleOld.replace(0, titleOld.find("/")+1, "");
      }
      histOldTH2F_p->SetTitle(("Old: " + titleOld).c_str());
	
      histOldTH2F_p->DrawCopy("COLZ"); 
      gPad->SetTicks();

      canv_p->cd();
      pads_p[1]->cd();

      histNewTH2F_p->GetYaxis()->SetTitleOffset(yOffset);

      std::string titleNew = goodObjectsNew[pos];
      while(titleNew.find("/") != std::string::npos){
	titleNew.replace(0, titleNew.find("/")+1, "");
      }
      histNewTH2F_p->SetTitle(("New: " + titleNew).c_str());
      histNewTH2F_p->DrawCopy("COLZ"); 
      gPad->SetTicks();
      
      if(histNewTH2F_p->GetSumw2()->fN == 0) histNewTH2F_p->Sumw2();
      if(histOldTH2F_p->GetSumw2()->fN == 0) histOldTH2F_p->Sumw2();

      bool allGood = true;
      bool allGoodBins = true;
      Int_t newNBinsX = histNewTH2F_p->GetXaxis()->GetNbins();
      Int_t oldNBinsX = histOldTH2F_p->GetXaxis()->GetNbins();
      Int_t newNBinsY = histNewTH2F_p->GetYaxis()->GetNbins();
      Int_t oldNBinsY = histOldTH2F_p->GetYaxis()->GetNbins();

      if(newNBinsX != oldNBinsX) allGoodBins = false;
      if(newNBinsY != oldNBinsY) allGoodBins = false;

      if(allGoodBins){
	for(Int_t bIX = 0; bIX < histOldTH2F_p->GetXaxis()->GetNbins()+2; ++bIX){
	  Float_t binLowEdgeOld = histOldTH2F_p->GetXaxis()->GetBinLowEdge(bIX+1);
	  Float_t binLowEdgeNew = histNewTH2F_p->GetXaxis()->GetBinLowEdge(bIX+1);

	  if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
	    std::cout << "Bin edges dont match. Precision test failure" << std::endl;
	    allGoodBins = false;
	    break;
	  }
	  
	  for(Int_t bIY = 0; bIY < histOldTH2F_p->GetYaxis()->GetNbins()+2; ++bIY){
	    if(bIX == 0){
	      binLowEdgeOld = histOldTH2F_p->GetYaxis()->GetBinLowEdge(bIY+1);
	      binLowEdgeNew = histNewTH2F_p->GetYaxis()->GetBinLowEdge(bIY+1);
	      
	      if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
		std::cout << "Bin edges dont match. Precision test failure" << std::endl;
		allGoodBins = false;
		break;
	      }
	    }	    

	    if(bIX < histOldTH2F_p->GetXaxis()->GetNbins()+1 && bIY < histOldTH2F_p->GetYaxis()->GetNbins()+1){           
	      double delta = TMath::Abs(histOldTH2F_p->GetBinContent(bIX+1, bIY+1) - histNewTH2F_p->GetBinContent(bIX+1, bIY+1));
	      if(delta > precision){allGood = false; break;}
	    }
	  }

	  if(!allGoodBins) break;
	  if(!allGood) break;
	}
      }
	
      if(!allGoodBins){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check of binning" << std::endl;
	histFailingBins.push_back(goodObjectsOld[gI]);
      }

      if(!allGood){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check at precision \'" << precision << "\'" << std::endl;
	histFailingPrecision.push_back(goodObjectsOld[gI]);
      }
      
      canv_p->cd();
      pads_p[2]->cd();

      if(allGoodBins){
	histNewTH2F_p->Divide(histOldTH2F_p);
	
	histNewTH2F_p->SetTitle("New/Old");
	histNewTH2F_p->SetMaximum(1.25);
	histNewTH2F_p->SetMinimum(0.75);
	
	histNewTH2F_p->GetXaxis()->SetTitleFont(titleFont);
	histNewTH2F_p->GetYaxis()->SetTitleFont(titleFont);
	histNewTH2F_p->GetXaxis()->SetLabelFont(titleFont);
	histNewTH2F_p->GetYaxis()->SetLabelFont(titleFont);
	
	histNewTH2F_p->GetYaxis()->SetNdivisions(505);
	
	histNewTH2F_p->DrawCopy("COLZ");
      }
      else{
	label_p->SetTextAlign(21);
	label_p->DrawLatex(0.5, 0.5, "Bin mismatch; No division performed");
	label_p->SetTextAlign(11);
      }
    }
    else if(isStrSame(objectClassOld[gI], "TH2D")){
      histOldTH2D_p = (TH2D*)oldFile_p->Get(goodObjectsOld[gI].c_str());
      histNewTH2D_p = (TH2D*)newFile_p->Get(goodObjectsNew[pos].c_str());

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
      
      if(histOldTH2D_p->GetEntries() < TMath::Power(10,-100) || histNewTH2D_p->GetEntries() < TMath::Power(10,-100)){
	std::cout << "GDJHISTDQM: Skipping \'" << goodObjectsOld[gI] << "\' as empty..." <<  std::endl;

	for(Int_t pI = 0; pI < nPadForThisHist; ++pI){
	  delete pads_p[pI];
	}    
	delete line_p;
	delete canv_p;
	delete leg_p;
	delete label_p;
	
	continue;
      }

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
    
      HIJet::Style::EquipHistogram(histOldTH2D_p, 0);
      HIJet::Style::EquipHistogram(histNewTH2D_p, 1);
     
      canv_p->cd();
      pads_p[0]->cd();

      histOldTH2D_p->GetYaxis()->SetTitleOffset(yOffset);
      histOldTH2D_p->SetTitle(("Old: " + goodObjectsOld[gI]).c_str());
      histOldTH2D_p->DrawCopy("COLZ"); 
      gPad->SetTicks();

      canv_p->cd();
      pads_p[1]->cd();

      histNewTH2D_p->GetYaxis()->SetTitleOffset(yOffset);
      histNewTH2D_p->SetTitle(("New: " + goodObjectsNew[pos]).c_str());
      histNewTH2D_p->DrawCopy("COLZ"); 
      gPad->SetTicks();
      
      if(histNewTH2D_p->GetSumw2()->fN == 0) histNewTH2D_p->Sumw2();
      if(histOldTH2D_p->GetSumw2()->fN == 0) histOldTH2D_p->Sumw2();

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;

      bool allGood = true;
      bool allGoodBins = true;
      Int_t newNBinsX = histNewTH2D_p->GetXaxis()->GetNbins();
      Int_t oldNBinsX = histOldTH2D_p->GetXaxis()->GetNbins();
      Int_t newNBinsY = histNewTH2D_p->GetYaxis()->GetNbins();
      Int_t oldNBinsY = histOldTH2D_p->GetYaxis()->GetNbins();

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
      
      if(newNBinsX != oldNBinsX) allGoodBins = false;
      if(newNBinsY != oldNBinsY) allGoodBins = false;

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
      
      if(allGoodBins){
	if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;

	for(Int_t bIX = 0; bIX < histOldTH2D_p->GetXaxis()->GetNbins()+2; ++bIX){
	  Float_t binLowEdgeOld = histOldTH2D_p->GetXaxis()->GetBinLowEdge(bIX+1);
	  Float_t binLowEdgeNew = histNewTH2D_p->GetXaxis()->GetBinLowEdge(bIX+1);

	  if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
	    std::cout << "Bin edges dont match. Precision test failure" << std::endl;
	    allGoodBins = false;
	    break;
	  }

	  for(Int_t bIY = 0; bIY < histOldTH2D_p->GetYaxis()->GetNbins()+2; ++bIY){

	    if(bIX == 0){
	      binLowEdgeOld = histOldTH2D_p->GetYaxis()->GetBinLowEdge(bIY+1);
	      binLowEdgeNew = histNewTH2D_p->GetYaxis()->GetBinLowEdge(bIY+1);
	      
	      if(TMath::Abs(binLowEdgeOld - binLowEdgeNew) > precision){
		std::cout << "Bin edges dont match. Precision test failure" << std::endl;
		allGoodBins = false;
		break;
	      }
	    }	    

	    if(bIX < histOldTH2D_p->GetXaxis()->GetNbins()+1 && bIY < histOldTH2D_p->GetYaxis()->GetNbins()+1){           	      
	      double delta = TMath::Abs(histOldTH2D_p->GetBinContent(bIX+1, bIY+1) - histNewTH2D_p->GetBinContent(bIX+1, bIY+1));
	      if(delta > precision){allGood = false; break;}
	    }
	  }

	  if(!allGoodBins) break;

	  if(!allGood) break;
	}

	if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
	
      }      
	
      if(!allGoodBins){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check of binning" << std::endl;
	histFailingBins.push_back(goodObjectsOld[gI]);
      }

      if(!allGood){
	std::cout << "Histogram \'" << goodObjectsOld[gI] << "\' fails check at precision \'" << precision << "\'" << std::endl;
	histFailingPrecision.push_back(goodObjectsOld[gI]);
      }

      if(doGlobalDebug) std::cout << "FILE, LINE, gI/nObjects, string: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << goodObjectsOld.size() << ", " << goodObjectsOld[gI] <<  std::endl;
      
      canv_p->cd();
      pads_p[2]->cd();

      if(allGoodBins){
	histNewTH2D_p->Divide(histOldTH2D_p);
	
	histNewTH2D_p->SetTitle("New/Old");
	histNewTH2D_p->SetMaximum(1.25);
	histNewTH2D_p->SetMinimum(0.75);
	
	histNewTH2D_p->GetXaxis()->SetTitleFont(titleFont);
	histNewTH2D_p->GetYaxis()->SetTitleFont(titleFont);
	histNewTH2D_p->GetXaxis()->SetLabelFont(titleFont);
	histNewTH2D_p->GetYaxis()->SetLabelFont(titleFont);
	
	histNewTH2D_p->GetYaxis()->SetNdivisions(505);
	
	histNewTH2D_p->DrawCopy("COLZ");    
      }
      else{
	label_p->SetTextAlign(21);
	label_p->DrawLatex(0.5, 0.5, "Bin mismatch; No division performed");
	label_p->SetTextAlign(11);
      }

    }


    
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    
    std::string saveName = goodObjectsOld[gI];
    while(saveName.find("/") != std::string::npos){saveName.replace(0, saveName.find("/")+1, "");}
    saveName = "pdfDir/" + dateStr + "/" + saveName + "_" + dateStr + ".png";
    //    if(isStrSame(objectClassOld[gI], "TH1F")){
      quietSaveAs(canv_p, saveName);
      //    }
    
    for(Int_t pI = 0; pI < nPadForThisHist; ++pI){
      delete pads_p[pI];
    }    
    delete line_p;
    delete canv_p;
    delete leg_p;
    delete label_p;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  oldFile_p->Close();
  delete oldFile_p;

  newFile_p->Close();
  delete newFile_p;
  
  delete config_p;

  std::cout << "JOB SUMMARY: " << std::endl;
  std::cout << "Skipped histograms in 'old' file: " << std::endl;
  for(unsigned int sI = 0; sI < oldHistSkipped.size(); ++sI){
    std::cout << " " << oldHistSkipped[sI] << std::endl;
  }
  std::cout << "Histograms failing bin comparisons: " << std::endl;
  for(unsigned int vI = 0; vI < histFailingBins.size(); ++vI){
    std::cout << " " << histFailingBins[vI] << std::endl;
  }  
  std::cout << "Histograms failing precision comparisons: " << std::endl;
  for(unsigned int vI = 0; vI < histFailingPrecision.size(); ++vI){
    std::cout << " " << histFailingPrecision[vI] << std::endl;
  }  
  std::cout << "'New' histograms with zero (0) matches: " << std::endl;
  for(auto const & val : newHistUsedCounter){
    if(val.second == 0) std::cout << " " << val.first << std::endl;
  }

  std::cout << "'New' histograms with multiple matches: " << std::endl;
  for(auto const & val : newHistUsedCounter){
    if(val.second >= 2) std::cout << " " << val.first << " (" << val.second << " matches)" << std::endl;
  }
  
  std::cout << "GDJHISTDQM COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHistDQM.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjHistDQM(argv[1]);
  return retVal;
}
