//Author: Chris McGinn (2021.02.18)
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
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/binFlattener.h"
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
#include "include/varUtil.h"

//Function to convert internal syst. names to pretty labeling strings
std::string systNameToLegName(std::string inStr)
{
  if(isStrSame("PURSIDEBANDLOOSE", inStr)) return "Loose Sideband";
  else if(isStrSame("PURSIDEBANDTIGHT", inStr)) return "Tight Sideband";
  else if(isStrSame("PURSIDEBANDISO", inStr)) return "Iso. Sideband";
  else if(isStrSame("PURBINORFIT", inStr)) return "Binned Purity";
  else if(isStrSame("PHOISOANDPUR", inStr)) return "#gamma Iso./Purity";
  else if(isStrSame("MIXING", inStr)) return "Mixing";
  else if(isStrSame("UNFOLDING", inStr)) return "Unfolding";
  else if(isStrSame("MIXNONCLOSURE", inStr)) return "MC Nonclosure";
  else if(isStrSame("JTPTCUT", inStr)) return "Vary Jet p_{T} Cut";
  else if(isStrSame("PRIOR", inStr)) return "Vary Prior";
  else if(isStrSame("JESRTRK", inStr)) return "JES R_{Trk}";
  else if(isStrSame("QGFRAC", inStr)) return "Frac._{q/g}";
  else if(isStrSame("QGRESP", inStr)) return "Resp._{q/g}";
  else if(isStrSame("MCSTAT", inStr)) return "MC Stat.";
  else if(isStrSame("GESGER", inStr)) return "#gamma ES/ER";
  else if(isStrSame("JESJER", inStr)) return "JES/JER";

  return inStr;
}

//function to return vector of syst, bin-by-bin
//As a start, assume symmetric syst
std::vector<std::vector<Double_t> > getSyst(TFile* histFile_p, TH1D* nominalHist_p, std::vector<std::string> systHistNames, std::vector<Int_t> yPos, std::vector<Double_t> scaledTotalSyst, std::string unfFileName, bool doDebug=false)
{
  std::vector<std::vector<Double_t> > systVals;

  if(doDebug) std::cout << __FILE__ << ", " << __LINE__ << ", systHistNames.size: " << systHistNames.size() << std::endl;

  //loop over syst. hist names
  for(unsigned int jI = 0; jI < systHistNames.size(); ++jI){
    if(doDebug) std::cout << __FILE__ << ", " << __LINE__ << ", jI/nSystHist = " << jI << "/" << systHistNames.size() << ": " << systHistNames[jI] << std::endl;

    TH2D* temp2D_p = nullptr;
    TH1D* temp1D_p = (TH1D*)nominalHist_p->Clone("temp1D_p");
    for(Int_t bIX = 0; bIX < temp1D_p->GetXaxis()->GetNbins(); ++bIX){
      temp1D_p->GetBinContent(bIX+1, 0.0);
      temp1D_p->GetBinError(bIX+1, 0.0);
    }


    if(systHistNames[jI].find("UNFNONCLOSURE") == std::string::npos){
      if(doDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " <<  std::endl;
      temp2D_p = (TH2D*)histFile_p->Get(systHistNames[jI].c_str());
      fineTH2ToCoarseTH1(temp2D_p, temp1D_p, yPos);
      binWidthAndScaleNorm(temp1D_p, scaledTotalSyst[jI]);
    }
    else{//Unfolding non-closure must be handled specially
      std::string gammaPtStr = nominalHist_p->GetName();
      if(gammaPtStr.find("GammaPt") != std::string::npos){
	gammaPtStr.replace(0, gammaPtStr.find("GammaPt"), "");

	while(gammaPtStr.find("_") != std::string::npos){
	  gammaPtStr.replace(gammaPtStr.rfind("_"), gammaPtStr.size(), "");
	}
      }

      TFile* tempFile_p = new TFile(unfFileName.c_str(), "READ");
      TH1D* corrHist_p = (TH1D*)tempFile_p->Get(("unfoldOverTruthMC_" + gammaPtStr + "_PP_MCNonClosure").c_str());
      fineHistToCoarseHist(corrHist_p, temp1D_p);

      for(Int_t bIX = 0; bIX < temp1D_p->GetXaxis()->GetNbins(); ++bIX){
	Double_t value = nominalHist_p->GetBinContent(bIX+1)/temp1D_p->GetBinContent(bIX+1);
	Double_t error = nominalHist_p->GetBinError(bIX+1)/temp1D_p->GetBinContent(bIX+1);

	temp1D_p->SetBinContent(bIX+1, value);
	temp1D_p->SetBinError(bIX+1, error);
      }

      tempFile_p->Close();
      delete tempFile_p;
    }


    systVals.push_back({});
    for(Int_t bIX = 0; bIX < temp1D_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t deltaVal = TMath::Abs(nominalHist_p->GetBinContent(bIX+1) - temp1D_p->GetBinContent(bIX+1));

      //MCStat must be handled special as the syst is stored in the errors not the delta w/ nominal
      if(systHistNames[jI].find("MCSTAT") != std::string::npos){
	//Check the deltaVal is 0; i.e. things are working correctly
	if(deltaVal > TMath::Power(10,-50)){
	  std::cout << "getSyst ERROR!!!! DeltaVal of Nominal and MCStat variation is greater than 0, return" << std::endl;
	  std::cout << deltaVal << std::endl;
	  std::cout << "CHECK 1: " << nominalHist_p->GetName() << std::endl;
	  std::cout << "CHECK 2: " << systHistNames[jI] << std::endl;
	  std::cout << "Print 1: " << std::endl;
	  nominalHist_p->Print("ALL");
	  std::cout << "Print 2: " << std::endl;
	  temp1D_p->Print("ALL");
	  return {};
	}

	//if things check out then just take the error bar as the deltaVal
	deltaVal = temp1D_p->GetBinError(bIX+1);
      }

      systVals[jI].push_back(deltaVal);
    }

    delete temp1D_p;
  }

  if(doDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return systVals;
}

std::vector<std::vector<Double_t> > getSystRat(TFile* histFilePbPb_p, TFile* histFilePP_p, TH1D* nominalHistPbPb_p, TH1D* nominalHistPP_p, std::vector<std::string> systHistNamesPbPb, std::vector<std::string> systHistNamesPP, std::vector<bool> isCorrelated, std::vector<Int_t> yPos, std::vector<Double_t> scaleTotalPbPbSyst, std::vector<Double_t> scaleTotalPPSyst, std::string unfFileName)
{
  std::vector<std::vector<Double_t> > systVals;

  for(unsigned int jI = 0; jI < systHistNamesPbPb.size(); ++jI){
    TH2D* temp2DPbPb_p = nullptr;
    TH1D* temp1DPbPb_p = (TH1D*)nominalHistPbPb_p->Clone("temp1DPbPb_p");

    TH2D* temp2DPP_p = nullptr;
    TH1D* temp1DPP_p = (TH1D*)nominalHistPP_p->Clone("temp1DPP_p");

    for(Int_t bIX = 0; bIX < temp1DPbPb_p->GetXaxis()->GetNbins(); ++bIX){
      temp1DPbPb_p->SetBinContent(bIX+1, 0.0);
      temp1DPbPb_p->SetBinError(bIX+1, 0.0);

      temp1DPP_p->SetBinContent(bIX+1, 0.0);
      temp1DPP_p->SetBinError(bIX+1, 0.0);
    }

    if(systHistNamesPP[jI].find("UNFNONCLOSURE") == std::string::npos){
      temp2DPbPb_p = (TH2D*)histFilePbPb_p->Get(systHistNamesPbPb[jI].c_str());
      temp2DPP_p = (TH2D*)histFilePP_p->Get(systHistNamesPP[jI].c_str());

      fineTH2ToCoarseTH1(temp2DPbPb_p, temp1DPbPb_p, yPos);
      fineTH2ToCoarseTH1(temp2DPP_p, temp1DPP_p, yPos);

      binWidthAndScaleNorm(temp1DPbPb_p, scaleTotalPbPbSyst[jI]);
      binWidthAndScaleNorm(temp1DPP_p, scaleTotalPPSyst[jI]);
    }
    else{
      std::string gammaPtStr = nominalHistPP_p->GetName();
      if(gammaPtStr.find("GammaPt") != std::string::npos){
        gammaPtStr.replace(0, gammaPtStr.find("GammaPt"), "");

        while(gammaPtStr.find("_") != std::string::npos){
          gammaPtStr.replace(gammaPtStr.rfind("_"), gammaPtStr.size(), "");
        }
      }

      TFile* tempFile_p = new TFile(unfFileName.c_str(), "READ");
      TH1D* corrHist_p = (TH1D*)tempFile_p->Get(("unfoldOverTruthMC_" + gammaPtStr + "_PP_MCNonClosure").c_str());
      fineHistToCoarseHist(corrHist_p, temp1DPP_p);
      fineHistToCoarseHist(corrHist_p, temp1DPbPb_p);

      for(Int_t bIX = 0; bIX < temp1DPP_p->GetXaxis()->GetNbins(); ++bIX){
        Double_t value = nominalHistPP_p->GetBinContent(bIX+1)/temp1DPP_p->GetBinContent(bIX+1);
        Double_t error = nominalHistPP_p->GetBinError(bIX+1)/temp1DPP_p->GetBinContent(bIX+1);

        temp1DPP_p->SetBinContent(bIX+1, value);
        temp1DPP_p->SetBinError(bIX+1, error);

        value = nominalHistPbPb_p->GetBinContent(bIX+1)/temp1DPbPb_p->GetBinContent(bIX+1);
        error = nominalHistPbPb_p->GetBinError(bIX+1)/temp1DPbPb_p->GetBinContent(bIX+1);

        temp1DPbPb_p->SetBinContent(bIX+1, value);
        temp1DPbPb_p->SetBinError(bIX+1, error);
      }

      tempFile_p->Close();
      delete tempFile_p;
    }

    systVals.push_back({});
    for(Int_t bIX = 0; bIX < temp1DPP_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t deltaVal = 0.0;
      Float_t ratNom = nominalHistPbPb_p->GetBinContent(bIX+1)/nominalHistPP_p->GetBinContent(bIX+1);

      if(isCorrelated[jI]){
	Float_t ratSyst = temp1DPbPb_p->GetBinContent(bIX+1)/temp1DPP_p->GetBinContent(bIX+1);
	deltaVal += TMath::Abs(ratSyst - ratNom);
      }
      else{
	Float_t ratSystPbPb = temp1DPbPb_p->GetBinContent(bIX+1)/nominalHistPP_p->GetBinContent(bIX+1);
	Float_t ratSystPP = nominalHistPbPb_p->GetBinContent(bIX+1)/temp1DPP_p->GetBinContent(bIX+1);

	Float_t deltaValPbPb = TMath::Abs(ratSystPbPb - ratNom);
	Float_t deltaValPP = TMath::Abs(ratSystPP - ratNom);

	if(systHistNamesPbPb[jI].find("MCSTAT") != std::string::npos){
	  //Check the deltaVal is 0; i.e. things are working correctly
	  if(deltaValPbPb > TMath::Power(10,-50)){
	    std::cout << "getSyst ERROR!!!! DeltaVal of Nominal and MCStat variation is greater than 0, " << deltaValPbPb << " for Pb+Pb, return" << std::endl;
	    return {};
	  }

	  if(deltaValPP > TMath::Power(10,-50)){
	    std::cout << "getSyst ERROR!!!! DeltaVal of Nominal and MCStat variation is greater than 0, "<< deltaValPbPb <<" for pp, return" << std::endl;
	    return {};
	  }

	  //if things check out then just take the error bar as the deltaVal
	  ratSystPbPb = (temp1DPbPb_p->GetBinContent(bIX+1) + temp1DPbPb_p->GetBinError(bIX+1))/nominalHistPP_p->GetBinContent(bIX+1);
	  ratSystPP = nominalHistPbPb_p->GetBinContent(bIX+1)/(temp1DPP_p->GetBinContent(bIX+1) + temp1DPP_p->GetBinError(bIX+1));

	  deltaValPbPb = TMath::Abs(ratSystPbPb - ratNom);
	  deltaValPP = TMath::Abs(ratSystPP - ratNom);
	}

	deltaVal += TMath::Sqrt(deltaValPbPb*deltaValPbPb + deltaValPP*deltaValPP);
      }

      systVals[jI].push_back(deltaVal);
    }

    delete temp1DPbPb_p;
    delete temp1DPP_p;
  }

  return systVals;
}


std::vector<Double_t> drawSyst(TPad* pad_p, TH1D* nominalHist_p, std::vector<std::vector<Double_t> > systVals)
{
  std::vector<Double_t> quadSumSystVect;

  pad_p->cd();
  Int_t color = nominalHist_p->GetMarkerColor();
  TBox* box_p = new TBox();

  for(Int_t bIX = 0; bIX < nominalHist_p->GetXaxis()->GetNbins(); ++bIX){
    Float_t systQuadSum = 0.0;
    for(unsigned int sI = 0; sI < systVals.size(); ++sI){
      systQuadSum += systVals[sI][bIX]*systVals[sI][bIX];
    }
    systQuadSum = TMath::Sqrt(systQuadSum);

    quadSumSystVect.push_back(systQuadSum);

    if(nominalHist_p->GetBinContent(bIX+1) < TMath::Power(10,-100)) continue;

    Float_t binVal = nominalHist_p->GetBinContent(bIX+1);
    Float_t binLow = nominalHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
    Float_t binHigh = nominalHist_p->GetXaxis()->GetBinLowEdge(bIX+2);

    box_p->SetFillColorAlpha(color, 0.3);

    Double_t minVal = TMath::Max((Double_t)(binVal - systQuadSum), (Double_t)nominalHist_p->GetMinimum());
    Double_t maxVal = TMath::Min((Double_t)(binVal + systQuadSum), (Double_t)nominalHist_p->GetMaximum());

    box_p->DrawBox(binLow, minVal, binHigh, maxVal);
  }

  delete box_p;

  return quadSumSystVect;
}


std::vector<Double_t> plotSyst(TH1D* nominalHist_p, std::vector<std::vector<Double_t> > systVals, std::vector<std::string> systNames, std::vector<std::string> labels, std::string yStr, std::string saveName, std::string totalStr = "Total", float systMin = -1.0, float systMax = -1.0)
{
  std::vector<Double_t> systTotal;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  //Series of consistency checks
  const Int_t nMaxHist = 36;
  if(systVals.size() > nMaxHist){
    std::cout << "Needed number of histograms '" << systVals.size()+1 << "' exceeds plotSyst nMaxHist="  << nMaxHist << ". return {}" << std::endl;
    return {};
  }

  if(systVals.size() != systNames.size()){
    std::cout << "SystVals and systNames have different sizes \'" << systVals.size() << "\', \'" << systNames.size() << "\', respectively. return {}" << std::endl;
    return {};
  }

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

  //Define some plotting parameters
  const int titleFont = 42;
  const double titleSize = 0.04;
  const double labelSize = titleSize*1.0;

  const Float_t height = 900.0;
  const Float_t width = 1200.0;

  const Double_t leftMargin = 0.14;
  const Double_t bottomAbsFracOfLeft = 0.8;
  const Double_t bottomMargin = bottomAbsFracOfLeft*leftMargin*width/height;

  //Define a bin array and grab the relevant bin edges; check it is within array size
  const Int_t nBinsMax = 100;
  Int_t nBinsX = nominalHist_p->GetXaxis()->GetNbins();
  if(nBinsX > nBinsMax){
    std::cout << "nBinsX needed in plotSyst \'" << nBinsX << "\' exceeds binsMax \'" << nBinsMax << "\'. return {}" << std::endl;
    return {};
  }
  Float_t binsX[nBinsMax];
  for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
    binsX[bIX] = nominalHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
  }

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

  //Define our histogram array and new the total syst hist
  std::string titleStr = ";" + std::string(nominalHist_p->GetXaxis()->GetTitle()) + ";" + yStr + " (Relative Systematic Uncertainty)";
  TH1D* systHist_p[nMaxHist+1];
  systHist_p[0] = new TH1D("systHist_Total_h", titleStr.c_str(), nBinsX, binsX);
  HIJet::Style::EquipHistogram(systHist_p[0], 0);

  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    systTotal.push_back(0.0);
  }

  for(unsigned int sI = 0; sI < systVals.size(); ++sI){
    systHist_p[sI+1] = new TH1D(("systHist_" + systNames[sI] + "_h").c_str(), titleStr.c_str(), nBinsX, binsX);

    for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
      if(nominalHist_p->GetBinContent(bIX+1) < TMath::Power(10,-100)) continue;

      systHist_p[sI+1]->SetBinContent(bIX+1, systVals[sI][bIX]/nominalHist_p->GetBinContent(bIX+1));
      systHist_p[sI+1]->SetBinError(bIX+1, 0);

      systTotal[bIX] += systVals[sI][bIX]*systVals[sI][bIX];
    }

    HIJet::Style::EquipHistogram(systHist_p[sI+1], sI+1);
    systHist_p[sI+1]->SetMarkerSize(1.2*systHist_p[sI+1]->GetMarkerSize());
    systHist_p[sI+1]->SetMarkerSize(1.2*systHist_p[sI+1]->GetMarkerSize());
  }

  Double_t maxVal = 0.0;
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    if(nominalHist_p->GetBinContent(bIX+1) < TMath::Power(10,-100)) continue;

    Float_t val = TMath::Sqrt(systTotal[bIX])/nominalHist_p->GetBinContent(bIX+1);
    if(val > maxVal) maxVal = val;

    systHist_p[0]->SetBinContent(bIX+1, val);
    systHist_p[0]->SetBinError(bIX+1, 0.0);
  }

  //  if(maxVal > 1.0) maxVal = 1.0;

  if(systMax < 0.0){
    systHist_p[0]->SetMaximum(maxVal*1.2);
    //systHist_p[0]->SetMinimum(-maxVal*1.2);
    //Changed for getting rid of the negative values given everything is being handled as symmetric
    systHist_p[0]->SetMinimum(-maxVal*0.2);
  }
  else{
    systHist_p[0]->SetMaximum(systMax);
    systHist_p[0]->SetMinimum(systMin);
  }

  TCanvas* canv_p = new TCanvas("canvSyst_p", "", width, height);
  Double_t topMargin = 0.20;
  Double_t rightMargin = 0.25;
  if(saveName.find("SystType") != std::string::npos){
    topMargin = 0.02;
    rightMargin = 0.02;
  }

  setMargins(canv_p, leftMargin, topMargin, rightMargin, 2.0*bottomMargin/3.0);

  canv_p->cd();

  systHist_p[0]->GetXaxis()->SetTitleSize(titleSize*1.25);
  systHist_p[0]->GetYaxis()->SetTitleSize(titleSize);
  systHist_p[0]->GetYaxis()->SetLabelSize(labelSize);
  systHist_p[0]->GetXaxis()->SetLabelSize(labelSize);
  systHist_p[0]->GetXaxis()->SetTitleOffset(0.9);

  for(unsigned int sI = 0; sI < systVals.size()+1; ++sI){
    if(sI == 0) systHist_p[sI]->DrawCopy("HIST");
    else systHist_p[sI]->DrawCopy("HIST SAME");

    systHist_p[sI]->DrawCopy("P SAME");

    //Create the mirrored version for negative values
    //You are symmetrizing everything rn so dont bother with this
    /*
    for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
      systHist_p[sI]->SetBinContent(bIX+1, -systHist_p[sI]->GetBinContent(bIX+1));
    }
    */
    //    systHist_p[sI]->DrawCopy("HIST P SAME");
  }

  systHist_p[0]->SetLineStyle(2);
  systHist_p[0]->DrawCopy("HIST SAME");
  systHist_p[0]->SetLineStyle(1);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSize);
  //  label_p->SetTextAlign(31);

  Float_t xPos = 0.01;
  Float_t yPos = 0.96;
  Float_t deltaY = 0.06;
  if(saveName.find("SystType") != std::string::npos){
    xPos = 0.185;
    yPos = 0.92;
  }


  Float_t maxLength = 0.0;
  for(unsigned int i = 0; i < labels.size(); ++i){
    if(i != 0 && i%3 == 0){
      xPos += maxLength + 0.03;
      maxLength = 0;
    }

    label_p->SetText(xPos, yPos - deltaY*(i%3), labels[i].c_str());
    if(label_p->GetXsize() > maxLength) maxLength = label_p->GetXsize();

    label_p->DrawLatex(xPos, yPos - deltaY*(i%3), labels[i].c_str());
  }
  delete label_p;

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

  std::vector<TLegend*> legs_p = {nullptr, nullptr};

  if(systVals.size() <= nMaxHist/2){
    if(saveName.find("SystType") != std::string::npos) legs_p[0] = new TLegend(0.46, 0.76 - (systVals.size()+1)*0.045, 0.60, 0.76);
    else legs_p[0] = new TLegend(0.76, 0.81 - (systVals.size()+1)*0.045, 0.90, 0.81);
  }
  else{
    legs_p[0] = new TLegend(0.76, 0.81 - (nMaxHist/2)*0.045, 0.87, 0.81);
    legs_p[1] = new TLegend(0.87, 0.81 - (systVals.size() - nMaxHist/2 + 1)*0.045 - 0.045, 0.99, 0.81 - 0.045);

    legs_p[1]->SetTextFont(titleFont);
    legs_p[1]->SetTextSize(titleSize);
    legs_p[1]->SetBorderSize(0);
    legs_p[1]->SetFillColor(0);
  }
  legs_p[0]->SetTextFont(titleFont);
  legs_p[0]->SetTextSize(titleSize);
  legs_p[0]->SetBorderSize(0);
  legs_p[0]->SetFillColor(0);
  legs_p[0]->SetFillStyle(0);

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << ", " << systVals.size() << ", " << systNames.size() << std::endl;

  legs_p[0]->AddEntry(systHist_p[0], totalStr.c_str(), "L");
  for(unsigned int sI = 0; sI < systVals.size(); ++sI){
    legs_p[(sI+1)/(nMaxHist/2)]->AddEntry(systHist_p[sI+1], systNameToLegName(systNames[sI]).c_str(), "L P");
  }

  legs_p[0]->Draw("SAME");
  if(systVals.size() > nMaxHist/2) legs_p[1]->Draw("SAME");

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(binsX[0], 0.0, binsX[nBinsX], 0.0);
  delete line_p;

  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;
  gPad->SetTicks();

  quietSaveAs(canv_p, saveName);
  delete canv_p;
  delete legs_p[0];
  if(systVals.size() > nMaxHist/2) delete legs_p[1];

  for(unsigned int sI = 0; sI < systVals.size()+1; ++sI){delete systHist_p[sI];}

  return systTotal;
}


void constructJETSCAPEHist(std::string fileName, TH1D** histPointer_p, int color)
{
  const Int_t nMaxPtBins = 100;
  std::ifstream inJETSCAPEFile(fileName.c_str());
  std::string binStr;

  Int_t nBinsJETSCAPE = 0;
  Double_t binsJETSCAPE[nMaxPtBins];
  std::vector<double> values;
  std::vector<double> errors;

  while(std::getline(inJETSCAPEFile, binStr)){
    while(binStr.find("  ") != std::string::npos) binStr.replace(binStr.find("  "), 1,"");
    while(binStr.substr(0,1).find(" ") != std::string::npos) binStr.replace(0,1,"");
    while(binStr.substr(binStr.size()-1,1).find(" ") != std::string::npos) binStr.replace(binStr.size()-1,1,"");
    while(binStr.find(" ") != std::string::npos) binStr.replace(binStr.find(" "), 1,",");

    std::vector<double> binVect = strToVectD(binStr);
    //       std::cout << "BIN STR: " << binStr << std::endl;
    //       for(unsigned int i = 0; i < binVect.size(); ++i){
    //	 std::cout << " " << i << ": " << binVect[i] << std::endl;
    //       }

    /*
    if(fileName.find("Jetscape") != std::string::npos){
      if(binVect[2] >= 1.3) break;
    }
    */

    if(nBinsJETSCAPE == 0) binsJETSCAPE[0] = binVect[1];
    ++nBinsJETSCAPE;
    binsJETSCAPE[nBinsJETSCAPE] = binVect[2];
    values.push_back(binVect[3]);
    errors.push_back(binVect[4]);
  }
  inJETSCAPEFile.close();

  std::string histName = fileName;
  while(histName.find("/") != std::string::npos){
    histName.replace(0, histName.find("/")+1, "");
  }

  (*histPointer_p) = new TH1D(histName.c_str(), ";;", nBinsJETSCAPE, binsJETSCAPE);
  for(unsigned int jI = 0; jI < values.size(); ++jI){
    (*histPointer_p)->SetBinContent(jI+1, values[jI]);
    (*histPointer_p)->SetBinError(jI+1, errors[jI]);
  }

  (*histPointer_p)->SetMarkerStyle(1);
  (*histPointer_p)->SetMarkerSize(0.1);
  (*histPointer_p)->SetLineWidth(4);
  (*histPointer_p)->SetMarkerColor(color);
  (*histPointer_p)->SetLineColor(color);
  if(fileName.find("Jetscape") != std::string::npos) (*histPointer_p)->SetFillColorAlpha(color, 0.3);

  (*histPointer_p)->Print("ALL");

  return;
}


int gdjPlotResults(std::string inConfigFileName)
{
  const int titleFont = 42;
  const double titleSize = 0.0325;
  const double labelSize = titleSize*1.0;
  const double yOffset = 1.2;

  cppWatch globalTimer;
  globalTimer.start();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  //Define necessary params for the input config file
  std::vector<std::string> necessaryParams = {"INPBPBFILENAME",
					      "INPPFILENAME",
					      "INPBPBTERMFILENAME",
					      "INPPTERMFILENAME",
					      "INUNFNONCLOSUREFILENAME",
					      "DOJTPTCUT",
					      "DOSIMPLESYST",
					      "SAVETAG",
					      "SAVEEXT"};

  //Define a second set once we've accessed the file and know what var we are working with
  std::vector<std::string> necessaryVarParams = {"VARMINY",
						 "VARMAXY",
						 "VARRATMINY",
						 "VARRATMAXY",
						 "VARPYTRATMINY",
						 "VARPYTRATMAXY",
						 "VARDOLOGX",
						 "VARDOLOGY",
						 "VARLEGX",
						 "VARLEGY",
  						 "VARLABELX",
						 "VARLABELY",
  						 "VARRATLEGX",
						 "VARRATLEGY",
  						 "VARRATLABELX",
						 "VARRATLABELY",
						 "VAR2JETXLABELX",
						 "VAR2JETXLABELY",
						 "VARSYSTMAX",
						 "VARSYSTMIN"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  const std::string inPbPbFileName = config_p->GetValue("INPBPBFILENAME", "");
  const std::string inPPFileName = config_p->GetValue("INPPFILENAME", "");

  const std::string inPbPbTermFileName = config_p->GetValue("INPBPBTERMFILENAME", "");
  const std::string inPPTermFileName = config_p->GetValue("INPPTERMFILENAME", "");

  const std::string inUnfoldNonClosureFileName = config_p->GetValue("INUNFNONCLOSUREFILENAME", "");

  const std::string inJETSCAPEPPFile = config_p->GetValue("JETSCAPEPP", "");
  const std::string inJETSCAPEPbPbFile = config_p->GetValue("JETSCAPEPBPB", "");

  const std::string inLBTPPFile = config_p->GetValue("LBTPP", "");
  const std::string inLBTPbPbFile = config_p->GetValue("LBTPBPB", "");

  //Check all files exist before continuing to grab params
  if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;
  if(!check.checkFileExt(inPPFileName, ".root")) return 1;
  if(!check.checkFileExt(inPbPbTermFileName, ".config")) return 1;
  if(!check.checkFileExt(inPPTermFileName, ".config")) return 1;
  if(!check.checkFileExt(inUnfoldNonClosureFileName, ".root")) return 1;

  const Bool_t doJtPtCut = config_p->GetValue("DOJTPTCUT", 0);
  const Bool_t doSimpleSyst = config_p->GetValue("DOSIMPLESYST", 0);

  TEnv* inPbPbTermFileConfig_p = new TEnv(inPbPbTermFileName.c_str());
  TEnv* inPPTermFileConfig_p = new TEnv(inPPTermFileName.c_str());

  std::vector<std::string> termNecessaryParams = {"UNFOLDFILENAME",
						  "VARNAME"};

  if(!checkEnvForParams(inPbPbTermFileConfig_p, termNecessaryParams)) return 1;
  if(!checkEnvForParams(inPPTermFileConfig_p, termNecessaryParams)) return 1;
  if(!compEnvParams(inPbPbTermFileConfig_p, inPPTermFileConfig_p, {"VARNAME"})) return 1;

  const std::string inPbPbTempFileName = inPbPbTermFileConfig_p->GetValue("UNFOLDFILENAME", "");
  const std::string inPPTempFileName = inPPTermFileConfig_p->GetValue("UNFOLDFILENAME", "");

  if(!isStrSame(inPbPbTempFileName, inPbPbFileName)){
    std::cout << "Given unfolding termination file for Pb+Pb \'" << inPbPbTermFileName << "\' is based on \'" << inPbPbTempFileName << "\', not matching given input \'" << inPbPbFileName << "\'. return 1" << std::endl;
    return 1;
  }
  if(!isStrSame(inPPTempFileName, inPPFileName)){
    std::cout << "Given unfolding termination file for p+p \'" << inPPTermFileName << "\' is based on \'" << inPPTempFileName << "\', not matching given input \'" << inPPFileName << "\'. return 1" << std::endl;
    return 1;
  }

  const std::string termVarName = returnAllLowercaseString(inPbPbTermFileConfig_p->GetValue("VARNAME", ""));

  const std::string saveTag = config_p->GetValue("SAVETAG", "");
  const std::string saveExt = config_p->GetValue("SAVEEXT", "");

  std::vector<std::string> validExt = {"pdf", "png"};

  if(!vectContainsStr(saveExt, &validExt)){
    std::cout << "Given SAVEEXT \'" << saveExt << "\' is not valid. return 1" << std::endl;
    return 1;
  }

 //We should check the file for possible other inputs, like JEWEL for overlay
 std::string inJEWELPPFileName = config_p->GetValue("JEWELPPFILENAME", "");
 std::vector<std::string> inJEWELPbPbFileNames;
 //Just cap it at ten to check for inputs - should never be more than that
 const Int_t nMaxPbPbJEWELFiles = 10;
 for(Int_t i = 0; i < nMaxPbPbJEWELFiles; ++i){
   std::string inJEWELPbPbFileName = config_p->GetValue(("JEWELPBPBFILENAME." + std::to_string(i)).c_str(), "");
   if(check.checkFileExt(inJEWELPbPbFileName, ".root")) inJEWELPbPbFileNames.push_back(inJEWELPbPbFileName);
 }

 //Scaling factor taking from HEP 13 TeV Photon+multijet
 // const Double_t scaleFactorPYT = 1.0;
 //Do jewel overlay only if pp and at least one valid Pb+Pb file exist
 bool doJEWEL = check.checkFileExt(inJEWELPPFileName, ".root") && inJEWELPbPbFileNames.size() > 0;

 bool doJETSCAPE = check.checkFileExt(inJETSCAPEPPFile, ".txt") && check.checkFileExt(inJETSCAPEPbPbFile, ".txt");

 bool doLBT = check.checkFileExt(inLBTPPFile, ".txt") && check.checkFileExt(inLBTPbPbFile, ".txt");

 //Check that the jewel inputs are matched to our input files (this guarantees things like same-binning
 TFile* ppJEWELFile_p = nullptr;
 TEnv* ppJEWELConfig_p = nullptr;
 TFile* pbpbJEWELFile_p[nMaxPbPbJEWELFiles];
 TEnv* pbpbJEWELConfig_p[nMaxPbPbJEWELFiles];
 if(doJEWEL){
   ppJEWELFile_p = new TFile(inJEWELPPFileName.c_str(), "READ");
   ppJEWELConfig_p = (TEnv*)ppJEWELFile_p->Get("config");

   std::string ppUnfoldFileFromJEWEL = ppJEWELConfig_p->GetValue("INUNFOLDFILENAME", "");
   if(!isStrSame(ppUnfoldFileFromJEWEL, inPPFileName)){
     std::cout << "Given JEWEL file was created from unfold file \'" << ppUnfoldFileFromJEWEL << "\' does not match input unfold file \'" << inPPFileName << "\'. Fix for JEWEL overlay. doJEWEL set to false"  << std::endl;
     doJEWEL = false;
   }

   if(doJEWEL){
     for(unsigned int jI = 0; jI < inJEWELPbPbFileNames.size(); ++jI){
 	pbpbJEWELFile_p[jI] = new TFile(inJEWELPbPbFileNames[jI].c_str(), "READ");
 	pbpbJEWELConfig_p[jI] = (TEnv*)pbpbJEWELFile_p[jI]->Get("config");

 	std::string pbpbUnfoldFileFromJEWEL = pbpbJEWELConfig_p[jI]->GetValue("INUNFOLDFILENAME", "");

 	if(!isStrSame(pbpbUnfoldFileFromJEWEL, inPbPbFileName)){
 	  std::cout << "Given JEWEL file was created from unfold file \'" << pbpbUnfoldFileFromJEWEL << "\' does not match input unfold file \'" << inPbPbFileName << "\'. Fix for JEWEL overlay.  doJEWEL set to false" << std::endl;
 	  doJEWEL = false;
 	}
     }
   }
 }



 //Get global labels
 std::vector<std::string> globalLabels;
 for(unsigned int i = 0; i < 100; ++i){
   std::string tempStr = config_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
   if(tempStr.size() != 0) globalLabels.push_back(tempStr);
 }

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;
 if(!check.checkFileExt(inPPFileName, ".root")) return 1;

 TFile* inPbPbFile_p = new TFile(inPbPbFileName.c_str(), "READ");
 TEnv* inPbPbFileConfig_p = (TEnv*)inPbPbFile_p->Get("config");
 TEnv* inPbPbFileLabel_p = (TEnv*)inPbPbFile_p->Get("label");

 TFile* inPPFile_p = new TFile(inPPFileName.c_str(), "READ");
 TEnv* inPPFileConfig_p = (TEnv*)inPPFile_p->Get("config");
 //  TEnv* inPPFileLabel_p = (TEnv*)inPPFile_p->Get("label");

 std::vector<std::string> paramsToCompare = {"ISMC",
 					      "NGAMMAPTBINS",
 					      "GAMMAPTBINSLOW",
 					      "GAMMAPTBINSHIGH",
 					      "GAMMAPTBINSDOLOG",
 					      "GAMMAPTBINSDOCUSTOM",
 					      "GAMMAPTBINSCUSTOM",
 					      "NSUBJTPTBINS",
 					      "SUBJTPTBINSLOW",
 					      "SUBJTPTBINSHIGH",
 					      "SUBJTPTBINSDOLOG",
 					      "SUBJTPTBINSDOCUSTOM",
 					      "SUBJTPTBINSCUSTOM",
 					      "SUBJTGAMMAPTMIN",
 					      "SUBJTGAMMAPTMAX",
 					      "JTPTBINSLOWRECO",
 					      "NITER",
 					      "JETR",
 					      "VARNAME"};

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 std::vector<std::string> necessaryFileParams = paramsToCompare;
 necessaryFileParams.push_back("ISPP");

 const std::string varName = returnAllLowercaseString(inPbPbFileConfig_p->GetValue("VARNAME", ""));
 if(!isStrSame(termVarName, varName)){
   std::cout << "Varname of input files \'" << varName << "\' does not match the varName of termination files \'" << termVarName << "\'. return 1" << std::endl;
   return 1;
 }
 const std::string varNameUpper = strLowerToUpper(varName);
 const std::string varNameLower = returnAllLowercaseString(varName);
 const std::string varNameLabel = varNameToLabel(varName);

 const Bool_t isMultijet = isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");


 //Now that we have the varName, fix necessaryVarParams and check our input config has them before getting the values
 for(unsigned int vI = 0; vI < necessaryVarParams.size(); ++vI){
   std::string repStr = necessaryVarParams[vI];
   repStr.replace(repStr.find("VAR"), 3, varNameUpper);
   necessaryVarParams[vI] = repStr;
 }
 if(!checkEnvForParams(config_p, necessaryVarParams)) return 1;

 //Only optional value truncating the x-axis range; hence the absurd default values
 float minValX = config_p->GetValue((varNameUpper + "MINX").c_str(), -1000.0);
 float maxValX = config_p->GetValue((varNameUpper + "MAXX").c_str(), -10000.0);
 bool doXTrunc = minValX < maxValX;


 //Min. necessary values for plotting
 float minVal = config_p->GetValue((varNameUpper + "MINY").c_str(), 0.0);
 float maxVal = config_p->GetValue((varNameUpper + "MAXY").c_str(), 0.0);
 float minRatVal = config_p->GetValue((varNameUpper + "RATMINY").c_str(), 0.0);
 float maxRatVal = config_p->GetValue((varNameUpper + "RATMAXY").c_str(), 0.0);
 float minPytRatVal = config_p->GetValue((varNameUpper + "PYTRATMINY").c_str(), 0.0);
 float maxPytRatVal = config_p->GetValue((varNameUpper + "PYTRATMAXY").c_str(), 0.0);
 bool doLogX = config_p->GetValue((varNameUpper + "DOLOGX").c_str(), 0);
 bool doLogY = config_p->GetValue((varNameUpper + "DOLOGY").c_str(), 0);
 float legX = config_p->GetValue((varNameUpper + "LEGX").c_str(), 0.1);
 float legY = config_p->GetValue((varNameUpper + "LEGY").c_str(), 0.1);
 float labelX = config_p->GetValue((varNameUpper + "LABELX").c_str(), 0.1);
 float labelY = config_p->GetValue((varNameUpper + "LABELY").c_str(), 0.1);

 float systMax = config_p->GetValue((varNameUpper + "SYSTMAX").c_str(), -1.0);
 float systMin = config_p->GetValue((varNameUpper + "SYSTMIN").c_str(), -1.0);

 float systPPMax = config_p->GetValue((varNameUpper + "SYSTPPMAX").c_str(), -1.0);
 float systPPMin = config_p->GetValue((varNameUpper + "SYSTPPMIN").c_str(), -1.0);

 float var2JetXLabelX = config_p->GetValue((varNameUpper + "2JETXLABELX").c_str(), 0.1);
 float var2JetXLabelY = config_p->GetValue((varNameUpper + "2JETXLABELY").c_str(), 0.1);

 float ratLegX = config_p->GetValue((varNameUpper + "RATLEGX").c_str(), 0.1);
 float ratLegY = config_p->GetValue((varNameUpper + "RATLEGY").c_str(), 0.1);
 float ratLabelX = config_p->GetValue((varNameUpper + "RATLABELX").c_str(), 0.1);
 float ratLabelY = config_p->GetValue((varNameUpper + "RATLABELY").c_str(), 0.1);

 float legXCent = config_p->GetValue((varNameUpper + "LEGXCENT").c_str(), legX);
 float legYCent = config_p->GetValue((varNameUpper + "LEGYCENT").c_str(), legY);

 float labelXCent = config_p->GetValue((varNameUpper + "LABELXCENT").c_str(), labelX);
 float labelYCent = config_p->GetValue((varNameUpper + "LABELYCENT").c_str(), labelY);
 float labelAlignRightCent = config_p->GetValue((varNameUpper + "LABELALIGNRIGHTCENT").c_str(), 1);

 std::string varPrefix = "";
 std::string varNameHist = varNameUpper;
 if(isStrSame(varNameUpper, "PT")){
   varNameHist = "Pt";
   varPrefix = "JT";
 }

 std::string binVarStr = varPrefix + varNameUpper;
 if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";
 if(isStrSame(binVarStr, "DRJJ")) binVarStr = "DR";

 const Int_t nMaxPtBins = 200;
 Int_t nVarBins = inPbPbFileConfig_p->GetValue(("N" + binVarStr + "BINS").c_str(), -1);
 Float_t varBinsLow = inPbPbFileConfig_p->GetValue((binVarStr + "BINSLOW").c_str(), -1);
 Float_t varBinsHigh = inPbPbFileConfig_p->GetValue((binVarStr + "BINSHIGH").c_str(), -1);
 Bool_t varBinsDoLog = inPbPbFileConfig_p->GetValue((binVarStr + "BINSDOLOG").c_str(), 0);
 Bool_t varBinsDoCustom = inPbPbFileConfig_p->GetValue((binVarStr + "BINSDOCUSTOM").c_str(), 0);
 std::string varBinsStr;
 Double_t varBins[nMaxPtBins+1];
 Int_t nVarBinsTrunc = -1;
 Double_t varBinsTrunc[nMaxPtBins+1];

 std::string gammaJtDPhiCutLabel = inPbPbFileConfig_p->GetValue("GAMMAJTDPHI", "");
 if(gammaJtDPhiCutLabel.size() == 0){
   std::cout << "GammaJtDPhi not found, check line L" << __LINE__ << ". return 1" << std::endl;
   return 1;
 }
 if(gammaJtDPhiCutLabel.find("pi") != std::string::npos) gammaJtDPhiCutLabel.replace(gammaJtDPhiCutLabel.find("pi"), 2, "#pi");
 gammaJtDPhiCutLabel = "|#Delta#phi_{#gamma,jet}|>" + gammaJtDPhiCutLabel;

 //Grab the bins and throw them in the varBins array
 if(varBinsDoLog) getLogBins(varBinsLow, varBinsHigh, nVarBins, varBins);
 else if(varBinsDoCustom){
   varBinsStr = inPbPbFileConfig_p->GetValue((binVarStr + "BINSCUSTOM").c_str(), "");
   if(varBinsStr.size() == 0){
     std::cout << "No custom bins given. return 1" << std::endl;
     return 1;
   }
   std::vector<float> varBinsTemp = strToVectF(varBinsStr);
   Int_t nVarBinsTemp = ((Int_t)varBinsTemp.size()) - 1;
   if(nVarBins != nVarBinsTemp){
     return 1;
   }

   for(unsigned int gI = 0; gI < varBinsTemp.size(); ++gI){
     varBins[gI] = varBinsTemp[gI];
   }
 }
 else getLinBins(varBinsLow, varBinsHigh, nVarBins, varBins);

 //if an x-axis truncation is requested, test that given values match array
 if(doXTrunc){
   if(!checkHistContainsBins({minValX, maxValX}, nVarBins, varBins, 0.0001)){
    std::cout << "gdjPlotResults ERROR: Requested plotting truncation in var \'" << varNameUpper << "\' to range " << minValX << "-" << maxValX << " is not found in bins \'" << varBinsStr << "\'. return 1" << std::endl;

     return 1;
   }

   nVarBinsTrunc = truncateBinRange(nVarBins, varBins, minValX, maxValX, 0.0001, varBinsTrunc);
   /*
   std::cout << "Original bins: " << std::endl;
   for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
     std::cout << " " << varBins[bIX] << ",";
   }
   std::cout << " " <<  varBins[nVarBins] << "." << std::endl;

   std::cout << std::endl;
   std::cout << "New bins: " << std::endl;
   for(Int_t bIX = 0; bIX < nVarBinsTrunc; ++bIX){
     std::cout << " " << varBinsTrunc[bIX] << ",";
   }
   std::cout << " " << varBinsTrunc[nVarBinsTrunc] << ".";

   return 1;
   */
 }

 const std::string gammaJtDPhiStr = "DPhi0";
 //  const std::string multiJtCutGlobalStr = "MultiJt0";
 const std::string jtPtBinsGlobalStr = "GlobalJtPt0";

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 std::vector<std::string> pbpbParams = {"CENTBINS"};
 if(!checkEnvForParams(inPbPbFileConfig_p, necessaryFileParams)) return 1;
 if(!checkEnvForParams(inPPFileConfig_p, necessaryFileParams)) return 1;
 if(!compEnvParams(inPbPbFileConfig_p, inPPFileConfig_p, paramsToCompare)) return 1;

 // const Int_t nMaxSyst = 100;
 //Grab all systs
 Int_t nSyst = inPPFileConfig_p->GetValue("NSYST", -1);
 std::vector<std::string> systStrVect = strToVect(inPPFileConfig_p->GetValue("SYST", ""));
 std::vector<std::string> systTypeVect = strToVect(inPPFileConfig_p->GetValue("SYSTTYPE", ""));
 std::map<std::string, int> systTypesMap;
 std::map<std::string, std::string> systStrToSystType;
 std::vector<std::string> uniqueSystTypes;

 //We need to first insert the unfolding syst
 ++nSyst;
 systStrVect.push_back("UNFNONCLOSURE");
 systTypeVect.push_back("UNFOLDING");

 //Do some prelim splitting of JESJER
 const Int_t nMaxJESJERFORLEG = 11;
 Int_t nJESJERCurr = 0;
 if(!doSimpleSyst){
   for(unsigned int sI = 0; sI < systTypeVect.size(); ++sI){
     if(isStrSame(systTypeVect[sI], "JESJER")){
       if(systStrVect[sI].substr(0,3).find("JES") != std::string::npos || systStrVect[sI].substr(0,2).find("QG") != std::string::npos){
	 systTypeVect[sI] = "JESSET" + std::to_string(nJESJERCurr/nMaxJESJERFORLEG);
	 ++nJESJERCurr;
       }
       else if(systStrVect[sI].substr(0,3).find("JER") != std::string::npos){
	 systTypeVect[sI] = "JER";
       }
     }
   }
 }

 for(int systI = 0; systI < nSyst; ++systI){
   if(systTypeVect[systI].find("NOMINAL") != std::string::npos) continue;//This is not a syst but a placeholder for convenience

   systStrToSystType[systStrVect[systI]] = systTypeVect[systI];
   ++(systTypesMap[systTypeVect[systI]]);
 }

 for(auto const & val : systTypesMap){
   uniqueSystTypes.push_back(val.first);
 }

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 //Add to global labels
 int jetR = inPbPbFileConfig_p->GetValue("JETR", 100);
 std::string jetRStr = "anti-k_{t} R=";
 if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
 else{
   jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
 }
 globalLabels.push_back(jetRStr);
 globalLabels.push_back(jtPtBinsGlobalStr);
 globalLabels.push_back(gammaJtDPhiStr);

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 configParser labels(inPbPbFileLabel_p);
 std::map<std::string, std::string> labelMap = labels.GetConfigMap();

 for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
   if(labelMap.count(globalLabels[gI]) != 0){
     globalLabels[gI] = labelMap[globalLabels[gI]];
   }
 }

 //Grab first set of parameters defined in necessaryInFileParams
 const bool isPP = inPPFileConfig_p->GetValue("ISPP", 0);
 const bool isPbPb = !(inPbPbFileConfig_p->GetValue("ISPP", 0));

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 if(!isPP){
   std::cout << "Given pp file is not pp return 1" << std::endl;
   return 1;
 }
 if(!isPbPb){
   std::cout << "Given pbpb file is not pbpb return 1" << std::endl;
   return 1;
 }

 std::vector<std::string> tempCentBinsStr = strToVect(inPbPbFileConfig_p->GetValue("CENTBINS", ""));
 std::vector<std::string> centBinsStr ;
 std::vector<std::string> centBinsLabel;
 std::vector<int> centColors;
 int ppColor = -1;
 //pp custom color object
 TColor* ppColor_p = new TColor();
 ppColor = ppColor_p->GetColor(0, 34, 85);

 int pbpbColor = -1;
 TColor* pbpbColor_p = new TColor();
 // pbpbColor = pbpbColor_p->GetColor(255, 51, 51);//kRed-4, GammaTest
 pbpbColor = pbpbColor_p->GetColor(191, 38, 38);//kRed-4, 3/4 values, DeltaTest

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 for(unsigned int cI = 0; cI < tempCentBinsStr.size()-1; ++cI){
   centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
   centBinsLabel.push_back("Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%");

   Int_t centPos = -1;
   if(centBinsStr[cI].find("Cent0to10") != std::string::npos) centPos = 1;
   else if(centBinsStr[cI].find("Cent10to30") != std::string::npos) centPos = 2;
   else if(centBinsStr[cI].find("Cent30to80") != std::string::npos) centPos = 3;

   TH1D* dummyHist_p = new TH1D("dummyHist_p", "", 10, 0, 1);
   HIJet::Style::EquipHistogram(dummyHist_p, centPos);

   centColors.push_back(dummyHist_p->GetMarkerColor());

   /*
   if(cI == 0){
     HIJet::Style::EquipHistogram(dummyHist_p, 2);
     ppColor = dummyHist_p->GetMarkerColor();
   }
   */

   delete dummyHist_p;
 }

 delete ppColor_p;

 //Define a bunch of parameters for the TCanvas
 const Int_t nPad = 2;
 Double_t padSplit = 0.45;
 const Double_t leftMargin = 0.17;
 const Double_t bottomMargin = 0.125/padSplit;
 const Double_t topMargin = 0.02;
 const Double_t rightMargin = 0.03;
 const Double_t noMargin = 0.001;
 Double_t height = 900.0;
 Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);
 height = 1200.0;

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 const Int_t nGammaPtBins = inPbPbFileConfig_p->GetValue("NGAMMAPTBINS", 0);
 const Float_t gammaPtBinsLow = inPbPbFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
 const Float_t gammaPtBinsHigh = inPbPbFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
 const Bool_t gammaPtBinsDoLog = (bool)inPbPbFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
 const Bool_t gammaPtBinsDoCustom = (bool)inPbPbFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
 Double_t gammaPtBins[nMaxPtBins+1];

 if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
 else if(gammaPtBinsDoCustom){
   std::string gammaPtBinsStr = inPbPbFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
   if(gammaPtBinsStr.size() == 0){
     std::cout << "No custom bins given. return 1" << std::endl;
     return 1;
   }
   std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
   Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
   if(nGammaPtBins != nGammaPtBinsTemp){
     std::cout << "GAMMAPTBINS DONT MATCH: " << std::endl;
     return 1;
   }

   for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
     gammaPtBins[gI] = gammaPtBinsTemp[gI];
   }
 }
 else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

 const Int_t nSubJtPtBins = inPbPbFileConfig_p->GetValue("NSUBJTPTBINS", 0);
 const Float_t subJtPtBinsLow = inPbPbFileConfig_p->GetValue("SUBJTPTBINSLOW", 80.0);
 const Float_t subJtPtBinsHigh = inPbPbFileConfig_p->GetValue("SUBJTPTBINSHIGH", 280.0);
 const Bool_t subJtPtBinsDoLog = (bool)inPbPbFileConfig_p->GetValue("SUBJTPTBINSDOLOG", 1);
 const Bool_t subJtPtBinsDoCustom = (bool)inPbPbFileConfig_p->GetValue("SUBJTPTBINSDOCUSTOM", 1);
 Double_t subJtPtBins[nMaxPtBins+1];

 if(subJtPtBinsDoLog) getLogBins(subJtPtBinsLow, subJtPtBinsHigh, nSubJtPtBins, subJtPtBins);
 else if(subJtPtBinsDoCustom){
   std::string subJtPtBinsStr = inPbPbFileConfig_p->GetValue("SUBJTPTBINSCUSTOM", "");
   if(subJtPtBinsStr.size() == 0){
     std::cout << "No custom bins given. return 1" << std::endl;
     return 1;
   }
   std::vector<float> subJtPtBinsTemp = strToVectF(subJtPtBinsStr);
   Int_t nSubJtPtBinsTemp = ((Int_t)subJtPtBinsTemp.size()) - 1;
   if(nSubJtPtBins != nSubJtPtBinsTemp){
     std::cout << "SUBJTPTBINS DONT MATCH: " << std::endl;
     return 1;
   }

   for(unsigned int gI = 0; gI < subJtPtBinsTemp.size(); ++gI){
     subJtPtBins[gI] = subJtPtBinsTemp[gI];
   }
 }
 else getLinBins(subJtPtBinsLow, subJtPtBinsHigh, nSubJtPtBins, subJtPtBins);

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 //Now that we have subjtpt bins, do the doJTPtCut check
 Double_t jtPtBinsLowReco = inPbPbFileConfig_p->GetValue("JTPTBINSLOWRECO", -10.0);
 if(doJtPtCut){
   bool jtPtCutFound = false;
   double deltaVal = 0.001;
   for(Int_t sI = 0; sI < nSubJtPtBins; ++sI){
     if(TMath::Abs(subJtPtBins[sI] - jtPtBinsLowReco) < deltaVal){
       jtPtCutFound = true;
       break;
     }
   }

   if(!jtPtCutFound){
     std::cout << "doJtPtCut is on but jtPtCut, " << jtPtBinsLowReco << " is not found in bins, ";
     for(Int_t sI = 0; sI < nSubJtPtBins; ++sI){
 	std::cout << subJtPtBins[sI] << ", ";
     }
     std::cout << subJtPtBins[nSubJtPtBins] << ". return 1" << std::endl;
     return 1;
   }
 }

 //We have subjt pt bins and gamma pt bins now setup the bin flattener
 const Float_t subJtGammaPtMin = inPbPbFileConfig_p->GetValue("SUBJTGAMMAPTBINSMIN", -1000.0);
 const Float_t subJtGammaPtMax = inPbPbFileConfig_p->GetValue("SUBJTGAMMAPTBINSMAX", -1000.0);

 binFlattener subJtGammaPtBinFlattener;
 if(!subJtGammaPtBinFlattener.Init("subJtGammaPtBinFlattener", nGammaPtBins, gammaPtBins, nSubJtPtBins, subJtPtBins)){
   std::cout << "Error gdjPlotResults: Failure to init binFlattener, L" << __LINE__ << ". return 1" << std::endl;
   return 1;
 }

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

 std::vector<double> subJtGammaPtBinsV = subJtGammaPtBinFlattener.GetFlattenedBins(subJtGammaPtMin, subJtGammaPtMax);
 Int_t nSubJtGammaPtBins = nGammaPtBins*nSubJtPtBins;

 gStyle->SetOptStat(0);


 const Int_t iterPPPhoPt = inPPTermFileConfig_p->GetValue("GAMMAPT_PP_Nominal", -1);

 const Int_t nCentIntBins = 5;
 Double_t centIntBins[nCentIntBins+1] = {0,10,30,80,100,120};

 for(Int_t gI = 1; gI < nGammaPtBins-1; ++gI){
   // Create a few summary plots
   TH1D* integralPerCentrality_h = new TH1D("integralPerCentrality_h", ";Centrality/p+p;R_{JJ#gamma}", nCentIntBins, centIntBins);
   TH1D* meanPerCentrality_h = new TH1D("meanPerCentrality_h", ";Centrality/p+p;#LT#vec{x}_{JJ#gamma}#GT", nCentIntBins, centIntBins);

   Int_t iterPP = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_PP_Nominal").c_str(), -1);

   if(iterPP < 0){
     std::cout << "Terminating iteration is \'" << iterPP << "\' in pp, for \'" << varNameUpper + "_GAMMAPT_PP_NOMINAL" << "\'. return 1" << std::endl;
     return 1;
   }

   TH1D* ppHistPhoPt_p = (TH1D*)inPPFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPPPhoPt) + "_PP_Nominal_PURCORR_COMBINED_h").c_str());
   Double_t scaledTotalPP = 0.0;
   for(Int_t bIX = 0; bIX < ppHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
     Float_t binCenter = ppHistPhoPt_p->GetBinCenter(bIX+1);
     if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

     scaledTotalPP += ppHistPhoPt_p->GetBinContent(bIX+1);
   }

   std::cout << "SCALED TOTAL PP: " << scaledTotalPP << std::endl;

   std::string ppHistName = "photonPtJet" + varNameHist + "Reco_Iter" + std::to_string(iterPP) + "_PP_Nominal_PURCORR_COMBINED_h";
   std::string pyt8TruthName = "photonPtJet" + varNameHist + "Truth_PreUnfold_PP_NOMINAL_TRUTH_COMBINED_h";

   std::vector<std::string> ppSystHistNames;
   std::vector<bool> isSystCorrelated;
   std::vector<std::string> systNames;
   std::vector<Double_t> scaledTotalPPSyst;
   for(unsigned int jI = 0; jI < systStrVect.size(); ++jI){
     if(isStrSame(systStrVect[jI], "Nominal")) continue;

     systNames.push_back(systStrVect[jI]);

     if(isStrSame(systStrVect[jI], "JESRTRK")) isSystCorrelated.push_back(false);
     else if(isStrSame(systStrVect[jI], "PRIOR")) isSystCorrelated.push_back(false);
     else if(isStrSame(systStrVect[jI], "MCSTAT")) isSystCorrelated.push_back(false);
     else isSystCorrelated.push_back(true);

     //     Int_t iterPPSyst = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_PP_" + systStrVect[jI]).c_str(), -1);
     Int_t iterPPSyst = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_PP_Nominal").c_str(), -1);
     ppSystHistNames.push_back("photonPtJet" + varNameHist + "Reco_Iter" + std::to_string(iterPPSyst) + "_PP_" + systStrVect[jI] + "_PURCORR_COMBINED_h");

     if(isStrSame(systTypeVect[jI], "GESGER") || isStrSame(systStrVect[jI], "ISO85") || isStrSame(systStrVect[jI], "ISO95")){
       // 	const Int_t iterPPPhoPtSyst = inPPTermFileConfig_p->GetValue(("GAMMAPT_PP_" + systStrVect[jI]).c_str(), -1);
 	const Int_t iterPPPhoPtSyst = inPPTermFileConfig_p->GetValue("GAMMAPT_PP_Nominal", -1);

 	TH1D* ppHistPhoPtSyst_p = (TH1D*)inPPFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPPPhoPtSyst) + "_PP_" + systStrVect[jI] + "_PURCORR_COMBINED_h").c_str());
 	Double_t scaledTotalPPTemp = 0.0;
 	for(Int_t bIX = 0; bIX < ppHistPhoPtSyst_p->GetXaxis()->GetNbins(); ++bIX){
 	  Float_t binCenter = ppHistPhoPtSyst_p->GetBinCenter(bIX+1);
 	  if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

 	  scaledTotalPPTemp += ppHistPhoPtSyst_p->GetBinContent(bIX+1);
 	}
 	scaledTotalPPSyst.push_back(scaledTotalPPTemp);
     }
     else scaledTotalPPSyst.push_back(scaledTotalPP);
   }

   TH1D* ratioJEWEL0to10Curve_p = nullptr;
   /*
   TH1D* ratioJEWEL0to10_p = nullptr;
   if(!doXTrunc) ratioJEWEL0to10_p = new TH1D("ratioJEWEL0to10_h", ";;", nVarBins, varBins);
   else ratioJEWEL0to10_p = new TH1D("ratioJEWEL0to10_h", ";;", nVarBinsTrunc, varBinsTrunc);
   */

   std::vector<TH1D*> ratioHistsPerCentrality;
   std::vector<std::vector<Double_t> > ratioSystsPerCentrality;
   std::vector<std::vector<std::string> > ratioLabelsPerCentrality;

   // We need to find all gamma bins matching our current bin for creating our final histogram

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   std::vector<int> gammaMatchedBins;
   if(isMultijet){
     for(Int_t sgI = 0; sgI < nSubJtGammaPtBins; ++sgI){
 	if(subJtGammaPtBinFlattener.GetBin1PosFromGlobal(sgI) == gI){
 	  if(doJtPtCut){
 	    int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(sgI);
 	    Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+1])/2.0;
 	    if(subJtPtBinCenter > jtPtBinsLowReco) gammaMatchedBins.push_back(sgI);
 	  }
 	  else gammaMatchedBins.push_back(sgI);
 	}
     }
   }
   else gammaMatchedBins.push_back(gI);

   std::cout << "PYT8NAME: " << pyt8TruthName << std::endl;
   TH1D* pyt8TruthPhoPt_p = (TH1D*)inPPFile_p->Get("photonPtTruth_PreUnfold_PP_TRUTH_COMBINED_h");
   TH2D* pyt8Truth2D_p = (TH2D*)inPPFile_p->Get(pyt8TruthName.c_str());
   TH2D* ppHist2D_p = (TH2D*)inPPFile_p->Get(ppHistName.c_str());

   if(gI == 2 && false){
     std::cout << "PP HIST PRINT ON GI == 2" << std::endl;
     ppHist2D_p->Print("ALL");
     std::cout << "END PRINT, gammaMatchedBins.size(): " << gammaMatchedBins.size() << std::endl;
     std::cout << "GammaMatchedBins: " << std::endl;
     for(unsigned int sgI = 0; sgI < gammaMatchedBins.size(); ++sgI){
 	std::cout << " " << sgI << ": " << gammaMatchedBins[sgI] << std::endl;
     }
   }

   TH1D* ppHist_p = nullptr;
   if(!doXTrunc) ppHist_p = new TH1D(("ppHist1D_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBins, varBins);
   else ppHist_p = new TH1D(("ppHist1D_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBinsTrunc, varBinsTrunc);
   fineTH2ToCoarseTH1(ppHist2D_p, ppHist_p, gammaMatchedBins);

   binWidthAndScaleNorm(ppHist_p, scaledTotalPP);
   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
   std::vector<std::vector<Double_t > > ppSyst = getSyst(inPPFile_p, ppHist_p, ppSystHistNames, gammaMatchedBins, scaledTotalPPSyst, inUnfoldNonClosureFileName, doGlobalDebug);
   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   // Now we gotta make a few subsets of systematics for clearer plotting
   std::map<std::string, std::vector<std::string> > systTypeToSystNames;
   std::map<std::string, std::vector<std::vector<Double_t> > > systTypeToSystValsPP;
   for(auto const &val : systTypesMap){
     std::vector<std::string> tempVectNames;
     std::vector<std::vector<Double_t> > tempVectVals;

     for(unsigned int sI = 0; sI < systNames.size(); ++sI){
 	if(!isStrSame(val.first, systStrToSystType[systNames[sI]])) continue;

 	tempVectNames.push_back(systNames[sI]);
 	tempVectVals.push_back(ppSyst[sI]);
     }

     systTypeToSystValsPP[val.first] = tempVectVals;
     systTypeToSystNames[val.first] = tempVectNames;
   }

 if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   HIJet::Style::EquipHistogram(ppHist_p, 0);
   ppHist_p->SetMarkerColor(ppColor);
   ppHist_p->SetMarkerStyle(72);
   ppHist_p->SetLineColor(ppColor);
   ppHist_p->SetLineWidth(2);
   ppHist_p->GetXaxis()->SetTitleFont(titleFont);
   ppHist_p->GetYaxis()->SetTitleFont(titleFont);
   ppHist_p->GetXaxis()->SetTitleSize(0.00001);
   ppHist_p->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));

   ppHist_p->GetXaxis()->SetLabelFont(titleFont);
   ppHist_p->GetYaxis()->SetLabelFont(titleFont);
   ppHist_p->GetXaxis()->SetLabelSize(0.00001);
   ppHist_p->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));

   ppHist_p->GetYaxis()->SetTitleOffset(yOffset);

   ppHist_p->GetYaxis()->SetNdivisions(505);

   if(isMultijet) ppHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{#gamma}} #frac{dN_{JJ#gamma}}{d" + varNameLabel + "}").c_str());
   else ppHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{#gamma}} #frac{dN_{J#gamma}}{d" + varNameLabel + "}}").c_str());
   ppHist_p->GetXaxis()->SetTitle(varNameLabel.c_str());

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   // JEWEL handling for pp
   TH1D* inJEWELPPHist_p = nullptr;
   TH1D* inJEWELPPHistCurve_p = nullptr;
   if(doJEWEL){
     inJEWELPPHist_p = (TH1D*)ppJEWELFile_p->Get((varNameLower + "_R" + std::to_string(jetR) + "_GammaPt" + std::to_string(gI) + "_h").c_str());
     inJEWELPPHistCurve_p = (TH1D*)ppJEWELFile_p->Get((varNameLower + "Curve_R" + std::to_string(jetR) + "_GammaPt" + std::to_string(gI) + "_h").c_str());
   }


   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   //JETSCAPE handling
   TH1D* inJETSCAPEPPHist_p = nullptr;
   if(doJETSCAPE){
     constructJETSCAPEHist(inJETSCAPEPPFile, &inJETSCAPEPPHist_p, kGreen+2);
     inJETSCAPEPPHist_p->SetLineStyle(2);
   }
   //END JETSCAPE HANDLING

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   //LBT handling
   TH1D* inLBTPPHist_p = nullptr;
   if(doLBT){
     constructJETSCAPEHist(inLBTPPFile, &inLBTPPHist_p, kAzure-3);//kMagenta);
     inLBTPPHist_p->SetLineStyle(2);
     inLBTPPHist_p->Print("ALL");
   }
   //END LBT HANDLING

   if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
     TH1D* pyt8Truth_p = nullptr;

     if(!doXTrunc) pyt8Truth_p = new TH1D("pyt8Truth1D_h", ";;", nVarBins, varBins);
     else pyt8Truth_p = new TH1D("pyt8Truth1D_h", ";;", nVarBinsTrunc, varBinsTrunc);
     fineTH2ToCoarseTH1(pyt8Truth2D_p, pyt8Truth_p, gammaMatchedBins);

     binWidthAndScaleNorm(pyt8Truth_p, pyt8TruthPhoPt_p->GetBinContent(gI+1));

     TH1D* jewelPPHist_p = nullptr;
     TH1D* jewelPPHistCurve_p = nullptr;
     if(doJEWEL){
       if(!doXTrunc){
	 jewelPPHist_p = new TH1D("jewelPPHist_h", ";;", nVarBins, varBins);
	 jewelPPHistCurve_p = new TH1D("jewelPPHistCurve_h", ";;", inJEWELPPHistCurve_p->GetNbinsX(), varBins[0], varBins[nVarBins]);
	 if(isStrSame(centBinsStr[cI], "Cent0to10")) ratioJEWEL0to10Curve_p = new TH1D("ratioJEWEL0to10Curve_h", ";;", inJEWELPPHistCurve_p->GetNbinsX(), varBins[0], varBins[nVarBins]);
       }
       else{
	 jewelPPHist_p = new TH1D("jewelPPHist_h", ";;", nVarBinsTrunc, varBinsTrunc);

	 Int_t nBinsCurve = 0;
	 for(Int_t bIX = 0; bIX < inJEWELPPHistCurve_p->GetNbinsX(); ++bIX){
	   Double_t binC = inJEWELPPHistCurve_p->GetXaxis()->GetBinCenter(bIX+1);
	   if(binC < varBinsTrunc[0]) continue;
	   if(binC > varBinsTrunc[nVarBinsTrunc]) continue;

	   ++nBinsCurve;
	 }

	 jewelPPHistCurve_p = new TH1D("jewelPPHistCurve_h", ";;", nBinsCurve, varBinsTrunc[0], varBinsTrunc[nVarBinsTrunc]);
	 if(isStrSame(centBinsStr[cI], "Cent0to10")) ratioJEWEL0to10Curve_p = new TH1D("ratioJEWEL0to10Curve_h", ";;", nBinsCurve, varBinsTrunc[0], varBinsTrunc[nVarBinsTrunc]);
       }

 	fineHistToCoarseHist(inJEWELPPHist_p, jewelPPHist_p);
 	fineHistToCoarseHist(inJEWELPPHistCurve_p, jewelPPHistCurve_p);

 	//For the tlegend
 	inJEWELPPHist_p->SetLineColor(1);
 	inJEWELPPHist_p->SetLineStyle(1);

  	inJEWELPPHistCurve_p->SetLineColor(1);
 	inJEWELPPHistCurve_p->SetLineStyle(1);
     }
     //END DO JEWEL

     TH1D* inJETSCAPEPbPbHist_p = nullptr;
     bool doJETSCAPEPbPb = doJETSCAPE && centBinsStr[cI].find("Cent0to10") != std::string::npos;
     if(doJETSCAPEPbPb) constructJETSCAPEHist(inJETSCAPEPbPbFile, &inJETSCAPEPbPbHist_p, kGreen+4);

     TH1D* inLBTPbPbHist_p = nullptr;
     bool doLBTPbPb = doLBT && centBinsStr[cI].find("Cent0to10") != std::string::npos;
     if(doLBTPbPb) constructJETSCAPEHist(inLBTPbPbFile, &inLBTPbPbHist_p, kAzure-3);//kMagenta);

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     const Int_t iterPbPbPhoPt = inPbPbTermFileConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI] + "_Nominal").c_str(), -1);
     TH1D* pbpbHistPhoPt_p = nullptr;
     if(iterPbPbPhoPt > 0) pbpbHistPhoPt_p = (TH1D*)inPbPbFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPbPbPhoPt) + "_" + centBinsStr[cI] + "_Nominal_PURCORR_COMBINED_h").c_str());
     else{
 	std::cout << "No good termination point found. return 1" << std::endl;
 	return 1;
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     Double_t scaledTotalPbPb = 0.0;
     for(Int_t bIX = 0; bIX < pbpbHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
 	Float_t binCenter = pbpbHistPhoPt_p->GetBinCenter(bIX+1);
 	if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

 	scaledTotalPbPb += pbpbHistPhoPt_p->GetBinContent(bIX+1);
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     Int_t iterPbPb = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_" + centBinsStr[cI] + "_Nominal").c_str(), -1);
     if(iterPbPb < 0){
 	std::cout << "Terminating iteration is \'" << iterPbPb << "\' in pbpb, \'" << centBinsStr[cI] << "\' for \'" << varNameUpper + "_GAMMAPT_" + centBinsStr[cI] << "\'. return 1" << std::endl;
 	return 1;
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     std::string centBins2Str = centBinsStr[cI];
     if(centBins2Str.find("to") != std::string::npos){
 	centBins2Str.replace(0, 4, "");
 	centBins2Str.replace(centBins2Str.find("to"), 2, "-");
 	centBins2Str = centBins2Str + "%";
     }

     std::string pbpbHistName = "photonPtJet" + varNameHist + "Reco_Iter" + std::to_string(iterPbPb) + "_" + centBinsStr[cI] + "_Nominal_PURCORR_COMBINED_h";
     std::vector<std::string> pbpbSystHistNames;


     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     std::vector<Double_t> scaledTotalPbPbSyst;
     for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
 	if(isStrSame(systStrVect[systI], "Nominal")) continue;

	//6/06 - you want to keep the termination point fixed for systematics per arc feedback
	// 	Int_t iterPbPbSyst = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_" + centBinsStr[cI] + "_" + systStrVect[systI]).c_str(), -1);
 	Int_t iterPbPbSyst = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPT_" + centBinsStr[cI] + "_Nominal").c_str(), -1);
 	pbpbSystHistNames.push_back("photonPtJet" + varNameHist + "Reco_Iter" + std::to_string(iterPbPbSyst) + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");

 	if(isStrSame(systTypeVect[systI], "GESGER") /* || isStrSame(systStrVect[systI], "ISO85") || isStrSame(systStrVect[systI], "ISO95")*/){
	  // 	  const Int_t iterPbPbPhoPtSyst = inPbPbTermFileConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI] + "_" + systStrVect[systI]).c_str(), -1);
 	  const Int_t iterPbPbPhoPtSyst = inPbPbTermFileConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI] + "_Nominal").c_str(), -1);

 	  TH1D* pbpbHistPhoPtSyst_p = (TH1D*)inPbPbFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPbPbPhoPtSyst) + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h").c_str());
 	  Double_t scaledTotalPbPbTemp = 0.0;

 	  for(Int_t bIX = 0; bIX < pbpbHistPhoPtSyst_p->GetXaxis()->GetNbins(); ++bIX){
 	    Float_t binCenter = pbpbHistPhoPtSyst_p->GetBinCenter(bIX+1);
 	    if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

 	    scaledTotalPbPbTemp += pbpbHistPhoPtSyst_p->GetBinContent(bIX+1);
 	  }
 	  scaledTotalPbPbSyst.push_back(scaledTotalPbPbTemp);
 	}
 	else scaledTotalPbPbSyst.push_back(scaledTotalPbPb);
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     TH2D* pbpbHist2D_p = (TH2D*)inPbPbFile_p->Get(pbpbHistName.c_str());
     TH1D* pbpbHist_p = nullptr;
     if(!doXTrunc) pbpbHist_p = new TH1D(("pbpbHist1D_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBins, varBins);
     else pbpbHist_p = new TH1D(("pbpbHist1D_GammaPt" + std::to_string(gI) + "_h").c_str(), ";;", nVarBinsTrunc, varBinsTrunc);
     fineTH2ToCoarseTH1(pbpbHist2D_p, pbpbHist_p, gammaMatchedBins);

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     binWidthAndScaleNorm(pbpbHist_p, scaledTotalPbPb);
     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
     std::vector<std::vector<Double_t > > pbpbSyst = getSyst(inPbPbFile_p, pbpbHist_p, pbpbSystHistNames, gammaMatchedBins, scaledTotalPbPbSyst, inUnfoldNonClosureFileName, doGlobalDebug);
     //pbpbHist_p->Scale(1./pbpbHist_p->Integral());

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
     std::map<std::string, std::vector<std::vector<Double_t> > > systTypeToSystValsPbPb;
     for(auto const &val : systTypesMap){
 	std::vector<std::vector<Double_t> > tempVectVals;

	if(doGlobalDebug) std::cout << val.first << std::endl;
 	for(unsigned int sI = 0; sI < systNames.size(); ++sI){
	  if(doGlobalDebug) std::cout << " " << sI << "/" << systNames.size() << ": " << systNames[sI] << std::endl;

 	  if(!isStrSame(val.first, systStrToSystType[systNames[sI]])) continue;
 	  tempVectVals.push_back(pbpbSyst[sI]);
 	}

 	systTypeToSystValsPbPb[val.first] = tempVectVals;
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     //   JEWEL handling for pbpb
     TH1D* inJEWELPbPbHist_p = nullptr;
     TH1D* inJEWELPbPbHistCurve_p = nullptr;
     TH1D* jewelPbPbHist_p = nullptr;
     TH1D* jewelPbPbHistCurve_p = nullptr;
     TLegend* jewelLeg_p = nullptr;
     if(doJEWEL){
       if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

       Int_t jewelCentPos = -1;
 	for(unsigned int jewelI = 0; jewelI < inJEWELPbPbFileNames.size(); ++jewelI){
 	  std::string jewelCent = pbpbJEWELConfig_p[jewelI]->GetValue("CENT", "");

	  std::cout << "JEWELCENT: " << jewelCent << std::endl;

 	  while(jewelCent.find("/") != std::string::npos){
 	    jewelCent.replace(jewelCent.rfind("/"), jewelCent.size(), "");
 	  }

 	  if(isStrSame(centBinsStr[cI], jewelCent)){
 	    jewelCentPos = jewelI;
 	    break;
 	  }
 	}

 	//Grab legend parameters for jewel
 	Double_t jewelLegX = config_p->GetValue("JEWELLEGX", -1.0);
 	Double_t jewelLegY = config_p->GetValue("JEWELLEGY", -1.0);

 	if(jewelCentPos == -1){
 	  std::cout << "JEWEL CENT STR MATCHING \'" << centBinsStr[cI] << "\' not found. Pb+Pb in this bin will not be plotted, please fix for overlay" << std::endl;
	  doJEWEL = false;
 	}
 	else{
 	  inJEWELPbPbHist_p = (TH1D*)pbpbJEWELFile_p[jewelCentPos]->Get((varNameLower + "_R" + std::to_string(jetR) + "_GammaPt" + std::to_string(gI) + "_h").c_str());
 	  inJEWELPbPbHistCurve_p = (TH1D*)pbpbJEWELFile_p[jewelCentPos]->Get((varNameLower + "Curve_R" + std::to_string(jetR) + "_GammaPt" + std::to_string(gI) + "_h").c_str());

 	  if(!doXTrunc){
	    jewelPbPbHist_p = new TH1D("jewelPbPbHist_h", ";;", nVarBins, varBins);
	    jewelPbPbHistCurve_p = new TH1D("jewelPbPbHistCurve_h", ";;", nVarBins, varBins);
	  }
 	  else{
	    jewelPbPbHist_p = new TH1D("jewelPbPbHist_h", ";;", nVarBinsTrunc, varBinsTrunc);
	    Int_t tempNBins = jewelPPHistCurve_p->GetXaxis()->GetNbins();
	    jewelPbPbHistCurve_p = new TH1D("jewelPbPbHistCurve_h", ";;", tempNBins, varBinsTrunc[0], varBinsTrunc[nVarBinsTrunc]);
	  }

 	  fineHistToCoarseHist(inJEWELPbPbHist_p, jewelPbPbHist_p);
 	  fineHistToCoarseHist(inJEWELPbPbHistCurve_p, jewelPbPbHistCurve_p);
 	}

 	jewelLeg_p = new TLegend(jewelLegX, jewelLegY - 0.065*2, jewelLegX+0.25, jewelLegY);
 	jewelLeg_p->SetTextFont(42);
 	jewelLeg_p->SetTextSize(0.035/(1.0 - padSplit));
 	jewelLeg_p->SetBorderSize(0);
 	jewelLeg_p->SetFillStyle(0);
     }

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     TCanvas* canvPP_p = nullptr;
     TPad* padsPP_p[nPad];
     if(cI == 0){
 	canvPP_p = new TCanvas("canvPP_p", "", width, height);
 	padsPP_p[0] = new TPad("padPP0", "", 0.0, padSplit, 1.0, 1.0);
 	setMargins(padsPP_p[0], leftMargin, topMargin, rightMargin, noMargin);
 	canvPP_p->cd();
 	padsPP_p[0]->Draw("SAME");

 	canvPP_p->cd();
 	padsPP_p[1] = new TPad("padPP1", "", 0.0, 0.0, 1.0, padSplit);
 	setMargins(padsPP_p[1], leftMargin, noMargin, rightMargin, bottomMargin);
 	canvPP_p->cd();
 	padsPP_p[1]->Draw("SAME");
     }

     TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
     TPad* pads_p[nPad];
     setMargins(canv_p, noMargin, noMargin, noMargin, noMargin);
     canv_p->cd();

     pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
     setMargins(pads_p[0], leftMargin, topMargin, rightMargin, noMargin);
     canv_p->cd();
     pads_p[0]->Draw("SAME");

     canv_p->cd();
     pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
     setMargins(pads_p[1], leftMargin, noMargin, rightMargin, bottomMargin);
     canv_p->cd();
     pads_p[1]->Draw("SAME");

     pads_p[0]->cd();
     HIJet::Style::EquipHistogram(pbpbHist_p, 1);
     pbpbHist_p->SetMarkerColor(1);
     pbpbHist_p->SetLineColor(1);
     pbpbHist_p->SetMarkerStyle(20);

     ppHist_p->SetMinimum(minVal);
     ppHist_p->SetMaximum(maxVal);

     //   Need to set this for both so the systematics dont overflow the plotting space but still draw correctly
     pbpbHist_p->SetMinimum(minVal);
     pbpbHist_p->SetMaximum(maxVal);

     if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

     ppHist_p->GetXaxis()->SetTitleSize(ppHist_p->GetXaxis()->GetTitleSize()*1.3);

     ppHist_p->DrawCopy("HIST E X0 P");
     pbpbHist_p->DrawCopy("HIST E X0 P SAME");

     if(doJEWEL){
 	HIJet::Style::EquipHistogram(jewelPPHistCurve_p, 1);
 	jewelPPHistCurve_p->SetLineStyle(2);
 	jewelPPHistCurve_p->SetLineWidth(8);
 	jewelPPHistCurve_p->SetMarkerSize(0.0001);
 	jewelPPHist_p->SetMarkerSize(0.0001);
 	jewelPPHist_p->SetMarkerColorAlpha(0, 0.0);
 	jewelPPHist_p->SetMarkerStyle(1);

 	if(jewelPbPbHistCurve_p != nullptr){
 	  HIJet::Style::EquipHistogram(jewelPbPbHistCurve_p, 1);
	  // 	  jewelPbPbHistCurve_p->SetLineStyle(9);
 	  jewelPbPbHistCurve_p->SetLineWidth(8);
 	  jewelPbPbHistCurve_p->SetMarkerSize(0.0001);
 	  jewelPbPbHistCurve_p->SetMarkerColorAlpha(0, 0.0);
 	  jewelPbPbHistCurve_p->SetMarkerStyle(1);
 	}

	//Getting rid of scale factor!!!
	// 	jewelPPHistCurve_p->Scale(scaleFactorPYT);
	// 	jewelPPHist_p->Scale(scaleFactorPYT);

 	pads_p[0]->cd();
 	jewelPPHistCurve_p->DrawCopy("E X0 HIST SAME");

 	//Change the line color to black for the tlegend and add both jewels
	// 	jewelPPHistCurve_p->SetLineColor(1);
	//	jewelLeg_p->AddEntry(inJEWELPPHist_p, "Data", "L");
	//	jewelLeg_p->AddEntry(jewelPPHistCurve_p, ("JEWEL (x " + prettyString(scaleFactorPYT,1,false) + ")").c_str(), "L F");
	jewelLeg_p->AddEntry(jewelPPHistCurve_p, "JEWEL", "L F");

	//	jewelLeg_p->Draw("SAME");

	if(jewelPbPbHistCurve_p != nullptr){
	  //	  jewelPbPbHistCurve_p->Scale(scaleFactorPYT);
	  jewelPbPbHistCurve_p->DrawCopy("E X0 HIST SAME");
	}
     }

     if(doJETSCAPE){
       inJETSCAPEPPHist_p->DrawCopy("SAME C E3");
       //       inJETSCAPEPPHist_p->DrawCopy("SAME E P HIST");
       if(doJETSCAPEPbPb){
	 inJETSCAPEPbPbHist_p->DrawCopy("SAME E3 C");
	 //	 inJETSCAPEPbPbHist_p->DrawCopy("SAME E P HIST ");
       }
     }

     if(doLBT){
       inLBTPPHist_p->DrawCopy("HIST C SAME");
       if(doLBTPbPb) inLBTPbPbHist_p->DrawCopy("HIST C SAME");
     }


      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      pads_p[0]->cd();


      //      HIJet::Style::EquipHistogram(ppHist_p, 2);
      drawSyst(pads_p[0], ppHist_p, ppSyst);
      //      HIJet::Style::EquipHistogram(pbpbHist_p, 1);
      pbpbHist_p->SetMarkerColor(pbpbColor);
      drawSyst(pads_p[0], pbpbHist_p, pbpbSyst);
      pbpbHist_p->SetMarkerColor(1);
      //      pbpbHist_p->SetMarkerColor(1);
      //      pbpbHist_p->SetLineColor(1);
      //      pbpbHist_p->SetMarkerStyle(72);

      ppHist_p->DrawCopy("HIST E X0 P SAME");
      pbpbHist_p->DrawCopy("HIST E X0 P SAME");

      //Do the same but just pp
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      if(cI == 0){
	padsPP_p[0]->cd();

	ppHist_p->DrawCopy("HIST E X0 P");
	drawSyst(padsPP_p[0], ppHist_p, ppSyst);
	ppHist_p->DrawCopy("HIST E X0 P SAME");

	HIJet::Style::EquipHistogram(pyt8Truth_p, 5);
	pyt8Truth_p->SetLineWidth(8);
	//	pyt8Truth_p->Scale(scaleFactorPYT);
	pyt8Truth_p->SetMarkerSize(0.00001);
	pyt8Truth_p->SetMarkerStyle(1);
	pyt8Truth_p->DrawCopy("E X0 HIST SAME");
	//	pyt8Truth_p->Scale(1.0/scaleFactorPYT);


	if(doJEWEL){
	  HIJet::Style::EquipHistogram(jewelPPHist_p, 2);
	  jewelPPHist_p->SetMarkerSize(0.00001);
	  jewelPPHist_p->SetMarkerColorAlpha(0, 0.0);
	  jewelPPHist_p->SetMarkerStyle(1);
	  jewelPPHist_p->SetLineStyle(2);
	  jewelPPHist_p->SetLineWidth(8);

	  jewelPPHist_p->DrawCopy("E X0 HIST SAME");
	  //	  jewelPPHist_p->Scale(1.0/scaleFactorPP);

	  jewelPPHist_p->SetLineColor(1);
	}

	ppHist_p->DrawCopy("HIST E X0 P SAME");

	gPad->SetTicks();
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      Float_t integral = 0.0;
      Float_t integralError = 0.0;
      for(Int_t bIX = 0; bIX < pbpbHist_p->GetXaxis()->GetNbins(); ++bIX){
	integral += pbpbHist_p->GetBinContent(bIX+1)*pbpbHist_p->GetXaxis()->GetBinWidth(bIX+1);
	Float_t err = pbpbHist_p->GetBinError(bIX+1)*pbpbHist_p->GetXaxis()->GetBinWidth(bIX+1);
	integralError += err*err;
      }
      integralPerCentrality_h->SetBinContent(cI+1, integral);
      integralPerCentrality_h->SetBinError(cI+1, TMath::Sqrt(integralError));

      meanPerCentrality_h->SetBinContent(cI+1, pbpbHist_p->GetMean());
      meanPerCentrality_h->SetBinError(cI+1, pbpbHist_p->GetMeanError());

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      std::vector<std::string> labelsAlt;
      std::vector<std::string> labelsTemp = getLabels(inPbPbFileConfig_p, pbpbHist_p, &labelMap, &labelsAlt);
      std::vector<std::string> labels = {"#bf{#it{ATLAS Internal}}", "2018 Pb+Pb 1.72 nb^{-1}", "2017 #it{pp} 260 pb^{-1}", "#sqrt{s} = 5.02 TeV"};

      for(unsigned int lI = 0; lI < labelsTemp.size(); ++lI){
	if(isStrSame(labelsTemp[lI], "h")) continue;

	if(doJtPtCut){
	  if(labelsTemp[lI].find("p_{T,Jet}^{Truth}") != std::string::npos){
	    //	  std::cout << "LABEL FOUND: " << labelsTemp[lI] << std::endl;
	    labelsTemp[lI].replace(0,labelsTemp[lI].find("<"), "");
	    labelsTemp[lI] = std::to_string(((int)jtPtBinsLowReco)) + labelsTemp[lI];
	  }
	}

	if(labelsTemp[lI].find("p_{T,Jet}^{Reco.}") != std::string::npos) continue;
	else if(labelsTemp[lI].find("^{Truth}") != std::string::npos) labelsTemp[lI].replace(labelsTemp[lI].find("^{Truth}"), 8, "");
	else if(labelsTemp[lI].find("#eta_{Jet}") != std::string::npos){
	  std::string frontStr = labelsTemp[lI].substr(0, labelsTemp[lI].find("<"));
	  std::string backStr = labelsTemp[lI].substr(labelsTemp[lI].rfind("<")+1, labelsTemp[lI].size());

	  if(isStrSame(frontStr, "-" + backStr)){
	    labelsTemp[lI] = "|#eta_{Jet}|<" + backStr;
	  }
	}


	if(labelsTemp[lI].find("p_{T,#gamma}") != std::string::npos){
	  while(labelsTemp[lI].find(".0") != std::string::npos){
	    labelsTemp[lI].replace(labelsTemp[lI].find(".0"), 2, "");
	  }
	}

	if(labelsTemp[lI].find("#eta_{Jet}") != std::string::npos) continue;
	if(labelsTemp[lI].find("#DeltaR_{JJ}") != std::string::npos) continue;

	if(labelsTemp[lI].find("p_{T") != std::string::npos) labelsTemp[lI] = labelsTemp[lI] + " GeV";


	if(labelsTemp[lI].find("#it{R}") != std::string::npos){
	  labelsTemp[lI].replace(labelsTemp[lI].find("#it{R}"), 6, "R");
	}

	labels.push_back(labelsTemp[lI]);
      }


      //Put the R=0.X label at the end
      std::string rLabel = "";
      for(unsigned int lI = 0; lI < labels.size(); ++lI){
	if(labels[lI].find("R=") != std::string::npos){
	  rLabel = labels[lI];
	  labels.erase(labels.begin()+lI);
	  break;
	}
      }
      labels.push_back(rLabel);


      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      unsigned int pos = 0;
      while(pos < labels.size()){
	if(labels[pos].find("PURCORR") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("COMBINED") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("Iter") != std::string::npos) labels.erase(pos + labels.begin());
	else ++pos;
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      //Create some labels and names for systematics
      std::string plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_PP_" + saveTag + "_" + dateStr + "." + saveExt;
      std::vector<std::string> systLabels = labels;
      //      systLabels[1] = "#it{pp} at #sqrt{s_{NN}}=5.02 TeV";

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      //Create a bunch of systematics plots
      if(cI == 0){
	systLabels.erase(systLabels.begin()+1);
	plotSyst(ppHist_p, ppSyst, systNames, systLabels, "pp", plotSaveName);

	std::vector<std::vector<Double_t > > systTypeQuadSumVals;
	for(unsigned int sysI = 0; sysI < uniqueSystTypes.size(); ++sysI){
	  plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_PP_" + uniqueSystTypes[sysI] + "_" + saveTag + "_" + dateStr + "." + saveExt;

	  std::vector<std::vector<Double_t > > tempSystValVect = systTypeToSystValsPP[uniqueSystTypes[sysI]];
	  std::vector<std::string> tempSystNameVect = systTypeToSystNames[uniqueSystTypes[sysI]];

	  plotSyst(ppHist_p, tempSystValVect, tempSystNameVect, systLabels, "pp", plotSaveName, systNameToLegName(uniqueSystTypes[sysI]) + " Total");

	  Int_t nBinsTemp = nVarBins;
	  if(doXTrunc) nBinsTemp = nVarBinsTrunc;

	  std::vector<Double_t> totalSysts;
	  for(Int_t bI = 0; bI < nBinsTemp; ++bI){

	    Float_t total = 0.0;
	    for(unsigned int sI = 0; sI < tempSystValVect.size(); ++sI){
	      total += tempSystValVect[sI][bI]*tempSystValVect[sI][bI];
	    }

	    totalSysts.push_back(TMath::Sqrt(total));
	  }

	  systTypeQuadSumVals.push_back(totalSysts);
	}

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

	plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_PP_SystType_" + saveTag + "_" + dateStr + "." + saveExt;
	plotSyst(ppHist_p, systTypeQuadSumVals, uniqueSystTypes, systLabels, "pp", plotSaveName, "Total", systPPMin, systPPMax);


	systLabels = labels;
      }

      pbpbHist_p->GetXaxis()->SetTitle(ppHist_p->GetXaxis()->GetTitle());
      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + "." + saveExt;
      //      systLabels[1] = centBinsLabel[cI] + " at #sqrt{s_{NN}}=5.02 TeV";
      plotSyst(pbpbHist_p, pbpbSyst, systNames, systLabels, centBinsLabel[cI], plotSaveName);

      std::vector<std::vector<Double_t > > systTypeQuadSumVals;
      for(unsigned int sysI = 0; sysI < uniqueSystTypes.size(); ++sysI){
	plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + uniqueSystTypes[sysI] + "_" + saveTag + "_" + dateStr + "." + saveExt;

	std::vector<std::vector<Double_t > > tempSystValVect = systTypeToSystValsPbPb[uniqueSystTypes[sysI]];
	std::vector<std::string> tempSystNameVect = systTypeToSystNames[uniqueSystTypes[sysI]];

	plotSyst(pbpbHist_p, tempSystValVect, tempSystNameVect, systLabels, centBinsLabel[cI], plotSaveName, systNameToLegName(uniqueSystTypes[sysI]) + " Total");

	Int_t nBinsTemp = nVarBins;
	if(doXTrunc) nBinsTemp = nVarBinsTrunc;

	std::vector<Double_t> totalSysts;
	for(Int_t bI = 0; bI < nBinsTemp; ++bI){

	  Float_t total = 0.0;
	  for(unsigned int sI = 0; sI < tempSystValVect.size(); ++sI){
	    total += tempSystValVect[sI][bI]*tempSystValVect[sI][bI];
	  }

	  totalSysts.push_back(TMath::Sqrt(total));
	}

	systTypeQuadSumVals.push_back(totalSysts);
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_SystType_" + saveTag + "_" + dateStr + "." + saveExt;
      plotSyst(pbpbHist_p, systTypeQuadSumVals, uniqueSystTypes, systLabels, centBinsLabel[cI], plotSaveName, "Total", systMin, systMax);


      canv_p->cd();
      pads_p[0]->cd();

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(titleSize/(1.0 - padSplit));
      label_p->SetTextAlign(31);

      std::vector<std::string> tempLabels;
      for(unsigned int i = 0; i < labels.size(); ++i){
	Float_t yLow = 0.063;
	if(i > 0){
	  if(labels[i-1].find("_{T,Jet}") != std::string::npos) yLow += 0.003;
	}

	label_p->DrawLatex(labelX, labelY - yLow*i, labels[i].c_str());
	tempLabels.push_back(labels[i]);
      }
      ratioLabelsPerCentrality.push_back(tempLabels);
      delete label_p;

      Int_t nLegY = 2;
      if(doJEWEL){
	++nLegY;
	if(jewelPbPbHistCurve_p != nullptr) ++nLegY;
      }
      if(doJETSCAPE){
	++nLegY;
	if(doJETSCAPEPbPb) ++nLegY;
      }
      if(doLBT){
	++nLegY;
	if(doLBTPbPb) ++nLegY;
      }

      TLegend* leg_p = new TLegend(legX, legY - 0.065*nLegY, legX+0.25, legY);
      leg_p->SetTextFont(42);
      leg_p->SetTextSize(titleSize/(1.0 - padSplit));
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);

      leg_p->AddEntry(pbpbHist_p, ("Pb+Pb #bf{" + centBins2Str + "}").c_str(), "P L F");
      leg_p->AddEntry(ppHist_p, "#it{pp}", "P L F");
      if(doJEWEL){
	//	leg_p->AddEntry(jewelPPHistCurve_p, ("JEWEL #it{pp} (x " + prettyString(scaleFactorPYT, 1, false) + ")").c_str(), "L");
	//	if(jewelPbPbHistCurve_p != nullptr) leg_p->AddEntry(jewelPbPbHistCurve_p, ("JEWEL Pb+Pb (x " + prettyString(scaleFactorPYT, 1, false) + ")").c_str(), "L");

	leg_p->AddEntry(jewelPPHistCurve_p, "JEWEL #it{pp}", "L");
	if(jewelPbPbHistCurve_p != nullptr) leg_p->AddEntry(jewelPbPbHistCurve_p, "JEWEL Pb+Pb", "L");
      }

      if(doJETSCAPE){
	inJETSCAPEPPHist_p->SetLineWidth(0.0);

	leg_p->AddEntry(inJETSCAPEPPHist_p, "JETSCAPE #it{pp}", "F");
	if(doJETSCAPEPbPb){
	  inJETSCAPEPbPbHist_p->SetLineWidth(0.0);
	  leg_p->AddEntry(inJETSCAPEPbPbHist_p, "JETSCAPE Pb+Pb", "F");
	}
      }

      if(doLBT){
	leg_p->AddEntry(inLBTPPHist_p, "LBT #it{pp}", "L");
	if(doLBTPbPb){
	  leg_p->AddEntry(inLBTPbPbHist_p, "LBT Pb+Pb", "L");
	}
      }

      leg_p->Draw("SAME");

      gPad->SetTicks();
      if(doLogX) pads_p[0]->SetLogx();
      if(doLogY) pads_p[0]->SetLogy();

      canv_p->cd();
      pads_p[1]->cd();

      if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;
      std::vector<std::vector<Double_t > > ratSyst = getSystRat(inPbPbFile_p, inPPFile_p, pbpbHist_p, ppHist_p, pbpbSystHistNames, ppSystHistNames, isSystCorrelated, gammaMatchedBins, scaledTotalPbPbSyst, scaledTotalPPSyst, inUnfoldNonClosureFileName);
  if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

      std::map<std::string, std::vector<std::vector<Double_t> > > systTypeToSystValsRatio;
      for(auto const &val : systTypesMap){
	std::vector<std::vector<Double_t> > tempVectVals;

	for(unsigned int sI = 0; sI < systNames.size(); ++sI){
	  if(!isStrSame(val.first, systStrToSystType[systNames[sI]])) continue;

	  tempVectVals.push_back(ratSyst[sI]);
	}

	systTypeToSystValsRatio[val.first] = tempVectVals;
      }

      if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

      pbpbHist_p->Divide(ppHist_p);

      if(doGlobalDebug) std::cout << "LINE, FILE: " << __LINE__ << ", " << __FILE__ << std::endl;

      //We need to save the ratio histograms for plotting all centrality together
      std::string ratioHistName = "ratioHist_" + centBinsStr[cI];
      std::string ratioHistTitle = ";" + varNameLabel + ";Pb+Pb/#it{pp}";
      if(!doXTrunc) ratioHistsPerCentrality.push_back(new TH1D(ratioHistName.c_str(), ratioHistTitle.c_str(), nVarBins, varBins));
      else ratioHistsPerCentrality.push_back(new TH1D(ratioHistName.c_str(), ratioHistTitle.c_str(), nVarBinsTrunc, varBinsTrunc));

      for(Int_t bIX = 0; bIX < pbpbHist_p->GetXaxis()->GetNbins(); ++bIX){
	ratioHistsPerCentrality[cI]->SetBinContent(bIX+1, pbpbHist_p->GetBinContent(bIX+1));
	ratioHistsPerCentrality[cI]->SetBinError(bIX+1, pbpbHist_p->GetBinError(bIX+1));
      }

      pbpbHist_p->GetYaxis()->SetNdivisions(505);
      pbpbHist_p->GetXaxis()->SetTitleOffset(1.2);

      pbpbHist_p->GetXaxis()->SetTitleFont(titleFont);
      pbpbHist_p->GetYaxis()->SetTitleFont(titleFont);
      pbpbHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit));
      pbpbHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));

      pbpbHist_p->GetXaxis()->SetLabelFont(titleFont);
      pbpbHist_p->GetYaxis()->SetLabelFont(titleFont);
      pbpbHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit));
      pbpbHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));

      pbpbHist_p->GetYaxis()->SetTitle("Pb+Pb/#it{pp}");
      pbpbHist_p->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0 - padSplit));

      pbpbHist_p->SetMinimum(minRatVal);
      pbpbHist_p->SetMaximum(maxRatVal);

      gPad->SetTicks();

      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength()*(1.0 - padSplit)/padSplit);
      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength("Y")*(1.0 - padSplit)/padSplit, "Y");

      pbpbHist_p->GetXaxis()->SetTitle(varNameLabel.c_str());
      pbpbHist_p->GetXaxis()->SetTitleSize(pbpbHist_p->GetXaxis()->GetTitleSize()*1.3);
      pbpbHist_p->DrawCopy("HIST E X0 P");

      if(doJEWEL && jewelPbPbHistCurve_p != nullptr){
	jewelPbPbHistCurve_p->Divide(jewelPPHistCurve_p);

	if(isStrSame("Cent0to10", centBinsStr[cI])) fineHistToCoarseHist(jewelPbPbHistCurve_p, ratioJEWEL0to10Curve_p);

	jewelPbPbHistCurve_p->DrawCopy("E X0 HIST SAME");
      }

      if(doJETSCAPEPbPb){
	inJETSCAPEPbPbHist_p->Divide(inJETSCAPEPPHist_p);
	inJETSCAPEPbPbHist_p->DrawCopy("E3 C SAME");
	//	inJETSCAPEPbPbHist_p->DrawCopy("E HIST P SAME");
      }

      if(doLBTPbPb){
	inLBTPbPbHist_p->Divide(inLBTPPHist_p);
	inLBTPbPbHist_p->DrawCopy("HIST C SAME");
      }

      //drawSyst here returns the quadrature sum of the syst error for convenience just put it right in the push_back
      //     HIJet::Style::EquipHistogram(pbpbHist_p, 1);
      pbpbHist_p->SetMarkerColor(pbpbColor);
      ratioSystsPerCentrality.push_back(drawSyst(pads_p[1], pbpbHist_p, ratSyst));
      pbpbHist_p->SetMarkerColor(1);
      //     pbpbHist_p->SetMarkerColor(1);
      //     pbpbHist_p->SetLineColor(1);
      //     pbpbHist_p->SetMarkerStyle(72);
      pbpbHist_p->DrawCopy("HIST E X0 P SAME");

      if(cI == 0){
	padsPP_p[1]->cd();
	if(doJEWEL){
	  //	  jewelPPHist_p->Scale(scaleFactorPP);
	  //	  pyt8Truth_p->Scale(scaleFactorPYT);

	  jewelPPHist_p->Divide(ppHist_p);
	  pyt8Truth_p->Divide(ppHist_p);

	  HIJet::Style::EquipHistogram(jewelPPHist_p, 2);
	  jewelPPHist_p->SetLineStyle(2);
	  jewelPPHist_p->SetLineWidth(8);
	  jewelPPHist_p->SetMarkerSize(0.0001);

	  jewelPPHist_p->GetYaxis()->SetTitle("Scaled Monte Carlo/Data");

	  jewelPPHist_p->GetYaxis()->SetNdivisions(505);
	  jewelPPHist_p->GetXaxis()->SetTitleOffset(1.2);

	  jewelPPHist_p->GetXaxis()->SetTitleFont(titleFont);
	  jewelPPHist_p->GetYaxis()->SetTitleFont(titleFont);
	  jewelPPHist_p->GetXaxis()->SetTitleSize(titleSize/(padSplit));
	  jewelPPHist_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));

	  jewelPPHist_p->GetXaxis()->SetLabelFont(titleFont);
	  jewelPPHist_p->GetYaxis()->SetLabelFont(titleFont);
	  jewelPPHist_p->GetXaxis()->SetLabelSize(labelSize/(padSplit));
	  jewelPPHist_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));

	  jewelPPHist_p->GetYaxis()->SetTitle("Scaled Monte Carlo/Data");
	  jewelPPHist_p->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0 - padSplit));

	  jewelPPHist_p->GetXaxis()->SetTitle(ppHist_p->GetXaxis()->GetTitle());

	  /*
	  Double_t max = TMath::Max(jewelPPHist_p->GetMaximum(), pyt8Truth_p->GetMaximum());
	  Double_t min = TMath::Min(jewelPPHist_p->GetMinimum(), pyt8Truth_p->GetMinimum());

	  if(max < 1.0) max = 1.0;
	  */
	  jewelPPHist_p->SetMaximum(maxPytRatVal);
	  jewelPPHist_p->SetMinimum(minPytRatVal);

	  jewelPPHist_p->GetXaxis()->SetTitleSize(jewelPPHist_p->GetXaxis()->GetTitleSize()*1.3);
	  jewelPPHist_p->DrawCopy("E X0 HIST");
	  pyt8Truth_p->DrawCopy("E X0 HIST SAME");

	  TLine* line_p = new TLine();
	  line_p->SetLineStyle(2);
	  line_p->DrawLine(jewelPPHist_p->GetBinLowEdge(1), 1.0, jewelPPHist_p->GetBinLowEdge(jewelPPHist_p->GetXaxis()->GetNbins()+1), 1.0);
	  delete line_p;

	  gPad->SetTicks();

	  jewelPPHist_p->SetLineColor(1);
	}
	pads_p[1]->cd();
      }

      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "OverPP_" + saveTag + "_" + dateStr + "." + saveExt;
      plotSyst(pbpbHist_p, ratSyst, systNames, systLabels, "#frac{" + centBinsLabel[cI] + "}{pp}", plotSaveName);

      systTypeQuadSumVals.clear();
      for(unsigned int sysI = 0; sysI < uniqueSystTypes.size(); ++sysI){
	plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "OverPP_" + uniqueSystTypes[sysI] + "_" + saveTag + "_" + dateStr + "." + saveExt;

        std::vector<std::vector<Double_t > > tempSystValVect = systTypeToSystValsRatio[uniqueSystTypes[sysI]];
        std::vector<std::string> tempSystNameVect = systTypeToSystNames[uniqueSystTypes[sysI]];

	plotSyst(pbpbHist_p, tempSystValVect, tempSystNameVect, systLabels, "#frac{" + centBinsLabel[cI] + "}{pp}", plotSaveName, systNameToLegName(uniqueSystTypes[sysI]) + " Total");

	Int_t nBinsTemp = nVarBins;
        if(doXTrunc) nBinsTemp = nVarBinsTrunc;

        std::vector<Double_t> totalSysts;
        for(Int_t bI = 0; bI < nBinsTemp; ++bI){

          Float_t total = 0.0;
          for(unsigned int sI = 0; sI < tempSystValVect.size(); ++sI){
            total += tempSystValVect[sI][bI]*tempSystValVect[sI][bI];
          }

          totalSysts.push_back(TMath::Sqrt(total));
        }

        systTypeQuadSumVals.push_back(totalSysts);
      }

      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "OverPP_SystType_" + saveTag + "_" + dateStr + "." + saveExt;
      plotSyst(pbpbHist_p, systTypeQuadSumVals, uniqueSystTypes, systLabels, "#frac{" + centBinsLabel[cI] + "}{pp}", plotSaveName, "Total", systMin, systMax);

      canv_p->cd();
      pads_p[1]->cd();

      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(pbpbHist_p->GetBinLowEdge(1), 1.0, pbpbHist_p->GetBinLowEdge(pbpbHist_p->GetXaxis()->GetNbins()+1), 1.0);
      delete line_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      canv_p->cd();

      TBox* box_p = new TBox();
      box_p->SetFillColor(0);
      box_p->DrawBox(0.11, 0.435, leftMargin - 0.0025, 0.46);
      delete box_p;

      pbpbHist_p->SetFillColorAlpha(pbpbColor, 0.3);
      ppHist_p->SetFillColorAlpha(ppColor, 0.3);

      //      gPad->Modified();
      //CMcGinn: Testing this is the results plot save
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + "." + saveExt);
      for(Int_t pI = 0; pI < nPad; ++pI){
	delete pads_p[pI];
      }
      delete canv_p;
      pbpbHist_p->SetFillColor(0);
      ppHist_p->SetFillColor(0);

      if(cI == 0){
	Int_t nLegEntriesPP = 2;
	if(doJEWEL) nLegEntriesPP = 3;
	padsPP_p[0]->cd();

	TLegend* legPP_p = new TLegend(legX, legY - 0.065*nLegEntriesPP - 0.07, legX+0.25, legY - 0.07);
	legPP_p->SetTextFont(42);
	legPP_p->SetTextSize(titleSize/(1.0 - padSplit));
	legPP_p->SetBorderSize(0);
	legPP_p->SetFillStyle(0);

	if(doJEWEL){
	  HIJet::Style::EquipHistogram(jewelPPHist_p, 2);
	  jewelPPHist_p->SetMarkerSize(0.00001);
	  jewelPPHist_p->SetLineWidth(8);
	}

	legPP_p->AddEntry(ppHist_p, "#it{pp}", "P L F");

	/*
	legPP_p->AddEntry(pyt8Truth_p, ("PYTHIA 8 (x " + prettyString(scaleFactorPYT, 1, false) + ")").c_str(), "L P");
	if(doJEWEL) legPP_p->AddEntry(jewelPPHist_p, ("JEWEL #it{pp} (x " + prettyString(scaleFactorPYT, 1, false) + ")").c_str(), "L P");
	*/

	legPP_p->AddEntry(pyt8Truth_p, "PYTHIA 8", "L P");
	if(doJEWEL) legPP_p->AddEntry(jewelPPHist_p, "JEWEL #it{pp}", "L P");

	legPP_p->Draw("SAME");

	TLatex* label_p = new TLatex();
	label_p->SetNDC();
	label_p->SetTextFont(42);
	label_p->SetTextSize(titleSize/(1.0 - padSplit));
	label_p->SetTextAlign(31);

	std::vector<std::string> tempLabels;
	Int_t currLabelPos=0;
	for(unsigned int i = 0; i < labels.size(); ++i){
	  std::string label = labels[i];
	  if(label.find("Pb+Pb") != std::string::npos) continue;

	  Float_t yLow = 0.063;
	  if(i > 0){
	    if(labels[i-1].find("_{T,Jet}") != std::string::npos) yLow += 0.003;
	  }

	  label_p->DrawLatex(labelX, labelY - yLow*currLabelPos, label.c_str());
	  ++currLabelPos;
	}
	label_p->SetTextAlign(11);
	label_p->DrawLatex(var2JetXLabelX, var2JetXLabelY, "#it{pp} #rightarrow #gamma + 2 jets + X");
	delete label_p;

	//Draw a white box to get rid of the zero
	canvPP_p->cd();
	TBox* box_p = new TBox();
	box_p->SetFillColor(0);
	box_p->DrawBox(0.13, .45, leftMargin-0.0025, 0.5);
	if(isStrSame(varNameLower, "drjj")){
	  box_p->DrawBox(0.1, .41, leftMargin-0.0025, 0.45);
	}
	delete box_p;


	ppHist_p->SetFillColorAlpha(ppColor, 0.3);
	quietSaveAs(canvPP_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI) + "_PPTheory_" + saveTag + "_" + dateStr + "." + saveExt);
	for(Int_t pI = 0; pI < nPad; ++pI){
	  delete padsPP_p[pI];
	}
	delete canvPP_p;

	ppHist_p->SetFillColor(0);
      }
      delete leg_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      delete pbpbHist_p;
      if(doJEWEL){
	delete jewelLeg_p;
	delete jewelPPHist_p;
	delete jewelPPHistCurve_p;
	if(jewelPbPbHist_p != nullptr){
	  delete jewelPbPbHist_p;
	  delete jewelPbPbHistCurve_p;
	}
      }

      if(doJETSCAPEPbPb){
	delete inJETSCAPEPbPbHist_p;
      }

      if(doLBTPbPb){
	delete inLBTPbPbHist_p;
      }

      delete pyt8Truth_p;
    }

    Float_t integral = 0.0;
    Float_t integralError = 0.0;
    for(Int_t bIX = 0; bIX < ppHist_p->GetXaxis()->GetNbins(); ++bIX){
      integral += ppHist_p->GetBinContent(bIX+1)*ppHist_p->GetXaxis()->GetBinWidth(bIX+1);
      Float_t err = ppHist_p->GetBinError(bIX+1)*ppHist_p->GetXaxis()->GetBinWidth(bIX+1);
      integralError += err*err;
    }
    integralPerCentrality_h->SetBinContent(nCentIntBins, integral);
    integralPerCentrality_h->SetBinError(nCentIntBins, integralError);

    meanPerCentrality_h->SetBinContent(nCentIntBins, ppHist_p->GetMean());
    meanPerCentrality_h->SetBinError(nCentIntBins, ppHist_p->GetMeanError());

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    //Plot integral per centrality, in brackets to just have a contained scope
    {
      TCanvas* canv_p = new TCanvas("canv_p", "", width, 900.0);
      setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin*padSplit);

      HIJet::Style::EquipHistogram(integralPerCentrality_h, 0);
      integralPerCentrality_h->SetMaximum(1.5);
      integralPerCentrality_h->SetMinimum(0.0);

      integralPerCentrality_h->DrawCopy("HIST E P");

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(titleSize);
      label_p->SetTextAlign(31);

      label_p->DrawLatex(labelX, labelY, std::string(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false)).c_str());
      delete label_p;

      canv_p->SaveAs(std::string("pdfDir/" + dateStr + "/integralPerCentrality_GammaPt" + gI + "_" + saveTag + "_" + dateStr + "." + saveExt).c_str());

      delete canv_p;
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    //Plot mean per centrality, in brackets to just have a contained scope
    {
      TCanvas* canv_p = new TCanvas("canv_p", "", width, 900.0);
      setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin*padSplit);

      HIJet::Style::EquipHistogram(meanPerCentrality_h, 0);
      meanPerCentrality_h->SetMaximum(1.5);
      meanPerCentrality_h->SetMinimum(0.0);

      meanPerCentrality_h->DrawCopy("HIST E P");

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(titleSize);
      label_p->SetTextAlign(31);

      label_p->DrawLatex(labelX, labelY, std::string(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false)).c_str());
      delete label_p;

      canv_p->SaveAs(std::string("pdfDir/" + dateStr + "/meanPerCentrality_GammaPt" + gI + "_" + saveExt + "_" + dateStr + "." + saveExt).c_str());

      delete canv_p;
    }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    //Plot centrality overlay
    {
      TCanvas* canv_p = new TCanvas("canv_p", "", width, 900.0);
      setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin*padSplit);

      TLegend* leg_p = new TLegend(legXCent, legYCent - 0.04*ratioHistsPerCentrality.size(), legXCent+0.25, legYCent);
      leg_p->SetTextFont(42);
      leg_p->SetTextSize(titleSize);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	HIJet::Style::EquipHistogram(ratioHistsPerCentrality[cI], cI+1);

	ratioHistsPerCentrality[cI]->SetMinimum(minRatVal);
	ratioHistsPerCentrality[cI]->SetMaximum(maxRatVal);

	if(cI == 0) ratioHistsPerCentrality[cI]->DrawCopy("HIST E P X0");
	else ratioHistsPerCentrality[cI]->DrawCopy("HIST E P SAME X0");

	drawSyst(canv_p, ratioHistsPerCentrality[cI], {ratioSystsPerCentrality[cI]});

	std::string legLabel = centBinsLabel[cI].substr(centBinsLabel[cI].find(",")+2, centBinsLabel[cI].size());
	if(legLabel.rfind("}") != std::string::npos) legLabel.replace(legLabel.rfind("}"), 1, "");
	leg_p->AddEntry(ratioHistsPerCentrality[cI], legLabel.c_str(), "P L F");
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(ratioHistsPerCentrality[0]->GetBinLowEdge(1), 1.0, ratioHistsPerCentrality[0]->GetBinLowEdge(ratioHistsPerCentrality[0]->GetXaxis()->GetNbins()+1), 1.0);
      delete line_p;

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(titleSize);
      if(labelAlignRightCent) label_p->SetTextAlign(31);

      for(unsigned int i = 0; i < ratioLabelsPerCentrality[0].size(); ++i){
	label_p->DrawLatex(labelXCent, labelYCent - 0.05*i, ratioLabelsPerCentrality[0][i].c_str());
      }
      delete label_p;

      leg_p->Draw("SAME");
      gPad->SetTicks();

      /*
      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	ratioHistsPerCentrality[cI]->SetFillColorAlpha(centColors[cI], 0.3);
	}*/


      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;
      //      gPad->Modified();
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI)	 + "_AllCentRatio_" + saveTag + "_" + dateStr + "." + saveExt);
      delete canv_p;
      delete leg_p;

      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	ratioHistsPerCentrality[cI]->SetFillColor(0);
      }
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    //Plot centrality 0-10% RAA w/ JEWEL
    if(ratioJEWEL0to10Curve_p != nullptr){
      Double_t newHeight = 900.0;
      TCanvas* canv_p = new TCanvas("canv_p", "", width, newHeight);
      setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin*padSplit);

      TLegend* leg_p = new TLegend(ratLegX, ratLegY - 0.04*ratioHistsPerCentrality.size(), ratLegX+0.25, ratLegY);
      leg_p->SetTextFont(42);
      leg_p->SetTextSize(titleSize*height/newHeight);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);

      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	std::string legLabel = centBinsLabel[cI].substr(centBinsLabel[cI].find(",")+2, centBinsLabel[cI].size());
	if(legLabel.rfind("}") != std::string::npos) legLabel.replace(legLabel.rfind("}"), 1, "");

	if(legLabel.find("0-10%") == std::string::npos) continue;

	HIJet::Style::EquipHistogram(ratioHistsPerCentrality[cI], cI+1);
	ratioHistsPerCentrality[cI]->SetMarkerColor(1);
	ratioHistsPerCentrality[cI]->SetLineColor(1);
	ratioHistsPerCentrality[cI]->SetMarkerStyle(20);

	ratioHistsPerCentrality[cI]->SetMinimum(minRatVal);
	ratioHistsPerCentrality[cI]->SetMaximum(maxRatVal);

	ratioHistsPerCentrality[cI]->GetXaxis()->SetTitleSize(titleSize*height/newHeight);
	ratioHistsPerCentrality[cI]->GetYaxis()->SetTitleSize(titleSize*height/newHeight);

	ratioHistsPerCentrality[cI]->GetXaxis()->SetTitleSize(ratioHistsPerCentrality[cI]->GetXaxis()->GetTitleSize()*1.3);
	if(cI == 0) ratioHistsPerCentrality[cI]->DrawCopy("HIST E X0 P");
	else ratioHistsPerCentrality[cI]->DrawCopy("HIST E P X0 SAME");

	//	HIJet::Style::EquipHistogram(ratioHistsPerCentrality[cI], cI+1);
	ratioHistsPerCentrality[cI]->SetFillColorAlpha(pbpbColor, 0.3);
	ratioHistsPerCentrality[cI]->SetMarkerColor(pbpbColor);
	drawSyst(canv_p, ratioHistsPerCentrality[cI], {ratioSystsPerCentrality[cI]});
	ratioHistsPerCentrality[cI]->SetMarkerColor(1);
	ratioHistsPerCentrality[cI]->SetMarkerColor(1);
	ratioHistsPerCentrality[cI]->SetLineColor(1);
	ratioHistsPerCentrality[cI]->SetMarkerStyle(20);

	leg_p->AddEntry(ratioHistsPerCentrality[cI], legLabel.c_str(), "P L F");
      }

      HIJet::Style::EquipHistogram(ratioJEWEL0to10Curve_p, 1);
      //      ratioJEWEL0to10Curve_p->SetLineStyle(9);
      ratioJEWEL0to10Curve_p->SetLineWidth(8);
      ratioJEWEL0to10Curve_p->SetMarkerSize(0.0000001);

            if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;


      ratioJEWEL0to10Curve_p->DrawCopy("HIST E X0 SAME");
      leg_p->AddEntry(ratioJEWEL0to10Curve_p, "JEWEL", "L");

      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(ratioHistsPerCentrality[0]->GetBinLowEdge(1), 1.0, ratioHistsPerCentrality[0]->GetBinLowEdge(ratioHistsPerCentrality[0]->GetXaxis()->GetNbins()+1), 1.0);
      delete line_p;

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(titleSize*height/newHeight);
      if(labelAlignRightCent) label_p->SetTextAlign(31);

      if(isStrSame(varNameLower, "drjj")){
	for(unsigned int i = 0; i < ratioLabelsPerCentrality[0].size(); ++i){
	  Float_t yLow = 0.063*newHeight/height;
	  if(i > 0){
	    if(ratioLabelsPerCentrality[0][i-1].find("_{T,Jet}") != std::string::npos) yLow += 0.003*newHeight/height;
	  }

	  if(i < 4){
	    label_p->DrawLatex(ratLabelX, ratLabelY - yLow*i, ratioLabelsPerCentrality[0][i].c_str());
	  }
	  else{
	    label_p->SetTextAlign(31);
	    label_p->DrawLatex(0.95, ratLabelY - yLow*(i-4), ratioLabelsPerCentrality[0][i].c_str());
	    if(!labelAlignRightCent) label_p->SetTextAlign(11);
	  }
	}
      }
      else{
	for(unsigned int i = 0; i < ratioLabelsPerCentrality[0].size(); ++i){
	  Float_t yLow = 0.063*newHeight/height;
	  if(i > 0){
	    if(ratioLabelsPerCentrality[0][i-1].find("_{T,Jet}") != std::string::npos) yLow += 0.003*newHeight/height;
	  }
	  std::cout << "LABEL, YLOW: " << ratioLabelsPerCentrality[0][i] << ", " << yLow << std::endl;
	  label_p->DrawLatex(ratLabelX, ratLabelY - yLow*i, ratioLabelsPerCentrality[0][i].c_str());
	}
      }
      delete label_p;

      leg_p->Draw("SAME");
      gPad->SetTicks();

      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	//	ratioHistsPerCentrality[cI]->SetFillColorAlpha(centColors[cI], 0.3);
	ratioHistsPerCentrality[cI]->SetFillColorAlpha(pbpbColor, 0.3);
      }

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;


      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI)	 + "_Cent0to10RatioWithJEWEL_" + saveTag + "_" + dateStr + "." + saveExt);
      delete canv_p;
      delete leg_p;

      for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
	ratioHistsPerCentrality[cI]->SetFillColor(0);
      }

    }
    delete ratioJEWEL0to10Curve_p;

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << std::endl;

    for(unsigned int cI = 0; cI < ratioHistsPerCentrality.size(); ++cI){
      delete ratioHistsPerCentrality[cI];
    }
    ratioHistsPerCentrality.clear();
    ratioSystsPerCentrality.clear();


    if(doJETSCAPE){
      delete inJETSCAPEPPHist_p;
    }

    if(doLBT){
      delete inLBTPPHist_p;
    }

    delete integralPerCentrality_h;
    delete meanPerCentrality_h;
    delete ppHist_p;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  inPbPbFile_p->Close();
  delete inPbPbFile_p;

  inPPFile_p->Close();
  delete inPPFile_p;


  std::cout << "GDJPLOTRESULTS COMPLETE. return 0." << std::endl;
  return 0;
}

//Main block; takes in config file and handles wrong inputs before calling gdjPlotResults
int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotResults.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjPlotResults(argv[1]);
  return retVal;
}
