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

//vector of syst, bin-by-bin
//As a start, assume symmetric syst
std::vector<std::vector<Double_t> > getSyst(TFile* histFile_p, TH1D* nominalHist_p, std::vector<std::string> systHistNames, std::vector<Int_t> yPos, Double_t scaleTotal)
{
  std::vector<std::vector<Double_t> > systVals;

  for(unsigned int jI = 0; jI < systHistNames.size(); ++jI){
    TH2D* temp2D_p = (TH2D*)histFile_p->Get(systHistNames[jI].c_str());
    TH1D* temp1D_p = (TH1D*)nominalHist_p->Clone("temp1D_p");
    for(Int_t bIX = 0; bIX < temp1D_p->GetXaxis()->GetNbins(); ++bIX){
      temp1D_p->GetBinContent(bIX+1, 0.0);
      temp1D_p->GetBinError(bIX+1, 0.0);
    }
    fineTH2ToCoarseTH1(temp2D_p, temp1D_p, yPos);
    binWidthAndScaleNorm(temp1D_p, scaleTotal);
    
    systVals.push_back({});   
    for(Int_t bIX = 0; bIX < temp1D_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t deltaVal = TMath::Abs(nominalHist_p->GetBinContent(bIX+1) - temp1D_p->GetBinContent(bIX+1));
      systVals[jI].push_back(deltaVal);
    }    

    delete temp1D_p;
  }
  
  return systVals;
}

std::vector<std::vector<Double_t> > getSystRat(TFile* histFilePbPb_p, TFile* histFilePP_p, TH1D* nominalHistPbPb_p, TH1D* nominalHistPP_p, std::vector<std::string> systHistNamesPbPb, std::vector<std::string> systHistNamesPP, std::vector<bool> isCorrelated, std::vector<Int_t> yPos, Double_t scaleTotalPbPb, Double_t scaleTotalPP)
{
  std::vector<std::vector<Double_t> > systVals;
  
  for(unsigned int jI = 0; jI < systHistNamesPbPb.size(); ++jI){
    TH2D* temp2DPbPb_p = (TH2D*)histFilePbPb_p->Get(systHistNamesPbPb[jI].c_str());
    TH1D* temp1DPbPb_p = (TH1D*)nominalHistPbPb_p->Clone("temp1DPbPb_p");

    TH2D* temp2DPP_p = (TH2D*)histFilePP_p->Get(systHistNamesPP[jI].c_str());
    TH1D* temp1DPP_p = (TH1D*)nominalHistPP_p->Clone("temp1DPP_p");

    for(Int_t bIX = 0; bIX < temp1DPbPb_p->GetXaxis()->GetNbins(); ++bIX){
      temp1DPbPb_p->SetBinContent(bIX+1, 0.0);
      temp1DPbPb_p->SetBinError(bIX+1, 0.0);

      temp1DPP_p->SetBinContent(bIX+1, 0.0);
      temp1DPP_p->SetBinError(bIX+1, 0.0);
    }
    fineTH2ToCoarseTH1(temp2DPbPb_p, temp1DPbPb_p, yPos);
    fineTH2ToCoarseTH1(temp2DPP_p, temp1DPP_p, yPos);
    
    binWidthAndScaleNorm(temp1DPbPb_p, scaleTotalPbPb);
    binWidthAndScaleNorm(temp1DPP_p, scaleTotalPP);
     
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

	deltaVal += TMath::Sqrt(deltaValPbPb*deltaValPbPb + deltaValPP*deltaValPP);
      }
      
      systVals[jI].push_back(deltaVal);
    }    

    delete temp1DPbPb_p;
    delete temp1DPP_p;    
  }
  
  return systVals;
}


void drawSyst(TPad* pad_p, TH1D* nominalHist_p, std::vector<std::vector<Double_t> > systVals)
{
  pad_p->cd();
  Int_t color = nominalHist_p->GetMarkerColor();
  TBox* box_p = new TBox();

  for(Int_t bIX = 0; bIX < nominalHist_p->GetXaxis()->GetNbins(); ++bIX){
    Float_t systQuadSum = 0.0;
    for(unsigned int sI = 0; sI < systVals.size(); ++sI){
      systQuadSum += systVals[sI][bIX]*systVals[sI][bIX];
    }    
    systQuadSum = TMath::Sqrt(systQuadSum);

    Float_t binVal = nominalHist_p->GetBinContent(bIX+1);
    Float_t binLow = nominalHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
    Float_t binHigh = nominalHist_p->GetXaxis()->GetBinLowEdge(bIX+2);

    box_p->SetFillColorAlpha(color, 0.3);

    Double_t minVal = TMath::Max(binVal - systQuadSum, nominalHist_p->GetMinimum());
    Double_t maxVal = TMath::Min(binVal + systQuadSum, nominalHist_p->GetMaximum());

    box_p->DrawBox(binLow, minVal, binHigh, maxVal);
  }

  delete box_p;
  
  return;
}


std::vector<Double_t> plotSyst(TH1D* nominalHist_p, std::vector<std::vector<Double_t> > systVals, std::vector<std::string> systNames, std::vector<std::string> labels, std::string yStr, std::string saveName)
{
  std::vector<Double_t> systTotal;
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  //Series of consistency checks
  const Int_t nMaxHist = 36;
  if(systVals.size() + 1 > nMaxHist){
    std::cout << "Needed number of histograms '" << systVals.size()+1 << "' exceeds plotSyst nMaxHist="  << nMaxHist << ". return {}" << std::endl;
    return {};
  }

  if(systVals.size() != systNames.size()){
    std::cout << "SystVals and systNames have different sizes \'" << systVals.size() << "\', \'" << systNames.size() << "\', respectively. return {}" << std::endl;
    return {};
  }

  //Define some plotting parameters
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;

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
  
  //Define our histogram array and new the total syst hist
  std::string titleStr = ";" + std::string(nominalHist_p->GetXaxis()->GetTitle()) + ";" + yStr + " Syst. (%)";
  TH1D* systHist_p[nMaxHist];
  systHist_p[0] = new TH1D("systHist_Total_h", titleStr.c_str(), nBinsX, binsX);
  HIJet::Style::EquipHistogram(systHist_p[0], 0);
  
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    systTotal.push_back(0.0);
  }
  
  for(unsigned int sI = 0; sI < systVals.size(); ++sI){
    systHist_p[sI+1] = new TH1D(("systHist_" + systNames[sI] + "_h").c_str(), titleStr.c_str(), nBinsX, binsX);

    for(Int_t bIX = 0; bIX < nBinsX; ++bIX){ 
      systHist_p[sI+1]->SetBinContent(bIX+1, systVals[sI][bIX]/nominalHist_p->GetBinContent(bIX+1));
      systHist_p[sI+1]->SetBinError(bIX+1, 0);

      systTotal[bIX] += systVals[sI][bIX]*systVals[sI][bIX];
    }    

    HIJet::Style::EquipHistogram(systHist_p[sI+1], sI+1);
  }

  Double_t maxVal = 0.0;
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    Float_t val = TMath::Sqrt(systTotal[bIX])/nominalHist_p->GetBinContent(bIX+1);
    if(val > maxVal) maxVal = val;

    systHist_p[0]->SetBinContent(bIX+1, val);
    systHist_p[0]->SetBinError(bIX+1, 0.0);
  }

  systHist_p[0]->SetMaximum(maxVal*1.2);
  systHist_p[0]->SetMinimum(-maxVal*1.2);
  
  TCanvas* canv_p = new TCanvas("canvSyst_p", "", width, height);
  canv_p->SetTopMargin(0.20);
  canv_p->SetRightMargin(0.25);
  canv_p->SetLeftMargin(leftMargin);
  canv_p->SetBottomMargin(bottomMargin);

  canv_p->cd();

  for(Int_t sI = 0; sI < systVals.size(); ++sI){
    if(sI == 0) systHist_p[sI]->DrawCopy("HIST");
    else systHist_p[sI]->DrawCopy("HIST SAME");

    //Create the mirrored version for negative values
    for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
      systHist_p[sI]->SetBinContent(bIX+1, -systHist_p[sI]->GetBinContent(bIX+1));
    }

    systHist_p[sI]->DrawCopy("HIST SAME"); 
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(titleFont);
  label_p->SetTextSize(titleSize);
  //  label_p->SetTextAlign(31);

  Float_t xPos = 0.01;
  Float_t maxLength = 0.0;
  for(unsigned int i = 0; i < labels.size(); ++i){
    if(i != 0 && i%3 == 0){
      xPos += maxLength + 0.05;
      maxLength = 0;
    }

    label_p->SetText(xPos, 0.96 - 0.06*(i%3), labels[i].c_str());
    if(label_p->GetXsize() > maxLength) maxLength = label_p->GetXsize();

    label_p->DrawLatex(xPos, 0.96 - 0.06*(i%3), labels[i].c_str());
  }
  delete label_p;

  std::vector<TLegend*> legs_p = {nullptr, nullptr};
  if(systVals.size() <= nMaxHist/2) legs_p[0] = new TLegend(0.81, 0.81 - systVals.size()*0.045, 0.90, 0.81);
  else{
    legs_p[0] = new TLegend(0.76, 0.81 - (nMaxHist/2)*0.045, 0.87, 0.81);
    legs_p[1] = new TLegend(0.87, 0.81 - (systVals.size() - nMaxHist/2)*0.045, 0.99, 0.81);

    legs_p[1]->SetTextFont(titleFont);
    legs_p[1]->SetTextSize(titleSize);
    legs_p[1]->SetBorderSize(0);
    legs_p[1]->SetFillColor(0);
  }   
  legs_p[0]->SetTextFont(titleFont);
  legs_p[0]->SetTextSize(titleSize);
  legs_p[0]->SetBorderSize(0);
  legs_p[0]->SetFillColor(0);
  
  legs_p[0]->AddEntry(systHist_p[0], "Total", "L");
  for(Int_t sI = 0; sI < systVals.size(); ++sI){
    legs_p[(sI+1)/(nMaxHist/2)]->AddEntry(systHist_p[sI+1], systNames[sI].c_str(), "L");
  }

  legs_p[0]->Draw("SAME");
  if(systVals.size() > nMaxHist/2) legs_p[1]->Draw("SAME");
  

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(binsX[0], 0.0, binsX[nBinsX], 0.0);  
  delete line_p;
  
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  delete legs_p[0];
  if(systVals.size() > nMaxHist/2) delete legs_p[1];

  for(unsigned int sI = 0; sI < systVals.size()+1; ++sI){delete systHist_p[sI];}
  
  return systTotal;
}


int gdjPlotResults(std::string inConfigFileName)
{
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 1.0;
  
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
					      "DOJTPTCUT",
					      "SAVETAG",
					      "SAVEEXT"};

  //Define a second set once we've accessed the file and know what var we are working with
  std::vector<std::string> necessaryVarParams = {"VARMINY",
						 "VARMAXY",
						 "VARRATMINY",
						 "VARRATMAXY",
						 "VARDOLOGX",
						 "VARDOLOGY",
						 "VARLEGX",
						 "VARLEGY",
  						 "VARLABELX",
						 "VARLABELY"};
  
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  const std::string inPbPbFileName = config_p->GetValue("INPBPBFILENAME", "");
  const std::string inPPFileName = config_p->GetValue("INPPFILENAME", "");

  const std::string inPbPbTermFileName = config_p->GetValue("INPBPBTERMFILENAME", "");
  const std::string inPPTermFileName = config_p->GetValue("INPPTERMFILENAME", "");

  //Check all files exist before continuing to grab params
  if(!check.checkFileExt(inPbPbFileName, ".root")) return 1;
  if(!check.checkFileExt(inPPFileName, ".root")) return 1;
  if(!check.checkFileExt(inPbPbTermFileName, ".config")) return 1;
  if(!check.checkFileExt(inPPTermFileName, ".config")) return 1;

  
  const Bool_t doJtPtCut = config_p->GetValue("DOJTPTCUT", 0);  
					        
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
  bool doLogX = config_p->GetValue((varNameUpper + "DOLOGX").c_str(), 0);
  bool doLogY = config_p->GetValue((varNameUpper + "DOLOGY").c_str(), 0);
  float legX = config_p->GetValue((varNameUpper + "LEGX").c_str(), 0.1);
  float legY = config_p->GetValue((varNameUpper + "LEGY").c_str(), 0.1);
  float labelX = config_p->GetValue((varNameUpper + "LABELX").c_str(), 0.1);
  float labelY = config_p->GetValue((varNameUpper + "LABELY").c_str(), 0.1);
  
  std::string varPrefix = "";
  if(isStrSame(varNameUpper, "PT")) varPrefix = "JT";

  std::string binVarStr = varPrefix + varNameUpper;
  if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";

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
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(!checkEnvForParams(inPbPbFileConfig_p, necessaryFileParams)) return 1;
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(!checkEnvForParams(inPPFileConfig_p, necessaryFileParams)) return 1;
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!compEnvParams(inPbPbFileConfig_p, inPPFileConfig_p, paramsToCompare)) return 1;

  //Grab the jes/jer systs
  const Int_t nMaxSyst = 100;
  Int_t nJESJER = inPPFileConfig_p->GetValue("NJETSYSANDNOM", -1);
  std::vector<std::string> jesJERStrVect = strToVect(inPPFileConfig_p->GetValue("JETSYSANDNOM", ""));

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Grab all systs
  Int_t nSyst = inPPFileConfig_p->GetValue("NSYST", -1);
  std::vector<std::string> systStrVect = strToVect(inPPFileConfig_p->GetValue("SYST", ""));

  //Temp for jtptcut
  for(unsigned int sI = 0; sI < systStrVect.size(); ++sI){
    if(systStrVect[sI].find("JTPTCUT") != std::string::npos){
      systStrVect.erase(systStrVect.begin() + sI);
      --nSyst;
      std::cout << "REMOVING JTPTCUT" << std::endl;
      break;
    }
  }
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  
  //Add to global labels
  int jetR = inPbPbFileConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
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

  for(unsigned int cI = 0; cI < tempCentBinsStr.size()-1; ++cI){
    centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
    centBinsLabel.push_back("#bf{Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%}");
  }  
  
  //Define a bunch of parameters for the TCanvas
  const Int_t nPad = 2;
  Double_t padSplit = 0.45;
  const Double_t leftMargin = 0.16;
  const Double_t bottomMargin = 0.125/padSplit;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;
  Double_t height = 900.0;
  Double_t width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);
  height = 1200.0;
  
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

  std::vector<double> subJtGammaPtBinsV = subJtGammaPtBinFlattener.GetFlattenedBins(subJtGammaPtMin, subJtGammaPtMax); 
  Int_t nSubJtGammaPtBins = nGammaPtBins*nSubJtPtBins;
  
  gStyle->SetOptStat(0);

  const Int_t iterPPPhoPt = inPPTermFileConfig_p->GetValue("GAMMAPT_PP", -1); 

  const Int_t nCentIntBins = 5;
  Double_t centIntBins[nCentIntBins+1] = {0,10,30,80,100,120};
  
  for(Int_t gI = 1; gI < nGammaPtBins-1; ++gI){
    //Create a few summary plots
    TH1D* integralPerCentrality_h = new TH1D("integralPerCentrality_h", ";Centrality/p+p;R_{JJ#gamma}", nCentIntBins, centIntBins);
    TH1D* meanPerCentrality_h = new TH1D("meanPerCentrality_h", ";Centrality/p+p;#LT#vec{x}_{JJ#gamma}#GT", nCentIntBins, centIntBins);

    Int_t iterPP = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_PP_Nominal").c_str(), -1);

    if(iterPP < 0){
      std::cout << "Terminating iteration is \'" << iterPP << "\' in pp, for \'" << varNameUpper + "_GAMMAPTALL_PP_NOMINAL" << "\'. return 1" << std::endl;
      return 1;
    }

    TH1D* ppHistPhoPt_p = (TH1D*)inPPFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPPPhoPt) + "_PP_PURCORR_COMBINED_h").c_str());
    Double_t scaledTotalPP = 0.0;
    for(Int_t bIX = 0; bIX < ppHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = ppHistPhoPt_p->GetBinCenter(bIX+1);
      if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;

      scaledTotalPP += ppHistPhoPt_p->GetBinContent(bIX+1);
    }
    
    std::string ppHistName = "photonPtJet" + varNameUpper + "Reco_Iter" + std::to_string(iterPP) + "_PP_Nominal_PURCORR_COMBINED_h";
    std::vector<std::string> ppSystHistNames;
    std::vector<bool> isSystCorrelated;
    std::vector<std::string> systNames;
    for(unsigned int jI = 0; jI < systStrVect.size(); ++jI){
      if(isStrSame(systStrVect[jI], "Nominal")) continue;

      systNames.push_back(systStrVect[jI]);
      
      if(isStrSame(systStrVect[jI], "JESRTRK")) isSystCorrelated.push_back(false);
      else if(isStrSame(systStrVect[jI], "PRIOR")) isSystCorrelated.push_back(false);
      else isSystCorrelated.push_back(true);
      
      //WE WILL NOT BE USING THIS CURRENTLY, INSTEAD TAKE THE SAME ITERATION AS NOMINAL FOR NOW
      Int_t iterPPSyst = inPPTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_PP_" + systStrVect[jI]).c_str(), -1);
      ppSystHistNames.push_back("photonPtJet" + varNameUpper + "Reco_Iter" + std::to_string(iterPPSyst) + "_PP_" + systStrVect[jI] + "_PURCORR_COMBINED_h");
    }

    //We need to find all gamma bins matching our current bin for creating our final histogram
    std::vector<int> gammaMatchedBins;
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
    
    TH2D* ppHist2D_p = (TH2D*)inPPFile_p->Get(ppHistName.c_str());
    TH1D* ppHist_p = nullptr;
    if(!doXTrunc) ppHist_p = new TH1D("ppHist1D_h", ";;", nVarBins, varBins);
    else ppHist_p = new TH1D("ppHist1D_h", ";;", nVarBinsTrunc, varBinsTrunc);
    fineTH2ToCoarseTH1(ppHist2D_p, ppHist_p, gammaMatchedBins);

    binWidthAndScaleNorm(ppHist_p, scaledTotalPP);
    std::vector<std::vector<Double_t > > ppSyst = getSyst(inPPFile_p, ppHist_p, ppSystHistNames, gammaMatchedBins, scaledTotalPP);
 
    //    ppHist_p->Scale(1./ppHist_p->Integral());
    HIJet::Style::EquipHistogram(ppHist_p, 2); 
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
    
    if(isMultijet) ppHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{#gamma}*}#frac{dN_{JJ#gamma}}{d" + varNameLabel + "}").c_str());
    else ppHist_p->GetYaxis()->SetTitle(("#frac{1}{N_{#gamma}*}#frac{dN_{J#gamma}}{d" + varNameLabel + "}}").c_str());
    ppHist_p->GetXaxis()->SetTitle(varNameLabel.c_str());

  
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){   
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
 
      const Int_t iterPbPbPhoPt = inPbPbTermFileConfig_p->GetValue(("GAMMAPT_" + centBinsStr[cI]).c_str(), -1);
   
      Int_t iterPbPb = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_" + centBinsStr[cI] + "_Nominal").c_str(), -1);
      if(iterPbPb < 0){
	std::cout << "Terminating iteration is \'" << iterPbPb << "\' in pbpb, \'" << centBinsStr[cI] << "\' for \'" << varNameUpper + "_GAMMAPTALL_" + centBinsStr[cI] << "\'. return 1" << std::endl;
	return 1;
      }

      std::string centBins2Str = centBinsStr[cI];
      if(centBins2Str.find("to") != std::string::npos){
	centBins2Str.replace(0, 4, "");
	centBins2Str.replace(centBins2Str.find("to"), 2, "-");
	centBins2Str = centBins2Str + "%";
      }
    
      std::string pbpbHistName = "photonPtJet" + varNameUpper + "Reco_Iter" + std::to_string(iterPbPb) + "_" + centBinsStr[cI] + "_Nominal_PURCORR_COMBINED_h";
      std::vector<std::string> pbpbSystHistNames;      
      for(unsigned int systI = 0; systI < systStrVect.size(); ++systI){
	if(isStrSame(systStrVect[systI], "Nominal")) continue;

	//WE WILL NOT BE USING THIS CURRENTLY, INSTEAD TAKE THE SAME ITERATION AS NOMINAL FOR NOW
	Int_t iterPbPbSyst = inPbPbTermFileConfig_p->GetValue((varNameUpper + "_GAMMAPTALL_" + centBinsStr[cI] + "_" + systStrVect[systI]).c_str(), -1);
	pbpbSystHistNames.push_back("photonPtJet" + varNameUpper + "Reco_Iter" + std::to_string(iterPbPbSyst) + "_" + centBinsStr[cI] + "_" + systStrVect[systI] + "_PURCORR_COMBINED_h");
      }

      TH2D* pbpbHist2D_p = (TH2D*)inPbPbFile_p->Get(pbpbHistName.c_str());
      TH1D* pbpbHist_p = nullptr;
      if(!doXTrunc) pbpbHist_p = new TH1D("pbpbHist1D_h", ";;", nVarBins, varBins);
      else pbpbHist_p = new TH1D("pbpbHist1D_h", ";;", nVarBinsTrunc, varBinsTrunc);
      fineTH2ToCoarseTH1(pbpbHist2D_p, pbpbHist_p, gammaMatchedBins);

      TH1D* pbpbHistPhoPt_p = nullptr;          
      if(iterPbPbPhoPt > 0) pbpbHistPhoPt_p = (TH1D*)inPbPbFile_p->Get(("photonPtReco_Iter" + std::to_string(iterPbPbPhoPt) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      else{
	std::cout << "No good termination point found. return 1" << std::endl;
	return 1;
      }
      
      Double_t scaledTotalPbPb = 0.0;
      for(Int_t bIX = 0; bIX < pbpbHistPhoPt_p->GetXaxis()->GetNbins(); ++bIX){
	Float_t binCenter = pbpbHistPhoPt_p->GetBinCenter(bIX+1);
	if(binCenter < gammaPtBins[gI] || binCenter >= gammaPtBins[gI+1]) continue;	

	scaledTotalPbPb += pbpbHistPhoPt_p->GetBinContent(bIX+1);
      }

          
      binWidthAndScaleNorm(pbpbHist_p, scaledTotalPbPb);      
      std::vector<std::vector<Double_t > > pbpbSyst = getSyst(inPbPbFile_p, pbpbHist_p, pbpbSystHistNames, gammaMatchedBins, scaledTotalPbPb);
      //      pbpbHist_p->Scale(1./pbpbHist_p->Integral());

      TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
      TPad* pads_p[nPad];
      canv_p->SetTopMargin(noMargin);
      canv_p->SetRightMargin(noMargin);
      canv_p->SetLeftMargin(noMargin);
      canv_p->SetBottomMargin(noMargin);
      canv_p->cd();
      
      pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
      pads_p[0]->SetTopMargin(topMargin);
      pads_p[0]->SetRightMargin(rightMargin);
      pads_p[0]->SetLeftMargin(leftMargin);
      pads_p[0]->SetBottomMargin(noMargin);
      
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      pads_p[0]->cd();
      canv_p->cd();
      
      canv_p->cd();
      pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
      pads_p[1]->SetTopMargin(noMargin);
      pads_p[1]->SetRightMargin(rightMargin);
      pads_p[1]->SetLeftMargin(leftMargin);
      pads_p[1]->SetBottomMargin(bottomMargin);
      
      canv_p->cd();
      pads_p[1]->Draw("SAME");
      pads_p[1]->cd();
      canv_p->cd();
      
      pads_p[0]->cd();      
      HIJet::Style::EquipHistogram(pbpbHist_p, 1); 

      ppHist_p->SetMinimum(minVal);
      ppHist_p->SetMaximum(maxVal);

      //Need to set this for both so the systematics dont overflow the plotting space but still draw correctly
      pbpbHist_p->SetMinimum(minVal);
      pbpbHist_p->SetMaximum(maxVal);
      
      ppHist_p->DrawCopy("HIST E1 P");
      pbpbHist_p->DrawCopy("HIST E1 P SAME");

      drawSyst(pads_p[0], ppHist_p, ppSyst);
      drawSyst(pads_p[0], pbpbHist_p, pbpbSyst);

      ppHist_p->DrawCopy("HIST E1 P SAME");
      pbpbHist_p->DrawCopy("HIST E1 P SAME");

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
      
      std::vector<std::string> labelsAlt;
      std::vector<std::string> labelsTemp = getLabels(inPbPbFileConfig_p, pbpbHist_p, &labelMap, &labelsAlt);
      std::vector<std::string> labels = {"#bf{#it{ATLAS Internal}}", "Pb+Pb, #it{pp} at #sqrt{s_{NN}}=5.02 TeV"};

      for(unsigned int lI = 0; lI < labelsTemp.size(); ++lI){
	if(isStrSame(labelsTemp[lI], "h")) continue;

	if(doJtPtCut){
	  if(labelsTemp[lI].find("p_{T,Jet}^{Truth}") != std::string::npos){
	    //	  std::cout << "LABEL FOUND: " << labelsTemp[lI] << std::endl;
	    labelsTemp[lI].replace(0,labelsTemp[lI].find("<"), "");
	    labelsTemp[lI] = std::to_string(((int)jtPtBinsLowReco)) + labelsTemp[lI];
	  }
	}
	
	labels.push_back(labelsTemp[lI]);
      }
      labels.push_back(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false));
      labels.push_back(gammaJtDPhiCutLabel);

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
      
      unsigned int pos = 0;
      while(pos < labels.size()){
	if(labels[pos].find("PURCORR") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("COMBINED") != std::string::npos) labels.erase(pos + labels.begin());
	else if(labels[pos].find("Iter") != std::string::npos) labels.erase(pos + labels.begin());
	else ++pos;
      }
 
      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
     
      std::string plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_PP_" + saveTag + "_" + dateStr + "." + saveExt;
      std::vector<std::string> systLabels = labels;
      systLabels[1] = "#it{pp} at #sqrt{s_{NN}}=5.02 TeV";
      if(cI == 0) plotSyst(ppHist_p, ppSyst, systNames, systLabels, "pp", plotSaveName);            

      pbpbHist_p->GetXaxis()->SetTitle(ppHist_p->GetXaxis()->GetTitle());      
      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + "." + saveExt;      
      systLabels[1] = centBinsLabel[cI] + " at #sqrt{s_{NN}}=5.02 TeV";
      plotSyst(pbpbHist_p, pbpbSyst, systNames, systLabels, centBinsLabel[cI], plotSaveName);            
      
      canv_p->cd();
      pads_p[0]->cd();           
      
      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(0.035/(1.0 - padSplit));
      label_p->SetTextAlign(31);
      
      for(unsigned int i = 0; i < labels.size(); ++i){
	label_p->DrawLatex(labelX, labelY - 0.083*i, labels[i].c_str());
      }
      delete label_p;
      
      TLegend* leg_p = new TLegend(legX, legY - 0.065*2, legX+0.25, legY);
      leg_p->SetTextFont(42);
      leg_p->SetTextSize(0.035/(1.0 - padSplit));
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
   
      leg_p->AddEntry(pbpbHist_p, ("Pb+Pb #bf{" + centBins2Str + "}").c_str(), "P L");
      leg_p->AddEntry(ppHist_p, "p+p", "P L");      
      leg_p->Draw("SAME");
      
      gPad->SetTicks();
      if(doLogX) pads_p[0]->SetLogx();
      if(doLogY) pads_p[0]->SetLogy();
      
      canv_p->cd();
      pads_p[1]->cd();
      
      std::vector<std::vector<Double_t > > ratSyst = getSystRat(inPbPbFile_p, inPPFile_p, pbpbHist_p, ppHist_p, pbpbSystHistNames, ppSystHistNames, isSystCorrelated, gammaMatchedBins, scaledTotalPbPb, scaledTotalPP);

      pbpbHist_p->Divide(ppHist_p);
      
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
      
      pbpbHist_p->GetYaxis()->SetTitle("Pb+Pb/p+p");
      pbpbHist_p->GetYaxis()->SetTitleOffset(yOffset*padSplit/(1.0 - padSplit));

      pbpbHist_p->SetMinimum(minRatVal);
      pbpbHist_p->SetMaximum(maxRatVal);

      gPad->SetTicks();
      
      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength()*(1.0 - padSplit)/padSplit);
      pbpbHist_p->SetTickLength(pbpbHist_p->GetTickLength("Y")*(1.0 - padSplit)/padSplit, "Y");

      pbpbHist_p->GetXaxis()->SetTitle(varNameLabel.c_str());
      pbpbHist_p->DrawCopy("HIST E1 P");

      drawSyst(pads_p[1], pbpbHist_p, ratSyst);
      pbpbHist_p->DrawCopy("HIST E1 P SAME");      

      plotSaveName = "pdfDir/" + dateStr + "/" + varName + "UnfoldedSyst_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "OverPP_" + saveTag + "_" + dateStr + "." + saveExt;            
      plotSyst(pbpbHist_p, ratSyst, systNames, systLabels, centBinsLabel[cI] + "/pp", plotSaveName);
      
      canv_p->cd();
      pads_p[1]->cd();
      
      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(pbpbHist_p->GetBinLowEdge(1), 1.0, pbpbHist_p->GetBinLowEdge(pbpbHist_p->GetXaxis()->GetNbins()+1), 1.0);
      delete line_p;

      canv_p->cd();
      
      TBox* box_p = new TBox();
      box_p->SetFillColor(0);
      box_p->DrawBox(0.11, 0.435, 0.1575, 0.46);
      delete box_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;
         
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + varName + "Unfolded_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + "." + saveExt);
      delete canv_p;
      delete leg_p;

      if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << gI << "/" << nGammaPtBins << ", " << cI << "/" << centBinsStr.size() << std::endl;

      delete pbpbHist_p;
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

    {
      TCanvas* canv_p = new TCanvas("canv_p", "", width, 900.0);
      canv_p->SetTopMargin(topMargin);
      canv_p->SetRightMargin(rightMargin);
      canv_p->SetLeftMargin(leftMargin);
      canv_p->SetBottomMargin(bottomMargin*padSplit);
      
      HIJet::Style::EquipHistogram(integralPerCentrality_h, 0);
      integralPerCentrality_h->SetMaximum(1.5);
      integralPerCentrality_h->SetMinimum(0.0);    
      
      integralPerCentrality_h->DrawCopy("HIST E1 P");
      
      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(0.035);
      label_p->SetTextAlign(31);
      
      label_p->DrawLatex(labelX, labelY, std::string(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false)).c_str());
      delete label_p;
      
      canv_p->SaveAs(std::string("pdfDir/" + dateStr + "/integralPerCentrality_GammaPt" + gI + "_" + saveTag + "_" + dateStr + "." + saveExt).c_str());
      
      delete canv_p;
    }

    {
      TCanvas* canv_p = new TCanvas("canv_p", "", width, 900.0);
      canv_p->SetTopMargin(topMargin);
      canv_p->SetRightMargin(rightMargin);
      canv_p->SetLeftMargin(leftMargin);
      canv_p->SetBottomMargin(bottomMargin*padSplit);
      
      HIJet::Style::EquipHistogram(meanPerCentrality_h, 0);
      meanPerCentrality_h->SetMaximum(1.5);
      meanPerCentrality_h->SetMinimum(0.0);    
      
      meanPerCentrality_h->DrawCopy("HIST E1 P");
      
      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(42);
      label_p->SetTextSize(0.035);
      label_p->SetTextAlign(31);
      
      label_p->DrawLatex(labelX, labelY, std::string(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false)).c_str());
      delete label_p;
      
      canv_p->SaveAs(std::string("pdfDir/" + dateStr + "/meanPerCentrality_GammaPt" + gI + "_" + saveExt + "_" + dateStr + "." + saveExt).c_str());
    
      delete canv_p;
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
