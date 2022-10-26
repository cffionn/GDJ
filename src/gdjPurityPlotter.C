//c and cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

//Local
#include "include/binFlattener.h"
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/varUtil.h"

int gdjPurityPlotter(const std::string inConfigFileName)
{
  //Init file/directory tool + check input file
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1;

  //Init debug tool + test
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;

  //Grab input config file + check parameters
  TEnv* plotConfig_p = new TEnv(inConfigFileName.c_str()); 
  std::vector<std::string> reqPlotConfigParams = {"INFILENAME"};
  if(!checkEnvForParams(plotConfig_p, reqPlotConfigParams)) return 1;

  //Check for global labels
  std::vector<std::string> globalLabels;
  for(unsigned int gI = 0; gI < 10; ++gI){
    std::string globalLabelStr = "GLOBALLABEL." + std::to_string(gI);
    if(checkEnvForParams(plotConfig_p, {globalLabelStr})) globalLabels.push_back(plotConfig_p->GetValue(globalLabelStr.c_str(), ""));
  }

  //Optionally also plot the purity curves
  const std::string purityFileName = plotConfig_p->GetValue("PURITYFILENAME", "");
  bool doPurityFits = check.checkFileExt(purityFileName, ".root");
  TFile* purityFile_p = nullptr;
  if(doPurityFits) purityFile_p = new TFile(purityFileName.c_str(), "READ");
  
  //check infile exists and grab input file
  const std::string inFileName = plotConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");

  std::vector<std::string> reqInConfigParams = {"ISPP",
						"NGAMMAPTBINS",
						"GAMMAPTBINSLOW",
						"GAMMAPTBINSHIGH",
						"GAMMAPTBINSLOWRECO",
						"GAMMAPTBINSHIGHRECO",
						"GAMMAPTBINSDOLOG",
						"GAMMAPTBINSDOCUSTOM",
                                                "NSUBJTPTBINS",
						"SUBJTPTBINSLOW",
						"SUBJTPTBINSHIGH",
						"SUBJTPTBINSLOWRECO",
						"SUBJTPTBINSHIGHRECO",
						"SUBJTPTBINSDOLOG",
						"SUBJTPTBINSDOCUSTOM",
						"SUBJTGAMMAPTMIN",
						"SUBJTGAMMAPTMAX"};
  if(!checkEnvForParams(inConfig_p, reqInConfigParams)) return 1;
  
  //If its Pb+Pb we need some additional params, check and grab
  std::vector<std::string> pbpbParams = {"CENTBINS"};    
  const Bool_t isPP = inConfig_p->GetValue("ISPP", 0);
  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!checkEnvForParams(inConfig_p, pbpbParams)) return 1;

    centBins = strToVect(inConfig_p->GetValue("CENTBINS", ""));
    nCentBins = centBins.size() - 1;
  }

  //need to define a bin cap
  const Int_t nMaxBins = 200;

  //Grab gammaPtBins and 2-D binning of gammaptbins x subjtptbins
  const Int_t nGammaPtBins = inConfig_p->GetValue("NGAMMAPTBINS", -1);
  const Float_t gammaPtBinsLow = inConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Float_t gammaPtBinsLowReco = inConfig_p->GetValue("GAMMAPTBINSLOWRECO", 1000.0);
  const Float_t gammaPtBinsHighReco = inConfig_p->GetValue("GAMMAPTBINSHIGHRECO", -1000.0);
  const Bool_t gammaPtBinsDoLog = (bool)inConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxBins+1];
  std::string gammaPtBinsStr = "";
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    gammaPtBinsStr = inConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "gammaPtBinsDoCustom is true but gammaPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBinsTemp != nGammaPtBins){
      std::cout << "Number of gamma pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  const Int_t nSubJtPtBins = inConfig_p->GetValue("NSUBJTPTBINS", -1);
  const Float_t subJtPtBinsLow = inConfig_p->GetValue("SUBJTPTBINSLOW", 80.0);
  const Float_t subJtPtBinsHigh = inConfig_p->GetValue("SUBJTPTBINSHIGH", 280.0);
  const Float_t subJtPtBinsLowReco = inConfig_p->GetValue("SUBJTPTBINSLOWRECO", 1000.0);
  const Float_t subJtPtBinsHighReco = inConfig_p->GetValue("SUBJTPTBINSHIGHRECO", -1000.0);
  const Bool_t subJtPtBinsDoLog = (bool)inConfig_p->GetValue("SUBJTPTBINSDOLOG", 1);
  const Bool_t subJtPtBinsDoCustom = (bool)inConfig_p->GetValue("SUBJTPTBINSDOCUSTOM", 1);
  Double_t subJtPtBins[nMaxBins+1];
  std::string subJtPtBinsStr = "";

  if(subJtPtBinsDoLog) getLogBins(subJtPtBinsLow, subJtPtBinsHigh, nGammaPtBins, subJtPtBins);
  else if(subJtPtBinsDoCustom){
    subJtPtBinsStr = inConfig_p->GetValue("SUBJTPTBINSCUSTOM", "");
    if(subJtPtBinsStr.size() == 0){
      std::cout << "subJtPtBinsDoCustom is true but subJtPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> subJtPtBinsTemp = strToVectF(subJtPtBinsStr);
    Int_t nSubJtPtBinsTemp = ((Int_t)subJtPtBinsTemp.size()) - 1;
    if(nSubJtPtBinsTemp != nSubJtPtBins){
      std::cout << "Number of subJt pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < subJtPtBinsTemp.size(); ++gI){
      subJtPtBins[gI] = subJtPtBinsTemp[gI];
    }
  }
  else getLinBins(subJtPtBinsLow, subJtPtBinsHigh, nGammaPtBins, subJtPtBins);


  //We need the flattened bin array to correctly construct the mixed event hists
  const Int_t nSubJtGammaPtBins = nGammaPtBins*nSubJtPtBins;
  const Float_t subJtGammaPtMin = inConfig_p->GetValue("SUBJTGAMMAPTMIN", -1.0);
  const Float_t subJtGammaPtMax = inConfig_p->GetValue("SUBJTGAMMAPTMAX", -1.0);
  Double_t subJtGammaPtBins[nMaxBins+1];

  //Now construct the flattened subjt pt-gamma pt bins
  binFlattener subJtGammaPtBinFlattener;
  //Now create your flattened bins using the binFlattener class  
  if(!subJtGammaPtBinFlattener.Init("subJtGammaPtBinFlattener", nGammaPtBins, gammaPtBins, nSubJtPtBins, subJtPtBins)){
    std::cout << "Error gdjNTupleToHist: Failure to init binFlattener, L" << __LINE__ << ". return 1" << std::endl;
    return 1;
  }
  std::vector<double> subJtGammaPtBinsV = subJtGammaPtBinFlattener.GetFlattenedBins(subJtGammaPtMin, subJtGammaPtMax);
  for(unsigned int sgI = 0; sgI < subJtGammaPtBinsV.size(); ++sgI){
    subJtGammaPtBins[sgI] = subJtGammaPtBinsV[sgI];
  }

  //check date-tagged output directory exist and if not create
  const std::string dateStr = getDateStr();
  std::string pdfStr = "pdfDir/" + dateStr;
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);

  std::vector<std::string> observables = {"pt", "xj", "dphi", "xjj", "drjj", "dphijjg", "dphijj", "ajj"};

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

  
  for(int cI = 0; cI < nCentBins; ++cI){
    inFile_p->cd();

    std::string centStr = centBins[cI];
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
   
    for(unsigned int oI = 0; oI < observables.size(); ++oI){
      std::string histVarName = varNameToHistName(observables[oI]);
      std::string varNameUpper = returnAllCapsString(observables[oI]);

      //Get Plotting Params
      bool doLogX = false;
      bool doLogY = false;
      Float_t legX = 0.7;
      Float_t legY = 0.9;
      Float_t labelX = 0.2;
      Float_t labelY = 0.4;
      Bool_t labelDoLeftAlign = false;
      
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "DOLOGX"})) doLogX = plotConfig_p->GetValue((varNameUpper + "DOLOGX").c_str(), 0);
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "DOLOGY"})) doLogY = plotConfig_p->GetValue((varNameUpper + "DOLOGY").c_str(), 0);
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "LEGX"})) legX = plotConfig_p->GetValue((varNameUpper + "LEGX").c_str(), 0.0);
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "LEGY"})) legY = plotConfig_p->GetValue((varNameUpper + "LEGY").c_str(), 0.0);
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "LABELX"})) labelX = plotConfig_p->GetValue((varNameUpper + "LABELX").c_str(), 0.0);
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "LABELY"})) labelY = plotConfig_p->GetValue((varNameUpper + "LABELY").c_str(), 0.0);      
      if(checkEnvForParams(plotConfig_p, {varNameUpper + "DOLABELLEFTALIGN"})) labelDoLeftAlign = plotConfig_p->GetValue((varNameUpper + "DOLABELLEFTALIGN").c_str(), 0);      

      std::cout << "ALIGNMENT: " << labelDoLeftAlign << std::endl;
      
      bool isMultijet = varNameToLabelIsMultijet(observables[oI]);
      std::string mixModeStr = "MIXMODE0";
      if(!isPP){
	if(isMultijet) mixModeStr = "MIXMODE2";
	else mixModeStr = "MIXMODE1";
      }
      
      std::string histNameSub = centStr + "/photonPtJt" + histVarName + "VCent_" + centStr + "_BarrelAndEC_NOMINAL_DPhi0_" + mixModeStr + "_SUB_h"; 
      std::string histNameSideband = histNameSub;
      std::string histNamePurCorr = histNameSub;
      histNameSideband.replace(histNameSideband.find(mixModeStr), mixModeStr.size(), "Sideband_" + mixModeStr);
      histNamePurCorr.replace(histNamePurCorr.find(mixModeStr), (mixModeStr + "_SUB").size(), "PURCORR");
      std::string histNameSidebandScaled = histNamePurCorr;
      histNameSidebandScaled.replace(histNameSidebandScaled.find("PURCORR"), std::string("PURCORR").size(), "PURCORRBkgd");

      TH2D* subHist_p = (TH2D*)inFile_p->Get(histNameSub.c_str());
      TH2D* sidebandHist_p = (TH2D*)inFile_p->Get(histNameSideband.c_str());
      TH2D* sidebandScaledHist_p = (TH2D*)inFile_p->Get(histNameSidebandScaled.c_str());
      TH2D* purCorrHist_p = (TH2D*)inFile_p->Get(histNamePurCorr.c_str());

      std::cout << "HistNamePurCorr: " << histNamePurCorr << std::endl;
      
      Int_t nBinsX = subHist_p->GetXaxis()->GetNbins();
      Double_t binsX[nMaxBins+1];
      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	binsX[bIX] = subHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }

      std::cout << "BINSX: " << std::endl;
      for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
	std::cout << " " << binsX[bIX] << ", ";	
      }
      std::cout << binsX[nBinsX] << std::endl;

      for(Int_t bIY = 0; bIY < subHist_p->GetYaxis()->GetNbins(); ++bIY){
	Int_t gammaPos = bIY;
	if(isMultijet){
	  Int_t subJtPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
	  Float_t subJtPtBinCenter = (subJtPtBins[subJtPos] + subJtPtBins[subJtPos+1])/2.0;
	  if(subJtPtBinCenter < subJtPtBinsLowReco) continue;

	  gammaPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
	}

	Float_t gammaPtBinCenter = (gammaPtBins[gammaPos] + gammaPtBins[gammaPos+1])/2.0;
	if(gammaPtBinCenter < gammaPtBinsLowReco) continue;
	if(gammaPtBinCenter > gammaPtBinsHighReco) continue;
	
	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;

	std::string titleStr = subHist_p->GetXaxis()->GetTitle();
	titleStr = ";" + titleStr + ";Counts (#color[2]{Sideband} is scaled)";

	TH1D* subHist1D_p = new TH1D("subHist1D_p", titleStr.c_str(), nBinsX, binsX);
	TH1D* sidebandHist1D_p = new TH1D("sidebandHist1D_p", titleStr.c_str(), nBinsX, binsX);
	TH1D* sidebandScaledHist1D_p = new TH1D("sidebandScaledHist1D_p", titleStr.c_str(), nBinsX, binsX);
	TH1D* purCorrHist1D_p = new TH1D("purCorrHist1D_p", titleStr.c_str(), nBinsX, binsX);

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;
	fineTH2ToCoarseTH1(subHist_p, subHist1D_p, bIY);
	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;
	fineTH2ToCoarseTH1(sidebandHist_p, sidebandHist1D_p, bIY);
	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;
	fineTH2ToCoarseTH1(sidebandScaledHist_p, sidebandScaledHist1D_p, bIY);
	fineTH2ToCoarseTH1(purCorrHist_p, purCorrHist1D_p, bIY);
	
	TCanvas* canv_p = new TCanvas("canv_p", "", 900, 900);
	setMargins(canv_p, leftMargin, 0.03, 0.03, bottomMargin);

	if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", L" << __LINE__ << std::endl;
	
	HIJet::Style::EquipHistogram(subHist1D_p, 2);
	HIJet::Style::EquipHistogram(sidebandScaledHist1D_p, 1); 
	HIJet::Style::EquipHistogram(purCorrHist1D_p, 0); 
	
	Double_t maxVal = TMath::Max(subHist1D_p->GetMaximum(), sidebandScaledHist1D_p->GetMaximum());
	maxVal = TMath::Max(maxVal, purCorrHist1D_p->GetMaximum());
	
	subHist1D_p->SetMaximum(maxVal*1.1);
	if(!doLogY) subHist1D_p->SetMinimum(0.0);
	else{
	  Double_t minVal = TMath::Min(getMinGTZero(subHist1D_p), getMinGTZero(sidebandScaledHist1D_p));
	  minVal = TMath::Min(minVal, getMinGTZero(purCorrHist1D_p));

	  subHist1D_p->SetMinimum(minVal/2.0);
	  subHist1D_p->SetMaximum(maxVal*2.0);
	}

	
	subHist1D_p->DrawCopy("HIST E1 P");
	sidebandScaledHist1D_p->DrawCopy("HIST E1 P SAME");
	purCorrHist1D_p->DrawCopy("HIST E1 P SAME");

	if(doLogX) gPad->SetLogx();
	if(doLogY) gPad->SetLogy();

	gPad->SetTicks();
	
	TLegend* leg_p = new TLegend(legX, TMath::Max(legY - 0.15, 0.01), TMath::Min(legX+0.2,0.99), legY);
	leg_p->SetBorderSize(0);
	leg_p->SetFillColor(0);
	leg_p->SetFillStyle(0);
	leg_p->SetTextFont(titleFont);
	leg_p->SetTextSize(titleSize);

	leg_p->AddEntry(subHist1D_p, "RAW", "P L");
	leg_p->AddEntry(sidebandScaledHist1D_p, "SIDEBAND", "P L");
	leg_p->AddEntry(purCorrHist1D_p, "CORRECTED", "P L");

	leg_p->Draw("SAME");

	std::vector<std::string> labels = globalLabels;
	if(isPP) labels.push_back("pp");
	else labels.push_back("Pb+Pb, " + centBins[cI] + "-" + centBins[cI+1] + "%");
	labels.push_back(prettyString(gammaPtBins[gammaPos],1,false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gammaPos+1],1,false));	
	
	TLatex* label_p = new TLatex();
	label_p->SetNDC();
	label_p->SetTextFont(titleFont);
	label_p->SetTextSize(titleSize);
	if(labelDoLeftAlign) label_p->SetTextAlign(31);
	
	for(unsigned int lI = 0; lI < labels.size(); ++lI){
	  label_p->DrawLatex(labelX, labelY - lI*0.05, labels[lI].c_str());
	}

	gStyle->SetOptStat(0);
	
	std::string saveName = pdfStr + "/purityJt" + histVarName + "_" + centStr + "_GammaPt" + std::to_string(gammaPos) + "_" + dateStr + ".png";
	quietSaveAs(canv_p, saveName);
	delete canv_p;

	delete leg_p;
	delete label_p;
	
	delete subHist1D_p;
	delete sidebandHist1D_p;
	delete sidebandScaledHist1D_p;
	delete purCorrHist1D_p;
      }
      //Insert a temp break for testing, only produce one set of observables plots      
      //      break;      
    }
  }  

  inFile_p->Close();
  delete inFile_p;


  if(doPurityFits){
    purityFile_p->cd();

    std::vector<TF1*> fits_p;
    
    TCanvas* canv_p = new TCanvas("canv_p", "", 900, 900);
    setMargins(canv_p, leftMargin, 0.03, 0.03, bottomMargin);
   
    TLatex* label_p = new TLatex();
    label_p->SetNDC();
    label_p->SetTextFont(titleFont);
    label_p->SetTextSize(titleSize);

    TLegend* leg_p = new TLegend(0.6,0.3,0.8,0.5);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    leg_p->SetTextFont(titleFont);
    leg_p->SetTextSize(titleSize);
	    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = centBins[cI];
      std::string centLabel = centBins[cI];
      if(!isPP){
	centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
	centLabel = centBins[cI] + "-" + centBins[cI+1] + "%";
      }

      fits_p.push_back((TF1*)purityFile_p->Get(("fit_purity_" + centStr + "_Eta0p00to2p37").c_str()));
      HIJet::Style::EquipTF1(fits_p[cI], nCentBins - cI);

      fits_p[cI]->SetMaximum(1.05);
      fits_p[cI]->SetMinimum(0.45);
      fits_p[cI]->SetTitle("");
      /*
      fits_p[cI]->GetXaxis()->SetTitleFont(titleFont);
      fits_p[cI]->GetYaxis()->SetTitleFont(titleFont);
      fits_p[cI]->GetXaxis()->SetTitleSize(titleSize);
      fits_p[cI]->GetYaxis()->SetTitleSize(titleSize);
      fits_p[cI]->GetXaxis()->SetTitleColor(1);
      fits_p[cI]->GetYaxis()->SetTitleColor(1);
      */
      fits_p[cI]->GetXaxis()->SetTitle("Photon p_{T} [GeV]");
      fits_p[cI]->GetYaxis()->SetTitle("Purity");
      
      fits_p[cI]->GetXaxis()->SetTitleOffset(2.0);
      fits_p[cI]->GetYaxis()->SetTitleOffset(1.0);
      
      if(cI == 0) fits_p[cI]->Draw("L");
      else fits_p[cI]->Draw("L SAME");

      leg_p->AddEntry(fits_p[cI], centLabel.c_str(), "L");
    }
    
    gPad->SetLogx();
    
    for(unsigned int lI = 0; lI < globalLabels.size(); ++lI){
      label_p->DrawLatex(0.6, 0.6 - lI*0.05, globalLabels[lI].c_str());
    }

    leg_p->Draw("SAME");

    TLine* line_p = new TLine();
    line_p->SetLineStyle(2);
    line_p->DrawLine(fits_p[0]->GetXmin(), 1.0, fits_p[0]->GetXmax(), 1.0);
    delete line_p;
    
    gStyle->SetOptStat(0);
    gPad->SetTicks();

    std::string ppStr = "PP";
    if(!isPP) ppStr = "PbPb";
    
    std::string saveName = pdfStr + "/purityFit_" + ppStr + "_"  + dateStr + ".png";
    quietSaveAs(canv_p, saveName);
    delete canv_p;

    delete leg_p;
  }

  
  if(doPurityFits){
    purityFile_p->Close();
    delete purityFile_p;
  }
  
  delete plotConfig_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPurityPlotter.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjPurityPlotter(argv[1]);
  return retVal;
}
