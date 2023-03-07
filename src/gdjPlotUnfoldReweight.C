//Author: Chris McGinn (2022.08.09)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c and cpp
#include <string>
#include <vector>

//Local
#include "include/binFlattener.h"
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"

void plotReweightHists(TH1D* reweight_p, TH1D* data_p, TH1D* reweighted_p, std::string saveName, std::vector<std::string> labels, bool doLogX = false)
{
  const int titleFont = 42;
  const double titleSize = 0.035;
  const double labelSize = titleSize*0.9;
  const double yOffset = 0.9;

  const Float_t height = 1200.0;
  const Float_t width = 942.831;
  const Double_t leftMargin = 0.16;
  const Double_t bottomAbsFracOfLeft = 0.8;

  const Double_t bottomMarginFracHeight = bottomAbsFracOfLeft*leftMargin*width/height;

  const Double_t topPanelFrac = 0.5;
  const Double_t padSplit = 1.0 - topPanelFrac;

  const Double_t bottomMargin = bottomMarginFracHeight/padSplit;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;
  const Double_t noMargin = 0.001;

  TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
  setMargins(canv_p, noMargin, noMargin, noMargin, noMargin);

  const Int_t nPad = 2;
  TPad* pads_p[2];
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  setMargins(pads_p[0], leftMargin, topMargin, rightMargin, noMargin);
  canv_p->cd();
  pads_p[0]->Draw("SAME");

  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  setMargins(pads_p[1], leftMargin, noMargin, rightMargin, bottomMargin);
  canv_p->cd();
  pads_p[1]->Draw("SAME");

  canv_p->cd();
  pads_p[0]->cd();

  pads_p[0]->SetTicks();

  HIJet::Style::EquipHistogram(reweight_p, 2);
  HIJet::Style::EquipHistogram(data_p, 0);
  HIJet::Style::EquipHistogram(reweighted_p, 1);

  binWidthAndSelfNorm(reweight_p);
  binWidthAndSelfNorm(data_p);
  binWidthAndSelfNorm(reweighted_p);

  Float_t max = TMath::Max(getMax(reweight_p), getMax(data_p));
  max = TMath::Max(getMax(reweighted_p), max);

  Float_t min = TMath::Min(getMinGTZero(reweight_p), getMinGTZero(data_p));
  min = TMath::Min(getMinGTZero(reweighted_p), min);

  reweight_p->SetMaximum(max*1.2);
  reweight_p->SetMinimum(min/1.2);

  reweight_p->GetYaxis()->SetTitleFont(titleFont);
  reweight_p->GetYaxis()->SetTitleSize(titleSize/(1.0-padSplit));
  reweight_p->GetYaxis()->SetLabelFont(titleFont);
  reweight_p->GetYaxis()->SetLabelSize(labelSize/(1.0-padSplit));
  reweight_p->GetYaxis()->SetTitleOffset(yOffset);
  reweight_p->GetYaxis()->SetNdivisions(505);
  
  reweight_p->DrawCopy("HIST E1 P");
  data_p->DrawCopy("HIST E1 P SAME");
  reweighted_p->DrawCopy("HIST E1 P SAME");

  if(doLogX) gPad->SetLogx();
  
  pads_p[1]->cd();

  pads_p[1]->SetTicks();


  reweight_p->Divide(data_p);
  reweighted_p->Divide(data_p);

  reweight_p->SetMaximum(1.6);
  reweight_p->SetMinimum(0.4);

  reweight_p->GetYaxis()->SetTitleFont(titleFont);
  reweight_p->GetYaxis()->SetTitleSize(titleSize/(padSplit));
  reweight_p->GetYaxis()->SetLabelFont(titleFont);
  reweight_p->GetYaxis()->SetLabelSize(labelSize/(padSplit));
  reweight_p->GetYaxis()->SetTitleOffset(yOffset);

  reweight_p->GetXaxis()->SetTitleFont(titleFont);
  reweight_p->GetXaxis()->SetTitleSize(titleSize/(padSplit));
  reweight_p->GetXaxis()->SetLabelFont(titleFont);
  reweight_p->GetXaxis()->SetLabelSize(labelSize/(padSplit));

  reweight_p->GetYaxis()->SetTitle("Ratio MC/Data");
  reweight_p->GetYaxis()->SetNdivisions(505);

  reweight_p->DrawCopy("HIST E1 P");
  reweighted_p->DrawCopy("HIST E1 P SAME");

  if(doLogX) gPad->SetLogx();
  
  TLine* line_p = new TLine();
  line_p->SetLineWidth(2);
  line_p->SetLineStyle(2);

  line_p->DrawLine(reweight_p->GetXaxis()->GetBinLowEdge(1), 1.0, reweight_p->GetXaxis()->GetBinLowEdge(reweight_p->GetXaxis()->GetNbins()+1), 1.0);

  delete line_p;
  pads_p[0]->cd();

  Float_t xLow = 0.5;
  if(saveName.find("DRJJ") != std::string::npos) xLow = 0.3;
  Float_t xHigh = xLow + 0.25;

  TLegend* leg_p = new TLegend(xLow, 0.7, xHigh, 0.9);
  leg_p->SetTextFont(titleFont);
   leg_p->SetTextSize(titleSize/(1.0 - padSplit));
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  leg_p->AddEntry(reweight_p, "MC (Unweighted)", "P L");
  leg_p->AddEntry(reweighted_p, "MC (Reweighted)", "P L");
  leg_p->AddEntry(data_p, "Data", "P L");

  leg_p->Draw("SAME");

  TLatex* label_p = new TLatex();
  xLow = 0.94;
  Float_t yLow = 0.65;

  if(saveName.find("DRJJ") != std::string::npos){
    xLow = 0.3;
    yLow = 0.2;
    initLabel(label_p, titleFont, titleSize/(1.0 - padSplit), 0);
  }
  else initLabel(label_p, titleFont, titleSize/(1.0 - padSplit), 1);

  for(unsigned int lI = 0; lI < labels.size(); ++lI){
    label_p->DrawLatex(xLow, yLow - lI*0.07, labels[lI].c_str());
  }

  gStyle->SetOptStat(0);

  quietSaveAs(canv_p, saveName);

  //Cleanup
  delete leg_p;
  delete label_p;
  for(Int_t pI = 0; pI < nPad; ++pI){delete pads_p[pI]; pads_p[pI] = nullptr;}
  delete canv_p;
  
  return;
}

int gdjPlotUnfoldReweight(std::string inConfigFileName)
{
  //Check config exists
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Create date-stamped and pdf dir
  const std::string dateStr = getDateStr();
  const std::string pdfStr = "pdfDir/" + dateStr;
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);
  
//Global debug handler - do 'export DOGLOBALDEBUGROOT=1' to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  //Grab our input config, define necessary internal params + test they exist
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {"INUNFOLDFILENAME"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  const std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;
  
  TFile* inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");
  TEnv* inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");

  std::vector<std::string> necessaryUnfoldFileParams = {
    "VARNAME",
    "JTPTBINSLOWRECO",
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
    "SUBJTPTBINSDOLOG",
    "SUBJTPTBINSDOCUSTOM",
    "SUBJTGAMMAPTMIN",
    "SUBJTGAMMAPTMAX",
    "ISPP",
    "DOREWEIGHTVAR",
    "SAVETAG"};
  
  std::vector<std::string> necessaryUnfoldFileParamsPbPb = {"CENTBINS"};
  if(!checkEnvForParams(inUnfoldFileConfig_p, necessaryUnfoldFileParams)) return 1;

  
  const Int_t nMaxBins = 200;
  const std::string varName = inUnfoldFileConfig_p->GetValue("VARNAME", "");
  const std::string varNameLower = returnAllLowercaseString(varName);

  const Bool_t isMultijet = isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");

  const Double_t jtPtBinsLowReco = inUnfoldFileConfig_p->GetValue("JTPTBINSLOWRECO", -1.0);
  
  
  //Gamma Pt Bin handling
  const Int_t nGammaPtBins = inUnfoldFileConfig_p->GetValue("NGAMMAPTBINS", -1);
  const Float_t gammaPtBinsLow = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOW", -1);
  const Float_t gammaPtBinsHigh = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGH", -1);
  const Float_t gammaPtBinsLowReco = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSLOWRECO", -1);
  const Float_t gammaPtBinsHighReco = inUnfoldFileConfig_p->GetValue("GAMMAPTBINSHIGHRECO", -1);
  const Bool_t isPP = inUnfoldFileConfig_p->GetValue("ISPP", 0);
  const Bool_t doReweightVar = inUnfoldFileConfig_p->GetValue("DOREWEIGHTVAR", 0);
  const std::string saveTag = inUnfoldFileConfig_p->GetValue("SAVETAG", "");
  Double_t gammaPtBins[nMaxBins+1];
  customBinFiller("GAMMAPT", inUnfoldFileConfig_p, nMaxBins, gammaPtBins);

  //SubJtPt Bin Handling
  const Int_t nSubJtPtBins = inUnfoldFileConfig_p->GetValue("NSUBJTPTBINS", 0);
  const Float_t subJtPtBinsLow = inUnfoldFileConfig_p->GetValue("SUBJTPTBINSLOW", -1.0);
  const Float_t subJtPtBinsHigh = inUnfoldFileConfig_p->GetValue("SUBJTPTBINSHIGH", -1.0);
  const Bool_t subJtPtBinsDoLog = (bool)inUnfoldFileConfig_p->GetValue("SUBJTPTBINSDOLOG", 1);
  const Bool_t subJtPtBinsDoCustom = (bool)inUnfoldFileConfig_p->GetValue("SUBJTPTBINSDOCUSTOM", 1);
  Double_t subJtPtBins[nMaxBins+1];
  customBinFiller("SUBJTPT", inUnfoldFileConfig_p, nMaxBins, subJtPtBins);

  //SubjtGammaPtBin handling
  Double_t subJtGammaPtMin = inUnfoldFileConfig_p->GetValue("SUBJTGAMMAPTMIN", -1000.0);
  Double_t subJtGammaPtMax = inUnfoldFileConfig_p->GetValue("SUBJTGAMMAPTMAX", -1000.0);
  Double_t subJtGammaPtBins[nMaxBins+1];

  //binflattener  
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
  
  
  if(!isPP){
    if(!checkEnvForParams(inUnfoldFileConfig_p, necessaryUnfoldFileParamsPbPb)) return 1;
  }
  
  if(!doReweightVar){
    std::cout << "DOREWEIGHTVAR=" << doReweightVar << "; nothing to plot. return 1" << std::endl;    
    return 1;
  }
  
  std::vector<std::string> centBinsStr = {"PP"};
  std::vector<std::string> centBinsLabel = {"#bf{p+p}"};
  Int_t nCentBins = 1;
  if(!isPP){
    std::vector<std::string> tempCentBinsStr = strToVect(inUnfoldFileConfig_p->GetValue("CENTBINS", ""));
    centBinsStr.clear();
    centBinsLabel.clear();

    nCentBins = tempCentBinsStr.size()-1;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
      centBinsLabel.push_back("#bf{Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%}");
    }
  }

  std::vector<std::string> varNamesForReweight = {varName};
  //You aren't reweighting in these
  //  if(!isStrSame("DRJJ", varName)) varNamesForReweight.push_back("DRJJ");
  //  if(!isStrSame("AJJ", varName)) varNamesForReweight.push_back("AJJ");

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Loop over nCentBins
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TH1D* reweight_p = (TH1D*)inUnfoldFile_p->Get(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_NOMINAL_PURCORR_COMBINED_ForReweight_h").c_str());
    TH1D* reweighted_p = (TH1D*)inUnfoldFile_p->Get(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_NOMINAL_PURCORR_COMBINED_Reweighted_h").c_str());
    TH1D* data_p = (TH1D*)inUnfoldFile_p->Get(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_NOMINAL_PURCORR_COMBINED_h").c_str());

    std::string saveName = pdfStr + "/photonPtReco_Reweight_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + ".png";
    std::vector<std::string> labels = {"#bf{#it{ATLAS Internal}}", centBinsLabel[cI]};

    plotReweightHists(reweight_p, data_p, reweighted_p, saveName, labels, true);
    
    for(unsigned int varI = 0; varI < varNamesForReweight.size(); ++varI){
      if(doGlobalDebug) std::cout << "FILE, LINE, cI/nCentBins, varI/nVar: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << ", " << varI << "/" << varNamesForReweight.size() << std::endl;

      TH2D* reweight2D_p = (TH2D*)inUnfoldFile_p->Get(("photonPtJet" + varNamesForReweight[varI] + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str());
      TH2D* reweighted2D_p = (TH2D*)inUnfoldFile_p->Get(("photonPtJet" + varNamesForReweight[varI] + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str());
      TH2D* data2D_p = (TH2D*)inUnfoldFile_p->Get(("photonPtJet" + varNamesForReweight[varI] + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      
      Int_t nBinsX = reweight2D_p->GetXaxis()->GetNbins();
      const Int_t nMaxBins = 200;
      if(nBinsX > nMaxBins){
	std::cout << "nBinsX: " << nBinsX << ", exceeds max, " << nMaxBins << ". return 1" << std::endl;
	return 1;
      }
      Float_t binsX[nMaxBins];
      
      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	binsX[bIX] = reweight2D_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }
      std::string xTitle = ";" + ((std::string)reweight2D_p->GetXaxis()->GetTitle()) + ";Shape (Bin-width norm.)";
      //Loop over nGammaPtBins
      for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
	Float_t gammaPtBinCenter = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;
	if(gammaPtBinCenter < gammaPtBinsLowReco) continue;
	if(gammaPtBinCenter > gammaPtBinsHighReco) continue;
	
	TH1D* reweight1D_p = new TH1D("reweight1D_h", xTitle.c_str(), nBinsX, binsX);
	TH1D* reweighted1D_p = new TH1D("reweighted1D_h", xTitle.c_str(), nBinsX, binsX);
	TH1D* data1D_p = new TH1D("data1D_h", xTitle.c_str(), nBinsX, binsX);

	std::vector<TH1D*> hists1D_p = {reweight1D_p, reweighted1D_p, data1D_p};
	std::vector<TH2D*> hists2D_p = {reweight2D_p, reweighted2D_p, data2D_p};
	
	for(Int_t bIX = 0; bIX < reweight2D_p->GetXaxis()->GetNbins(); ++bIX){
	  if(isMultijet){
	    for(Int_t bIY = 0; bIY < reweight2D_p->GetYaxis()->GetNbins(); ++bIY){
	      Int_t gammaPtBinPos = subJtGammaPtBinFlattener.GetBin1PosFromGlobal(bIY);
	      if(gammaPtBinPos != gI) continue;

	      int subJtPtBinPos = subJtGammaPtBinFlattener.GetBin2PosFromGlobal(bIY);
	      Float_t subJtPtBinCenter = (subJtPtBins[subJtPtBinPos] + subJtPtBins[subJtPtBinPos+ 1])/2.0;
	      if(subJtPtBinCenter < jtPtBinsLowReco) continue;	      

	      for(unsigned int hI = 0; hI < hists2D_p.size(); ++hI){
		Float_t value = hists1D_p[hI]->GetBinContent(bIX+1);
		Float_t error = hists1D_p[hI]->GetBinError(bIX+1);

		value += hists2D_p[hI]->GetBinContent(bIX+1, bIY+1);
		error = TMath::Sqrt(error*error + hists2D_p[hI]->GetBinError(bIX+1, bIY+1)*hists2D_p[hI]->GetBinError(bIX+1, bIY+1));

		hists1D_p[hI]->SetBinContent(bIX+1, value);
		hists1D_p[hI]->SetBinError(bIX+1, error);
	      }
	    }
	  }
	  else{
	    for(unsigned int hI = 0; hI < hists2D_p.size(); ++hI){
	      hists1D_p[hI]->SetBinContent(bIX+1, hists2D_p[hI]->GetBinContent(bIX+1, gI+1));
	      hists1D_p[hI]->SetBinError(bIX+1, hists2D_p[hI]->GetBinError(bIX+1, gI+1));
	    }	    
	  }
	}
	
	std::vector<std::string> tempLabels = labels;
	tempLabels.push_back(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false));
      
	saveName = pdfStr + "/photonPtRecoJt" + varNamesForReweight[varI] + "_Reweight_" + centBinsStr[cI] + "_GammaPt" + std::to_string(gI) + "_" + saveTag + "_" + dateStr + ".png";
	plotReweightHists(reweight1D_p, data1D_p, reweighted1D_p, saveName, tempLabels);
	
	delete reweight1D_p;
	delete reweighted1D_p;
	delete data1D_p;

	if(doGlobalDebug) std::cout << "FILE, LINE, cI/nCentBins, varI/nVar: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << nCentBins << ", " << varI << "/" << varNamesForReweight.size() << std::endl;
      }
    }
  }      
  
  inUnfoldFile_p->Close();
  delete inUnfoldFile_p;
  
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotUnfoldReweight.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjPlotUnfoldReweight(argv[1]);
  return retVal;
}
