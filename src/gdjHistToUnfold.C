//Author: Chris McGinn (2021.01.11)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TStyle.h"

//RooUnfold
//https://gitlab.cern.ch/RooUnfold/RooUnfold
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
//#include "include/configParser.h"
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
#include "include/treeUtil.h"


void makeReweightHist(TH1* weightHist_p, TH1* numHist_p, TH1* denomHist_p, std::string saveName)
{
  Double_t zeroThresh = TMath::Power(10, -100);
  
  std::string className = numHist_p->ClassName();
  bool isTH1 = false;
  if(isStrSame(className, "TH1D") || isStrSame(className, "TH1F")) isTH1 = true;
  else if(isStrSame(className, "TH2D") || isStrSame(className, "TH2F")) isTH1 = false;
  else{
    std::cout << "MAKEREWEIGHTHIST ERROR: Class \'" << className << "\' is not recognized. Please fix. return" << std::endl;
    return;
  }
  
  std::vector<double> weightVals;
  Int_t nYBins = 1;
  if(!isTH1) nYBins = numHist_p->GetYaxis()->GetNbins();
  
  for(Int_t bIX = 0; bIX < numHist_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < nYBins; ++bIY){

      Float_t ratioVal = 0.0;
      Float_t numeratorVal = 0.0;
      Float_t denominatorVal = 0.0;
      
      if(isTH1){
	numeratorVal = numHist_p->GetBinContent(bIX+1);
	denominatorVal = denomHist_p->GetBinContent(bIX+1);
      }
      else{
	numeratorVal = numHist_p->GetBinContent(bIX+1, bIY+1);
	denominatorVal = denomHist_p->GetBinContent(bIX+1, bIY+1);	
      }

      if(denominatorVal > zeroThresh){
	ratioVal = numeratorVal/denominatorVal;
	weightVals.push_back(ratioVal);
      }

      if(isTH1){
	weightHist_p->SetBinContent(bIX+1, ratioVal);
	weightHist_p->SetBinError(bIX+1, 0.0);
      }
      else{
	weightHist_p->SetBinContent(bIX+1, bIY+1, ratioVal);
	weightHist_p->SetBinError(bIX+1, bIY+1, 0.0);
      }
    }
  }

  std::sort(std::begin(weightVals), std::end(weightVals));
  Double_t scaleVal = 0.0;
  if(weightVals.size()%2 == 0) scaleVal = weightVals[(weightVals.size()-1)/2];
  else scaleVal = (weightVals[weightVals.size()/2] + weightVals[weightVals.size()/2 - 1])/2.0;
  
  weightHist_p->Scale(1.0/scaleVal);

  Int_t minXFill = -1;
  Int_t maxXFill = 100000;
  Int_t minYFill = -1;
  Int_t maxYFill = 100000;

  for(Int_t bIX = 0; bIX < weightHist_p->GetXaxis()->GetNbins(); ++bIX){
    bool allBinsZero = true;

    if(isTH1){
      if(weightHist_p->GetBinContent(bIX+1) > zeroThresh) allBinsZero = false;
    }
    else{
      for(Int_t bIY = 0; bIY < nYBins; ++bIY){
	if(weightHist_p->GetBinContent(bIX+1, bIY+1) > zeroThresh){allBinsZero = false; break;}
      }
    }

    if(!allBinsZero){minXFill = bIX; break;}
  }

  for(Int_t bIX = weightHist_p->GetXaxis()->GetNbins()-1; bIX >= 0; --bIX){
    bool allBinsZero = true;

    if(isTH1){
      if(weightHist_p->GetBinContent(bIX+1) > zeroThresh) allBinsZero = false;
    }
    else{
      for(Int_t bIY = 0; bIY < nYBins; ++bIY){
	if(weightHist_p->GetBinContent(bIX+1, bIY+1) > zeroThresh){allBinsZero = false; break;}
      }
    }

    if(!allBinsZero){maxXFill = bIX; break;}
  }

  if(!isTH1){
    for(Int_t bIY = 0; bIY < nYBins; ++bIY){
      bool allBinsZero = true;

      for(Int_t bIX = 0; bIX < weightHist_p->GetXaxis()->GetNbins(); ++bIX){
	if(weightHist_p->GetBinContent(bIX+1, bIY+1) > zeroThresh){allBinsZero = false; break;}
      }      
      
      if(!allBinsZero){minYFill = bIY; break;}
    }
    
    for(Int_t bIY = nYBins-1; bIY >= 0; --bIY){
      bool allBinsZero = true;
      
      for(Int_t bIX = 0; bIX < weightHist_p->GetXaxis()->GetNbins(); ++bIX){
	if(weightHist_p->GetBinContent(bIX+1, bIY+1) > zeroThresh){allBinsZero = false; break;}
      }    
      
      if(!allBinsZero){maxYFill = bIY; break;}
    }    
  }

  //Having found the range of fills on the weighting hist - need to reweight globally
  for(Int_t bIX = 0; bIX < weightHist_p->GetXaxis()->GetNbins(); ++bIX){
    if(bIX >= minXFill && bIX <= maxXFill) continue;

    for(Int_t bIY = 0; bIY < nYBins; ++bIY){                                                  
      if(bIY < minYFill || bIY > maxYFill) continue;

      if(bIX < minXFill) weightHist_p->SetBinContent(bIX+1, bIY+1, weightHist_p->GetBinContent(minXFill+1, bIY+1));
      else if(bIX > maxXFill) weightHist_p->SetBinContent(bIX+1, bIY+1, weightHist_p->GetBinContent(maxXFill+1, bIY+1));                                                         
    }                                                                                               
  }                                                                                                 
                                                                                                       
  for(Int_t bIY = 0; bIY < nYBins; ++bIY){                                                    
    if(bIY >= minYFill && bIY <= maxYFill) continue;
    
    for(Int_t bIX = 0; bIX < weightHist_p->GetXaxis()->GetNbins(); ++bIX){
      if(bIY < minYFill) weightHist_p->SetBinContent(bIX+1, bIY+1, weightHist_p->GetBinContent(bIX+1, minYFill+1));                                                              
      else if(bIY > maxYFill) weightHist_p->SetBinContent(bIX+1, bIY+1, weightHist_p->GetBinContent(bIX+1, maxYFill+1));                                                         
    }                                                                                               
  }   

  //Dump the reweighting matrix to a canvas
  const int titleFont = 42;
  const double titleSize = 0.03;
  const double labelSize = titleSize*0.9;
  const double yOffset = 0.9;

  const Float_t height = 1200.0;
  const Float_t width = 942.831;
  const Double_t leftMargin = 0.16;
  const Double_t bottomAbsFracOfLeft = 0.8;

  const Double_t bottomMarginFracHeight = bottomAbsFracOfLeft*leftMargin*width/height;

  const Double_t bottomMargin = bottomMarginFracHeight;
  const Double_t topMargin = 0.02;
  const Double_t rightMargin = 0.02;

  TCanvas* canv_p = new TCanvas("canv_p", "", width, height);
  setMargins(canv_p, leftMargin, topMargin, rightMargin, bottomMargin);

  if(isTH1){
    HIJet::Style::EquipHistogram(weightHist_p, 0);
    weightHist_p->DrawCopy("HIST E1 P");
  }
  else{
    canv_p->SetRightMargin(leftMargin);
    weightHist_p->DrawCopy("COLZ");
  }

  gStyle->SetOptStat(0);
  
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  
  return;
}

//Simple plotting function for checking the results of the reweighting
void plotReweightHists(TH1D* reweight_p, TH1D* data_p, TH1D* reweighted_p, std::string saveName, std::vector<std::string> labels)
{
  const int titleFont = 42;
  const double titleSize = 0.04;
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
  
  reweight_p->DrawCopy("HIST E1 P");
  data_p->DrawCopy("HIST E1 P SAME");
  reweighted_p->DrawCopy("HIST E1 P SAME");

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
  reweight_p->DrawCopy("HIST E1 P");
  reweighted_p->DrawCopy("HIST E1 P SAME");

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

Double_t getVar(std::string varName, TLorentzVector jet, TLorentzVector jet2, TLorentzVector photon)
{
  Double_t jetVar = -1000.0;
  if(isStrSame(varName, "dphi")) jetVar = TMath::Abs(getDPHI(jet.Phi(), photon.Phi()));
  else if(isStrSame(varName, "pt")) jetVar = jet.Pt();
  else if(isStrSame(varName, "xj")) jetVar = jet.Pt()/photon.Pt();
  else if(isStrSame(varName, "dphijj")) jetVar = TMath::Abs(getDPHI(jet.Phi(), jet2.Phi()));
  else if(isStrSame(varName, "drjj")) jetVar = getDR(jet.Eta(), jet.Phi(), jet2.Eta(), jet2.Phi());
  else if(isStrSame(varName, "ajj")) jetVar = (jet.Pt() - jet2.Pt())/photon.Pt();
  else{
    jet2 += jet;
    if(isStrSame(varName, "xjj")) jetVar = jet2.Pt()/photon.Pt();
    else if(isStrSame(varName, "dphijjg")) jetVar = TMath::Abs(getDPHI(jet2.Phi(), photon.Pt()));
    else{
      std::cout << "GETVAR ERROR: varName \'" << varName << "\' is not found. return -1000.0" << std::endl;
    }
  }
  return jetVar;
}

//Primary function called in main (should be the only one - support functions go elsewhere)
int gdjHistToUnfold(std::string inConfigFileName)
{
  //Watch just for timing tests if things look low
  cppWatch globalTimer;
  globalTimer.start();

  //Random number generator seed is fixed but determined randomly for debuggable results
  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101             
  TRandom3* randGen_p = new TRandom3(randSeed);

  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  const std::string pdfStr = "pdfDir/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);
  
  //Debugging is set as an environment variable - we use a simple class to pick up
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  //Global debug default spits out file and line as follows
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INRESPONSEFILENAME",
					      "INUNFOLDFILENAME",
					      "VARNAME",
                                              "OUTFILENAME",
					      "NITER",
					      "SAVETAG",
					      "DOREWEIGHTVAR",
					      "DOREWEIGHTDR",
					      "UNFOLDERRTYPE",
					      "NTOYS",
					      "DOREBIN",
					      "UNFOLD"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string varName = config_p->GetValue("VARNAME", "");
  std::string varNameLower = returnAllLowercaseString(varName);

  const std::string saveTag = config_p->GetValue("SAVETAG", "");
  
  //Multijet variables are a well defined list and we must handle differently, so define a bool
  const bool isMultijet = isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");  
  
  const int nIter = config_p->GetValue("NITER", -1);
  const bool doReweightVar = config_p->GetValue("DOREWEIGHTVAR", 0);
  //second conditional to avoid double weighting
  const bool doReweightDR = config_p->GetValue("DOREWEIGHTDR", 0) && !isStrSame("drjj", varNameLower);

  const std::string inUnfoldName = config_p->GetValue("UNFOLD", "");
  
  
  //Rebinning requires additional params i.e. what is the rebin in X and Y
  const bool doRebin = (bool)config_p->GetValue("DOREBIN", 0);
  std::vector<std::string> rebinParams = {"REBINX",
					  "REBINY"};

  //Since arrays require predefined sizes, just define a large cap on them lie 1000 - we will limit ourselves to the subset of the array we need when we get to it
  const int nMaxBins = 1000;
  Int_t nRebinX, nRebinY;
  Double_t rebinX[nMaxBins+1], rebinY[nMaxBins+1];  
  std::vector<float> rebinXVect, rebinYVect;

  //We need some param delta defining when two bin edges are the same - here I pick something small
  const Float_t deltaValue = 0.001;

  //If we are doing rebin, grab the new arrays of bins 
  if(doRebin){
    if(!checkEnvForParams(config_p, rebinParams)) return 1;

    std::string rebinXStr = config_p->GetValue("REBINX", "");
    std::string rebinYStr = config_p->GetValue("REBINY", "");

    //Currently we require both rebin arrays specified even if we only wanna rebin one dimension
    if(rebinXStr.size() == 0){
      std::cout << "ERROR: 'DOREBIN' is " << doRebin << ", but rebinXStr is empty. return 1" << std::endl;
      return 1;
    }
    if(rebinYStr.size() == 0){
      std::cout << "ERROR: 'DOREBIN' is " << doRebin << ", but rebinYStr is empty. return 1" << std::endl;
      return 1;
    }

    //Convert the string to a vector of floats
    rebinXVect = strToVectF(rebinXStr);
    rebinYVect = strToVectF(rebinYStr);

    //Since we are making a choice w/ delta for bin uniqueness we need to check there is nothing below that level
    //Pre-check the bins requested are unique at the level of 2xDeltaValue
    if(!checkMinBinWidth(rebinXVect, 2*deltaValue)) return 1;
    if(!checkMinBinWidth(rebinYVect, 2*deltaValue)) return 1;

    //Finally, fill our final bin arrays
    nRebinX = rebinXVect.size() - 1;
    for(Int_t bIX = 0; bIX < nRebinX+1; ++bIX){
      rebinX[bIX] = rebinXVect[bIX];
    }
    nRebinY = rebinYVect.size() - 1;
    for(Int_t bIY = 0; bIY < nRebinY+1; ++bIY){
      rebinY[bIY] = rebinYVect[bIY];
    }   
  }

  //RooUnfold allows for many different methods for error calculation
  //My preference is toys - simple to understand, precision well defined by number of toys
  //Downside: Not exact, significant CPU time increase 
  const int unfoldErrType = config_p->GetValue("UNFOLDERRTYPE", 0);
  const int nToys = config_p->GetValue("NTOYS", 0);
  if(unfoldErrType < 0 || unfoldErrType > 1){//Valid error types are 0: default kErrors bin-by-bin errors (diagonal covariance matrix), 1: toys
    std::cout << "Given UNFOLDERRTYPE, " << unfoldErrType << ", is invalid. Please give valid type, return 1" << std::endl;    
    return 1;
  }
  if(nToys <= 0 && unfoldErrType == 1){
    std::cout << "NTOYS given \'" << nToys << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  //Define a cap on the number of iterations you can do to not waste time on absurd high amounts
  const int nIterMax = 100;
  if(nIter > nIterMax){
    std::cout << "Requested number of iterations \'" << nIter << "\' exceeds maximum allowed \'" << nIterMax << "\'" << std::endl;
  }

  //Modify the output file name for missing bits (.root, output directory, etc)
  if(outFileName.find(".root") != std::string::npos){
    outFileName.replace(outFileName.rfind(".root"), 5, "");
    outFileName = outFileName + "_" + saveTag + "_" + dateStr + ".root";
  }
  if(outFileName.find("/") == std::string::npos){
    outFileName = outputStr + "/" + outFileName;
  }
    
  //Define some unfolding variables in arrays along with paired-by-position stylized versions
  //Previously i rebuilt everything from a stripped down unfolding TTree - ive re-worked this to be histogram based and temporarily comment out the *Tree variabls  
  std::vector<std::string> validUnfoldVar = {"Pt", "XJ", "DPhi", "XJJ", "AJJ", "DPHIJJ", "DPHIJJG", "DRJJ"};
  //  std::vector<std::string> validUnfoldVarTree = {"XJ", "XJ", "XJJ"};
  std::vector<std::string> validUnfoldVarStyle = {"Jet p_{T}", "x_{J}", "#Delta#phi_{J#gamma}", "#vec{x}_{JJ#gamma}", "A_{JJ}", "#Delta#phi_{JJ}", "#Delta#phi_{JJ#gamma}", "#DeltaR_{JJ}"};

  //Check that our requested unfolding Var matches what is in the vectors + grab position
  const int vectPos = vectContainsStrPos(varName, &validUnfoldVar);
  if(vectPos < 0){
    std::cout << "GIVEN VARNAME \'" << varName << "\' is invalid. Please pick from:" << std::endl;
    for(unsigned int i = 0; i < validUnfoldVar.size(); ++i){
      std::cout << " \'" << validUnfoldVar[i] << "\'" << std::endl;
    }
    std::cout << "return 1" << std::endl;
    return 1;
  }
  //  std::string varNameTree = validUnfoldVarTree[vectPos];
  std::string varNameStyle = validUnfoldVarStyle[vectPos];
  
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;

  //We  need to define a cap here too - given we are extremely unlikely to ever drop below 10% level 10 is sufficient
  const Int_t nMaxCentBins = 10;

  //Grab the response file - this will also be the unfold file while working with MC - with data it will differ
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");  
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  
  //declare as nullptr; actual file / configs will depend on if its data or MC
  TFile* inUnfoldFile_p = nullptr;
  TEnv* inUnfoldFileConfig_p = nullptr;
  TEnv* inUnfoldFileLabels_p = nullptr;

  //If MC - same file; if data - different files
  const bool isSameInputFile = isStrSame(inResponseFileName, inUnfoldFileName);

  //Either have two var point to the same file and handle it that way or separate files (in case of data)
  if(isSameInputFile){
    inUnfoldFile_p = inResponseFile_p;
    inUnfoldFileConfig_p = inResponseFileConfig_p;
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }
  else{
    inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");  
    inUnfoldFileConfig_p = (TEnv*)inUnfoldFile_p->Get("config");    
    inUnfoldFileLabels_p = (TEnv*)inUnfoldFile_p->Get("label");
  }

  //there is an additional set of parameters that we must make sure are defined from the file we are passed 
  std::vector<std::string> necessaryInFileParams = {"ISMC",
						    "DOUNIFIEDPURITY",
						    "NGAMMAPTBINS",
						    "GAMMAPTBINSLOW",
						    "GAMMAPTBINSHIGH",
						    "GAMMAPTBINSDOLOG",
						    "GAMMAPTBINSDOCUSTOM",
						    "NJTPTBINS",
						    "JTPTBINSLOW",
						    "JTPTBINSHIGH",
						    "JTPTBINSDOLOG",
						    "JTPTBINSDOCUSTOM",
						    "JTETABINSLOW",
						    "JTETABINSHIGH",
						    "GAMMAJTDPHI",
						    "ISPP"};

 
  std::vector<std::string> pbpbParams = {"CENTBINS"};
  //  std::vector<std::string> mcParams = {"ASSOCGENMINPT"};

  //Some of these params should be identical between data and MC else the unfold wont make any sense
  std::vector<std::string> paramsToCompare = {"DOUNIFIEDPURITY",
					       "JETR",
					      "ISPP",
					      "CENTBINS",
					      "GAMMAEXCLUSIONDR",
					      "MIXJETEXCLUSIONDR",
					      "NGAMMAPTBINS",
					      "GAMMAPTBINSLOW",
					      "GAMMAPTBINSHIGH",
					      "GAMMAPTBINSLOWRECO",
					      "GAMMAPTBINSHIGHRECO",
					      "GAMMAPTBINSDOLOG",
					      "NJTPTBINS",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "JTPTBINSDOLOG",
					      "GAMMAJTDPHI",
					      "NXJBINS",
					      "XJBINSLOW",
					      "XJBINSHIGH",
					      "XJBINSDOCUSTOM",
					      "XJBINSCUSTOM",
					      "NAJBINS",
					      "AJBINSLOW",
					      "AJBINSHIGH",
					      "AJBINSDOCUSTOM",
					      "AJBINSCUSTOM",
					      "NDRBINS",
					      "DRBINSLOW",
					      "DRBINSHIGH",
					      "DRBINSDOCUSTOM",
					      "DRBINSCUSTOM",
					      "NDPHIBINS",
					      "DPHIBINSLOW",
					      "DPHIBINSHIGH",
					      "DPHIBINSDOCUSTOM",
					      "DPHIBINSCUSTOM"};

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //These are almost always true - throw them in the names and use the label map later to make sure we get the cuts on the plot
  //DPhi0 - the dphi cut
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string gammaJtDPhiStr = "DPhi0";
  //Multijtcut - the minimum pT cut on the jets to construct the multijet object
  const std::string multiJtCutGlobalStr = "MultiJt0";
  //dont remember
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string jtPtBinsGlobalStr = "GlobalJtPt0";
 
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Return failed if we are missing required parameters
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  const bool isMCResponse = inResponseFileConfig_p->GetValue("ISMC", 1);
  //make sure response file given is MC
  if(!isMCResponse){
    std::cout << "Given responseFile \'" << inResponseFileName << "\' is not an MC file. return 1" << std::endl;
    return 1;
  }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Grab first set of parameters defined in necessaryInFileParams
  const bool isMC = inUnfoldFileConfig_p->GetValue("ISMC", 1);
  const bool isPP = inResponseFileConfig_p->GetValue("ISPP", 0);
  const bool doUnifiedPurity = inResponseFileConfig_p->GetValue("DOUNIFIEDPURITY", 0);

  //Unified purity refers to having purity defined barrel and endcap separately or together
  //Originally we were using a separately defined purity but have switched to a barrel + endcap combined definition
  std::cout << "UNIFIED PURITY?" << std::endl;
  std::cout << " " << doUnifiedPurity << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //MixMode - 0 is p+p (no mixing), 1 is normal inclusive jet mixing, 2 is multi-jet mixing
  std::string mixModeStr = "MIXMODE0";
  if(!isPP){
    if(isMultijet) mixModeStr = "MIXMODE2";
    else mixModeStr = "MIXMODE1";
  }

  //Grab the number of jet systematics
  const Int_t nMaxSyst = 100;
  Int_t nJESJER = inResponseFileConfig_p->GetValue("NJETSYSANDNOM", -1);
  std::vector<std::string> jesJERStrVect = strToVect(inResponseFileConfig_p->GetValue("JETSYSANDNOM", ""));

  //Declare set of unfolding syst
  const Int_t nUnfSyst = 2;
  std::vector<std::string> unfSystStrVect = {"PRIOR", "JTPTCUT"};

  if(!doReweightVar){
    std::cout << "WARNING: doReweightVar=false; Unfolding error associated w/ prior variation will be zero by construction!" << std::endl;    
  } 

  //Create systematics vector combining all relevant systematics
  //Note that nJESJER includes nominal
  Int_t nSyst = nJESJER + nUnfSyst;
  std::vector<std::string> systStrVect = jesJERStrVect;
  systStrVect.insert(systStrVect.end(), unfSystStrVect.begin(), unfSystStrVect.end());
  /*
  std::cout << "Nominal + Syst. List: " << std::endl;
  for(unsigned int sI = 0; sI < systStrVect.size(); ++sI){
    std::cout << " " << sI << "/" << systStrVect.size() << ": " << systStrVect[sI] << std::endl;
  }  

  std::cout << "JESJER: " << nJESJER << std::endl;
  std::cout << "JESJER: " << jesJERStrVect.size() << std::endl;
  std::cout << "UNFSYST: " << unfSystStrVect.size() << std::endl;
  //  return 1;
  */
  
  //Grab the in-file defined bins
  Int_t nGammaPtBins = inResponseFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  Float_t gammaPtBinsLow = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  Float_t gammaPtBinsHigh = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  Bool_t gammaPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  Bool_t gammaPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxBins+1];
  Float_t gammaPtBinsLowReco = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOWRECO", gammaPtBinsLow);
  Float_t gammaPtBinsHighReco = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGHRECO", gammaPtBinsHigh);
  std::string gammaPtBinsStr = "";
  
  const Int_t nGammaFinePtBins = inResponseFileConfig_p->GetValue("NGAMMAFINEPTBINS", 0);
  std::string gammaFinePtBinsStr = inResponseFileConfig_p->GetValue("GAMMAFINEPTBINS", "");
  std::vector<float> gammaFinePtBinsVect = strToVectF(gammaFinePtBinsStr);
  Double_t gammaFinePtBins[nMaxBins+1];

  for(unsigned int i = 0; i < gammaFinePtBinsVect.size(); ++i){gammaFinePtBins[i] = gammaFinePtBinsVect[i];}
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    gammaPtBinsStr = inResponseFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
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

  const Int_t nJtPtBins = inResponseFileConfig_p->GetValue("NJTPTBINS", 0);
  const Float_t jtPtBinsLow = inResponseFileConfig_p->GetValue("JTPTBINSLOW", -1.0);
  const Float_t jtPtBinsHigh = inResponseFileConfig_p->GetValue("JTPTBINSHIGH", -1.0);
  const Float_t jtPtBinsLowReco = inResponseFileConfig_p->GetValue("JTPTBINSLOWRECO", -1.0);
  const Float_t jtPtBinsHighReco = inResponseFileConfig_p->GetValue("JTPTBINSHIGHRECO", -1.0);
  const Bool_t jtPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOLOG", 1);
  const Bool_t jtPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("JTPTBINSDOCUSTOM", 1);
  Double_t jtPtBins[nMaxBins+1];

  //Get the eta cuts - default absurd to break things when not stored
  const Float_t jtEtaBinsLow = inResponseFileConfig_p->GetValue("JTETABINSLOW", -100000.0);
  const Float_t jtEtaBinsHigh = inResponseFileConfig_p->GetValue("JTETABINSHIGH", -100000.0);

  if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
  else if(jtPtBinsDoCustom){
    std::string jtPtBinsStr = inResponseFileConfig_p->GetValue("JTPTBINSCUSTOM", "");
    if(jtPtBinsStr.size() == 0){
      std::cout << "jtPtBinsDoCustom is true but jtPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> jtPtBinsTemp = strToVectF(jtPtBinsStr);
    Int_t nJtPtBinsTemp = ((Int_t)jtPtBinsTemp.size()) - 1;
    if(nJtPtBinsTemp != nJtPtBins){
      std::cout << "Number of jt pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < jtPtBinsTemp.size(); ++gI){
      jtPtBins[gI] = jtPtBinsTemp[gI];
    }
  }
  else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);

  const Double_t gammaJtDRExclusionCut = inResponseFileConfig_p->GetValue("GAMMAEXCLUSIONDR", 100.0);
  const Double_t mixJtDRExclusionCut = inResponseFileConfig_p->GetValue("MIXJETEXCLUSIONDR", 100.0);
  const Double_t gammaJtDPhiCut = mathStringToNum(inResponseFileConfig_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);
  const Double_t gammaMultiJtDPhiCut = mathStringToNum(inResponseFileConfig_p->GetValue("GAMMAMULTIJTDPHI", ""), doGlobalDebug);
  
  const Int_t nXJBins = inResponseFileConfig_p->GetValue("NXJBINS", 20);
  const Float_t xjBinsLow = inResponseFileConfig_p->GetValue("XJBINSLOW", 0.0);
  const Float_t xjBinsHigh = inResponseFileConfig_p->GetValue("XJBINSHIGH", 2.0);
  const Bool_t xjBinsDoCustom = inResponseFileConfig_p->GetValue("XJBINSDOCUSTOM", 0);
  Double_t xjBins[nMaxBins+1];
  const Float_t xjBinsLowReco = inResponseFileConfig_p->GetValue("XJBINSLOWRECO", xjBinsLow);
  const Float_t xjBinsHighReco = inResponseFileConfig_p->GetValue("XJBINSHIGHRECO", xjBinsHigh);

  if(!xjBinsDoCustom) getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);
  else{
    std::string xjBinsCustomStr = inResponseFileConfig_p->GetValue("XJBINSCUSTOM", "");
    std::vector<float> xjBinsCustomVect = strToVectF(xjBinsCustomStr);
    for(unsigned int i = 0; i < xjBinsCustomVect.size(); ++i){
      xjBins[i] = xjBinsCustomVect[i];
    }
  }

  const Int_t nXJJBins = inResponseFileConfig_p->GetValue("NXJJBINS", 20);
  const Float_t xjjBinsLow = inResponseFileConfig_p->GetValue("XJJBINSLOW", 0.0);
  const Float_t xjjBinsHigh = inResponseFileConfig_p->GetValue("XJJBINSHIGH", 2.0);
  const Bool_t xjjBinsDoCustom = inResponseFileConfig_p->GetValue("XJJBINSDOCUSTOM", 0);
  Double_t xjjBins[nMaxBins+1];
  const Float_t xjjBinsLowReco = inResponseFileConfig_p->GetValue("XJJBINSLOWRECO", xjjBinsLow);
  const Float_t xjjBinsHighReco = inResponseFileConfig_p->GetValue("XJJBINSHIGHRECO", xjjBinsHigh);


  if(!xjjBinsDoCustom) getLinBins(xjjBinsLow, xjjBinsHigh, nXJJBins, xjjBins);
  else{
    std::string xjjBinsCustomStr = inResponseFileConfig_p->GetValue("XJJBINSCUSTOM", "");
    std::vector<float> xjjBinsCustomVect = strToVectF(xjjBinsCustomStr);
    for(unsigned int i = 0; i < xjjBinsCustomVect.size(); ++i){
      xjjBins[i] = xjjBinsCustomVect[i];
    }
  }

  const Int_t nAJBins = inResponseFileConfig_p->GetValue("NAJBINS", 20);
  const Float_t ajBinsLow = inResponseFileConfig_p->GetValue("AJBINSLOW", 0.0);
  const Float_t ajBinsHigh = inResponseFileConfig_p->GetValue("AJBINSHIGH", 2.0);
  const Bool_t ajBinsDoCustom = inResponseFileConfig_p->GetValue("AJBINSDOCUSTOM", 0);
  Double_t ajBins[nMaxBins+1];
  const Float_t ajBinsLowReco = inResponseFileConfig_p->GetValue("AJBINSLOWRECO", ajBinsLow);
  const Float_t ajBinsHighReco = inResponseFileConfig_p->GetValue("AJBINSHIGHRECO", ajBinsHigh);

  if(!ajBinsDoCustom) getLinBins(ajBinsLow, ajBinsHigh, nAJBins, ajBins);
  else{
    std::string ajBinsCustomStr = inResponseFileConfig_p->GetValue("AJBINSCUSTOM", "");
    std::vector<float> ajBinsCustomVect = strToVectF(ajBinsCustomStr);
    for(unsigned int i = 0; i < ajBinsCustomVect.size(); ++i){
      ajBins[i] = ajBinsCustomVect[i];
    }
  }
  
  const Int_t nDPhiBins = inResponseFileConfig_p->GetValue("NDPHIBINS", 20);
  const Float_t dphiBinsLow = inResponseFileConfig_p->GetValue("DPHIBINSLOW", 0.0);
  const Float_t dphiBinsHigh = inResponseFileConfig_p->GetValue("DPHIBINSHIGH", TMath::Pi());
  const Bool_t dphiBinsDoCustom = inResponseFileConfig_p->GetValue("DPHIBINSDOCUSTOM", 0);
  Double_t dphiBins[nMaxBins+1];
  if(!dphiBinsDoCustom) getLinBins(dphiBinsLow, dphiBinsHigh, nDPhiBins, dphiBins);
  else{
    std::string dphiBinsCustomStr = inResponseFileConfig_p->GetValue("DPHIBINSCUSTOM", "");
    std::vector<float> dphiBinsCustomVect = strToVectF(dphiBinsCustomStr);
    for(unsigned int i = 0; i < dphiBinsCustomVect.size(); ++i){
      dphiBins[i] = dphiBinsCustomVect[i];
    }
  }

  const Int_t nDRBins = inResponseFileConfig_p->GetValue("NDRBINS", 20);
  const Float_t drBinsLow = inResponseFileConfig_p->GetValue("DRBINSLOW", 0.0);
  const Float_t drBinsHigh = inResponseFileConfig_p->GetValue("DRBINSHIGH", TMath::Pi());
  const Bool_t drBinsDoCustom = inResponseFileConfig_p->GetValue("DRBINSDOCUSTOM", 0);
  Double_t drBins[nMaxBins+1];
  if(!drBinsDoCustom) getLinBins(drBinsLow, drBinsHigh, nDRBins, drBins);
  else{
    std::string drBinsCustomStr = inResponseFileConfig_p->GetValue("DRBINSCUSTOM", "");
    std::vector<float> drBinsCustomVect = strToVectF(drBinsCustomStr);
    for(unsigned int i = 0; i < drBinsCustomVect.size(); ++i){
      drBins[i] = drBinsCustomVect[i];
    }
  }
  
  //Enforce additional requirements on config depending on first set
  //  if(isMC){
  //    if(!checkEnvForParams(inResponseFileConfig_p, mcParams)) return 1;
  //  }

  int nCentBins = 1;
  //  Float_t centBins[nMaxCentBins];
  std::vector<std::string> centBinsStr = {"PP"};
  std::vector<float> centBinsF = {};
  if(!isPP){
    if(!checkEnvForParams(inResponseFileConfig_p, pbpbParams)) return 1;

    centBinsStr.clear();
    centBinsF = strToVectF(inResponseFileConfig_p->GetValue("CENTBINS", ""));
    nCentBins = centBinsF.size()-1;
    std::vector<std::string> tempCentBinsStr = strToVect(inResponseFileConfig_p->GetValue("CENTBINS", ""));
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
    }

    nCentBins = centBinsStr.size();
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  //Enforce matching over specified parameters
  if(!isSameInputFile){
    if(!compEnvParams(inResponseFileConfig_p, inUnfoldFileConfig_p, paramsToCompare)) return 1;
  }

  //Gotta do this bit manually
  Int_t nVarBins = nXJBins;
  Double_t varBins[nMaxBins+1];

  Float_t varBinsLow = xjBinsLow;
  Float_t varBinsHigh = xjBinsHigh;
  Float_t varBinsLowReco = xjBinsLowReco;
  Float_t varBinsHighReco = xjBinsHighReco;

  std::string varBinsStr = "";
  
  if(isStrSame(varNameLower, "pt")){
    nVarBins = nJtPtBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = jtPtBins[vI];
    }

    varBinsLow = jtPtBinsLow;
    varBinsHigh = jtPtBinsHigh;
    varBinsLowReco = jtPtBinsLowReco;
    varBinsHighReco = jtPtBinsHighReco;
  }
  else if(isStrSame(varNameLower, "xj")){
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = xjBins[vI];
    }
  }
  else if(isStrSame(varNameLower, "xjj")){
    nVarBins = nXJJBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = xjjBins[vI];
    }

    varBinsLow = xjjBinsLow;
    varBinsHigh = xjjBinsHigh;
    varBinsLowReco = xjjBinsLowReco;
    varBinsHighReco = xjjBinsHighReco;
  }
  else if(isStrSame(varNameLower, "ajj")){
    nVarBins = nAJBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = ajBins[vI];
    }

    varBinsLow = ajBinsLow;
    varBinsHigh = ajBinsHigh;
    varBinsLowReco = ajBinsLowReco;
    varBinsHighReco = ajBinsHighReco;
  }
  else if(isStrSame(varNameLower, "dphi") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "dphijjg")){
    nVarBins = nDPhiBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = dphiBins[vI];
    }

    varBinsLow = dphiBinsLow;
    varBinsHigh = dphiBinsHigh;
  }
  else if(isStrSame(varNameLower, "drjj")){
    nVarBins = nDRBins;
    for(Int_t vI = 0; vI < nVarBins+1; ++vI){
      varBins[vI] = drBins[vI];
    }
    varBinsLow = drBinsLow;
    varBinsHigh = drBinsHigh;
  }

  
  //Construct matrices
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  

  TH1D* photonPtReco_p[nMaxCentBins];
  TH1D* photonPtTruth_p[nMaxCentBins];  

  TH2D* photonPtJetVarReco_p[nMaxCentBins];
  TH2D* photonPtJetVarTruth_p[nMaxCentBins]; 
  
  //  const Int_t nBarrelAndEC = 2;
  //  const std::string barrelAndECStr[nBarrelAndEC] = {"Barrel", "EC"};

  //Define some necessary variables for handling Barrel/EC separately or unified like labels
  const Int_t nBarrelAndECMax = 2;
  Int_t nBarrelAndEC = nBarrelAndECMax;
  std::string barrelAndECStr[nBarrelAndECMax] = {"Barrel", "EC"};
  if(doUnifiedPurity){
    nBarrelAndEC = 1;
    barrelAndECStr[0] = "BarrelAndEC";

    std::cout << "SWITCHING TO UNIFIED PURITY MODE" << std::endl;
  }
  
  TH1D* photonPtReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];  
  TH1D* photonPtReco_PURCORR_COMBINED_p[nMaxCentBins];  
  TH1D* photonPtReco_PURCORR_ForReweight_p[nMaxCentBins][nBarrelAndECMax];  
  TH1D* photonPtReco_PURCORR_COMBINED_ForReweight_p[nMaxCentBins];  
  TH1D* photonPtReco_PURCORR_COMBINED_Reweighted_p[nMaxCentBins];  
  TH1D* photonPtReco_TRUTH_p[nMaxCentBins][nBarrelAndECMax];
  TH1D* photonPtReco_TRUTH_COMBINED_p[nMaxCentBins];

  TH2D* photonPtJetVarReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtJetVarReco_PURCORR_COMBINED_p[nMaxCentBins];
  TH2D* photonPtJetVarReco_PURCORR_ForReweight_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[nMaxCentBins];
  TH2D* photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[nMaxCentBins];
  TH2D* photonPtJetVarReco_TRUTH_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtJetVarReco_TRUTH_COMBINED_p[nMaxCentBins];
  
  RooUnfoldResponse* rooResGamma_p[nMaxCentBins];  
  TH2D* rooResGammaMatrix_p[nMaxCentBins];
  TH1D* rooResGammaMisses_p[nMaxCentBins];

  RooUnfoldResponse* rooResGammaJetVar_p[nMaxCentBins][nMaxSyst]; 
  TH2D* rooResGammaJetVarMatrixReco_p[nMaxCentBins][nMaxSyst];
  TH2D* rooResGammaJetVarMatrixTruth_p[nMaxCentBins][nMaxSyst];
  TH2D* rooResGammaJetVarMisses_p[nMaxCentBins][nMaxSyst];
  

  //Declare some additional histograms for reweighting
  TH2D* photonPtDRJJReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtDRJJReco_PURCORR_COMBINED_p[nMaxCentBins];
  TH2D* photonPtDRJJReco_PURCORR_ForReweight_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[nMaxCentBins];  
  TH2D* photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[nMaxCentBins];  

  TH2D* photonPtAJJReco_PURCORR_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtAJJReco_PURCORR_COMBINED_p[nMaxCentBins];
  TH2D* photonPtAJJReco_PURCORR_ForReweight_p[nMaxCentBins][nBarrelAndECMax];
  TH2D* photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[nMaxCentBins];  
  TH2D* photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[nMaxCentBins];  
  
  //Lets check if a rebin is requested that it is valid
  if(doRebin){
    std::string histName = centBinsStr[0] + "/photonPtJt" + varName + "VCent_" + centBinsStr[0] + "_" + barrelAndECStr[0] + "_" + gammaJtDPhiStr + "_" + mixModeStr + "_SUB_h";
    TH2D* tempTH2DorCheck_p = (TH2D*)inUnfoldFile_p->Get(histName.c_str());

    //Check that the macro histogram bins align with requested new binning                      
    if(!checkHistContainsBins(rebinXVect, tempTH2DorCheck_p, deltaValue, false)) return 1;
    if(!checkHistContainsBins(rebinYVect, tempTH2DorCheck_p, deltaValue, true)) return 1;    

    //Since we got out of previous check, reset the varbins!
    nVarBins = nRebinX;
    for(Int_t vI = 0; vI < nRebinX+1; ++vI){
      varBins[vI] = rebinXVect[vI];
    }
    varBinsLow = varBins[0];
    varBinsHigh = varBins[nVarBins];

    for(int bI = 0; bI < nVarBins+1; ++bI){
      std::string newStr = std::to_string(varBins[bI]);
      if(newStr.find(".") != std::string::npos){
	while(newStr.rfind("0") == newStr.size()-1){
	  if(newStr.rfind(".") == newStr.size()-2) break;
	  newStr = newStr.substr(0, newStr.size()-1);
	}
      }
      
      varBinsStr = varBinsStr + newStr + ",";
    }    
        
    nGammaPtBins = nRebinY;
    for(Int_t vI = 0; vI < nRebinY+1; ++vI){
      gammaPtBins[vI] = rebinYVect[vI];
    }

    gammaPtBinsLow = gammaPtBins[0];
    gammaPtBinsHigh = gammaPtBins[nGammaPtBins];
    gammaPtBinsDoLog = false;
    gammaPtBinsDoCustom = true;
    gammaPtBinsStr = "";
    for(int bI = 0; bI < nGammaPtBins+1; ++bI){
      std::string newStr = std::to_string(gammaPtBins[bI]);
      if(newStr.find(".") != std::string::npos){
	while(newStr.rfind("0") == newStr.size()-1){
	  if(newStr.rfind(".") == newStr.size()-2) break;
	  newStr = newStr.substr(0, newStr.size()-1);
	}
      }
      
      gammaPtBinsStr = gammaPtBinsStr + newStr + ",";
    }    
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  
  for(Int_t cI = 0; cI < nCentBins; ++ cI){    
    photonPtReco_p[cI] = new TH1D(("photonPtReco_" + centBinsStr[cI] + "_h").c_str(), ";Reco. Photon p_{T};Counts", nGammaPtBins, gammaPtBins);
    photonPtTruth_p[cI] = new TH1D(("photonPtTruth_" + centBinsStr[cI] + "_h").c_str(), ";Truth Photon p_{T};Counts", nGammaPtBins, gammaPtBins);

    photonPtJetVarReco_p[cI] = new TH2D(("photonPtJet" + varName + "Reco_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetVarTruth_p[cI] = new TH2D(("photonPtJet" + varName + "Truth_" + centBinsStr[cI] + "_h").c_str(), (";Truth " + varNameStyle + ";Truth Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    
    for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
      std::string baseName = centBinsStr[cI] + "/photonPt";
      std::string jtName = "Jt" + varName;
      std::string midName =  "VCent_" + centBinsStr[cI] + "_" + barrelAndECStr[eI];
      std::string photonPtPurCorrName = baseName + midName + "_PURCORR_h";
      std::string photonPtJtVarPurCorrName = baseName + jtName + midName + "_" + gammaJtDPhiStr + "_PURCORR_h";
      std::string photonPtTruthName = baseName + midName + "_TRUTH_h";
      std::string photonPtJtVarTruthName = baseName + jtName + midName + "_" + gammaJtDPhiStr + "_" + mixModeStr + "_TRUTH_h";
      
      TFile* tempFilePointer_p = inUnfoldFile_p;
      if(!isMC) tempFilePointer_p = inResponseFile_p;
      
      if(doRebin){
	TH1D* tempHistTH1D_p = (TH1D*)inUnfoldFile_p->Get(photonPtPurCorrName.c_str());
	strReplace(&photonPtPurCorrName, "PURCORR_h", "PURCORR_REBIN_h");
	std::string xTitle = tempHistTH1D_p->GetXaxis()->GetTitle();
	std::string yTitle = tempHistTH1D_p->GetYaxis()->GetTitle();
	std::string titleStr = ";" + xTitle + ";" + yTitle;
	photonPtReco_PURCORR_p[cI][eI] = new TH1D(photonPtPurCorrName.c_str(), titleStr.c_str(), nGammaPtBins, gammaPtBins);
	fineHistToCoarseHist(tempHistTH1D_p, photonPtReco_PURCORR_p[cI][eI]);
	
	TH2D* tempHistTH2D_p = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());
	strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_h");		  
	xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;
	photonPtJetVarReco_PURCORR_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);	
	fineHistToCoarseHist(tempHistTH2D_p, photonPtJetVarReco_PURCORR_p[cI][eI]);

	//For reweighting
	strReplace(&photonPtPurCorrName, "PURCORR_REBIN_h", "PURCORR_h");
	tempHistTH1D_p = (TH1D*)tempFilePointer_p->Get(photonPtPurCorrName.c_str());
	strReplace(&photonPtPurCorrName, "PURCORR_h", "PURCORR_REBIN_ForReweight_h");
	xTitle = tempHistTH1D_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH1D_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;
	photonPtReco_PURCORR_ForReweight_p[cI][eI] = new TH1D(photonPtPurCorrName.c_str(), titleStr.c_str(), nGammaPtBins, gammaPtBins);
	fineHistToCoarseHist(tempHistTH1D_p, photonPtReco_PURCORR_ForReweight_p[cI][eI]);

	strReplace(&photonPtJtVarPurCorrName, "PURCORR_REBIN_h", "PURCORR_h");
	tempHistTH2D_p = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());
	strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_ForReweight_h"); 
	xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;
	photonPtJetVarReco_PURCORR_ForReweight_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);	
	fineHistToCoarseHist(tempHistTH2D_p, photonPtJetVarReco_PURCORR_ForReweight_p[cI][eI]);

	if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
     	//For reweighting DR - note we want to check this even when the doReweightDR = false;
	if(doReweightVar){
	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_REBIN_ForReweight_h", "PURCORR_h");
	  strReplace(&photonPtJtVarPurCorrName, jtName, "JtDRJJ");

	  tempHistTH2D_p = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());
	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_h"); 
	  xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	  yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	  titleStr = ";" + xTitle + ";" + yTitle;
	  photonPtDRJJReco_PURCORR_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nDRBins, drBins, nGammaPtBins, gammaPtBins);	
	  fineHistToCoarseHist(tempHistTH2D_p, photonPtDRJJReco_PURCORR_p[cI][eI]);

	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_REBIN_h", "PURCORR_h");
	  tempHistTH2D_p = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());
	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_ForReweight_h"); 
	  xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	  yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	  titleStr = ";" + xTitle + ";" + yTitle;
	  photonPtDRJJReco_PURCORR_ForReweight_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nDRBins, drBins, nGammaPtBins, gammaPtBins);	
	  fineHistToCoarseHist(tempHistTH2D_p, photonPtDRJJReco_PURCORR_ForReweight_p[cI][eI]);

	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_REBIN_ForReweight_h", "PURCORR_h");
	  strReplace(&photonPtJtVarPurCorrName, "JtDRJJ", "JtAJJ");

	  tempHistTH2D_p = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());
	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_h"); 
	  xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	  yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	  titleStr = ";" + xTitle + ";" + yTitle;
	  photonPtAJJReco_PURCORR_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nAJBins, ajBins, nGammaPtBins, gammaPtBins);	
	  fineHistToCoarseHist(tempHistTH2D_p, photonPtAJJReco_PURCORR_p[cI][eI]);

	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_REBIN_h", "PURCORR_h");
	  tempHistTH2D_p = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());
	  strReplace(&photonPtJtVarPurCorrName, "PURCORR_h", "PURCORR_REBIN_ForReweight_h"); 
	  xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	  yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	  titleStr = ";" + xTitle + ";" + yTitle;
	  photonPtAJJReco_PURCORR_ForReweight_p[cI][eI] = new TH2D(photonPtJtVarPurCorrName.c_str(), titleStr.c_str(), nAJBins, ajBins, nGammaPtBins, gammaPtBins);	
	  fineHistToCoarseHist(tempHistTH2D_p, photonPtAJJReco_PURCORR_ForReweight_p[cI][eI]);
	}
	
	//Truth histograms
	tempHistTH1D_p = (TH1D*)tempFilePointer_p->Get(photonPtTruthName.c_str());
	strReplace(&photonPtTruthName, "TRUTH_h", "TRUTH_REBIN_h");
	xTitle = tempHistTH1D_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH1D_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;	
	photonPtReco_TRUTH_p[cI][eI] = new TH1D(photonPtTruthName.c_str(), titleStr.c_str(), nGammaPtBins, gammaPtBins);
	fineHistToCoarseHist(tempHistTH1D_p, photonPtReco_TRUTH_p[cI][eI]);
	
	tempHistTH2D_p = (TH2D*)tempFilePointer_p->Get(photonPtJtVarTruthName.c_str());
	strReplace(&photonPtJtVarTruthName, "TRUTH_h", "TRUTH_REBIN_h");		  
	xTitle = tempHistTH2D_p->GetXaxis()->GetTitle();
	yTitle = tempHistTH2D_p->GetYaxis()->GetTitle();
	titleStr = ";" + xTitle + ";" + yTitle;
	photonPtJetVarReco_TRUTH_p[cI][eI] = new TH2D(photonPtJtVarTruthName.c_str(), titleStr.c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
	fineHistToCoarseHist(tempHistTH2D_p, photonPtJetVarReco_TRUTH_p[cI][eI]);
      }
      else{
	photonPtReco_PURCORR_p[cI][eI] = (TH1D*)inUnfoldFile_p->Get(photonPtPurCorrName.c_str());     
	photonPtJetVarReco_PURCORR_p[cI][eI] = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());     
	
	photonPtReco_TRUTH_p[cI][eI] = (TH1D*)tempFilePointer_p->Get(photonPtTruthName.c_str());
	photonPtJetVarReco_TRUTH_p[cI][eI] = (TH2D*)tempFilePointer_p->Get(photonPtJtVarTruthName.c_str());

      	photonPtReco_PURCORR_ForReweight_p[cI][eI] = (TH1D*)tempFilePointer_p->Get(photonPtPurCorrName.c_str());     
	photonPtJetVarReco_PURCORR_ForReweight_p[cI][eI] = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());     
      
	//For dR based reweighitng
	photonPtJtVarPurCorrName.replace(photonPtJtVarPurCorrName.find(jtName), jtName.size(), "JtDRJJ");
	photonPtDRJJReco_PURCORR_p[cI][eI] = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());     
	photonPtDRJJReco_PURCORR_ForReweight_p[cI][eI] = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());     

	//For AJJ based reweighitng
	photonPtJtVarPurCorrName.replace(photonPtJtVarPurCorrName.find("JtDRJJ"), 6, "JtAJJ");
	photonPtAJJReco_PURCORR_p[cI][eI] = (TH2D*)inUnfoldFile_p->Get(photonPtJtVarPurCorrName.c_str());     
	photonPtAJJReco_PURCORR_ForReweight_p[cI][eI] = (TH2D*)tempFilePointer_p->Get(photonPtJtVarPurCorrName.c_str());     
      }
    }

    photonPtReco_PURCORR_COMBINED_p[cI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), ";Reco. Photon p_{T};Counts (Purity Corrected)", nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_PURCORR_COMBINED_p[cI] = new TH2D(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtDRJJReco_PURCORR_COMBINED_p[cI] = new TH2D(("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_h").c_str(), ";Reco. #DeltaR_{JJ};Reco. Photon p_{T}", nDRBins, drBins, nGammaPtBins, gammaPtBins);
    photonPtAJJReco_PURCORR_COMBINED_p[cI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_h").c_str(), ";Reco. A_{JJ};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
    
    photonPtReco_TRUTH_COMBINED_p[cI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_TRUTH_COMBINED_p[cI] = new TH2D(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_TRUTH_COMBINED_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);   

    photonPtReco_PURCORR_COMBINED_ForReweight_p[cI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str(), ";Reco. Photon p_{T};Counts (Purity Corrected)", nGammaPtBins, gammaPtBins);
    photonPtReco_PURCORR_COMBINED_Reweighted_p[cI] = new TH1D(("photonPtVCent_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str(), ";Reco. Photon p_{T};Counts (Purity Corrected)", nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI] = new TH2D(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_ForReweight_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
    photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[cI] = new TH2D(("photonPtJt" + varName + "VCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_Reweighted_h").c_str(), (";Reco. " + varNameStyle + ";Reco. Photon p_{T}").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);

    photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI] = new TH2D(("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_ForReweight_h").c_str(), ";Reco. #DeltaR_{JJ};Reco. Photon p_{T}", nDRBins, drBins, nGammaPtBins, gammaPtBins);
    photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[cI] = new TH2D(("photonPtJtDRJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_Reweighted_h").c_str(), ";Reco. #DeltaR_{JJ};Reco. Photon p_{T}", nDRBins, drBins, nGammaPtBins, gammaPtBins);

    photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_ForReweight_h").c_str(), ";Reco. A_{JJ};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
    photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[cI] = new TH2D(("photonPtJtAJJVCent_" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_PURCORR_COMBINED_Reweighted_h").c_str(), ";Reco. A_{JJ};Reco. Photon p_{T}", nAJBins, ajBins, nGammaPtBins, gammaPtBins);
     
    for(Int_t bIX = 0; bIX < photonPtReco_TRUTH_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      Float_t content = 0.0;
      Float_t error = 0.0;
    
      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	Float_t tempContent = photonPtReco_TRUTH_p[cI][eI]->GetBinContent(bIX+1);
	Float_t tempError = photonPtReco_TRUTH_p[cI][eI]->GetBinError(bIX+1);
	
	content += tempContent;
	error = TMath::Sqrt(error*error + tempError*tempError);
      }

      photonPtReco_TRUTH_COMBINED_p[cI]->SetBinContent(bIX+1, content);
      photonPtReco_TRUTH_COMBINED_p[cI]->SetBinError(bIX+1, error);
    }
    
    for(Int_t bIX = 0; bIX < photonPtReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      Float_t content = 0.0;
      Float_t error = 0.0;

      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	Float_t tempContent = photonPtReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1);
	Float_t tempError = photonPtReco_PURCORR_p[cI][eI]->GetBinError(bIX+1);
	
	content += tempContent;
	error = TMath::Sqrt(error*error + tempError*tempError);
      }

      photonPtReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, content);
      photonPtReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, error);

      content = 0.0;
      error = 0.0;

      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	Float_t tempContent = photonPtReco_PURCORR_ForReweight_p[cI][eI]->GetBinContent(bIX+1);
	Float_t tempError = photonPtReco_PURCORR_ForReweight_p[cI][eI]->GetBinError(bIX+1);
	
	content += tempContent;
	error = TMath::Sqrt(error*error + tempError*tempError);
      }

      photonPtReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinContent(bIX+1, content);
      photonPtReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinError(bIX+1, error);
    }
  
    for(Int_t bIX = 0; bIX < photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
      for(Int_t bIY = 0; bIY < photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins(); ++bIY){

	Float_t content = 0.0;
	Float_t error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetVarReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetVarReco_PURCORR_p[cI][eI]->GetBinError(bIX+1, bIY+1);

	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetVarReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetVarReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);	

	content = 0.0;
	error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetVarReco_TRUTH_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetVarReco_TRUTH_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	  
	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetVarReco_TRUTH_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetVarReco_TRUTH_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);		

	content = 0.0;
	error = 0.0;

	for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	  Float_t tempContent = photonPtJetVarReco_PURCORR_ForReweight_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	  Float_t tempError = photonPtJetVarReco_PURCORR_ForReweight_p[cI][eI]->GetBinError(bIX+1, bIY+1);

	  content += tempContent;
	  error = TMath::Sqrt(error*error + tempError*tempError);
	}

	photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinError(bIX+1, bIY+1, error);	
      }
    }
  
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    if(doReweightVar){
      for(Int_t bIX = 0; bIX < photonPtDRJJReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
	for(Int_t bIY = 0; bIY < photonPtDRJJReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins(); ++bIY){
	  
	  Float_t content = 0.0;
	  Float_t error = 0.0;

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	    Float_t tempContent = photonPtDRJJReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	    Float_t tempError = photonPtDRJJReco_PURCORR_p[cI][eI]->GetBinError(bIX+1, bIY+1);

	    content += tempContent;
	    error = TMath::Sqrt(error*error + tempError*tempError);
	  }

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  photonPtDRJJReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	  photonPtDRJJReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);	
	  
	  content = 0.0;
	  error = 0.0;

	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	    Float_t tempContent = photonPtDRJJReco_PURCORR_ForReweight_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	    Float_t tempError = photonPtDRJJReco_PURCORR_ForReweight_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	    
	    content += tempContent;
	    error = TMath::Sqrt(error*error + tempError*tempError);
	  }

	  photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	  photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinError(bIX+1, bIY+1, error);
	}
      }


      //Repeat for AJJ / consider rewriting as function
      for(Int_t bIX = 0; bIX < photonPtAJJReco_PURCORR_COMBINED_p[cI]->GetXaxis()->GetNbins(); ++bIX){
	for(Int_t bIY = 0; bIY < photonPtAJJReco_PURCORR_COMBINED_p[cI]->GetYaxis()->GetNbins(); ++bIY){
	  
	  Float_t content = 0.0;
	  Float_t error = 0.0;
	  
	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	    Float_t tempContent = photonPtAJJReco_PURCORR_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	    Float_t tempError = photonPtAJJReco_PURCORR_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	    
	    content += tempContent;
	    error = TMath::Sqrt(error*error + tempError*tempError);
	  }
	  
	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  photonPtAJJReco_PURCORR_COMBINED_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	  photonPtAJJReco_PURCORR_COMBINED_p[cI]->SetBinError(bIX+1, bIY+1, error);	
	  
	  content = 0.0;
	  error = 0.0;

	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	    Float_t tempContent = photonPtAJJReco_PURCORR_ForReweight_p[cI][eI]->GetBinContent(bIX+1, bIY+1);
	    Float_t tempError = photonPtAJJReco_PURCORR_ForReweight_p[cI][eI]->GetBinError(bIX+1, bIY+1);
	    
	    content += tempContent;
	    error = TMath::Sqrt(error*error + tempError*tempError);
	  }

	      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  
	  photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinContent(bIX+1, bIY+1, content);
	  photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->SetBinError(bIX+1, bIY+1, error);
	}
      }
    }
   
    std::vector<TH1*> hists_p = {photonPtJetVarReco_p[cI], photonPtJetVarTruth_p[cI]};
    hists_p.push_back(photonPtJetVarReco_PURCORR_COMBINED_p[cI]);
    hists_p.push_back(photonPtJetVarReco_TRUTH_COMBINED_p[cI]);
    hists_p.push_back(photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]);

    if(doReweightDR){     
      hists_p.push_back(photonPtDRJJReco_PURCORR_COMBINED_p[cI]);
      hists_p.push_back(photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]);

      hists_p.push_back(photonPtAJJReco_PURCORR_COMBINED_p[cI]);
      hists_p.push_back(photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]);
    }
    
    setSumW2(hists_p);
  }

  //Lets immediately contruct the reweighting histogram
  TH1D* reweightPhoPt_p[nMaxCentBins];
  TH2D* reweightPhoPtJetVar_p[nMaxCentBins];
  TH2D* reweightPhoPtDRJJ_p[nMaxCentBins];
  TH2D* reweightPhoPtAJJ_p[nMaxCentBins];

  if(!isMC){
    if(doReweightVar){
      for(Int_t cI = 0; cI < nCentBins; ++cI){    
	reweightPhoPt_p[cI] = new TH1D(("reweightPhoPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma p_{T} [GeV];Weight Factor", nGammaPtBins, gammaPtBins);	
	reweightPhoPtJetVar_p[cI] = new TH2D(("reweightPhoPtJetVar_" + centBinsStr[cI] + "_h").c_str(), (";Reco. " + varNameStyle + ";#gamma p_{T} [GeV]").c_str(), nVarBins, varBins, nGammaPtBins, gammaPtBins);
	
	std::string saveName = pdfStr + "/reweightPhoPt_" + centBinsStr[cI] + "_" + dateStr + ".png";
	makeReweightHist(reweightPhoPt_p[cI], photonPtReco_PURCORR_COMBINED_p[cI], photonPtReco_PURCORR_COMBINED_ForReweight_p[cI], saveName);

	saveName = pdfStr + "/reweightPhoPtJt" + varName + "_" + centBinsStr[cI] + "_" + dateStr + ".png";	      
	makeReweightHist(reweightPhoPtJetVar_p[cI], photonPtJetVarReco_PURCORR_COMBINED_p[cI], photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI], saveName); 
      }     
    }

   if(doReweightVar){
      for(Int_t cI = 0; cI < nCentBins; ++cI){    
	reweightPhoPtDRJJ_p[cI] = new TH2D(("reweightPhoPtDRJJ_" + centBinsStr[cI] + "_h").c_str(), ";Reco. #DeltaR_{JJ};#gamma p_{T}", nDRBins, drBins, nGammaPtBins, gammaPtBins);
	
	std::string saveName = pdfStr + "/reweightPhoPtDRJJ_" + centBinsStr[cI] + "_" + dateStr + ".png";
	makeReweightHist(reweightPhoPtDRJJ_p[cI], photonPtDRJJReco_PURCORR_COMBINED_p[cI], photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI], saveName); 

	reweightPhoPtAJJ_p[cI] = new TH2D(("reweightPhoPtAJJ_" + centBinsStr[cI] + "_h").c_str(), ";Reco. #DeltaR_{JJ};#gamma p_{T}", nDRBins, drBins, nGammaPtBins, gammaPtBins);	
	saveName = pdfStr + "/reweightPhoPtAJJ_" + centBinsStr[cI] + "_" + dateStr + ".png";
	makeReweightHist(reweightPhoPtAJJ_p[cI], photonPtAJJReco_PURCORR_COMBINED_p[cI], photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI], saveName); 

      }
    }
  }

  //Grab the unfolding ttree for matrix filling
  TTree* unfoldTree_p = (TTree*)inResponseFile_p->Get("unfoldTree_p");

  //Before declaring variables define array size for JES/JER syst. + check it matches input
  const Int_t nJESSys = 18;
  const Int_t nJERSys = 9;
  //1 (nominal) + nJESSys + 1 (centrality dependent rtrk JESSys manual) + nJERSys                
  const Int_t nJetSysAndNom = 1 + nJESSys + 1 + nJERSys;
  const Int_t nJetSysAndNomInput = inResponseFileConfig_p->GetValue("NJETSYSANDNOM", -1);

  if(nJetSysAndNom != nJetSysAndNomInput){
    std::cout << "Mismatch between number of jet systematics expected and input from file \'" << inResponseFileName << "\'. return 1" << std::endl;
    return 1;
  }

  //Run,lumi,event check for duplication
  UInt_t runNumber;
  UInt_t lumiBlock;
  ULong64_t eventNumber;
  
  //Declare unfolding TTree variables
  const Int_t nMaxJets = 100;

  Float_t recoGammaPt_;
  Float_t recoGammaPhi_;
  Float_t recoGammaEta_;  

  Float_t truthGammaPt_;
  Float_t truthGammaPhi_;
  Float_t truthGammaEta_;

  Int_t nRecoJt_;
  Float_t recoJtPt_[nMaxJets][nJetSysAndNom];
  Float_t recoJtPhi_[nMaxJets];
  Float_t recoJtEta_[nMaxJets];

  Float_t truthJtPt_[nMaxJets];
  Float_t truthJtPhi_[nMaxJets];
  Float_t truthJtEta_[nMaxJets];

  Int_t nTruthJtUnmatched_;
  Float_t truthJtUnmatchedPt_[nMaxJets];
  Float_t truthJtUnmatchedPhi_[nMaxJets];
  Float_t truthJtUnmatchedEta_[nMaxJets];

  Double_t unfoldWeight_;
  Float_t unfoldCent_;

  unfoldTree_p->SetBranchAddress("runNumber", &runNumber);
  unfoldTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  unfoldTree_p->SetBranchAddress("eventNumber", &eventNumber);
  
  unfoldTree_p->SetBranchAddress("recoGammaPt", &recoGammaPt_);
  unfoldTree_p->SetBranchAddress("recoGammaPhi", &recoGammaPhi_);
  unfoldTree_p->SetBranchAddress("recoGammaEta", &recoGammaEta_);

  unfoldTree_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldTree_p->SetBranchAddress("truthGammaPhi", &truthGammaPhi_);
  unfoldTree_p->SetBranchAddress("truthGammaEta", &truthGammaEta_);
  
  unfoldTree_p->SetBranchAddress("nRecoJt", &nRecoJt_);
  unfoldTree_p->SetBranchAddress("recoJtPt", recoJtPt_);
  unfoldTree_p->SetBranchAddress("recoJtPhi", recoJtPhi_);
  unfoldTree_p->SetBranchAddress("recoJtEta", recoJtEta_);
  
  unfoldTree_p->SetBranchAddress("truthJtPt", truthJtPt_);
  unfoldTree_p->SetBranchAddress("truthJtPhi", truthJtPhi_);
  unfoldTree_p->SetBranchAddress("truthJtEta", truthJtEta_);
  
  unfoldTree_p->SetBranchAddress("nTruthJtUnmatched", &nTruthJtUnmatched_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPt", truthJtUnmatchedPt_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPhi", truthJtUnmatchedPhi_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedEta", truthJtUnmatchedEta_);
  
  unfoldTree_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  unfoldTree_p->SetBranchAddress("unfoldCent", &unfoldCent_);    

  bool varIsPhi = varName.find("DPhiJJ") != std::string::npos;
  
  //Starting out we need to calc our gamma normalization
  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    rooResGamma_p[cI] = new RooUnfoldResponse(photonPtReco_p[cI], photonPtTruth_p[cI], ("rooResGamma_" + centBinsStr[cI]).c_str(), "");
    rooResGammaMatrix_p[cI] = new TH2D(("rooResGammaMatrix_" + centBinsStr[cI] + "_h").c_str(), ";Reco #gamma p_{T} [GeV]; Truth #gamma p_{T} [GeV]", nGammaPtBins, gammaPtBins, nGammaPtBins, gammaPtBins);
    rooResGammaMisses_p[cI] = new TH1D(("rooResGammaMisses_" + centBinsStr[cI] + "_h").c_str(), ";Truth #gamma p_{T} [GeV]; Counts (weighted)", nGammaPtBins, gammaPtBins);

    for(Int_t sI = 0; sI < nSyst; ++sI){            
      //If you wish to only check one syst, skip others - else it will say "ALL"
      if(!isStrSame(inUnfoldName, systStrVect[sI])){
	if(!isStrSame(inUnfoldName, "ALL")) continue;
      }      
    
      rooResGammaJetVar_p[cI][sI] = new RooUnfoldResponse(photonPtJetVarReco_p[cI], photonPtJetVarTruth_p[cI], ("rooResGammaJetVar_" + centBinsStr[cI] + "_" + systStrVect[sI]).c_str(), ("rooResGammaJet" + varName + "_" + centBinsStr[cI] + "_" + systStrVect[sI]).c_str());
      rooResGammaJetVarMatrixReco_p[cI][sI] = new TH2D(("rooResJetVarGammaMatrixReco_" + centBinsStr[cI] + "_" + systStrVect[sI] + "_h").c_str(), ";Reco. Jet p_{T} [GeV]; Reco. #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      rooResGammaJetVarMatrixTruth_p[cI][sI] = new TH2D(("rooResJetVarGammaMatrixTruth_" + centBinsStr[cI] + "_" + systStrVect[sI] + "_h").c_str(), ";Truth Jet p_{T} [GeV]; Truth #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
      rooResGammaJetVarMisses_p[cI][sI] = new TH2D(("rooResJetVarGammaMisses_" + centBinsStr[cI] + "_" + systStrVect[sI] + "_h").c_str(), ";Truth Misses Jet p_{T} [GeV]; Truth Misses #gamma p_{T} [GeV]", nJtPtBins, jtPtBins, nGammaPtBins, gammaPtBins);
    }
  }

  std::cout << "Constructing response matrices..." << std::endl;
  const ULong64_t nEntriesUnfold = (ULong64_t)unfoldTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntriesUnfold/20);

  //Check we don't have duplicate events
  std::vector<Int_t> runNumbers, eventNumbers;
  std::vector<UInt_t> lumiBlocks;
  
  for(ULong64_t entry = 0; entry < nEntriesUnfold; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntriesUnfold << "..." << std::endl;
    unfoldTree_p->GetEntry(entry);

    //pushbacks
    runNumbers.push_back(runNumber);
    lumiBlocks.push_back(lumiBlock);
    eventNumbers.push_back(eventNumber);
    
    Int_t centPos = -1;
    if(isPP) centPos = 0;
    else centPos = ghostPos(centBinsF, unfoldCent_, true, doGlobalDebug);
    if(centPos < 0) continue;

    bool truthGammaOutOfBounds = truthGammaPt_ < gammaPtBinsLow || truthGammaPt_ >= gammaPtBinsHigh;
    if(truthGammaOutOfBounds) continue;

    Double_t phoWeight = 1.0;
    if(doReweightVar) phoWeight = reweightPhoPt_p[centPos]->GetBinContent(reweightPhoPt_p[centPos]->FindBin(truthGammaPt_));

    TLorentzVector truthGammaTL;
    truthGammaTL.SetPtEtaPhiM(truthGammaPt_, truthGammaEta_, truthGammaPhi_, 0.0);
    
    //Construct 1-D unfold response matrix for gamma-pt (required in normalization)
    bool recoGammaOutOfBounds = recoGammaPt_ < gammaPtBinsLowReco || recoGammaPt_ >= gammaPtBinsHighReco;   

    if(recoGammaOutOfBounds){
      rooResGamma_p[centPos]->Miss(truthGammaPt_, unfoldWeight_*phoWeight);
      rooResGammaMisses_p[centPos]->Fill(truthGammaPt_, unfoldWeight_*phoWeight);
    }
    else{
      rooResGamma_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_*phoWeight);
      rooResGammaMatrix_p[centPos]->Fill(recoGammaPt_, truthGammaPt_, unfoldWeight_*phoWeight);

      //For checking the reweighting is effective
      if(doReweightVar) photonPtReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(recoGammaPt_, unfoldWeight_*phoWeight);
    }
  
    //Now construct 2-D unfold response matrix for gammapt-jetvariable    
    //Start w/ truth jets, no reco    
    TLorentzVector tL;

    //We have to do this for all syst.
    for(Int_t sysI = 0; sysI < nSyst; ++sysI){
      if(!isStrSame(inUnfoldName, systStrVect[sysI])){
	if(!isStrSame(inUnfoldName, "ALL")) continue;
      }
      std::vector<TLorentzVector> goodUnmatchedTruthJets;

      //Process unmatched truth jets
      for(Int_t tI = 0; tI < nTruthJtUnmatched_; ++tI){
	bool truthJetOutOfBounds = truthJtUnmatchedPt_[tI] < jtPtBinsLow || truthJtUnmatchedPt_[tI] >= jtPtBinsHigh;
	truthJetOutOfBounds = truthJetOutOfBounds || truthJtUnmatchedEta_[tI] < jtEtaBinsLow || truthJtUnmatchedEta_[tI] > jtEtaBinsHigh;
	if(truthJetOutOfBounds) continue;

	TLorentzVector truthJetTL;
	truthJetTL.SetPtEtaPhiM(truthJtUnmatchedPt_[tI], truthJtUnmatchedEta_[tI], truthJtUnmatchedPhi_[tI], 0.0);
      
	Double_t truthJetVar = getVar(varNameLower, truthJetTL, truthJetTL, truthGammaTL);
	Double_t phoJetWeight = 1.0;

	//Reweight if doReweightVar AND its not the prior variation systematic
	if(doReweightVar && !isStrSame(systStrVect[sysI], "PRIOR")) phoJetWeight = reweightPhoPtJetVar_p[centPos]->GetBinContent(reweightPhoPtJetVar_p[centPos]->GetXaxis()->FindBin(truthJetVar), reweightPhoPtJetVar_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));
	
	Float_t gammaJtDRTruth = getDR(truthJtUnmatchedEta_[tI], truthJtUnmatchedPhi_[tI], truthGammaEta_, truthGammaPhi_);
	bool gammaJtPassesDRTruth = gammaJtDRTruth >= gammaJtDRExclusionCut;
	if(!gammaJtPassesDRTruth) continue;	
	
	Float_t gammaJtDPhiTruth = TMath::Abs(getDPHI(truthJtUnmatchedPhi_[tI], truthGammaPhi_));
	if(isStrSame(varNameLower, "dphi")){
	  rooResGammaJetVar_p[centPos][sysI]->Miss(gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	  rooResGammaJetVarMisses_p[centPos][sysI]->Fill(gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	}      
	
	bool gammaJtPassesDPhiTruth = gammaJtDPhiTruth >= gammaJtDPhiCut;
	if(!gammaJtPassesDPhiTruth) continue;
	
	if(isStrSame(varNameLower, "pt")){
	  rooResGammaJetVar_p[centPos][sysI]->Miss(truthJtUnmatchedPt_[tI], truthGammaPt_, unfoldWeight_*phoJetWeight);
	  rooResGammaJetVarMisses_p[centPos][sysI]->Fill(truthJtUnmatchedPt_[tI], truthGammaPt_, unfoldWeight_*phoJetWeight);
	}
	else if(isStrSame(varNameLower, "xj")){
	  Float_t varVal = truthJtUnmatchedPt_[tI]/truthGammaPt_;
	  Bool_t varValGood = varVal >= varBinsLow && varVal < varBinsHigh;
	  
	  if(varValGood){
	    rooResGammaJetVar_p[centPos][sysI]->Miss(varVal, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMisses_p[centPos][sysI]->Fill(varVal, truthGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	}
	else if(isMultijet){
	  tL.SetPtEtaPhiM(truthJtUnmatchedPt_[tI], truthJtUnmatchedEta_[tI], truthJtUnmatchedPhi_[tI], 0.0);
	  goodUnmatchedTruthJets.push_back(tL);
	}
      }
     
      std::vector<TLorentzVector> goodRecoJets, goodRecoTruthMatchJets;
      //Now process reco-truth jet matched pairs
      for(Int_t jI = 0; jI < nRecoJt_; ++jI){
	//The jet vector position will be zero unless we are doing jet related systematics
	Int_t jtVPos = sysI;
	if(sysI >= nJESJER) jtVPos = 0;
	
	bool recoJtOutOfBounds = recoJtPt_[jI][jtVPos] < jtPtBinsLowReco || recoJtPt_[jI][jtVPos] >= jtPtBinsHighReco;
	recoJtOutOfBounds = recoJtOutOfBounds || recoJtEta_[jI] < jtEtaBinsLow || recoJtEta_[jI] >= jtEtaBinsHigh;
	
	bool truthJetOutOfBounds = truthJtPt_[jI] < jtPtBinsLow || truthJtPt_[jI] >= jtPtBinsHigh;
	truthJetOutOfBounds = truthJetOutOfBounds || truthJtEta_[jI] < jtEtaBinsLow || truthJtEta_[jI] > jtEtaBinsHigh;
	
	if(truthJetOutOfBounds) continue;

	Float_t gammaJtDRTruth = getDR(truthJtEta_[jI], truthJtPhi_[jI], truthGammaEta_, truthGammaPhi_);
	bool gammaJtPassesDRTruth = gammaJtDRTruth >= gammaJtDRExclusionCut;
	if(!gammaJtPassesDRTruth) continue;       
     	       	
	Float_t gammaJtDRReco = -999;
	bool gammaJtPassesDRReco = false;
	if(!recoGammaOutOfBounds){
	  gammaJtDRReco = getDR(recoJtEta_[jI], recoJtPhi_[jI], recoGammaEta_, recoGammaPhi_);
	  gammaJtPassesDRReco = gammaJtDRReco >= gammaJtDRExclusionCut;
	}            
	
	Float_t gammaJtDPhiTruth = TMath::Abs(getDPHI(truthJtPhi_[jI], truthGammaPhi_));
	Float_t gammaJtDPhiReco = -999;
	bool gammaJtPassesDPhiReco = false;
	if(!recoGammaOutOfBounds){
	  gammaJtDPhiReco = TMath::Abs(getDPHI(recoJtPhi_[jI], recoGammaPhi_));
	  gammaJtPassesDPhiReco = gammaJtDPhiReco >= gammaJtDPhiCut;
	}
      
	//Find the reweighting bin, first define your observable
	TLorentzVector truthJetTL;
	truthJetTL.SetPtEtaPhiM(truthJtPt_[jI], truthJtEta_[jI], truthJtPhi_[jI], 0.0);

	Double_t truthJetVar = getVar(varNameLower, truthJetTL, truthJetTL, truthGammaTL);
           
	Double_t phoJetWeight = 1.0;
	//Reweight if doReweightVar AND its not the prior variation systematic
	if(doReweightVar && !isStrSame(systStrVect[sysI], "PRIOR")) phoJetWeight = reweightPhoPtJetVar_p[centPos]->GetBinContent(reweightPhoPtJetVar_p[centPos]->GetXaxis()->FindBin(truthJetVar), reweightPhoPtJetVar_p[centPos]->GetYaxis()->FindBin(truthGammaPt_)); 
	
	if(isStrSame(varNameLower, "dphi")){
	  if(recoJtOutOfBounds || recoGammaOutOfBounds || !gammaJtPassesDRReco){
	    rooResGammaJetVar_p[centPos][sysI]->Miss(gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMisses_p[centPos][sysI]->Fill(gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	  else{
	    rooResGammaJetVar_p[centPos][sysI]->Fill(gammaJtDPhiReco, recoGammaPt_, gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMatrixReco_p[centPos][sysI]->Fill(gammaJtDPhiReco, recoGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMatrixTruth_p[centPos][sysI]->Fill(gammaJtDPhiTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);

	    if(doReweightVar && sysI == 0) photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(gammaJtDPhiReco, recoGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	}
        
	bool gammaJtPassesDPhiTruth = gammaJtDPhiTruth >= gammaJtDPhiCut;
	if(!gammaJtPassesDPhiTruth) continue;
	
	//We will add a condition on recoJt being good - passes dphi cut
	if(!gammaJtPassesDPhiReco) recoJtOutOfBounds = true;
      
	//Now that the truth has passed all cuts - fill out the vectors goodRecoJets OR Fill the inclusive jet variables
	if(isStrSame(varNameLower, "pt")){
	  if(recoJtOutOfBounds || recoGammaOutOfBounds){
	    rooResGammaJetVar_p[centPos][sysI]->Miss(truthJtPt_[jI], truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMisses_p[centPos][sysI]->Fill(truthJtPt_[jI], truthGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	  else{
	    rooResGammaJetVar_p[centPos][sysI]->Fill(recoJtPt_[jI][jtVPos], recoGammaPt_, truthJtPt_[jI], truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMatrixReco_p[centPos][sysI]->Fill(recoJtPt_[jI][jtVPos], recoGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMatrixTruth_p[centPos][sysI]->Fill(truthJtPt_[jI], truthGammaPt_, unfoldWeight_*phoJetWeight);

	    if(doReweightVar && sysI == 0) photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(recoJtPt_[jI][jtVPos], recoGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	}
	else if(isStrSame(varNameLower, "xj")){
	  Float_t varValTruth = truthJtPt_[jI]/truthGammaPt_;
	  Float_t varValTruthGood = varValTruth >= varBins[0] && varValTruth < varBins[nVarBins];
	  
	  if(varValTruthGood){
	    Float_t varValReco = recoJtPt_[jI][jtVPos]/recoGammaPt_;
	    Float_t varValRecoGood = varValReco >= varBinsLowReco && varValReco < varBinsHighReco;
	    
	    if(recoJtOutOfBounds || recoGammaOutOfBounds || !varValRecoGood){
	      rooResGammaJetVar_p[centPos][sysI]->Miss(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	      rooResGammaJetVarMisses_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    }
	    else{
	      rooResGammaJetVar_p[centPos][sysI]->Fill(varValReco, recoGammaPt_, varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	      rooResGammaJetVarMatrixReco_p[centPos][sysI]->Fill(varValReco, recoGammaPt_, unfoldWeight_*phoJetWeight);
	      rooResGammaJetVarMatrixTruth_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);

	    if(doReweightVar && sysI == 0) photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(varValReco, recoGammaPt_, unfoldWeight_*phoJetWeight);	      
	    }
 	  }     	
	}
	else if(isMultijet){
	  tL.SetPtEtaPhiM(truthJtPt_[jI], truthJtEta_[jI], truthJtPhi_[jI], 0.0);
	  
	  if(recoJtOutOfBounds) goodUnmatchedTruthJets.push_back(tL);
	  else{
	    goodRecoTruthMatchJets.push_back(tL);
	    
	    tL.SetPtEtaPhiM(recoJtPt_[jI][jtVPos], recoJtEta_[jI], recoJtPhi_[jI], 0.0);
	    goodRecoJets.push_back(tL);
	  }
	}      
      }//end for(Int_t jI = 0; jI < nRecoJt_; ....
          
      if(isMultijet){
	//First do it for unmatched truth jets
	for(unsigned int tI = 0; tI < goodUnmatchedTruthJets.size(); ++tI){
	  TLorentzVector goodTruthJet1 = goodUnmatchedTruthJets[tI];
	  
	  for(unsigned int tI2 = tI+1; tI2 < goodUnmatchedTruthJets.size(); ++tI2){
	    TLorentzVector goodTruthJet2 = goodUnmatchedTruthJets[tI2];	   

	    Double_t varValTruth = getVar(varNameLower, goodTruthJet1, goodTruthJet2, truthGammaTL);
	    Double_t phoJetWeight = 1.0;
	    //Reweight if doReweightVar AND its not the prior variation systematic
	    if(doReweightVar && !isStrSame(systStrVect[sysI], "PRIOR")) phoJetWeight = reweightPhoPtJetVar_p[centPos]->GetBinContent(reweightPhoPtJetVar_p[centPos]->GetXaxis()->FindBin(varValTruth), reweightPhoPtJetVar_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));       
	    
	    Float_t multiJtTruthDR = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodTruthJet2.Eta(), goodTruthJet2.Phi());
	    //Reweight if doReweightDR AND its not the prior variation systematic
	    if(doReweightDR && !isStrSame(systStrVect[sysI], "PRIOR")){
	      phoJetWeight *= reweightPhoPtDRJJ_p[centPos]->GetBinContent(reweightPhoPtDRJJ_p[centPos]->GetXaxis()->FindBin(multiJtTruthDR), reweightPhoPtDRJJ_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));
	    }

	    goodTruthJet2 += goodTruthJet1;
	    
	    //Test it passes the only additional multijet cut - dphi
	    Float_t multiJtTruthDPhi = TMath::Abs(getDPHI(truthGammaPhi_, goodTruthJet2.Phi()));

	    //Phi truth cut, continue
	    if(multiJtTruthDPhi < gammaMultiJtDPhiCut) continue;
	    if(multiJtTruthDR < mixJtDRExclusionCut) continue;
	    if(varValTruth < varBinsLow || varValTruth >= varBinsHigh) continue;	  
	    
	    rooResGammaJetVar_p[centPos][sysI]->Miss(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMisses_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	  }
	}
      
	//Now process matched jets...
	if(sysI == 0 && false) std::cout << "N GOOD RECOJETS: " << goodRecoJets.size() << std::endl;

	for(unsigned int jI = 0; jI < goodRecoJets.size(); ++jI){
	  TLorentzVector goodTruthJet1 = goodRecoTruthMatchJets[jI];
	  TLorentzVector goodRecoJet1 = goodRecoJets[jI];
	  
	  //First w/ unmatched jets
	  for(unsigned int tI = 0; tI < goodUnmatchedTruthJets.size(); ++tI){
	    TLorentzVector goodTruthJet2 = goodUnmatchedTruthJets[tI];

	    Double_t varValTruth = getVar(varNameLower, goodTruthJet1, goodTruthJet2, truthGammaTL);
	    //Reweight if doReweightVar AND its not the prior variation systematic
	    Double_t phoJetWeight = 1.0;
	    if(doReweightVar && !isStrSame(systStrVect[sysI], "PRIOR")) phoJetWeight = reweightPhoPtJetVar_p[centPos]->GetBinContent(reweightPhoPtJetVar_p[centPos]->GetXaxis()->FindBin(varValTruth), reweightPhoPtJetVar_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));       
	    
	    Float_t multiJtTruthDR = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodTruthJet2.Eta(), goodTruthJet2.Phi());

	    //Reweight if doReweightDR AND its not the prior variation systematic
	    if(doReweightDR && !isStrSame(systStrVect[sysI], "PRIOR")){
	      phoJetWeight *= reweightPhoPtDRJJ_p[centPos]->GetBinContent(reweightPhoPtDRJJ_p[centPos]->GetXaxis()->FindBin(multiJtTruthDR), reweightPhoPtDRJJ_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));
	    }
	    
	    goodTruthJet2 += goodTruthJet1;
	    
	    Float_t multiJtTruthDPhi = TMath::Abs(getDPHI(truthGammaPhi_, goodTruthJet2.Phi()));
	    
	    if(multiJtTruthDPhi < gammaMultiJtDPhiCut) continue;
	    if(multiJtTruthDR < mixJtDRExclusionCut) continue;
	    if(varValTruth < varBinsLow || varValTruth >= varBinsHigh) continue;		    
	    
	    rooResGammaJetVar_p[centPos][sysI]->Miss(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	    rooResGammaJetVarMisses_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);	  	   
	  }

	  //Then with matched jets
	  for(unsigned int jI2 = jI+1; jI2 < goodRecoJets.size(); ++jI2){
	    Float_t varValReco = -9999.0;
	    
	    TLorentzVector goodRecoJet2 = goodRecoJets[jI2];
	    TLorentzVector goodTruthJet2 = goodRecoTruthMatchJets[jI2];

	    Double_t varValTruth = getVar(varNameLower, goodTruthJet1, goodTruthJet2, truthGammaTL);
	    Double_t phoJetWeight = 1.0;
	    //Reweight if doReweightVar AND its not the prior variation systematic
	    if(doReweightVar && !isStrSame(systStrVect[sysI], "PRIOR")) phoJetWeight = reweightPhoPtJetVar_p[centPos]->GetBinContent(reweightPhoPtJetVar_p[centPos]->GetXaxis()->FindBin(varValTruth), reweightPhoPtJetVar_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));       
	    
	    if(isStrSame(varNameLower, "ajj")) varValReco = TMath::Abs(goodRecoJet1.Pt() - goodRecoJet2.Pt())/recoGammaPt_;
	    else if(isStrSame(varNameLower, "dphijj")) varValReco = TMath::Abs(getDPHI(goodRecoJet1.Phi(), goodRecoJet2.Phi()));
	    else if(isStrSame(varNameLower, "drjj")) varValReco = getDR(goodRecoJet1.Eta(), goodRecoJet1.Phi(), goodRecoJet2.Eta(), goodRecoJet2.Phi());
	    
	    Float_t multiJtTruthDR = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodTruthJet2.Eta(), goodTruthJet2.Phi());
	    //Reweight if doReweightDR AND its not the prior variation systematic
	    if(doReweightDR && !isStrSame(systStrVect[sysI], "PRIOR")){
	      phoJetWeight *= reweightPhoPtDRJJ_p[centPos]->GetBinContent(reweightPhoPtDRJJ_p[centPos]->GetXaxis()->FindBin(multiJtTruthDR), reweightPhoPtDRJJ_p[centPos]->GetYaxis()->FindBin(truthGammaPt_));
	    }

	    Float_t multiJtRecoDR = getDR(goodRecoJet1.Eta(), goodRecoJet1.Phi(), goodRecoJet2.Eta(), goodRecoJet2.Phi());	  
	    Float_t multiJtRecoAJJ = TMath::Abs(goodRecoJet1.Pt() - goodRecoJet2.Pt())/recoGammaPt_;
	    goodTruthJet2 += goodTruthJet1;	  
	    goodRecoJet2 += goodRecoJet1;
	    
	    Float_t multiJtRecoDPhi = TMath::Abs(getDPHI(goodRecoJet2.Phi(), recoGammaPhi_));
	    Float_t multiJtTruthDPhi = TMath::Abs(getDPHI(truthGammaPhi_, goodTruthJet2.Phi()));
	    
	    if(sysI == 0 && false){
	      std::cout << "XJJ FILL entry " << entry << std::endl;
	      std::cout << " Reco Photon pt, eta, phi: " << recoGammaPt_ << ", " << recoGammaEta_ << ", " << recoGammaPhi_ << std::endl;
	      std::cout << " Good Reco Jets 1 pt, eta, phi: " << goodRecoJets[jI].Pt() << ", " << goodRecoJets[jI].Eta() << ", " << goodRecoJets[jI].Phi() << std::endl;
	      std::cout << " Good Reco Jets 2 pt, eta, phi: " << goodRecoJets[jI2].Pt() << ", " << goodRecoJets[jI2].Eta() << ", " << goodRecoJets[jI2].Phi() << std::endl;
	      std::cout << " composite jet pt, eta, phi: " << goodRecoJet2.Pt() << ", " << goodRecoJet2.Eta() << ", " << goodRecoJet2.Phi() << std::endl;
	      std::cout << " DPhi: " << multiJtRecoDPhi << std::endl;
	      
	      std::cout << std::endl;
	      
	      std::cout << " Truth Photon pt, eta, phi: " << truthGammaPt_ << ", " << truthGammaEta_ << ", " << truthGammaPhi_ << std::endl;
	      std::cout << " Good Truth Jets 1 pt, eta, phi: " << goodRecoTruthMatchJets[jI].Pt() << ", " << goodRecoTruthMatchJets[jI].Eta() << ", " << goodRecoTruthMatchJets[jI].Phi() << std::endl;
	      std::cout << " Good Truth Jets 2 pt, eta, phi: " << goodRecoTruthMatchJets[jI2].Pt() << ", " << goodRecoTruthMatchJets[jI2].Eta() << ", " << goodRecoTruthMatchJets[jI2].Phi() << std::endl;	      
	      std::cout << " DPhi: " << multiJtTruthDPhi << std::endl;
	    }	    

	    
	    if(isStrSame(varNameLower, "xjj")) varValReco = goodRecoJet2.Pt()/recoGammaPt_;
	    else if(isStrSame(varNameLower, "dphijjg")) varValReco = TMath::Abs(getDPHI(goodRecoJet2.Phi(), recoGammaPhi_));

	    if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    if(varValTruth < varBinsLow || varValTruth >= varBinsHigh) continue;
	  
	    if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    //Continue if the truth is invalid
	    if(multiJtTruthDPhi < gammaMultiJtDPhiCut) continue;
	    if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    if(multiJtTruthDR < mixJtDRExclusionCut) continue;
	    if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    //	  Float_t varValReco = goodRecoJet2.Pt()/recoGammaPt_;
	    Bool_t varValRecoGood = varValReco >= varBinsLowReco && varValReco < varBinsHighReco;
	    if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    if(multiJtRecoDR < mixJtDRExclusionCut || multiJtRecoDPhi < gammaMultiJtDPhiCut || recoGammaOutOfBounds || !varValRecoGood){
	      if(sysI == 0 && false) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      if(sysI == 0 && false) std::cout << " FILLING AS MISS" << std::endl;
	      rooResGammaJetVar_p[centPos][sysI]->Miss(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	      rooResGammaJetVarMisses_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	      
	      if(sysI == 0 && false){
		Float_t tempRecoDR = getDR(goodRecoJet1.Eta(), goodRecoJet1.Phi(), goodRecoJets[jI2].Eta(), goodRecoJets[jI2].Phi());
		Float_t tempTruthDR = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodRecoTruthMatchJets[jI2].Eta(), goodRecoTruthMatchJets[jI2].Phi());
	      }
	    }
	    else{
	      if(sysI == 0 && false) std::cout << " FILLING AS HIT" << std::endl;
	      
	      Float_t tempDRReco = getDR(goodRecoJets[jI].Eta(), goodRecoJets[jI].Phi(), goodRecoJets[jI2].Eta(), goodRecoJets[jI2].Phi());
	      Float_t tempDRTruth = getDR(goodTruthJet1.Eta(), goodTruthJet1.Phi(), goodTruthJet2.Eta(), goodTruthJet2.Phi());
	      
	      rooResGammaJetVar_p[centPos][sysI]->Fill(varValReco, recoGammaPt_, varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);	    
	      
	      rooResGammaJetVarMatrixTruth_p[centPos][sysI]->Fill(varValTruth, truthGammaPt_, unfoldWeight_*phoJetWeight);
	      rooResGammaJetVarMatrixReco_p[centPos][sysI]->Fill(varValReco, recoGammaPt_, unfoldWeight_*phoJetWeight);

	      if(doReweightVar && sysI == 0){
		photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(varValReco, recoGammaPt_, unfoldWeight_*phoJetWeight);
		photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(multiJtRecoDR, recoGammaPt_, unfoldWeight_*phoJetWeight);	      
		photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[centPos]->Fill(multiJtRecoAJJ, recoGammaPt_, unfoldWeight_*phoJetWeight);	      
	      }
	    }		  
	  }
	}//end goodreco jets loop      
      }//if(isStrSame("xjj", varNameLower){ ending
    }    
  }//end nEntries loop


  if(doReweightVar){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      //1-D case is easy - just declare a saveName plug into function
      std::string saveName = pdfStr + "/photonPtReco_Reweight_" + centBinsStr[cI] + "_" + saveTag + "_" + dateStr + ".png";
      std::vector<std::string> labels = {"#bf{#it{ATLAS Internal}}", centBinsStr[cI]};
      
      //      plotReweightHists(photonPtReco_PURCORR_COMBINED_ForReweight_p[cI], photonPtReco_PURCORR_COMBINED_p[cI], photonPtReco_PURCORR_COMBINED_Reweighted_p[cI], saveName, labels);

      //2-D case, do projections along gamma-pt (y-axis)
      //Iterate over photon pTBins
      for(Int_t gI = 0; gI < nGammaPtBins; ++gI){	
	Float_t binCenter = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;
	if(binCenter < gammaPtBinsLowReco) continue;
	if(binCenter > gammaPtBinsHighReco) continue;

	std::string xTitle = photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->GetXaxis()->GetTitle();
	std::string titleStr = ";" + xTitle + ";Shape (Bin-width norm.)";

	TH1D* histMCForReweight_p = new TH1D("histMCForReweight_h", titleStr.c_str(), nVarBins, varBins);
	TH1D* histData_p = new TH1D("histData_h", titleStr.c_str(), nVarBins, varBins);
	TH1D* histMCReweighted_p = new TH1D("histMCReweighted_h", titleStr.c_str(), nVarBins, varBins);

	std::vector<std::string> tempLabels = labels;
	tempLabels.push_back("#bf{" + prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBins[gI+1], 1, false) + "}");
	
	//Set content of each 1D histogram
	for(Int_t bIX = 0; bIX < nVarBins; ++bIX){
	  Float_t content = photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinContent(bIX+1, gI+1);
	  Float_t error = photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCForReweight_p->SetBinContent(bIX+1, content);
	  histMCForReweight_p->SetBinError(bIX+1, error);

	  content = photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histData_p->SetBinContent(bIX+1, content);
	  histData_p->SetBinError(bIX+1, error);

	  content = photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCReweighted_p->SetBinContent(bIX+1, content);
	  histMCReweighted_p->SetBinError(bIX+1, error);
	}
	//Run reweight plotter and clean
	saveName = pdfStr + "/photonPtJt" + varName + "Reco_Reweight_" + centBinsStr[cI] + "_GammaPt" + std::to_string(gI) + "_" + saveTag + "_" + dateStr + ".png";
	//	plotReweightHists(histMCForReweight_p, histData_p, histMCReweighted_p, saveName, tempLabels);
	
	delete histMCForReweight_p;
	delete histData_p;
	delete histMCReweighted_p;

	//repeat for drjj
	xTitle = photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetXaxis()->GetTitle();
	titleStr = ";" + xTitle + ";Shape (Bin-width norm.)";

	histMCForReweight_p = new TH1D("histMCForReweight_h", titleStr.c_str(), nDRBins, drBins);
	histData_p = new TH1D("histData_h", titleStr.c_str(), nDRBins, drBins);
	histMCReweighted_p = new TH1D("histMCReweighted_h", titleStr.c_str(), nDRBins, drBins);

	//Set content of each 1D histogram
	for(Int_t bIX = 0; bIX < nDRBins; ++bIX){
	  Float_t content = photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinContent(bIX+1, gI+1);
	  Float_t error = photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCForReweight_p->SetBinContent(bIX+1, content);
	  histMCForReweight_p->SetBinError(bIX+1, error);

	  content = photonPtDRJJReco_PURCORR_COMBINED_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtDRJJReco_PURCORR_COMBINED_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histData_p->SetBinContent(bIX+1, content);
	  histData_p->SetBinError(bIX+1, error);

	  content = photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCReweighted_p->SetBinContent(bIX+1, content);
	  histMCReweighted_p->SetBinError(bIX+1, error);
	}
	//Run reweight plotter and clean
	saveName = pdfStr + "/photonPtJtDRJJReco_Reweight_" + centBinsStr[cI] + "_GammaPt" + std::to_string(gI) + "_" + saveTag + "_" + dateStr + ".png";
	//	plotReweightHists(histMCForReweight_p, histData_p, histMCReweighted_p, saveName, tempLabels);
	
	delete histMCForReweight_p;
	delete histData_p;
	delete histMCReweighted_p;

	//repeat for ajj
	xTitle = photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetXaxis()->GetTitle();
	titleStr = ";" + xTitle + ";Shape (Bin-width norm.)";

	histMCForReweight_p = new TH1D("histMCForReweight_h", titleStr.c_str(), nAJBins, ajBins);
	histData_p = new TH1D("histData_h", titleStr.c_str(), nAJBins, ajBins);
	histMCReweighted_p = new TH1D("histMCReweighted_h", titleStr.c_str(), nAJBins, ajBins);

	//Set content of each 1D histogram
	for(Int_t bIX = 0; bIX < nAJBins; ++bIX){
	  Float_t content = photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinContent(bIX+1, gI+1);
	  Float_t error = photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCForReweight_p->SetBinContent(bIX+1, content);
	  histMCForReweight_p->SetBinError(bIX+1, error);

	  content = photonPtAJJReco_PURCORR_COMBINED_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtAJJReco_PURCORR_COMBINED_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histData_p->SetBinContent(bIX+1, content);
	  histData_p->SetBinError(bIX+1, error);

	  content = photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinContent(bIX+1, gI+1);
	  error = photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[cI]->GetBinError(bIX+1, gI+1);
	    
	  histMCReweighted_p->SetBinContent(bIX+1, content);
	  histMCReweighted_p->SetBinError(bIX+1, error);
	}
	//Run reweight plotter and clean
	saveName = pdfStr + "/photonPtJtAJJReco_Reweight_" + centBinsStr[cI] + "_GammaPt" + std::to_string(gI) + "_" + saveTag + "_" + dateStr + ".png";
	//	plotReweightHists(histMCForReweight_p, histData_p, histMCReweighted_p, saveName, tempLabels);      
	
	delete histMCForReweight_p;
	delete histData_p;
	delete histMCReweighted_p;
      }
    }    
  }
  
  std::cout << "Construction of response matrices is complete."   << std::endl;
  
  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtReco_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);
    
    for(int gI = 0; gI < nGammaPtBins; ++gI){
      TH1D* xjRecoProjection_PURCORR_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);
      TH1D* xjRecoProjection_TRUTH_COMBINED_p = new TH1D((varNameLower + "Reco_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), (";" + varNameStyle + ";Arbitrary").c_str(), nVarBins, varBins);

      for(int xI = 0; xI < nVarBins; ++xI){
	xjRecoProjection_PURCORR_COMBINED_p->SetBinContent(xI+1, photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_PURCORR_COMBINED_p->SetBinError(xI+1, photonPtJetVarReco_PURCORR_COMBINED_p[cI]->GetBinError(xI+1, gI+1));

	xjRecoProjection_TRUTH_COMBINED_p->SetBinContent(xI+1, photonPtJetVarReco_TRUTH_COMBINED_p[cI]->GetBinContent(xI+1, gI+1));
	xjRecoProjection_TRUTH_COMBINED_p->SetBinError(xI+1, photonPtJetVarReco_TRUTH_COMBINED_p[cI]->GetBinError(xI+1, gI+1));
      }
      
      xjRecoProjection_PURCORR_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_PURCORR_COMBINED_p;

      xjRecoProjection_TRUTH_COMBINED_p->Write("", TObject::kOverwrite);
      delete xjRecoProjection_TRUTH_COMBINED_p;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
      
    //Write the full histograms at reco level to file
    photonPtReco_PURCORR_COMBINED_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    if(doReweightVar){
      photonPtReco_PURCORR_COMBINED_ForReweight_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str(), TObject::kOverwrite);
      photonPtReco_PURCORR_COMBINED_Reweighted_p[cI]->Write(("photonPtReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str(), TObject::kOverwrite);
    }

    photonPtReco_TRUTH_COMBINED_p[cI]->Write(("photonPtTruth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    photonPtJetVarReco_PURCORR_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
    if(doReweightVar){
      photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI]->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str(), TObject::kOverwrite);
      photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[cI]->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str(), TObject::kOverwrite);

      photonPtDRJJReco_PURCORR_COMBINED_ForReweight_p[cI]->Write(("photonPtJetDRJJReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str(), TObject::kOverwrite);
      photonPtDRJJReco_PURCORR_COMBINED_Reweighted_p[cI]->Write(("photonPtJetDRJJReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str(), TObject::kOverwrite);

      photonPtAJJReco_PURCORR_COMBINED_ForReweight_p[cI]->Write(("photonPtJetAJJReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_ForReweight_h").c_str(), TObject::kOverwrite);
      photonPtAJJReco_PURCORR_COMBINED_Reweighted_p[cI]->Write(("photonPtJetAJJReco_PreUnfold_" + centBinsStr[cI] + "_PURCORR_COMBINED_Reweighted_h").c_str(), TObject::kOverwrite);
    }
    
    photonPtJetVarReco_TRUTH_COMBINED_p[cI]->Write(("photonPtJet" + varName + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH_COMBINED_h").c_str(), TObject::kOverwrite);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << cI << "/" << centBinsStr.size() << std::endl;  
    
    TH2D* tempTruthHist_p = (TH2D*)rooResGammaJetVar_p[cI][0]->Htruth();
    tempTruthHist_p->Write(("photonPtJet" + varName + "Truth_PreUnfold_" + centBinsStr[cI] + "_TRUTH2_COMBINED_h").c_str(), TObject::kOverwrite); 

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
    
    TH2D* tempRecoHist_p = (TH2D*)rooResGammaJetVar_p[cI][0]->Hmeasured();
    tempRecoHist_p->Write(("photonPtJet" + varName + "Reco_PreUnfold_" + centBinsStr[cI] + "_PURCORR2_COMBINED_h").c_str(), TObject::kOverwrite);     

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
    
    rooResGammaMatrix_p[cI]->Write("", TObject::kOverwrite);
    rooResGammaMisses_p[cI]->Write("", TObject::kOverwrite);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
    for(Int_t sysI = 0; sysI < nSyst; ++sysI){
      if(!isStrSame(inUnfoldName, systStrVect[sysI])){
	if(!isStrSame(inUnfoldName, "ALL")) continue;
      }

      rooResGammaJetVarMatrixReco_p[cI][sysI]->Write("", TObject::kOverwrite);
      rooResGammaJetVarMatrixTruth_p[cI][sysI]->Write("", TObject::kOverwrite);
      rooResGammaJetVarMisses_p[cI][sysI]->Write("", TObject::kOverwrite);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
    
    //Full unfold gamma only
    std::cout << "Photon unfold, " << centBinsStr[cI] << "..." << std::endl;
    for(int i = 1; i <= nIter; ++ i){
      std::cout << " Iter " << i << "/" << nIter << "..." << std::endl;
      RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGamma_p[cI], photonPtReco_PURCORR_COMBINED_p[cI], i);
      rooBayes_p->SetVerbose(0);
      
      TH1D* unfolded_p = nullptr;
      if(unfoldErrType == 0) unfolded_p = (TH1D*)rooBayes_p->Hreco()->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      else if(unfoldErrType == 1){
	rooBayes_p->SetNToys(nToys);
	unfolded_p = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str());
      }
      
      TH1D* refolded_p = (TH1D*)rooResGamma_p[cI]->ApplyToTruth(unfolded_p)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str());

      unfolded_p->Write(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
      refolded_p->Write(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_PURCORR_COMBINED_Refolded_h").c_str(), TObject::kOverwrite);

      delete unfolded_p;
      delete refolded_p;
      delete rooBayes_p;
    }
    std::cout << "Photon unfold, " << centBinsStr[cI] << ", complete." << std::endl;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
    
    //Gamma-jet unfold
    std::cout << "Photon-jet unfold, " << centBinsStr[cI] << "..." << std::endl;
    for(Int_t sysI = 0; sysI < nSyst; ++sysI){
      if(!isStrSame(inUnfoldName, systStrVect[sysI])){
	if(!isStrSame(inUnfoldName, "ALL")) continue;
      }

      std::cout << " " << systStrVect[sysI] << std::endl;
      for(int i = 1; i <= nIter; ++ i){
	std::cout << "  Iter " << i << "/" << nIter << std::endl;

	RooUnfoldBayes* rooBayes_p = new RooUnfoldBayes(rooResGammaJetVar_p[cI][sysI], photonPtJetVarReco_PURCORR_COMBINED_p[cI], i);
	rooBayes_p->SetVerbose(0);
	
	TH2D* unfolded_p = nullptr;
	if(unfoldErrType == 0) unfolded_p = (TH2D*)rooBayes_p->Hreco()->Clone(("photonPtJetVarReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_" + systStrVect[sysI] + "_PURCORR_COMBINED_h").c_str());
	else if(unfoldErrType == 1){
	  rooBayes_p->SetNToys(nToys);
	  unfolded_p = (TH2D*)rooBayes_p->Hreco(RooUnfold::kCovToy)->Clone(("photonPtJetVarReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_" + systStrVect[sysI] + "_PURCORR_COMBINED_h").c_str());
	}
		
	TH2D* refolded_p = (TH2D*)rooResGammaJetVar_p[cI][sysI]->ApplyToTruth(unfolded_p)->Clone(("photonPtReco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_" + systStrVect[sysI] + "_PURCORR_COMBINED_Refolded_h").c_str());
	
	unfolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_" + systStrVect[sysI] + "_PURCORR_COMBINED_h").c_str(), TObject::kOverwrite);
	refolded_p->Write(("photonPtJet" + varName + "Reco_Iter" + std::to_string(i) + "_" + centBinsStr[cI] + "_" + systStrVect[sysI] + "_PURCORR_COMBINED_Refolded_h").c_str(), TObject::kOverwrite);
	
	delete unfolded_p;
	delete refolded_p;
	delete rooBayes_p;
      }
    }
    std::cout << "Photon-jet unfold, " << centBinsStr[cI] << ", complete." << std::endl;
  }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    delete photonPtReco_p[cI];
    delete photonPtTruth_p[cI];

    delete photonPtJetVarReco_p[cI];
    delete photonPtJetVarTruth_p[cI];

    if(doRebin){
      for(Int_t eI = 0; eI < nBarrelAndEC; ++eI){
	delete photonPtReco_PURCORR_p[cI][eI];
	delete photonPtJetVarReco_PURCORR_p[cI][eI];
	delete photonPtReco_TRUTH_p[cI][eI];
	delete photonPtJetVarReco_TRUTH_p[cI][eI];
      }
    }

    delete photonPtReco_PURCORR_COMBINED_p[cI];
    delete photonPtJetVarReco_PURCORR_COMBINED_p[cI];
    if(doReweightVar){
      delete photonPtReco_PURCORR_COMBINED_ForReweight_p[cI];
      delete photonPtReco_PURCORR_COMBINED_Reweighted_p[cI];

      delete photonPtJetVarReco_PURCORR_COMBINED_ForReweight_p[cI];
      delete photonPtJetVarReco_PURCORR_COMBINED_Reweighted_p[cI];
    }
    
    delete photonPtReco_TRUTH_COMBINED_p[cI];
    delete photonPtJetVarReco_TRUTH_COMBINED_p[cI];

    delete rooResGamma_p[cI];
    delete rooResGammaMatrix_p[cI];
    delete rooResGammaMisses_p[cI];

    for(Int_t sysI = 0; sysI < nSyst; ++sysI){
      if(!isStrSame(inUnfoldName, systStrVect[sysI])){
	if(!isStrSame(inUnfoldName, "ALL")) continue;
      }

      delete rooResGammaJetVar_p[cI][sysI];
      delete rooResGammaJetVarMatrixReco_p[cI][sysI];
      delete rooResGammaJetVarMatrixTruth_p[cI][sysI];
      delete rooResGammaJetVarMisses_p[cI][sysI];
    }
    
    if(!isSameInputFile && doReweightVar){
      delete reweightPhoPtJetVar_p[cI];
    }
  }

  if(doRebin){
    inUnfoldFileConfig_p->SetValue("NGAMMAPTBINS", nGammaPtBins);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSLOW", gammaPtBinsLow);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSHIGH", gammaPtBinsHigh);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSDOLOG", gammaPtBinsDoLog);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSDOCUSTOM", gammaPtBinsDoCustom);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSCUSTOM", gammaPtBinsStr.c_str());
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSLOWRECO", gammaPtBinsLowReco);
    inUnfoldFileConfig_p->SetValue("GAMMAPTBINSHIGHRECO", gammaPtBinsHighReco);

    std::string varNameUpper = returnAllCapsString(varName);
    std::string varPrefix = "";
    if(isStrSame(varNameUpper, "PT")) varPrefix = "JT";

    std::string binVarStr = varPrefix + varNameUpper;
    if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";

    inUnfoldFileConfig_p->SetValue(("N" + binVarStr + "BINS").c_str(), nVarBins);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSLOW").c_str(), varBinsLow);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSHIGH").c_str(), varBinsHigh);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSDOLOG").c_str(), 0);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSDOCUSTOM").c_str(), 1);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSCUSTOM").c_str(), varBinsStr.c_str());
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSLOWRECO").c_str(), varBinsLowReco);
    inUnfoldFileConfig_p->SetValue((binVarStr + "BINSHIGHRECO").c_str(), varBinsHighReco);
  }
  
  inUnfoldFileConfig_p->SetValue("VARNAME", varName.c_str());
  inUnfoldFileConfig_p->SetValue("NITER", nIter);
  inUnfoldFileConfig_p->SetValue("UNFOLD", inUnfoldName.c_str());
  inUnfoldFileConfig_p->SetValue("DOREWEIGHTVAR", doReweightVar);
  inUnfoldFileConfig_p->SetValue("SAVETAG", saveTag.c_str());

  inUnfoldFileConfig_p->Write("config", TObject::kOverwrite);
  inUnfoldFileLabels_p->Write("label", TObject::kOverwrite);

  inResponseFile_p->Close();
  delete inResponseFile_p;  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  if(!isSameInputFile){
    inUnfoldFile_p->Close();
    delete inUnfoldFile_p;
  }

  outFile_p->Close();
  delete outFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  
  
  std::cout << "GDJHISTTOUNFOLD COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHistToUnfold.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjHistToUnfold(argv[1]);
  return retVal;
}
