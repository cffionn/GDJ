#ifndef BINUTILS_H
#define BINUTILS_H

//c+cpp
#include <cmath>
#include <iostream>
#include <string>

//Local
#include "include/stringUtil.h"

#define _USE_MATH_DEFINES

inline bool goodBinning(std::string config, int nMaxBins, int nBins, std::string binType)
{
  if(nBins > nMaxBins){
    std::cout << "BINUTILS ERROR IN GOODBINNING - config \'" << config << "\' contains " << binType << "  of size \'" << nBins << "', greater than max val \'" << nMaxBins << "\'. Please fix. return false" << std::endl;
    return false;
  }
  return true;
}

inline double mathStringToNum(std::string inStr, const bool doGlobalDebug=false)//currently just handles fractions and pi
{
  double outVal = 1.0;
  inStr = returnAllCapsString(inStr);
  if(inStr.find("PI") != std::string::npos){
    outVal *= M_PI;
    inStr.replace(inStr.find("PI"), 2, "");
  }
  if(inStr.find("/") != std::string::npos){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << inStr << std::endl;

    double num = 1.0;
    if(inStr.find("/") != 0) num = std::stod(inStr.substr(0, inStr.find("/")));

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << inStr << std::endl;

    double denom = std::stod(inStr.substr(inStr.find("/")+1, inStr.size()));
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    outVal *= num/denom;
  }
  else outVal *= std::stod(inStr);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
  return outVal;
}

inline std::string mathStringToNameStr(std::string inStr, const bool doGlobalDebug=false)//currently just handles fractions and pi
{
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inStr = returnAllCapsString(inStr);
  while(inStr.find("/") != std::string::npos){
    inStr.replace(inStr.find("/"), 1, "Over");
  }
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  while(inStr.find("*") != std::string::npos){
    inStr.replace(inStr.find("*"), 1, "");
  }
    
  return inStr;
}

inline void binWidthNorm(TH1* inHist_p)
{
  for(Int_t i = 0; i < inHist_p->GetNbinsX(); ++i){
    Double_t binContent = inHist_p->GetBinContent(i+1);
    Double_t binError = inHist_p->GetBinError(i+1);
    Double_t binWidth = inHist_p->GetBinWidth(i+1);

    inHist_p->SetBinContent(i+1, binContent/binWidth);
    inHist_p->SetBinError(i+1, binError/binWidth);
  }
  return;
}

inline void binWidthAndScaleNorm(TH1* inHist_p, double scale)
{
  binWidthNorm(inHist_p);
  inHist_p->Scale(1./scale);

  return;
}

inline void binWidthAndSelfNorm(TH1* inHist_p)
{
  binWidthAndScaleNorm(inHist_p, inHist_p->Integral());

  return;
}

void fineHistToCoarseHist(TH2F* fineHist_p, TH2F* coarseHist_p)
{
  //Initialize coarse histogram to 0                                                                  
  for(Int_t bIX = 0; bIX < coarseHist_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < coarseHist_p->GetYaxis()->GetNbins(); ++bIY){
      coarseHist_p->SetBinContent(bIX+1, bIY+1, 0.0);
      coarseHist_p->SetBinError(bIX+1, bIY+1, 0.0);
    }
  }

  //Declare our bin maps
  std::map<Int_t, Int_t> fineToCoarseXs, fineToCoarseYs;

  //Build out our x-map                                                                               
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Float_t binCenterFine = fineHist_p->GetXaxis()->GetBinCenter(bIX+1);
    Int_t fineToCoarseXPos = -1;

    for(Int_t bIX2 = 0; bIX2 < coarseHist_p->GetXaxis()->GetNbins(); ++bIX2){
      Float_t binLowCoarse = coarseHist_p->GetXaxis()->GetBinLowEdge(bIX2+1);
      Float_t binHighCoarse = coarseHist_p->GetXaxis()->GetBinLowEdge(bIX2+2);

      if(binCenterFine >= binLowCoarse && binCenterFine < binHighCoarse){
        fineToCoarseXPos = bIX2;
        break;
      }
    }

    if(fineToCoarseXPos < 0){
      std::cout << "FINEHISTTOCOARSEHIST WARNING: binCenterFine '" << binCenterFine << "' has no matc\
h in coarser binning. returning"  << std::endl;
      return;
    }

    fineToCoarseXs[bIX] = fineToCoarseXPos;
  }

  //Build out our y-map
  for(Int_t bIY = 0; bIY < fineHist_p->GetYaxis()->GetNbins(); ++bIY){
    Float_t binCenterFine = fineHist_p->GetYaxis()->GetBinCenter(bIY+1);
    Int_t fineToCoarseYPos = -1;

    for(Int_t bIY2 = 0; bIY2 < coarseHist_p->GetYaxis()->GetNbins(); ++bIY2){
      Float_t binLowCoarse = coarseHist_p->GetYaxis()->GetBinLowEdge(bIY2+1);
      Float_t binHighCoarse = coarseHist_p->GetYaxis()->GetBinLowEdge(bIY2+2);

      if(binCenterFine >= binLowCoarse && binCenterFine < binHighCoarse){
        fineToCoarseYPos = bIY2;
        break;
      }
    }

    if(fineToCoarseYPos < 0){
      std::cout << "FINEHISTTOCOARSEHIST WARNING: binCenterFine '" << binCenterFine << "' has no matc\
h in coarser binning. returning"  << std::endl;
      return;
    }

    fineToCoarseYs[bIY] = fineToCoarseYPos;
  }

  //Now fill out coarse from fine
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Int_t coarsePosX = fineToCoarseXs[bIX];

    for(Int_t bIY = 0; bIY < fineHist_p->GetYaxis()->GetNbins(); ++bIY){
      Int_t coarsePosY = fineToCoarseYs[bIY];

      Float_t newValue = coarseHist_p->GetBinContent(coarsePosX+1, coarsePosY+1) + fineHist_p->GetBinContent(bIX+1, bIY+1);

      Float_t newError = coarseHist_p->GetBinError(coarsePosX+1, coarsePosY+1);
      newError = TMath::Sqrt(newError*newError + fineHist_p->GetBinError(bIX+1, bIY+1)*fineHist_p->GetBinError(bIX+1, bIY+1));

      coarseHist_p->SetBinContent(coarsePosX+1, coarsePosY+1, newValue);
      coarseHist_p->SetBinError(coarsePosX+1, coarsePosY+1, newError);
    }
  }

  return;
}

bool checkMinBinWidth(std::vector<float> bins, float minBinWidth)
{
  bool allBinsGood = true;
  for(int bI = 0; bI < ((int)bins.size())-1; ++bI){
    if(TMath::Abs(bins[bI] - bins[bI+1]) < minBinWidth){
      std::cout << "CHECKMINBINWIDTH ERROR: Adjacent bins given '" << bins[bI] << "', '" << bins[bI+1] << "' have minimum distance '" << TMath::Abs(bins[bI] - bins[bI+1]) << "' less than required '" << minBinWidth << "'. return false" << std::endl;
      allBinsGood = false;
    }
  }
  return allBinsGood;
}

bool checkHistContainsBins(std::vector<float> bins, TH2F* hist_p, float deltaValue, bool isYAxis)
{
  bool allBinsFound = true;
  std::string xyStr = "X";
  TAxis* axis_p = (TAxis*)hist_p->GetXaxis();
  if(isYAxis){
    xyStr = "Y";
    axis_p = (TAxis*)hist_p->GetYaxis();
  }

  for(unsigned int xI = 0; xI < bins.size(); ++xI){
    bool binFound = false;
    for(Int_t bIX = 0; bIX < axis_p->GetNbins()+1; ++bIX){
      Float_t binLowEdge = axis_p->GetBinLowEdge(bIX+1);

      if(TMath::Abs(binLowEdge - bins[xI]) < deltaValue){
        binFound = true;
        break;
      }
    }

    if(!binFound){
      std::cout << "CHECKHISTCONTAINSBINS ERROR: Requested bin '" << bins[xI] << "' is not found in histogram, " << hist_p->GetName() << ", " << xyStr << "-axis. return false" << std::endl;
      allBinsFound = false;
    }
  }

  if(!allBinsFound){
    std::cout << "CHECKHISTCONTAINSBINS ERROR: Bin not found" << std::endl;
    std::cout << " Bins requested:" << std::endl;
    std::string reqBinsStr = "";
    for(unsigned int bI = 0; bI < bins.size()-1; ++bI){
      reqBinsStr = reqBinsStr + bins[bI] + ", ";
    }
    reqBinsStr = reqBinsStr + bins[bins.size()-1] + ".";
    std::cout << "  " << reqBinsStr << std::endl;
    std::cout << " Hist bins:" << std::endl;
    reqBinsStr = "";
    for(unsigned int bI = 0; bI < axis_p->GetNbins(); ++bI){
      reqBinsStr = reqBinsStr + bins[bI] + ", ";
    }
    reqBinsStr = reqBinsStr + bins[axis_p->GetNbins()] + ".";
    std::cout << "  " << reqBinsStr << std::endl;

    std::cout << "return false" << std::endl;
  }

  return allBinsFound;
}


#endif
