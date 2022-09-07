#ifndef BINUTILS_H
#define BINUTILS_H

//c+cpp
#include <cmath>
#include <iostream>
#include <string>

//ROOT
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

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

template <typename T>
void fineHistToCoarseHist(T* fineHist_p, T* coarseHist_p)
{
  std::string histClassName = fineHist_p->ClassName();
  bool isTH2 = histClassName.find("TH2") != std::string::npos;
  
  //Initialize coarse histogram to 0                                                                  
  for(Int_t bIX = 0; bIX < coarseHist_p->GetXaxis()->GetNbins(); ++bIX){

    if(isTH2){
      for(Int_t bIY = 0; bIY < coarseHist_p->GetYaxis()->GetNbins(); ++bIY){
	coarseHist_p->SetBinContent(bIX+1, bIY+1, 0.0);
	coarseHist_p->SetBinError(bIX+1, bIY+1, 0.0);
      }
    }
    else{
      coarseHist_p->SetBinContent(bIX+1, 0.0);
      coarseHist_p->SetBinError(bIX+1, 0.0);
    }
  }

  //Declare our bin maps
  std::map<Int_t, Int_t> fineToCoarseXs, fineToCoarseYs;

  //Build out our x-map                                                                               
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Float_t binCenterFine = fineHist_p->GetXaxis()->GetBinCenter(bIX+1);
    Int_t fineToCoarseXPos = -1;

    bool belowCoarse = binCenterFine < coarseHist_p->GetXaxis()->GetBinLowEdge(1);
    bool aboveCoarse = binCenterFine >= coarseHist_p->GetXaxis()->GetBinLowEdge(coarseHist_p->GetXaxis()->GetNbins()+1);

    if(belowCoarse || aboveCoarse){
      fineToCoarseXs[bIX] = -1;
      continue;
    }
    
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
  if(isTH2){
    for(Int_t bIY = 0; bIY < fineHist_p->GetYaxis()->GetNbins(); ++bIY){
      Float_t binCenterFine = fineHist_p->GetYaxis()->GetBinCenter(bIY+1);
      Int_t fineToCoarseYPos = -1;
      
      bool belowCoarse = binCenterFine < coarseHist_p->GetYaxis()->GetBinLowEdge(1);
      bool aboveCoarse = binCenterFine >= coarseHist_p->GetYaxis()->GetBinLowEdge(coarseHist_p->GetYaxis()->GetNbins()+1);
      if(belowCoarse || aboveCoarse){
	fineToCoarseYs[bIY] = -1;
	continue;
      }
      
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
  }

  for(auto const & fineToCoarseY : fineToCoarseYs){
    std::cout << " " << fineToCoarseY.first << ", " << fineToCoarseY.second << std::endl;
  }
  
  //Now fill out coarse from fine
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Int_t coarsePosX = fineToCoarseXs[bIX];

    if(coarsePosX < 0) continue;

    if(isTH2){
      for(Int_t bIY = 0; bIY < fineHist_p->GetYaxis()->GetNbins(); ++bIY){
	Int_t coarsePosY = fineToCoarseYs[bIY];

	if(coarsePosY < 0) continue;
	
	Float_t newValue = coarseHist_p->GetBinContent(coarsePosX+1, coarsePosY+1) + fineHist_p->GetBinContent(bIX+1, bIY+1);
	
	Float_t newError = coarseHist_p->GetBinError(coarsePosX+1, coarsePosY+1);
	newError = TMath::Sqrt(newError*newError + fineHist_p->GetBinError(bIX+1, bIY+1)*fineHist_p->GetBinError(bIX+1, bIY+1));
	
	coarseHist_p->SetBinContent(coarsePosX+1, coarsePosY+1, newValue);
	coarseHist_p->SetBinError(coarsePosX+1, coarsePosY+1, newError);
      }
    }
    else{
      Float_t newValue = coarseHist_p->GetBinContent(coarsePosX+1) + fineHist_p->GetBinContent(bIX+1);
	
      Float_t newError = coarseHist_p->GetBinError(coarsePosX+1);
      newError = TMath::Sqrt(newError*newError + fineHist_p->GetBinError(bIX+1)*fineHist_p->GetBinError(bIX+1));
	
      coarseHist_p->SetBinContent(coarsePosX+1, newValue);
      coarseHist_p->SetBinError(coarsePosX+1, newError);      
    }
  }
  
  return;
}

inline void fineTH2ToCoarseTH1(TH2* fineHist_p, TH1* coarseHist_p, Int_t yPos)
{
  //Initialize coarse histogram to 0                                                                  
  for(Int_t bIX = 0; bIX < coarseHist_p->GetXaxis()->GetNbins(); ++bIX){
    coarseHist_p->SetBinContent(bIX+1, 0.0);
    coarseHist_p->SetBinError(bIX+1, 0.0);
  }

  //Declare our bin maps
  std::map<Int_t, Int_t> fineToCoarseXs;

  //Build out our x-map                                                                               
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Float_t binCenterFine = fineHist_p->GetXaxis()->GetBinCenter(bIX+1);
    Int_t fineToCoarseXPos = -1;

    bool belowCoarse = binCenterFine < coarseHist_p->GetXaxis()->GetBinLowEdge(1);
    bool aboveCoarse = binCenterFine >= coarseHist_p->GetXaxis()->GetBinLowEdge(coarseHist_p->GetXaxis()->GetNbins()+1);

    if(belowCoarse || aboveCoarse){
      fineToCoarseXs[bIX] = -1;
      continue;
    }
    
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
    
  //Now fill out coarse from fine
  for(Int_t bIX = 0; bIX < fineHist_p->GetXaxis()->GetNbins(); ++bIX){
    Int_t coarsePosX = fineToCoarseXs[bIX];

    if(coarsePosX < 0) continue;

    Float_t newValue = coarseHist_p->GetBinContent(coarsePosX+1) + fineHist_p->GetBinContent(bIX+1, yPos+1);
	
    Float_t newError = coarseHist_p->GetBinError(coarsePosX+1);
    newError = TMath::Sqrt(newError*newError + fineHist_p->GetBinError(bIX+1, yPos+1)*fineHist_p->GetBinError(bIX+1, yPos+1));
    
    coarseHist_p->SetBinContent(coarsePosX+1, newValue);
    coarseHist_p->SetBinError(coarsePosX+1, newError);          
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

template <typename T>
bool checkHistContainsBins(std::vector<float> bins, T* hist_p, float deltaValue, bool isYAxis)
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
    for(Int_t bI = 0; bI < axis_p->GetNbins(); ++bI){
      reqBinsStr = reqBinsStr + std::to_string(axis_p->GetBinLowEdge(bI+1)) + ", ";
    }
    reqBinsStr = reqBinsStr + axis_p->GetBinLowEdge(axis_p->GetNbins()+1) + ".";
    std::cout << "  " << reqBinsStr << std::endl;

    std::cout << "return false" << std::endl;
  }

  return allBinsFound;
}

bool checkHistContainsBins(std::vector<float> bins, TH1F* hist_p, float deltaValue)
{
  bool allBinsFound = true;
  std::string xyStr = "X";
  TAxis* axis_p = (TAxis*)hist_p->GetXaxis();

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
    for(Int_t bI = 0; bI < axis_p->GetNbins(); ++bI){
      reqBinsStr = reqBinsStr + std::to_string(axis_p->GetBinLowEdge(bI+1)) + ", ";
    }
    reqBinsStr = reqBinsStr + axis_p->GetBinLowEdge(axis_p->GetNbins()) + ".";
    std::cout << "  " << reqBinsStr << std::endl;

    std::cout << "return false" << std::endl;
  }

  return allBinsFound;
}

bool checkHistContainsBins(std::vector<float> bins, Int_t nHistBins, Double_t histBins[], float deltaValue)
{
  std::string binsStr = "";
  for(Int_t hI = 0; hI < nHistBins+1; ++hI){
    binsStr = binsStr + std::to_string(histBins[hI]) + ",";
  }
  
  bool allBinsFound = true;  
  for(unsigned int xI = 0; xI < bins.size(); ++xI){
    bool binFound = false;

    for(Int_t bIX = 0; bIX < nHistBins+1; ++bIX){
      Float_t binLowEdge = histBins[bIX];

      if(TMath::Abs(binLowEdge - bins[xI]) < deltaValue){
        binFound = true;
        break;
      }
    }

    if(!binFound){           
      std::cout << "CHECKHISTCONTAINSBINS ERROR: Requested bin '" << bins[xI] << "' is not found in histogram defined by bins \'" << binsStr << "\'. return false" << std::endl;
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
    for(Int_t bI = 0; bI < nHistBins; ++bI){
      reqBinsStr = reqBinsStr + std::to_string(histBins[bI]) + ", ";
    }
    reqBinsStr = reqBinsStr + histBins[nHistBins] + ".";
    std::cout << "  " << reqBinsStr << std::endl;

    std::cout << "return false" << std::endl;
  }

  return allBinsFound;
}

bool isFloatSame(float float1, float float2, float deltaValue)
{ 
  return TMath::Abs(float1 - float2) < deltaValue;
}


Int_t truncateBinRange(Int_t inNBins, Double_t inBins[], Double_t low, Double_t high, Double_t deltaVal, Double_t outBins[])
{
  Int_t outNBins = -1;
  bool hitLow = false;
  for(Int_t bI = 0; bI < inNBins; ++bI){
    if(!hitLow){
      if(TMath::Abs(low - inBins[bI]) < deltaVal){
	hitLow = true;
	++outNBins;
	outBins[outNBins] = inBins[bI];
      }
    }
    else{
      ++outNBins;
      outBins[outNBins] = inBins[bI];
      if(TMath::Abs(high - inBins[bI]) < deltaVal){
	break;
      }      
    }    
  }

  return outNBins;
}


template <typename T>
int hist1BinToHist2(T* inHist_p, bool isX, Int_t inBinPos, Float_t deltaVal, TH1F* outHist_p)
{
  int binPos = -1;

  TAxis* axis_p = nullptr;
  if(isX) axis_p = (TAxis*)inHist_p->GetXaxis();
  else axis_p = (TAxis*)inHist_p->GetYaxis();

  Float_t lowEdgeIn = axis_p->GetBinLowEdge(inBinPos + 1);
  Float_t highEdgeIn = axis_p->GetBinLowEdge(inBinPos + 2);


  for(Int_t bI = 0; bI < outHist_p->GetXaxis()->GetNbins(); ++bI){
    Float_t lowEdgeOut = outHist_p->GetXaxis()->GetBinLowEdge(bI+1);
    Float_t highEdgeOut = outHist_p->GetXaxis()->GetBinLowEdge(bI+2);

    //First check if the bin is fully contained
    if(lowEdgeIn >= lowEdgeOut && highEdgeIn < highEdgeOut){
      binPos = bI;
      break;
    }

    //Next check if its just the same bin w/in deltaVal
    if(TMath::Abs(lowEdgeOut - lowEdgeIn) > deltaVal) continue;
    if(TMath::Abs(highEdgeOut - highEdgeIn) > deltaVal) continue;

    binPos = bI;
    break;
  }
  
  return binPos;
}


#endif
