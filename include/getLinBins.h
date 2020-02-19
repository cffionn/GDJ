#ifndef getLinBins_h
#define getLinBins_h

#include "TMath.h"

void getLinBins(const Float_t lower, const Float_t higher, const Int_t nBins, Double_t bins[])
{
  bins[0] = lower;
  bins[nBins] = higher;

  Float_t interval = (bins[nBins] - bins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    bins[iter] = bins[0] + iter*interval;
  }

  return;
}

#endif
