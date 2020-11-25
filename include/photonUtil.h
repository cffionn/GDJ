//Original isGoodPhoton function by Chris McGinn 
//Addition of photon isolation functions by Yeonju Go
//As taken from fork here: https://www.github.com/YeonjuGo/GDJ

#ifndef PHOTONUTIL_H
#define PHOTONUTIL_H

//ROOT
#include "TMath.h"

inline bool isGoodPhoton(bool isPP, bool phoTight, float phoIso, float phoEta)
{
  phoEta = TMath::Abs(phoEta);
  if(phoEta >= 1.37 && phoEta < 1.52) return false;
  if(!phoTight) return false;
  if(isPP){
    if(phoIso > 3.0) return false;
  }
  else{
    if(phoIso > 8.0) return false;
  }

  return true;
}

double getCorrectedPhotonIsolation(bool isPP, float phoIso, float phoPt, float phoEta, float cent)
{
  phoEta = TMath::Abs(phoEta);

  ////////// These correction is optimized for the correctedIso < 3 GeV cut for both pp and PbPb /////////
  double correctedIso = phoIso;
  if(phoEta < 1.37) correctedIso = correctedIso - 0.021049*phoPt - 2.855968 + 3.0;
  else if(phoEta >= 1.52 && phoEta < 2.37) correctedIso = correctedIso - 0.038043*phoPt - 2.860857 + 3.0;
  else return 999; 

  if(!isPP){
      if(phoEta < 1.37) correctedIso = correctedIso + 0.163862*cent - 0.000797*cent*cent - 10.268763 + 3.0;
      else correctedIso = correctedIso + 0.155623*cent - 0.000789*cent*cent - 9.301528 + 3.0;
  }

  return correctedIso;
}

double getPtCorrectedPhotonIsolation(float phoIso, float phoPt, float phoEta)
{
  phoEta = TMath::Abs(phoEta);

  ////////// These correction is optimized for the correctedIso < 3 GeV cut for both pp and PbPb /////////
  double correctedIso = phoIso;
  if(phoEta < 1.37) correctedIso = correctedIso - 0.021049*phoPt - 2.855968 + 3.0;
  else if(phoEta >= 1.52 && phoEta < 2.37) correctedIso = correctedIso - 0.038043*phoPt - 2.860857 + 3.0;
  else return 999; 

  return correctedIso;
}

double getCentCorrectedPhotonIsolation(float phoIso, float phoEta, float cent)
{
  phoEta = TMath::Abs(phoEta);

  ////////// These correction is optimized for the correctedIso < 3 GeV cut for both pp and PbPb /////////
  double correctedIso = phoIso;
  if(phoEta < 1.37) correctedIso = correctedIso + 0.163862*cent - 0.000797*cent*cent - 10.268763 + 3.0;
  else if(phoEta >= 1.52 && phoEta < 2.37) correctedIso = correctedIso + 0.155623*cent - 0.000789*cent*cent - 9.301528 + 3.0;
  else return 999; 

  return correctedIso;
}

#endif
