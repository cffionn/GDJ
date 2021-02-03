//Original isGoodPhoton function by Chris McGinn 
//Addition of photon isolation functions by Yeonju Go
//As taken from fork here: https://www.github.com/YeonjuGo/GDJ

#ifndef PHOTONUTIL_H
#define PHOTONUTIL_H

//ROOT
#include "TMath.h"

inline bool isIsolatedPhoton(bool isPP, bool isCorrected, float phoIso)
{
  const Float_t isoCutPPAndPbPbCorr = 3.0;
  const Float_t isoCutPbPbUncorr = 8.0;

  Float_t isoCut = isoCutPPAndPbPbCorr;
  if(!isPP && !isCorrected) isoCut = isoCutPbPbUncorr;

  bool isIsolated = phoIso < isoCut;
  
  return isIsolated;
}

//Adding a 'sideband' function for quick determination
//SIDEBAND TYPE 1: tight and non-isolated
//SIDEBAND TYPE 2: non-tight and isolated
//SIDEBAND TYPE 3: non-tight and non-isolated
//SIDEBAND TYPE 4: non-isolated
//SIDEBAND TYPE 5: non-tight

inline bool isSidebandPhoton(bool isIsolated, bool isTight, int sidebandType)
{
  if(sidebandType < 1 || sidebandType > 5){
    std::cout << "ERROR IN ISSIDEBANDPHOTON: sidebandType \'" << sidebandType << "\' is not valid please choose 1-5. return false" << std::endl;
    return false;
  }

  bool isSB = false;
  if(sidebandType == 1 && isTight && !isIsolated) isSB = true;
  else if(sidebandType == 2 && !isTight && isIsolated) isSB = true;
  else if(sidebandType == 3 && !isTight && !isIsolated) isSB = true;
  else if(sidebandType == 4 && !isIsolated) isSB = true;
  else if(sidebandType == 5 && !isTight) isSB = true;
  
  return isSB;
}

inline bool isSidebandPhoton(bool isPP, bool isCorrected, int sidebandType, bool phoTight, float phoIso)
{
  bool isIsolated = isIsolatedPhoton(isPP, isCorrected, phoIso);
  return isSidebandPhoton(isIsolated, phoTight, sidebandType);
}

//NOTE THAT THIS IS RE-OPTIMIZED BY YEONJU
//ISOLATION CUT SHOULD BE 3 GeV for Pb+Pb and p+p after correction, R=0.3
inline bool isGoodPhoton(bool isPP, bool isCorrected, int sidebandType, bool phoTight, float phoIso, float phoEta)
{
  phoEta = TMath::Abs(phoEta);
  //Always
  if(phoEta >= 1.37 && phoEta < 1.52) return false;

  bool isIsolated = isIsolatedPhoton(isPP, isCorrected, phoIso);
  bool isSideband = isSidebandPhoton(isIsolated, phoTight, sidebandType);
  bool retVal = (isIsolated && phoTight) || isSideband;

  return retVal;
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
