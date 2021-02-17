//Original isGoodPhoton function by Chris McGinn 
//Addition of photon isolation functions by Yeonju Go
//As taken from fork here: https://www.github.com/YeonjuGo/GDJ

#ifndef PHOTONUTIL_H
#define PHOTONUTIL_H

//cpp
#include <string>
#include <vector>

//ROOT
#include "TMath.h"

//SIDEBAND TYPE 1: tight and non-isolated
//SIDEBAND TYPE 2: non-tight and isolated
//SIDEBAND TYPE 3: non-tight and non-isolated
//SIDEBAND TYPE 4: non-isolated
//SIDEBAND TYPE 5: non-tight
enum photonType{TIGHT_ISO = 0,
		TIGHT_NONISO = 1,
		NONTIGHT_ISO = 2,
		NONTIGHT_NONISO = 3,
		NONISO = 4,
		NONTIGHT = 5,
		OTHER = 6};

std::string getPhotonTypeString(photonType inPhotonType)
{
  if(inPhotonType == TIGHT_ISO) return "TIGHT_ISO";
  else if(inPhotonType == TIGHT_NONISO) return "TIGHT_NONISO";
  else if(inPhotonType == NONTIGHT_ISO) return "NONTIGHT_ISO";
  else if(inPhotonType == NONTIGHT_NONISO) return "NONTIGHT_NONISO";
  else if(inPhotonType == NONISO) return "NONISO";
  else if(inPhotonType == NONTIGHT) return "NONTIGHT";
  else if(inPhotonType == OTHER) return "OTHER";

  std::cout << "getPhotonTypeString: Given photonType \'" << inPhotonType << "\' is not valid. return empty string" << std::endl;

  return "";
}

inline bool isIsolatedPhoton(bool isPP, bool isCorrected, float phoIso)
{
  //If less than these values it is isolated
  const Float_t isoCutPPAndPbPbCorr = 3.0;
  const Float_t isoCutPbPbUncorr = 8.0;

  Float_t isoCut = isoCutPPAndPbPbCorr;
  if(!isPP && !isCorrected) isoCut = isoCutPbPbUncorr;

  bool isIsolated = phoIso < isoCut;
  
  return isIsolated;
}

inline bool isNonIsolatedPhoton(bool isPP, bool isCorrected, float phoIso)
{
  //If greater than these values it is non-isolated
  const Float_t nonIsoCutPPAndPbPbCorr = 5.0;
  const Float_t nonIsoCutPbPbUncorr = 11.0; //This is an adhoc value

  Float_t nonIsoCut = nonIsoCutPPAndPbPbCorr;
  if(!isPP && !isCorrected) nonIsoCut = nonIsoCutPbPbUncorr;

  bool isNonIsolated = phoIso > nonIsoCut;
  
  return isNonIsolated;
}

inline bool isSidebandPhoton(bool isIsolated, bool isNonIsolated, bool isTight, photonType sidebandType)
{
  if(sidebandType < TIGHT_NONISO || sidebandType > OTHER){
    std::cout << "ERROR IN ISSIDEBANDPHOTON: sidebandType \'" << sidebandType << "\' is not valid please choose a valid type: " << std::endl;
    for(int i = TIGHT_ISO; i <= OTHER; ++i){
      photonType pType = static_cast<photonType>(i);
      std::cout << " " << i << ": \'" << getPhotonTypeString(pType) << "\'" << std::endl;
    }

    std::cout << "return false" << std::endl;
    return false;
  }

  bool isSB = false;
  if(sidebandType == 1 && isTight && isNonIsolated) isSB = true;
  else if(sidebandType == 2 && !isTight && isIsolated) isSB = true;
  else if(sidebandType == 3 && !isTight && isNonIsolated) isSB = true;
  else if(sidebandType == 4 && isNonIsolated) isSB = true;
  else if(sidebandType == 5 && !isTight) isSB = true;
  
  return isSB;
}

inline bool isSidebandPhoton(bool isPP, bool isCorrected, photonType sidebandType, bool phoTight, float phoIso)
{
  bool isIsolated = isIsolatedPhoton(isPP, isCorrected, phoIso);
  bool isNonIsolated = isNonIsolatedPhoton(isPP, isCorrected, phoIso);
  return isSidebandPhoton(isIsolated, isNonIsolated, phoTight, sidebandType);
}

//NOTE THAT THIS IS RE-OPTIMIZED BY YEONJU
//ISOLATION CUT SHOULD BE 3 GeV for Pb+Pb and p+p after correction, R=0.3
inline bool isGoodPhoton(bool isPP, bool isCorrected, photonType sidebandType, bool phoTight, float phoIso, float phoEta)
{
  phoEta = TMath::Abs(phoEta);
  //Always
  if(phoEta >= 1.37 && phoEta < 1.52) return false;

  bool isIsolated = isIsolatedPhoton(isPP, isCorrected, phoIso);
  bool isNonIsolated = isNonIsolatedPhoton(isPP, isCorrected, phoIso);
  bool isSideband = isSidebandPhoton(isIsolated, isNonIsolated, phoTight, sidebandType);
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
