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

#endif
