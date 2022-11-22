#ifndef etaPhiFunc_h
#define etaPhiFunc_h

#include "TMath.h"
#include <iostream>

inline Float_t getDPHI(Float_t phi1, Float_t phi2, std::string debug = ""){
  Float_t dphi = phi1 - phi2;

  if(dphi > TMath::Pi())
    dphi = dphi - 2. * TMath::Pi();
  if(dphi <= -TMath::Pi())
    dphi = dphi + 2. * TMath::Pi();

  if(TMath::Abs(dphi) > TMath::Pi()) {
    std::cout << "commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
    if(debug.size() != 0) std::cout << "DEBUGSTR: " << debug << std::endl;
    std::cout << " " << phi1 << ", " << phi2 << ", " << dphi << std::endl;
  }

  return dphi;
}


inline Float_t getDR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2, std::string debug = ""){
  Float_t theDphi = getDPHI(phi1, phi2, debug);
  Float_t theDeta = eta1 - eta2;
  return TMath::Sqrt(theDphi*theDphi + theDeta*theDeta);
}

/*
inline Bool_t photonEtaIsBarrel(Float_t photonEta)
{
  return TMath::Abs(photonEta) < 1.37;
}

inline Bool_t photonEtaIsEC(Float_t photonEta)
{
  return TMath::Abs(photonEta) > 1.52 && TMath::Abs(photonEta) < 2.37;
}

inline Bool_t photonEtaIsGood(Float_t photonEta)
{
  return photonEtaIsBarrel(photonEta) || photonEtaIsEC(photonEta);
}
*/

#endif
