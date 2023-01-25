#ifndef etaPhiFunc_h
#define etaPhiFunc_h

#include "TMath.h"
#include <iostream>

//Basic phi handler w/ wraparound accounted for
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

//dEta dPhi quadrature sum between objects
inline Float_t getDR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2, std::string debug = ""){
  Float_t theDphi = getDPHI(phi1, phi2, debug);
  Float_t theDeta = eta1 - eta2;
  return TMath::Sqrt(theDphi*theDphi + theDeta*theDeta);
}

#endif
