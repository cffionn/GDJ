#ifndef VARUTIL_H
#define VARUTIL_H

//c+cpp
#include <string>

//ROOT
#include "TLorentzVector.h"
#include "TMath.h"

//Local
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/stringUtil.h"

inline bool varNameToLabelIsMultijet(std::string varName)
{
  std::string varNameLower = returnAllLowercaseString(varName);
  bool retBool = false;
  
  if(isStrSame(varNameLower, "xj")) retBool = false;
  else if(isStrSame(varNameLower, "xjj")) retBool = true;
  else if(isStrSame(varNameLower, "ajj")) retBool = true;
  else if(isStrSame(varNameLower, "pt")) retBool = false;
  else if(isStrSame(varNameLower, "dphi")) retBool = false;
  else if(isStrSame(varNameLower, "drjj")) retBool = true;
  else if(isStrSame(varNameLower, "dphijj")) retBool = true;
  else if(isStrSame(varNameLower, "dphijjg")) retBool = true;
  else{
    std::cout << "VARUTIL::varNameToLabelIsMultijet() - Var \'" << varName << "\' not found. please add. return false" << std::endl;
  }
  
  return retBool;
}

inline std::string varNameToHistName(std::string varName)
{
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string retStr = "";

  if(isStrSame(varNameLower, "xj")) retStr = "XJ";
  else if(isStrSame(varNameLower, "xjj")) retStr = "XJJ";
  else if(isStrSame(varNameLower, "ajj")) retStr = "AJJ";
  else if(isStrSame(varNameLower, "pt")) retStr = "Pt";
  else if(isStrSame(varNameLower, "dphi")) retStr = "DPhi";
  else if(isStrSame(varNameLower, "drjj")) retStr = "DRJJ";
  else if(isStrSame(varNameLower, "dphijj")) retStr = "DPhiJJ";
  else if(isStrSame(varNameLower, "dphijjg")) retStr = "DPhiJJG";
  else{
    std::cout << "VARUTIL::varNameToHistName() - Var \'" << varName << "\' not found. please add. return empty string" << std::endl;
  }
  
  return retStr;
}

inline std::string varNameToLabel(std::string varName)
{
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string retStr = "";

  if(isStrSame(varNameLower, "xj")) retStr = "x_{J#gamma}";
  else if(isStrSame(varNameLower, "xjj")) retStr = "#vec{x}_{JJ#gamma}";
  else if(isStrSame(varNameLower, "ajj")) retStr = "A_{JJ#gamma}";
  else if(isStrSame(varNameLower, "pt")) retStr = "p_{T}";
  else if(isStrSame(varNameLower, "dphi")) retStr = "#Delta#phi_{J#gamma}";
  else if(isStrSame(varNameLower, "drjj")) retStr = "#DeltaR_{JJ}";
  else if(isStrSame(varNameLower, "dphijj")) retStr = "#Delta#phi_{JJ}";
  else if(isStrSame(varNameLower, "dphijjg")) retStr = "#Delta#phi_{JJ#gamma}";
  else{
    std::cout << "VARUTIL::varNameToLabel() - Var \'" << varName << "\' not found. please add. return empty string" << std::endl;
  }
  
  return retStr;
}

inline std::string varNameToNLabel(std::string varName)
{
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string retStr = "";

  if(isStrSame(varNameLower, "xj")) retStr = "N_{J#gamma}";
  else if(isStrSame(varNameLower, "xjj")) retStr = "N_{JJ#gamma}";
  else if(isStrSame(varNameLower, "ajj")) retStr = "N_{JJ#gamma}";
  else if(isStrSame(varNameLower, "pt")) retStr = "N_{J}";
  else if(isStrSame(varNameLower, "dphi")) retStr = "N_{J#gamma}";
  else if(isStrSame(varNameLower, "drjj")) retStr = "N_{JJ}";
  else if(isStrSame(varNameLower, "dphijj")) retStr = "N_{JJ}";
  else if(isStrSame(varNameLower, "dphijjg")) retStr = "N_{JJ#gamma}";
  else{
    std::cout << "VARUTIL::varNameToNLabel() - Var \'" << varName << "\' not found. please add. return empty string" << std::endl;
  }
  
  return retStr;
}

inline Double_t getVar(std::string varName, TLorentzVector jet, TLorentzVector jet2, TLorentzVector photon)
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


#endif
