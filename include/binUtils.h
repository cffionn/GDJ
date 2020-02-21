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

inline double mathStringToNum(std::string inStr)//currently just handles fractions and pi
{
  double outVal = 1.0;
  inStr = returnAllCapsString(inStr);
  if(inStr.find("PI") != std::string::npos){
    outVal *= M_PI;
    inStr.replace(inStr.find("PI"), 2, "");
  }
  if(inStr.find("/") != std::string::npos){
    double num = std::stod(inStr.substr(0, inStr.find("/")));
    double denom = std::stod(inStr.substr(inStr.find("/")+1, inStr.size()));
    outVal *= num/denom;
  }
  else outVal *= std::stod(inStr);
    
  return outVal;
}

#endif
