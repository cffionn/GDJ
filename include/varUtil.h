#ifndef VARUTIL_H
#define VARUTIL_H

//c+cpp
#include <string>

//Local
#include "include/stringUtil.h"

std::string varNameToLabel(std::string varName)
{
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string retStr = "";

  if(isStrSame(varNameLower, "xj")) retStr = "x_{J#gamma}";
  else if(isStrSame(varNameLower, "xjj")) retStr = "#vec{x}_{JJ#gamma}";
  else if(isStrSame(varNameLower, "ajj")) retStr = "A_{JJ}";
  else if(isStrSame(varNameLower, "pt")) retStr = "p_{T}";
  else if(isStrSame(varNameLower, "dphi")) retStr = "#Delta#phi_{J#gamma}";
  else if(isStrSame(varNameLower, "drjj")) retStr = "#DeltaR_{JJ}";
  else{
    std::cout << "VARUTIL::varNameToLabel() - Var \'" << varName << "\' not found. please add. return empty string" << std::endl;
  }
  
  return retStr;
}

#endif
