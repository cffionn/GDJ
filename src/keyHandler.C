//Author: Chris McGinn (2020.03.04)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>

//Local
#include "include/globalDebugHandler.h"
#include "include/keyHandler.h"

keyHandler::keyHandler(std::string in_handlerName)
{
  m_handlerName = in_handlerName;
  return;
}

keyHandler::keyHandler(std::string in_handlerName, std::vector<unsigned long long> in_valMaxes)
{
  Init(in_handlerName, in_valMaxes);
  return;
}

keyHandler::keyHandler(std::vector<std::string> in_valNames, std::vector<unsigned long long> in_valMaxes)
{

  Init(in_valNames, in_valMaxes);
  return;
}

keyHandler::keyHandler(std::string in_handlerName, std::vector<std::string> in_valNames, std::vector<unsigned long long> in_valMaxes)
{
  Init(in_handlerName, in_valNames, in_valMaxes);
  return;
}

keyHandler::~keyHandler()
{
  Clean();
  return;
}

bool keyHandler::Init(std::string in_handlerName, std::vector<unsigned long long> in_valMaxes)
{
  m_handlerName = in_handlerName;
  return Init(in_valMaxes);
}

bool keyHandler::Init(std::vector<std::string> in_valNames, std::vector<unsigned long long> in_valMaxes)
{
  m_valNames = in_valNames;
  return Init(in_valMaxes);
}

bool keyHandler::Init(std::string in_handlerName, std::vector<std::string> in_valNames, std::vector<unsigned long long> in_valMaxes)
{
  m_handlerName = in_handlerName;
  m_valNames = in_valNames;
  return Init(in_valMaxes);
}

bool keyHandler::Init(std::vector<unsigned long long> in_valMaxes)
{
  globalDebugHandler gDebug;
  m_doDebug = gDebug.GetDoGlobalDebug();
  m_valMaxes = in_valMaxes;
  m_multipliers.push_back(1);

  for(unsigned int mI = 1; mI < m_valMaxes.size(); ++mI){
    unsigned long long nearestTen = GetNearestTen(m_valMaxes[mI-1]);

    for(unsigned int mI2 = 0; mI2 < m_multipliers.size(); ++mI2){
      nearestTen *= m_multipliers[mI2];
    }

    m_multipliers.push_back(nearestTen);
  }

  for(unsigned int mI = 0; mI < m_valMaxes.size(); ++mI){
    std::cout << "MaxVal, multiplier for " << m_handlerName << ": " << m_valMaxes[mI] << ", " << m_multipliers[mI] << std::endl;
  }
  
  if(false && m_doDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  //Test that we will not hit our cap
  float val = 0;
  for(unsigned int mI = 0; mI < m_valMaxes.size(); ++mI){
    val += m_valMaxes[mI]*m_multipliers[mI];
  }

  if(val > __LONG_MAX__){
    std::cout << "keyHandler error - Given maxes would exceed max unsigned long long value \'" << __LONG_MAX__ << "\'. Initialization failed." << std::endl;
    Clean();
    return false;    
  }
  
  return true;
}

unsigned long long keyHandler::GetKey(std::vector<unsigned long long> in_vals)
{
  //test number of values and that they are below caps
  if(in_vals.size() != m_valMaxes.size()){
    std::cout << "keyHandler::GetKey() error - Given number of values \'" << in_vals.size() << "\' does not match number of initialized maxes, \'" << m_valMaxes.size() << "\'. return max unsigned long long val, " << __LONG_MAX__ << "." << std::endl;
    return __LONG_MAX__;
  }

  bool allValsBelowMax = true;
  for(unsigned int mI = 0; mI < m_valMaxes.size(); ++mI){
    if(in_vals[mI] > m_valMaxes[mI]){
      std::cout << "Given value \'" << in_vals[mI] << "\' exceeds corresponding max \'" << m_valMaxes[mI] << "\'. return max unsigned long long val, " << __LONG_MAX__ << "." << std::endl;
      allValsBelowMax = false;
    }
  }
  if(!allValsBelowMax) return __LONG_MAX__;

  if(false && m_doDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  unsigned long long retVal = 0;
  for(unsigned int mI = 0; mI < m_valMaxes.size(); ++mI){
    retVal += in_vals[mI]*m_multipliers[mI];
  }
  if(false && m_doDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  return retVal;
}


std::string keyHandler::GetKeyStr(unsigned long long inVal)
{
  std::string retStr = "";
  std::vector<unsigned long long> vals = InvertKey(inVal);

  for(unsigned int val = 0; val < vals.size(); ++val){
    retStr = retStr + "Bin \'" + m_valNames[val] + "\' value=" + std::to_string(vals[val]) + ",";
  }
  if(retStr.size() != 0) retStr.replace(retStr.size()-1, 1, ".");

  return retStr;
}


std::vector<unsigned long long> keyHandler::InvertKey(unsigned long long inVal)
{
  std::vector<unsigned long long> retVals;
  for(unsigned int mI = 1; mI < m_valMaxes.size(); ++mI){
    if(mI == 1){
      retVals.push_back(inVal%m_multipliers[mI]);
    }
    
    unsigned long long outVal = inVal/m_multipliers[mI];
    if(mI < m_valMaxes.size()-1){
      //      std::cout << "OUTVAL : " << outVal << ", " <<  m_multipliers[mI+1] << ", " << outVal%m_multipliers[mI+1] << std::endl;
      outVal = outVal%(m_multipliers[mI+1]/m_multipliers[mI]);
    }

    retVals.push_back(outVal);
  }
  
  return retVals;
}

void keyHandler::Clean()
{
  m_valMaxes.clear();
  m_multipliers.clear();
  
  return;
}

unsigned long long keyHandler::GetNearestTen(unsigned long long inVal)
{
  float tempVal = inVal;
  unsigned long long retVal = 10;
  while(tempVal >= 10){
    tempVal /= 10.;
    retVal *= 10;
  }
  
  return retVal;
}
