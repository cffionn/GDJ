//Author: Chris McGinn (2020.03.04)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>

//Local
#include "include/globalDebugHandler.h"
#include "include/keyHandler.h"

keyHandler::keyHandler(std::vector<unsigned long long> in_valMaxes)
{
  Init(in_valMaxes);
  return;
}
 
keyHandler::~keyHandler()
{
  Clean();
  return;
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
    std::cout << "MaxVal, multiplier: " << m_valMaxes[mI] << ", " << m_multipliers[mI] << std::endl;
  }
  
  if(m_doDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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

  if(m_doDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  unsigned long long retVal = 0;
  for(unsigned int mI = 0; mI < m_valMaxes.size(); ++mI){
    retVal += in_vals[mI]*m_multipliers[mI];
  }
  return retVal;
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
