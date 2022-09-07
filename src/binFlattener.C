//Local dependencies
#include "include/binFlattener.h"
#include "include/binUtils.h"
#include "include/getLinBins.h"

binFlattener::binFlattener(std::string inBinFlattenerName, std::vector<double> inBinSet1, std::vector<double> inBinSet2)
{
  if(!Init(inBinFlattenerName, inBinSet1, inBinSet2)) std::cout << "BINFLATTENER: Initialization failure for binFlattener \'" << inBinFlattenerName << "\'. return " << std::endl;
  return;
}

bool binFlattener::Init(std::string inBinFlattenerName, std::vector<double> inBinSet1, std::vector<double> inBinSet2)
{
  m_binFlattenerName = inBinFlattenerName;

  m_binSet1 = inBinSet1;
  m_binSet2 = inBinSet2;

  m_nBins1 = m_binSet1.size();
  m_nBins2 = m_binSet2.size();

  if(m_nBins1 > nMaxBins || m_nBins1 <= 0){
    if(m_nBins1 > nMaxBins) std::cout << "binFlattener Init Error: m_nBins1 \'" << m_nBins1 << "\' is greater than nMaxBins \'" << nMaxBins << "\'. Clean and return false" << std::endl;
    else if(m_nBins1 <= 0) std::cout << "binFlattener Init Error: m_nBins1 \'" << m_nBins1 << "\' is less-than-or-equal-to zero. Clean and return false" << std::endl;
      
    Clean();
    return false;
  }
  if(m_nBins2 > nMaxBins){
    if(m_nBins2 > nMaxBins) std::cout << "binFlattener Init Error: m_nBins2 \'" << m_nBins2 << "\' is greater than nMaxBins \'" << nMaxBins << "\'. Clean and return false" << std::endl;
    else if(m_nBins2 <= 0) std::cout << "binFlattener Init Error: m_nBins2 \'" << m_nBins2 << "\' is less-than-or-equal-to zero. Clean and return false" << std::endl;
    Clean();
    return false;
  }

  //Default values
  m_lowVal = 0.0;
  m_hiVal = 1.0;
  
  m_isInit = true;
  return m_isInit;
}

std::vector<double> binFlattener::GetFlattenedBins(double inLowVal, double inHiVal)
{
  if(!m_isInit){
    std::cout << "binFlattener::GetFlattenedBins() Error: Called w/o initializing. return empty vect {}" << std::endl;
    return {};
  }

  m_lowVal = inLowVal;
  m_hiVal = inHiVal;
  
  std::vector<double> flattenedBinArray;  
  getLinBins(m_lowVal, m_hiVal, m_nBins1*m_nBins2, &flattenedBinArray);
  
  //Create bin maps
  unsigned int globalBinPos = 0;
  for(unsigned int bI1 = 0; bI1 < m_nBins1; ++bI1){

    std::map<unsigned int, unsigned int> tempMap;
    for(unsigned int bI2 = 0; bI2 < m_nBins2; ++bI2){
      m_globalToBins1[globalBinPos] = bI1;
      m_globalToBins2[globalBinPos] = bI2;

      tempMap[bI2] = globalBinPos;
      
      ++globalBinPos;
    }   

    m_bin1And2ToGlobal[bI1] = tempMap;
  }
  m_nBinsGlobal = globalBinPos;
  
  return flattenedBinArray;
}

int binFlattener::GetBin1PosFromGlobal(int globalPos)
{
  if(!m_isInit){
    std::cout << "binFlattener::GetBin1PosFromGlobal() Error: Called w/o initializing. return bin1Pos -1" << std::endl;
    return -1;
  }

  if(globalPos < 0 || globalPos >= m_nBinsGlobal){
    std::cout << "binFlattener::GetBin1PosFromGlobal() Error: Given globalPos \'" << globalPos << "\' is outside accepted range (inclusive) \'0-" << m_nBinsGlobal-1 << "\'. return binPos -1" << std::endl;
    return -1;
  }

  return m_globalToBins1[globalPos];
}

int binFlattener::GetBin2PosFromGlobal(int globalPos)
{
  if(!m_isInit){
    std::cout << "binFlattener::GetBin2PosFromGlobal() Error: Called w/o initializing. return bin2Pos -1" << std::endl;
    return -1;
  }

  if(globalPos < 0 || globalPos >= m_nBinsGlobal){
    std::cout << "binFlattener::GetBin2PosFromGlobal() Error: Given globalPos \'" << globalPos << "\' is outside accepted range (inclusive) \'0-" << m_nBinsGlobal-1 << "\'. return binPos -1" << std::endl;
    return -1;
  }

  return m_globalToBins2[globalPos];
}

int binFlattener::GetGlobalFromBin12Pos(int bin1Pos, int bin2Pos)
{ 
  if(!m_isInit){
    std::cout << "binFlattener::GetGlobalFromBin12Pos() Error: Called w/o initializing. return bin2Pos -1" << std::endl;
    return -1;
  }

  if(bin1Pos < 0 || bin1Pos >= m_nBins1){
    std::cout << "binFlattener::GetGlobalFromBin12Pos() Error: Given bin1Pos \'" << bin1Pos << "\' is outside accepted range (inclusive) \'0-" << m_nBins1-1 << "\'. return bin1Pos -1" << std::endl;
    return -1;
  }
  if(bin2Pos < 0 || bin2Pos >= m_nBins2){
    std::cout << "binFlattener::GetGlobalFromBin12Pos() Error: Given bin2Pos \'" << bin2Pos << "\' is outside accepted range (inclusive) \'0-" << m_nBins2-1 << "\'. return bin2Pos -1" << std::endl;
    return -1;
  }

  return (m_bin1And2ToGlobal[bin1Pos])[bin2Pos];
}

void binFlattener::Clean()
{ 
  m_isInit = false;
  m_binFlattenerName = "";
 
  m_binSet1.clear();
  m_binSet2.clear();
  m_binSetGlobal.clear();
  
  m_globalToBins1.clear();
  m_globalToBins2.clear();
  m_bin1And2ToGlobal.clear();
  
  m_nBins1 = -1;
  m_nBins2 = -1;
  m_nBinsGlobal = -1;
  
  return;
}
