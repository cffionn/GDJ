#ifndef BINFLATTENER_H
#define BINFLATTENER_H

//cpp dependencies
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

class binFlattener{
 public:
  binFlattener(){m_isInit=false;};
  ~binFlattener(){};

  binFlattener(std::string inBinFlattenerName, std::vector<double> inBinSet1, std::vector<double> inBinSet2);
  bool Init(std::string inBinFlattenerName, std::vector<double> inBinSet1, std::vector<double> inBinSet2);
  std::vector<double> GetFlattenedBins(double inLowVal=0.0, double inHiVal=1.0);
  int GetBin1PosFromGlobal(int globalPos);
  int GetBin2PosFromGlobal(int globalPos);
  int GetGlobalFromBin12Pos(int bin1Pos, int bin2Pos);
  void Clean();

 private: 
  bool m_isInit;
  std::string m_binFlattenerName;
  std::vector<double> m_binSet1, m_binSet2, m_binSetGlobal;
  unsigned int m_nBins1, m_nBins2, m_nBinsGlobal;
  double m_lowVal, m_hiVal;
  std::map<unsigned int, unsigned int> m_globalToBins1, m_globalToBins2;
  std::map<unsigned int, std::map<unsigned int, unsigned int> > m_bin1And2ToGlobal;
  
  static const int nMaxBins = 1000;
};

#endif

