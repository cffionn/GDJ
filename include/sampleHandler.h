//Author: Chris McGinn (2020.02.27)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

#ifndef SAMPLEHANDLER_H
#define SAMPLEHANDLER_H

#include <map>
#include <string>
#include <vector>

class sampleHandler{
 public:
  sampleHandler(){Clean(); PreInit();}
  sampleHandler(bool in_isPP, bool in_isMC, int in_year, int in_minPthat);
  ~sampleHandler();

  bool Init(std::string sampleString);
  bool Init(bool in_isPP, bool in_isMC, int in_year, int in_minPthat);
  int GetTag();
  double GetXSection();
  double GetFilterEff();
  void Clean();
  void PrintTags();  
  
 private:
  int CreateTag();
  int CreateTag(bool in_isPP, bool in_isMC, int in_year, int in_minPthat);
  void PreInit();
  
  bool m_isInit;
  bool m_isPP;
  bool m_isMC;
  int m_year;
  int m_minPthat;

  int m_tagVal;

  std::vector<int> validYears = {2015, 2017, 2018};
  std::map<int, std::vector<int>> validMinPthatsByYear;
  std::map<int, double> tagToXSec;
  std::map<int, double> tagToFilterEff;

  std::map<std::string, bool> dataSetNameToIsPP;
  std::map<std::string, bool> dataSetNameToIsMC;
  std::map<std::string, int> dataSetNameToYear;
  std::map<std::string, int> dataSetNameToMinPthat;
};

#endif
