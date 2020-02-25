#ifndef CENTRALITYFROMINPUT_H
#define CENTRALITYFROMINPUT_H

//cpp
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"

class centralityFromInput
{
 public:
  centralityFromInput(){return;}
  centralityFromInput(std::string inTableFile);
  ~centralityFromInput(){};
  void SetTable(std::string inTableFile);
  double GetCent(double inVal);
  void PrintTableTex();
  
 private:
  std::string m_tableFileName;
  checkMakeDir m_check;
  bool m_isInit;
  bool m_isDescending;
  std::vector<double> m_centVals;
};

#endif
