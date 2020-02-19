#ifndef CENTRALITYFROMINPUT_H
#define CENTRALITYFROMINPUT_H

//cpp
#include <fstream>
#include <iostream>
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

 private:
  checkMakeDir check;
  bool isInit;
  bool isDescending;
  std::vector<double> centVals;
};

#endif
