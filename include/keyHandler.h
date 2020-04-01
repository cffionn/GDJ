//Author: Chris McGinn (2020.03.04)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

#ifndef KEYHANDLER_H
#define KEYHANDLER_H

//c+cpp
#include <vector>

class keyHandler{
 public:
  keyHandler(){};
  keyHandler(std::vector<unsigned long long> in_valMaxes);
  ~keyHandler();

  bool Init(std::vector<unsigned long long> in_valMaxes);
  unsigned long long GetKey(std::vector<unsigned long long> in_vals);
  void Clean();
  
 private:
  bool m_doDebug;
  std::vector<unsigned long long> m_valMaxes;
  std::vector<unsigned long long> m_multipliers;

  unsigned long long GetNearestTen(unsigned long long inVal);
};

#endif
