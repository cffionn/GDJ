#ifndef FILLUTILS_H
#define FILLUTILS_H

//c+cpp
#include <cmath>
#include <iostream>
#include <string>

//Local
#include "include/stringUtil.h"

#define _USE_MATH_DEFINES

inline void fillTH1(TH1F* inHist_p, Float_t fillVal, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal, weight);
  }
  return;
}

inline void fillTH2(TH2F* inHist_p, Float_t fillVal1, Float_t fillVal2, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal1, fillVal2);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal1, fillVal2, weight);
  }

  return;                                                                                              
}

#endif
