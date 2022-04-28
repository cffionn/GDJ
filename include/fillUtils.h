#ifndef FILLUTILS_H
#define FILLUTILS_H

//c+cpp
#include <cmath>
#include <iostream>
#include <string>

//ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"

//Local
#include "include/stringUtil.h"

#define _USE_MATH_DEFINES

template <typename T>
inline void fillTH1(T* inHist_p, Float_t fillVal, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal, weight);
  }
  return;
}

template <typename T>
inline void fillTH2(T* inHist_p, Float_t fillVal1, Float_t fillVal2, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal1, fillVal2);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal1, fillVal2, weight);
  }

  return;                                                                                              
}

#endif
