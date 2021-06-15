#ifndef HISTDEFUTILITY_H
#define HISTDEFUTILITY_H

#include <vector>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

void centerTitles(TH1* hist_p)
{
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();
  
  return;
}

void centerTitles(std::vector<TH1*> hists_)
{
  for(unsigned int pI = 0; pI < hists_.size(); ++pI){
    centerTitles(hists_.at(pI));
  }
  
  return;
}

void setSumW2(TH1* hist_p){hist_p->Sumw2(); return;}
void setSumW2(std::vector<TH1*> hists_)
{
  for(unsigned int pI = 0; pI < hists_.size(); ++pI){
    setSumW2(hists_.at(pI));
  }
  
  return;
}


Double_t getMinGTZero(TH1* inHist_p)
{
  Double_t min = inHist_p->GetMaximum();

  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) <= 0) continue;
    if(inHist_p->GetBinContent(bIX+1) < min) min = inHist_p->GetBinContent(bIX+1);
  }

  return min;
}

Double_t getMin(TH1* inHist_p)
{
  Double_t min = inHist_p->GetMaximum();

  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) < min) min = inHist_p->GetBinContent(bIX+1);
  }

  return min;
}

Double_t getMax(TH1* inHist_p)
{
  Double_t max = inHist_p->GetBinContent(1);

  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) > max) max = inHist_p->GetBinContent(bIX+1);
  }

  return max;
}

bool histToHistOverride(TH1F* inHist_p, TH1F* outHist_p)
{
  if(inHist_p == nullptr){
    std::cout << "histToHistOverride: 'inHist' is a nullptr. return" << std::endl;
    return false;
  }
  else if(outHist_p == nullptr){
    std::cout << "histToHistOverride: 'outHist' is a nullptr. return" << std::endl;
    return false;
  }

  Int_t inNBinsX = inHist_p->GetXaxis()->GetNbins();
  Int_t inNBinsY = inHist_p->GetYaxis()->GetNbins();

  if(inNBinsX != inNBinsY){
    
    return false;
  }

  for(Int_t bIX = 0; bIX < inNBinsX; ++bIX){
    
  }

  return true;
}


#endif
