#ifndef HISTDEFUTILITY_H
#define HISTDEFUTILITY_H

#include <vector>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCanvas.h"
#include "TPad.h"

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


template <typename T>
void setMargins(T* input, Float_t marginL, Float_t marginT, Float_t marginR, Float_t marginB)
{
  input->SetLeftMargin(marginL);
  input->SetTopMargin(marginT);
  input->SetRightMargin(marginR);
  input->SetBottomMargin(marginB);

  return;
}

bool prepTH1(TH1* hist_p, int font, float titleSize, float labelSize, int color, int style, float size, int width, float xOffset, float yOffset)
{
  hist_p->GetXaxis()->SetTitleFont(font);
  hist_p->GetXaxis()->SetTitleSize(titleSize);
  hist_p->GetYaxis()->SetTitleFont(font);
  hist_p->GetYaxis()->SetTitleSize(titleSize);

  hist_p->GetXaxis()->SetLabelFont(font);
  hist_p->GetXaxis()->SetLabelSize(labelSize);
  hist_p->GetYaxis()->SetLabelFont(font);
  hist_p->GetYaxis()->SetLabelSize(labelSize);

  hist_p->SetLineColor(color);
  hist_p->SetLineWidth(width);

  hist_p->SetMarkerColor(color);
  hist_p->SetMarkerSize(size);
  hist_p->SetMarkerStyle(style);

  hist_p->GetXaxis()->SetTitleOffset(xOffset);
  hist_p->GetYaxis()->SetTitleOffset(yOffset);

  return true;
}

#endif
