#ifndef UNFOLDINGUTIL_H
#define UNFOLDINGUTIL_H

//c and c++
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TH1F.h"
#include "TH2F.h"

void getIterativeHists(std::vector<TH1F*> inUnfoldedHists_p, TH1F* statDelta_p, TH1F* iterDelta_p, TH1F* totalDelta_p)
{
  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter - 1;

  if(nBins == 0){
    std::cout << "getIterativeHists Error: Given unfolding histogram vector has size \'" << nIter << "\', must have at least 2 iterations to eval. return" << std::endl;
    return;
  }

  std::vector<TH1F*> inHistsForChecks_p = {statDelta_p, iterDelta_p, totalDelta_p};
  std::vector<std::string> inHistNames = {"statDelta_p", "iterDelta_p", "totalDelta_p"};

  for(unsigned int i = 0; i < inHistsForChecks_p.size(); ++i){
    if(inHistsForChecks_p[i] == nullptr){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram is nullptr. return" << std::endl;
      return;
    }
    else if(inHistsForChecks_p[i]->GetXaxis()->GetNbins() != nBins){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram has nBins \'" << inHistsForChecks_p[i]->GetXaxis()->GetNbins() << "\' not equal to number needed for unfolding iterations \'" << nIter << "\'" << ", minus 1 for \'" << nBins << "\'. return 1" << std::endl;
      return;
    }
  }

  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      //Throw away the zeroes
      if(inUnfoldedHists_p[i]->GetBinContent(bIX+1) < TMath::Power(10, -50)) continue;

      Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1))/inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaIter += tempDeltaIter*tempDeltaIter;

      double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1)/inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaStat += tempDeltaStat*tempDeltaStat;    
    }

    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i-1, 0.0);
    
    iterDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i-1, 0.0);
    
    totalDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i-1, 0.0);
  }
    
  return;
}


void getIterativeHists2D(std::vector<TH2F*> inUnfoldedHists_p, TH1F* statDelta_p, TH1F* iterDelta_p, TH1F* totalDelta_p)
{
  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter - 1;

  if(nBins == 0){
    std::cout << "getIterativeHists Error: Given unfolding histogram vector has size \'" << nIter << "\', must have at least 2 iterations to eval. return" << std::endl;
    return;
  }

  std::vector<TH1F*> inHistsForChecks_p = {statDelta_p, iterDelta_p, totalDelta_p};
  std::vector<std::string> inHistNames = {"statDelta_p", "iterDelta_p", "totalDelta_p"};

  for(unsigned int i = 0; i < inHistsForChecks_p.size(); ++i){
    if(inHistsForChecks_p[i] == nullptr){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram is nullptr. return" << std::endl;
      return;
    }
    else if(inHistsForChecks_p[i]->GetXaxis()->GetNbins() != nBins){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram has nBins \'" << inHistsForChecks_p[i]->GetXaxis()->GetNbins() << "\' not equal to number needed for unfolding iterations \'" << nIter << "\'" << ", minus 1 for \'" << nBins << "\'. return 1" << std::endl;
      return;
    }
  }

  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      for(Int_t bIY = 0; bIY < inUnfoldedHists_p[i]->GetYaxis()->GetNbins(); ++bIY){
	//Throw away the zeroes
	if(inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) < TMath::Power(10, -50)) continue;

	Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1, bIY+1))/inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaIter += tempDeltaIter*tempDeltaIter;
	
	double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1, bIY+1)/inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaStat += tempDeltaStat*tempDeltaStat;    
      }
    }
    
    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i-1, 0.0);
    
    iterDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i-1, 0.0);
    
    totalDelta_p->SetBinContent(i-1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i-1, 0.0);
  }
    
  return;
}


#endif
