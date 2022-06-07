#ifndef UNFOLDINGUTIL_H
#define UNFOLDINGUTIL_H

//c and c++
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

void getIterativeHists(TH1F* reco_p, std::vector<TH1F*> inUnfoldedHists_p, TH1F* statDelta_p, TH1F* iterDelta_p, TH1F* totalDelta_p, bool doRelative, float xLow, float xHigh, bool doPrint = false)
{
  const Float_t deltaCheck = 0.0001;

  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter;

  if(nBins == 0){
    std::cout << "getIterativeHists Error: Given unfolding histogram vector has size \'" << nIter << "\', must have at least 2 iterations to eval. return" << std::endl;
    return;
  }

  if(TMath::Abs(statDelta_p->GetBinCenter(1) - 1.0) > deltaCheck){
    std::cout << "getIterativeHists ERROR: First bin should be centered on 1 (for iteration 1). return" << std::endl;
    return;
  }

  //Check our x-bounds are bins
  bool xLowFound = false;
  bool xHighFound = false;
  for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1; ++bIX){
    Float_t binEdge = inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1);

    if(TMath::Abs(binEdge - xLow) < deltaCheck) xLowFound = true;
    if(TMath::Abs(binEdge - xHigh) < deltaCheck) xHighFound = true;
  }

  if(!xLowFound || !xHighFound){
    std::cout << "getIterativeHists ERROR: x-bounds not found!" << std::endl;
    std::cout << " xLow, isFound: " << xLow << ", " << xLowFound << std::endl;
    std::cout << " xHigh, isFound: " << xHigh << ", " << xHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
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
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram has nBins \'" << inHistsForChecks_p[i]->GetXaxis()->GetNbins() << "\' not equal to number needed for unfolding iterations \'" << nIter << "\'" << ", inpuy number of unfolded hists \'" << nBins << "\'. return 1" << std::endl;
      return;
    }
  }

  //We start be doing first iteration, which is a special iteration because we have to do it w/ respect to pre-unfold reco. rather than a previous iteration
  Float_t deltaIterReco = 0.0;
  Float_t deltaStatReco = 0.0;

  if(reco_p != nullptr){
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);

      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenter < xLow) continue;
      if(binCenter >= xHigh) continue;
      
      if(doRelative){
	if(inUnfoldedHists_p[0]->GetBinContent(bIX+1) < TMath::Power(10, -50)) continue;
      }
      
      Float_t tempDeltaIter = (inUnfoldedHists_p[0]->GetBinContent(bIX+1) - reco_p->GetBinContent(bIX+1));
      if(doRelative) tempDeltaIter /= inUnfoldedHists_p[0]->GetBinContent(bIX+1);
      deltaIterReco += tempDeltaIter*tempDeltaIter;
      
      double tempDeltaStat = inUnfoldedHists_p[0]->GetBinError(bIX+1);
      if(doRelative) tempDeltaStat /= inUnfoldedHists_p[0]->GetBinContent(bIX+1);
      deltaStatReco += tempDeltaStat*tempDeltaStat;
    }
  }

  statDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco));
  statDelta_p->SetBinError(1, 0.0);

  iterDelta_p->SetBinContent(1, TMath::Sqrt(deltaIterReco));
  iterDelta_p->SetBinError(1, 0.0);

  totalDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco + deltaIterReco));
  totalDelta_p->SetBinError(1, 0.0);
      
  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);      

      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenter < xLow) continue;
      if(binCenter >= xHigh) continue;

      //Throw away the zeroes
      if(doRelative){
	if(inUnfoldedHists_p[i]->GetBinContent(bIX+1) < TMath::Power(10, -50)) continue;
      }
      
      //Temporarily modified 2021.02.22 to test hypothesis on termination criteria (currently using relative differences a la Martin R but older analyses used absolute differences

      Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1));///inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      if(doRelative) tempDeltaIter /= inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaIter += tempDeltaIter*tempDeltaIter;

      double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1);///inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      if(doRelative) tempDeltaStat /= inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaStat += tempDeltaStat*tempDeltaStat;    
      
      if(doPrint) std::cout << "TEMPDELTASTAT: " << bIX << ": " << tempDeltaStat << std::endl;
    }

    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i+1, 0.0);
    
    iterDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i+1, 0.0);
    
    totalDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i+1, 0.0);
  }
    
  return;
}


void getIterativeHists2D(TH2F* reco_p, std::vector<TH2F*> inUnfoldedHists_p, TH1F* statDelta_p, TH1F* iterDelta_p, TH1F* totalDelta_p, bool doRelative, float xLow, float xHigh, float yLow, float yHigh, Int_t excludeYPos = -1)
{
  const Float_t deltaCheck = 0.0001;
  
  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter;

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

  //Check our x-bounds are bins
  bool xLowFound = false;
  bool xHighFound = false;
  for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1; ++bIX){
    Float_t binEdge = inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1);

    if(TMath::Abs(binEdge - xLow) < deltaCheck) xLowFound = true;
    if(TMath::Abs(binEdge - xHigh) < deltaCheck) xHighFound = true;
  }

  if(!xLowFound || !xHighFound){
    std::cout << "getIterativeHists ERROR: x-bounds not found!" << std::endl;
    std::cout << " xLow, isFound: " << xLow << ", " << xLowFound << std::endl;
    std::cout << " xHigh, isFound: " << xHigh << ", " << xHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
    return;
  }

  //Check our y-bounds are bins
  bool yLowFound = false;
  bool yHighFound = false;
  for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins()+1; ++bIY){
    Float_t binEdge = inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(bIY+1);

    if(TMath::Abs(binEdge - yLow) < deltaCheck) yLowFound = true;
    if(TMath::Abs(binEdge - yHigh) < deltaCheck) yHighFound = true;
  }

  if(!yLowFound || !yHighFound){
    std::cout << "getIterativeHists ERROR: y-bounds not found!" << std::endl;
    std::cout << " yLow, isFound: " << yLow << ", " << yLowFound << std::endl;
    std::cout << " yHigh, isFound: " << yHigh << ", " << yHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins(); ++bIY){
      std::cout << inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(bIY+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetYaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
    return;
  }
  
  //We start be doing first iteration, which is a special iteration because we have to do it w/ respect to pre-unfold reco. rather than a previous iteration
  Float_t deltaIterReco = 0.0;
  Float_t deltaStatReco = 0.0;

  if(reco_p != nullptr){
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenterX = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenterX < xLow) continue;
      if(binCenterX >= xHigh) continue;

      for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins(); ++bIY){
	Float_t binCenterY = inUnfoldedHists_p[0]->GetYaxis()->GetBinCenter(bIY+1);
	//Only calculate the unfolding termination in the region of your acceptance
	if(binCenterY < yLow) continue;
	if(binCenterY >= yHigh) continue;
	
	if(doRelative){
	  if(inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1) < TMath::Power(10, -50)) continue;
	}
	
	Float_t tempDeltaIter = (inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1) - reco_p->GetBinContent(bIX+1, bIY+1));
	if(doRelative) tempDeltaIter /= inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1);
	deltaIterReco += tempDeltaIter*tempDeltaIter;
	
	double tempDeltaStat = inUnfoldedHists_p[0]->GetBinError(bIX+1, bIY+1);
	if(doRelative) tempDeltaStat /= inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1);
	deltaStatReco += tempDeltaStat*tempDeltaStat;
      }
    }
  }

  statDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco));
  statDelta_p->SetBinError(1, 0.0);

  iterDelta_p->SetBinContent(1, TMath::Sqrt(deltaIterReco));
  iterDelta_p->SetBinError(1, 0.0);

  totalDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco + deltaIterReco));
  totalDelta_p->SetBinError(1, 0.0);
  
  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenterX = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenterX < xLow) continue;
      if(binCenterX >= xHigh) continue;

      for(Int_t bIY = 0; bIY < inUnfoldedHists_p[i]->GetYaxis()->GetNbins(); ++bIY){
	Float_t binCenterY = inUnfoldedHists_p[0]->GetYaxis()->GetBinCenter(bIY+1);
	//Only calculate the unfolding termination in the region of your acceptance
	if(binCenterY < yLow) continue;
	if(binCenterY >= yHigh) continue;

	//Throw away the zeroes
	if(doRelative){
	  if(inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) < TMath::Power(10, -50)) continue;
	}
	
      //Temporarily modified 2021.02.22 to test hypothesis on termination criteria (currently using relative differences a la Martin R but older analyses used absolute differences
	Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1, bIY+1));//inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	if(doRelative) tempDeltaIter /= inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaIter += tempDeltaIter*tempDeltaIter;
	
	double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1, bIY+1);///inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	if(doRelative) tempDeltaStat /= inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaStat += tempDeltaStat*tempDeltaStat;    
      }
    }
    
    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i+1, 0.0);
    
    iterDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i+1, 0.0);
    
    totalDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i+1, 0.0);
  }
    
  return;
}


void getIterativeHists(TH1D* reco_p, std::vector<TH1D*> inUnfoldedHists_p, TH1D* statDelta_p, TH1D* iterDelta_p, TH1D* totalDelta_p, bool doRelative, float xLow, float xHigh, bool doPrint = false)
{
  const Float_t deltaCheck = 0.0001;

  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter;

  if(nBins == 0){
    std::cout << "getIterativeHists Error: Given unfolding histogram vector has size \'" << nIter << "\', must have at least 2 iterations to eval. return" << std::endl;
    return;
  }

  if(TMath::Abs(statDelta_p->GetBinCenter(1) - 1.0) > deltaCheck){
    std::cout << "getIterativeHists ERROR: First bin should be centered on 1 (for iteration 1). return" << std::endl;
    return;
  }

  //Check our x-bounds are bins
  bool xLowFound = false;
  bool xHighFound = false;
  for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1; ++bIX){
    Float_t binEdge = inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1);

    if(TMath::Abs(binEdge - xLow) < deltaCheck) xLowFound = true;
    if(TMath::Abs(binEdge - xHigh) < deltaCheck) xHighFound = true;
  }

  if(!xLowFound || !xHighFound){
    std::cout << "getIterativeHists ERROR: x-bounds not found!" << std::endl;
    std::cout << " xLow, isFound: " << xLow << ", " << xLowFound << std::endl;
    std::cout << " xHigh, isFound: " << xHigh << ", " << xHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
    return;
  }


  std::vector<TH1D*> inHistsForChecks_p = {statDelta_p, iterDelta_p, totalDelta_p};
  std::vector<std::string> inHistNames = {"statDelta_p", "iterDelta_p", "totalDelta_p"};

  for(unsigned int i = 0; i < inHistsForChecks_p.size(); ++i){
    if(inHistsForChecks_p[i] == nullptr){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram is nullptr. return" << std::endl;
      return;
    }
    else if(inHistsForChecks_p[i]->GetXaxis()->GetNbins() != nBins){
      std::cout << "getIterativeHists Error: " << inHistNames[i] << " histogram has nBins \'" << inHistsForChecks_p[i]->GetXaxis()->GetNbins() << "\' not equal to number needed for unfolding iterations \'" << nIter << "\'" << ", inpuy number of unfolded hists \'" << nBins << "\'. return 1" << std::endl;
      return;
    }
  }

  //We start be doing first iteration, which is a special iteration because we have to do it w/ respect to pre-unfold reco. rather than a previous iteration
  Float_t deltaIterReco = 0.0;
  Float_t deltaStatReco = 0.0;

  if(reco_p != nullptr){
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenter < xLow) continue;
      if(binCenter >= xHigh) continue;

      if(doRelative){
	if(inUnfoldedHists_p[0]->GetBinContent(bIX+1) < TMath::Power(10, -50)) continue;
      }
      
      Float_t tempDeltaIter = (inUnfoldedHists_p[0]->GetBinContent(bIX+1) - reco_p->GetBinContent(bIX+1));
      if(doRelative) tempDeltaIter /= inUnfoldedHists_p[0]->GetBinContent(bIX+1);
      deltaIterReco += tempDeltaIter*tempDeltaIter;
      
      double tempDeltaStat = inUnfoldedHists_p[0]->GetBinError(bIX+1);
      if(doRelative) tempDeltaStat /= inUnfoldedHists_p[0]->GetBinContent(bIX+1);
      deltaStatReco += tempDeltaStat*tempDeltaStat;
    }
  }

  statDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco));
  statDelta_p->SetBinError(1, 0.0);

  iterDelta_p->SetBinContent(1, TMath::Sqrt(deltaIterReco));
  iterDelta_p->SetBinError(1, 0.0);

  totalDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco + deltaIterReco));
  totalDelta_p->SetBinError(1, 0.0);
      
  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenter = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenter < xLow) continue;
      if(binCenter >= xHigh) continue;

      //Throw away the zeroes
      if(doRelative){
	if(inUnfoldedHists_p[i]->GetBinContent(bIX+1) < TMath::Power(10, -50)) continue;
      }
      
      //Temporarily modified 2021.02.22 to test hypothesis on termination criteria (currently using relative differences a la Martin R but older analyses used absolute differences

      Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1));///inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      if(doRelative) tempDeltaIter /= inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaIter += tempDeltaIter*tempDeltaIter;

      double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1);///inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      if(doRelative) tempDeltaStat /= inUnfoldedHists_p[i]->GetBinContent(bIX+1);
      deltaStat += tempDeltaStat*tempDeltaStat;    
      
      if(doPrint) std::cout << "TEMPDELTASTAT: " << bIX << ": " << tempDeltaStat << std::endl;
    }

    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i+1, 0.0);
    
    iterDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i+1, 0.0);
    
    totalDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i+1, 0.0);
  }
    
  return;
}


void getIterativeHists2D(TH2D* reco_p, std::vector<TH2D*> inUnfoldedHists_p, TH1D* statDelta_p, TH1D* iterDelta_p, TH1D* totalDelta_p, bool doRelative, float xLow, float xHigh, float yLow, float yHigh)
{
  const Float_t deltaCheck = 0.0001;
  const Int_t nIter = inUnfoldedHists_p.size();
  const Int_t nBins = nIter;

  if(nBins == 0){
    std::cout << "getIterativeHists Error: Given unfolding histogram vector has size \'" << nIter << "\', must have at least 2 iterations to eval. return" << std::endl;
    return;
  }

  //Check our x-bounds are bins
  bool xLowFound = false;
  bool xHighFound = false;
  for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1; ++bIX){
    Float_t binEdge = inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1);

    if(TMath::Abs(binEdge - xLow) < deltaCheck) xLowFound = true;
    if(TMath::Abs(binEdge - xHigh) < deltaCheck) xHighFound = true;
  }

  if(!xLowFound || !xHighFound){
    std::cout << "getIterativeHists ERROR: x-bounds not found!" << std::endl;
    std::cout << " xLow, isFound: " << xLow << ", " << xLowFound << std::endl;
    std::cout << " xHigh, isFound: " << xHigh << ", " << xHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(bIX+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetXaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetXaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
    return;
  }

  //Check our y-bounds are bins
  bool yLowFound = false;
  bool yHighFound = false;
  for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins()+1; ++bIY){
    Float_t binEdge = inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(bIY+1);

    if(TMath::Abs(binEdge - yLow) < deltaCheck) yLowFound = true;
    if(TMath::Abs(binEdge - yHigh) < deltaCheck) yHighFound = true;
  }

  if(!yLowFound || !yHighFound){
    std::cout << "getIterativeHists ERROR: y-bounds not found!" << std::endl;
    std::cout << " yLow, isFound: " << yLow << ", " << yLowFound << std::endl;
    std::cout << " yHigh, isFound: " << yHigh << ", " << yHighFound << std::endl;
    std::cout << " Bins: ";
    for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins(); ++bIY){
      std::cout << inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(bIY+1) << ", ";
    }
    std::cout << inUnfoldedHists_p[0]->GetYaxis()->GetBinLowEdge(inUnfoldedHists_p[0]->GetYaxis()->GetNbins()+1) << ".";
    std::cout << " return w/o creating iterative hists" << std::endl;
    return;
  }
  
  std::vector<TH1D*> inHistsForChecks_p = {statDelta_p, iterDelta_p, totalDelta_p};
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

  //We start be doing first iteration, which is a special iteration because we have to do it w/ respect to pre-unfold reco. rather than a previous iteration
  Float_t deltaIterReco = 0.0;
  Float_t deltaStatReco = 0.0;

  if(reco_p != nullptr){
    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[0]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenterX = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenterX < xLow) continue;
      if(binCenterX >= xHigh) continue;

      for(Int_t bIY = 0; bIY < inUnfoldedHists_p[0]->GetYaxis()->GetNbins(); ++bIY){
	Float_t binCenterY = inUnfoldedHists_p[0]->GetYaxis()->GetBinCenter(bIY+1);
        //Only calculate the unfolding termination in the region of your acceptance
        if(binCenterY < yLow) continue;
        if(binCenterY >= yHigh) continue;

	if(doRelative){
	  if(inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1) < TMath::Power(10, -50)) continue;
	}
	
	Float_t tempDeltaIter = (inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1) - reco_p->GetBinContent(bIX+1, bIY+1));
	if(doRelative) tempDeltaIter /= inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1);
	deltaIterReco += tempDeltaIter*tempDeltaIter;
	
	double tempDeltaStat = inUnfoldedHists_p[0]->GetBinError(bIX+1, bIY+1);
	if(doRelative) tempDeltaStat /= inUnfoldedHists_p[0]->GetBinContent(bIX+1, bIY+1);
	deltaStatReco += tempDeltaStat*tempDeltaStat;
      }
    }
  }
  
  statDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco));
  statDelta_p->SetBinError(1, 0.0);

  iterDelta_p->SetBinContent(1, TMath::Sqrt(deltaIterReco));
  iterDelta_p->SetBinError(1, 0.0);

  totalDelta_p->SetBinContent(1, TMath::Sqrt(deltaStatReco + deltaIterReco));
  totalDelta_p->SetBinError(1, 0.0);


  for(unsigned int i = 1; i < inUnfoldedHists_p.size(); ++i){
    Float_t deltaIter = 0.0;
    Float_t deltaStat = 0.0;

    for(Int_t bIX = 0; bIX < inUnfoldedHists_p[i]->GetXaxis()->GetNbins(); ++bIX){
      Float_t binCenterX = inUnfoldedHists_p[0]->GetXaxis()->GetBinCenter(bIX+1);
      //Only calculate the unfolding termination in the region of your acceptance
      if(binCenterX < xLow) continue;
      if(binCenterX >= xHigh) continue;

      for(Int_t bIY = 0; bIY < inUnfoldedHists_p[i]->GetYaxis()->GetNbins(); ++bIY){
        Float_t binCenterY = inUnfoldedHists_p[0]->GetYaxis()->GetBinCenter(bIY+1);
        //Only calculate the unfolding termination in the region of your acceptance
        if(binCenterY < yLow) continue;
        if(binCenterY >= yHigh) continue;       

	//Throw away the zeroes
	if(doRelative){
	  if(inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) < TMath::Power(10, -50)) continue;
	}
	
      //Temporarily modified 2021.02.22 to test hypothesis on termination criteria (currently using relative differences a la Martin R but older analyses used absolute differences
	Float_t tempDeltaIter = (inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1) - inUnfoldedHists_p[i-1]->GetBinContent(bIX+1, bIY+1));//inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	if(doRelative) tempDeltaIter /= inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaIter += tempDeltaIter*tempDeltaIter;
	
	double tempDeltaStat = inUnfoldedHists_p[i]->GetBinError(bIX+1, bIY+1);///inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	if(doRelative) tempDeltaStat /= inUnfoldedHists_p[i]->GetBinContent(bIX+1, bIY+1);
	deltaStat += tempDeltaStat*tempDeltaStat;    
      }
    }
    
    Float_t deltaTot = deltaIter + deltaStat;

    statDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaStat));
    statDelta_p->SetBinError(i+1, 0.0);
    
    iterDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaIter));
    iterDelta_p->SetBinError(i+1, 0.0);
    
    totalDelta_p->SetBinContent(i+1, TMath::Sqrt(deltaTot));
    totalDelta_p->SetBinError(i+1, 0.0);
  }
    
  return;
}


#endif
