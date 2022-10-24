#ifndef PURITYUTIL_H
#define PURITYUTIL_H

//c and cpp
#include <vector>
#include <string>

//ROOT
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"

//Local
#include "include/binFlattener.h"

inline void doPurityCorr(TH1D* phoHist_p, TH1D* phoMeanHist_p, TH1D* phoSidebandHist_p, TH2D* subHist_p, TH2D* sidebandHist_p, TH2D* purCorrHist_p, binFlattener* inBinFlattener, TF1* purityFit_p, TH2D* backgroundHist_p)
{
  std::vector<double> meanValPerGammaBin, rawPerGammaBin, sidebandPerGammaBin;

  for(Int_t bIX = 0; bIX < phoHist_p->GetXaxis()->GetNbins(); ++bIX){
    double numVal = phoMeanHist_p->GetBinContent(bIX+1);
    double denomVal = phoHist_p->GetBinContent(bIX+1);

    if(numVal < TMath::Power(10,-100)) meanValPerGammaBin.push_back(0.0);
    else meanValPerGammaBin.push_back(numVal/denomVal);

    rawPerGammaBin.push_back(denomVal);

    double sidebandVal = phoSidebandHist_p->GetBinContent(bIX+1);
    if(sidebandVal < TMath::Power(10,-100)) sidebandPerGammaBin.push_back(0.0);
    else sidebandPerGammaBin.push_back(sidebandVal);
  }
  
  for(Int_t bIX = 0; bIX < subHist_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < subHist_p->GetYaxis()->GetNbins(); ++bIY){
      Float_t binContent = subHist_p->GetBinContent(bIX+1, bIY+1);
      Float_t binError = subHist_p->GetBinError(bIX+1, bIY+1);

      Float_t binSidebandContent = sidebandHist_p->GetBinContent(bIX+1, bIY+1);
      Float_t binSidebandError = sidebandHist_p->GetBinError(bIX+1, bIY+1);

      if(binContent < TMath::Power(10,-50) && binSidebandContent < TMath::Power(10,-50)){
	purCorrHist_p->SetBinContent(bIX+1, bIY+1, 0.0);
	purCorrHist_p->SetBinError(bIX+1, bIY+1, 0.0);
	continue;
      }

      Int_t gammaPtBin = bIY;
      if(inBinFlattener != nullptr){gammaPtBin = inBinFlattener->GetBin1PosFromGlobal(bIY);}
   
      Float_t purity = purityFit_p->Eval(meanValPerGammaBin[gammaPtBin]);
      binSidebandContent = (1.0 - purity)*binSidebandContent*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];
      binSidebandError = (1.0 - purity)*binSidebandError*rawPerGammaBin[gammaPtBin]/sidebandPerGammaBin[gammaPtBin];

      if(backgroundHist_p != nullptr){
	backgroundHist_p->SetBinContent(bIX+1, bIY+1, binSidebandContent);
	backgroundHist_p->SetBinError(bIX+1, bIY+1, binSidebandError);
      }

      binContent -= binSidebandContent;
      binError = TMath::Sqrt(binError*binError + binSidebandError*binSidebandError);

      purCorrHist_p->SetBinContent(bIX+1, bIY+1, binContent);
      purCorrHist_p->SetBinError(bIX+1, bIY+1, binError);
    }
  }
 

  return;  
}



#endif
