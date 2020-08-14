//Author: Chris McGinn (2020.08.11)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/stringUtil.h"

int gdjToyMultiMix(std::string inConfigFileName)
{
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const std::string dateStr = getDateStr();
  checkMakeDir check;
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"NEVT",
					"OUTFILENAME"};
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const ULong64_t nEvt = inConfig_p->GetValue("NEVT", 1);

  std::string outFileName = inConfig_p->GetValue("OUTFILENAME", "");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";


  TRandom3* randGen_p = new TRandom3(0);

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");  
  
  TH1F* gammaHist_p = new TH1F("gammaHist_h", ";#gamma p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{#gamma}}{dp_{T}^{#gamma}}", 60, 60, 120);
  TH1F* jetSignalHist_p = new TH1F("jetSignalHist_h", ";Signal Jet p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{Signal}}{dp_{T}^{Signal}}", 120, 0, 120);
  TH1F* jet1Hist_p = new TH1F("jet1Hist_h", ";Leading Jet p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{Leading}}{dp_{T}^{Leading}}", 120, 0, 120);
  TH1F* jet2Hist_p = new TH1F("jet2Hist_h", ";Subleading Jet p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{Subleading}}{dp_{T}^{Subleading}}", 120, 0, 120);
  TH1F* jetBkgdHist_p = new TH1F("jetBkgdHist_h", ";Background Jet p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{Fake Jet}}{dp_{T}^{Fake Jet}}", 120, 0, 120);

  TH1F* jetTotalHist_p = new TH1F("jetTotalHist_h", ";All Jet p_{T} [GeV];#frac{1}{N_{#gamma}} #frac{dN_{Jet}}{dp_{T}^{Jet}}", 120, 0, 120);
  
  TH1F* signalHist_p = new TH1F("signalHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* bkgdHist_p = new TH1F("bkgdHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* pureBkgdHist_p = new TH1F("pureBkgdHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* mixedBkgdHist_p = new TH1F("mixedBkgdHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* signalAndBkgdHist_p = new TH1F("signalAndBkgdHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);

  TH1F* mixedHistTrue_p = new TH1F("mixedHistTrue_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* mixedHist_p = new TH1F("mixedHist_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);
  TH1F* mixedHistCorrection_p = new TH1F("mixedHistCorrection_h", ";#vec{x}_{JJ#gamma}=(#Sigma #vec{p}_{Dijet})_{T}/#gamma p_{T};#frac{1}{N_{#gamma}} #frac{dN}{d#vec{x}_{JJ#gamma}}", 40, 0, 2.0);

  std::vector<TH1F*> hists_p = {gammaHist_p, jetSignalHist_p, jet1Hist_p, jet2Hist_p, jetBkgdHist_p, jetTotalHist_p, signalHist_p, bkgdHist_p, pureBkgdHist_p, mixedBkgdHist_p, signalAndBkgdHist_p, mixedHistTrue_p, mixedHist_p, mixedHistCorrection_p};  

  const Double_t absEtaMax = 2.8;
  const Int_t nDraws = 2.0*absEtaMax*2.0/(0.4*0.4);
  const double minPt = 20.0;
  
  ULong64_t nEvts = 0;
  //  for(Int_t eI = 0; eI < nEvt; ++eI){
  while(nEvts < nEvt){
    Double_t leadPt = randGen_p->Uniform(80.0, 100.0);

    double prob = randGen_p->Uniform(1.0/leadPt, 1.0);

    double pt2X = 10.0/prob;
    double pt1X = leadPt - pt2X;

    double pt2Y = pt2X*randGen_p->Gaus(1.0, 0.2);
    double pt1Y = -pt2Y;

    double e2 = TMath::Sqrt(pt2X*pt2X + pt2Y*pt2Y);
    double e1 = TMath::Sqrt(pt1X*pt1X + pt1Y*pt1Y);

    if(e1 < minPt) continue;
    if(e2 < minPt) continue;
    
    std::vector<TLorentzVector> jets, jetsSignal, jetsBkgd, jetsBkgd2, jetsBkgd3;
    jetsSignal.push_back(TLorentzVector(pt1X, pt1Y, 0, e1));
    jetsSignal.push_back(TLorentzVector(pt2X, pt2Y, 0, e2));

    for(Int_t bI = 0; bI < 3; ++bI){
      for(Int_t dI = 0; dI < nDraws; ++dI){
	double pT = randGen_p->Gaus(0.0, 20.0);
	if(pT < minPt) continue;
	double phi = randGen_p->Uniform(-TMath::Pi(), TMath::Pi());
	
	//Following cut reflects choice of gamma phi = 0
	if(phi < -TMath::Pi()/2. || phi >= TMath::Pi()/2.) continue;
	
	if(bI == 0) jetsBkgd.push_back(TLorentzVector(pT*TMath::Cos(phi), pT*TMath::Sin(phi), 0.0, pT));
	else if(bI == 1) jetsBkgd2.push_back(TLorentzVector(pT*TMath::Cos(phi), pT*TMath::Sin(phi), 0.0, pT));
	else jetsBkgd3.push_back(TLorentzVector(pT*TMath::Cos(phi), pT*TMath::Sin(phi), 0.0, pT));
      }
    }
    
    gammaHist_p->Fill(leadPt);//, 1.0/(double)nEvt);
    jet1Hist_p->Fill(jetsSignal[0].Pt());//, 1.0/(double)nEvt);
    jet2Hist_p->Fill(jetsSignal[1].Pt());//, 1.0/(double)nEvt);

    jetSignalHist_p->Fill(jetsSignal[0].Pt());
    jetSignalHist_p->Fill(jetsSignal[1].Pt());   
    
    jetTotalHist_p->Fill(jetsSignal[0].Pt());
    jetTotalHist_p->Fill(jetsSignal[1].Pt());
    for(unsigned int gI = 0; gI < jetsBkgd.size(); ++gI){
      jetBkgdHist_p->Fill(jetsBkgd[gI].Pt());//, 1.0/(double)nEvt);
      jetTotalHist_p->Fill(jetsBkgd[gI].Pt());
    }
  
    jets = jetsSignal;
    jets.insert(std::end(jets), std::begin(jetsBkgd), std::end(jetsBkgd));     
    
    for(unsigned int jI = 0; jI < jets.size(); ++jI){
      for(unsigned int jI2 = jI+1; jI2 < jets.size(); ++jI2){

	TLorentzVector tL = jets[jI] + jets[jI2];

	if(jI2 == 1) signalHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
	else{
	  if(jI == 0 || jI == 1) mixedBkgdHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
	  else pureBkgdHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);

	  bkgdHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
	}

	signalAndBkgdHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
      }
    }

    for(unsigned int jI = 0; jI < jetsBkgd2.size(); ++jI){
      //Mixing step 1 - in event associations
      for(unsigned int jI2 = jI+1; jI2 < jetsBkgd2.size(); ++jI2){
	TLorentzVector tL = jetsBkgd2[jI] + jetsBkgd2[jI2];
	mixedHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
	mixedHistTrue_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
      }

      //mixing step 2 - mixed associations      
      for(unsigned int jI2 = 0; jI2 < jets.size(); ++jI2){	
	TLorentzVector tL = jetsBkgd2[jI] + jets[jI2];
	mixedHist_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
	if(jI2 == 0 || jI2 == 1) mixedHistTrue_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);
      }

      for(unsigned int jI2 = 0; jI2 < jetsBkgd3.size(); ++jI2){
	TLorentzVector tL = jetsBkgd2[jI] + jetsBkgd3[jI2];
        mixedHistCorrection_p->Fill(tL.Px()/leadPt);//, 1.0/(double)nEvt);	
      }
    }
        
    ++nEvts;
  }
    
  outFile_p->cd();

  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    hists_p[hI]->Scale(1.0/(double)nEvt);

    
    for(Int_t bIX = 0; bIX < hists_p[hI]->GetXaxis()->GetNbins(); ++bIX){
      double width = hists_p[hI]->GetBinWidth(bIX+1);
      hists_p[hI]->SetBinContent(bIX+1, hists_p[hI]->GetBinContent(bIX+1)/width);
      hists_p[hI]->SetBinError(bIX+1, hists_p[hI]->GetBinError(bIX+1)/width);
    }
    
    hists_p[hI]->Write("", TObject::kOverwrite);
    delete hists_p[hI];
  }

  inConfig_p->Write("config", TObject::kOverwrite);
  delete inConfig_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  delete randGen_p;
  
  std::cout << "GDJTOYMULTIMIX COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjToyMultiMix.exe <nEvt>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjToyMultiMix(argv[1]);
  return retVal;
}
