//Author: Chris McGinn (2020.08.06)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs
//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjNTupleToSignalHist(std::string inConfigFileName)
{
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"OUTFILENAME",
					"JETR",
                                        "NGAMMAPTBINS",
                                        "GAMMAPTBINSLOW",
                                        "GAMMAPTBINSHIGH",
                                        "GAMMAETABINSLOW",
                                        "GAMMAETABINSHIGH",
                                        "GAMMAETABINSDOABS",
                                        "NJETPTBINS",
                                        "JETPTBINSLOW",
                                        "JETPTBINSHIGH",
                                        "NTWOJETPTBINS",
                                        "TWOJETPTBINSLOW",
                                        "TWOJETPTBINSHIGH",
                                        "JETETABINSLOW",
                                        "JETETABINSHIGH",
                                        "JETETABINSDOABS"};
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  std::string outFileName = inConfig_p->GetValue("OUTFILENAME", "");
  if(outFileName.rfind(".") != std::string::npos) outFileName.replace(outFileName.rfind("."), outFileName.size(), "");
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  const Int_t nGammaPtBins = inConfig_p->GetValue("NGAMMAPTBINS", 20);
  const Float_t gammaPtBinsLow = inConfig_p->GetValue("GAMMAPTBINSLOW", 50);
  const Float_t gammaPtBinsHigh = inConfig_p->GetValue("GAMMAPTBINSHIGH", 70);
  const Float_t gammaPtMax = inConfig_p->GetValue("GAMMAPTMAX", 250);

  const Float_t gammaEtaBinsLow = inConfig_p->GetValue("GAMMAETABINSLOW", -2.8);
  const Float_t gammaEtaBinsHigh = inConfig_p->GetValue("GAMMAETABINSHIGH", 2.8);
  const Bool_t gammaEtaBinsDoAbs = inConfig_p->GetValue("GAMMAETABINSDOABS", 0);

  const Int_t nJetPtBins = inConfig_p->GetValue("NJETPTBINS", 20);
  const Float_t jetPtBinsLow = inConfig_p->GetValue("JETPTBINSLOW", 30);
  const Float_t jetPtBinsHigh = inConfig_p->GetValue("JETPTBINSHIGH", 50);

  const Int_t nTwoJetPtBins = inConfig_p->GetValue("NTWOJETPTBINS", 20);
  const Float_t twoJetPtBinsLow = inConfig_p->GetValue("TWOJETPTBINSLOW", 30);
  const Float_t twoJetPtBinsHigh = inConfig_p->GetValue("TWOJETPTBINSHIGH", 50); 
  const Float_t jetEtaBinsLow = inConfig_p->GetValue("JETETABINSLOW", -2.8);
  const Float_t jetEtaBinsHigh = inConfig_p->GetValue("JETETABINSHIGH", 2.8);
  const Bool_t jetEtaBinsDoAbs = inConfig_p->GetValue("JETETABINSDOABS", 0);

  const int jetR = inConfig_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);

  const Int_t nMaxBins = 500;
  const Int_t nMaxBinsHist = 50;

  if(nGammaPtBins > nMaxBins){
    std::cout << "REQUESTED NGAMMAPTBINS \'" << nGammaPtBins << "\' EXCEEDS MAX BINS \'" << nMaxBins << "\'. Either expand max bins or reduce request. return 1" << std::endl;
    return 1;
  }
  
  if(nJetPtBins > nMaxBins){
    std::cout << "REQUESTED NJETPTBINS \'" << nJetPtBins << "\' EXCEEDS MAX BINS \'" << nMaxBins << "\'. Either expand max bins or reduce request. return 1" << std::endl;
    return 1;
  }

  std::vector<double> dphiCuts;
  std::vector<std::string> dphiCutLabels;
  
  for(unsigned int dI = 0; dI < 20; ++dI){
    double dphiCut = inConfig_p->GetValue(("DPHICUT." + std::to_string(dI)).c_str(), -1.0);
    if(dphiCut < -0.0001) continue;
    
    std::string dphiCutLabel = inConfig_p->GetValue(("DPHICUTLABEL." + std::to_string(dI)).c_str(), "");

    if(dphiCutLabel.size() == 0){
      std::cout << "DPHICUT." << dI << " GIVEN \'" << dphiCut << "\' must have label DPHICUTLABEL." << dI << ". Please edit inConfigFile \'" << inConfigFileName << ". Skipping" << std::endl;
      continue;
    }

    dphiCuts.push_back(dphiCut);
    dphiCutLabels.push_back(dphiCutLabel);
  }

  Double_t gammaPtBins[nMaxBins+1];
  getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  Double_t jetPtBins[nMaxBins+1];
  getLinBins(jetPtBinsLow, jetPtBinsHigh, nJetPtBins, jetPtBins);

  Double_t twoJetPtBins[nMaxBins+1];
  getLinBins(twoJetPtBinsLow, twoJetPtBinsHigh, nTwoJetPtBins, twoJetPtBins);
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  std::vector<unsigned long long> nGammaPerMinPt;
  TH1F* jetPtPerGammaPtDPhi_h[nMaxBinsHist][nMaxBinsHist];
  TH1F* jetPtPerGammaPtDPhi_TwoJets_h[nMaxBinsHist][nMaxBinsHist][nMaxBinsHist];
  TH1F* jetVectPtPerGammaPtDPhi_TwoJets_h[nMaxBinsHist][nMaxBinsHist][nMaxBinsHist];

  if(nMaxBinsHist < nGammaPtBins) return 1;
  if(nMaxBinsHist < (int)dphiCuts.size()) return 1;
  if(nMaxBinsHist < nTwoJetPtBins) return 1;

  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    const std::string gammaStr = "GammaPt" + std::to_string(gI);
    nGammaPerMinPt.push_back(0);
    
    for(unsigned int dI = 0; dI < dphiCuts.size(); ++dI){
      const std::string dphiStr = "DPhi" + std::to_string(dI);
      jetPtPerGammaPtDPhi_h[gI][dI] = new TH1F(("jetPtPerGammaPtDPhi_" + gammaStr + "_" + dphiStr + "_h").c_str(), (";Reco. Jet p_{T} [GeV];N_{Jet,R=" + jetRStr+ "}/N_{#gamma}").c_str(), nJetPtBins, jetPtBins);

      jetPtPerGammaPtDPhi_h[gI][dI]->Sumw2();
      
      for(int jI = 0; jI < nTwoJetPtBins; ++jI){
	const std::string jStr = "JetPt" + std::to_string(jI);
	jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI] = new TH1F(("jetPtPerGammaPtDPhi_" + gammaStr + "_" + dphiStr + "_" + jStr + "_TwoJets_h").c_str(), (";Reco. Jet p_{T} [GeV];N_{Jet,R=" + jetRStr+ "}/N_{#gamma}").c_str(), nJetPtBins, jetPtBins);
	jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI] = new TH1F(("jetVectPtPerGammaPtDPhi_" + gammaStr + "_" + dphiStr + "_" + jStr + "_TwoJets_h").c_str(), (";Reco. Jet #Sigma #vec{p}_{T} [GeV];N_{Jet,R=" + jetRStr+ "}/N_{#gamma}").c_str(), nJetPtBins, jetPtBins);

	setSumW2({jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI], jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]});
      }
    }
  }
    
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* inFileConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> inFileReqParams = {"ISPP",
                                              "ISCOLLISIONS"};
  if(!checkEnvForParams(inFileConfig_p, inFileReqParams)) return 1;

  bool isPP = inFileConfig_p->GetValue("ISPP", 0);
  bool isCollisions = inFileConfig_p->GetValue("ISCOLLISIONS", 0);

  if(!isPP){
    std::cout << "Given input \'" << inFileName << "\' is not p+p. return 1" << std::endl;
    return 1;
  }
  if(!isCollisions){
    std::cout << "Given input \'" << inFileName << "\' is not data. return 1" << std::endl;
    return 1;
  }
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;
  std::vector<float>* photon_etcone30_p=nullptr;

  std::vector<float>* aktRhi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_phi_p=nullptr;

  inTree_p->SetBranchStatus("*", 0);
  
  inTree_p->SetBranchStatus("photon_pt", 1);
  inTree_p->SetBranchStatus("photon_eta", 1);
  inTree_p->SetBranchStatus("photon_phi", 1);
  inTree_p->SetBranchStatus("photon_tight", 1);
  inTree_p->SetBranchStatus("photon_etcone30", 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);  

  inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
  inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
  inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
  inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
  inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);                              

  const ULong64_t coutInterval = inConfig_p->GetValue("COUTINTERVAL", 20);
  ULong64_t nEvt = inConfig_p->GetValue("NEVT", -1);
  ULong64_t nEntries = (ULong64_t)inTree_p->GetEntries();
  if(nEvt > 0) nEntries = TMath::Min(nEntries, nEvt);
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/coutInterval);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    const bool doPrint = entry%nDiv == 0;
    if(doPrint) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;

    inTree_p->GetEntry(entry);

    std::vector<double> goodPhotonPts, goodPhotonPhis, goodPhotonEtas, goodJetPts, goodJetPhis, goodJetEtas;

    //Grab all the photons that pass cuts
    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      if(!isGoodPhoton(isPP, photon_tight_p->at(pI), photon_etcone30_p->at(pI), photon_eta_p->at(pI))) continue;
      
      if(photon_pt_p->at(pI) < gammaPtBinsLow) continue;
      if(photon_pt_p->at(pI) >= gammaPtMax) continue; //NOTE THIS IS NOT THE BIN MAX WE ARE NOT TESTING THAT

      double gammaEtaVal = photon_eta_p->at(pI);
      if(gammaEtaBinsDoAbs) gammaEtaVal = TMath::Abs(gammaEtaVal);

      if(gammaEtaVal < gammaEtaBinsLow) continue;
      if(gammaEtaVal >= gammaEtaBinsHigh) continue;

      goodPhotonPts.push_back(photon_pt_p->at(pI));
      goodPhotonPhis.push_back(photon_phi_p->at(pI));
      goodPhotonEtas.push_back(photon_eta_p->at(pI));
    }

    //Grab all the jets that pass cuts
    for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
      if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jetPtBinsLow) continue;
      if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jetPtBinsHigh) continue;
      
      double jetEtaVal = aktRhi_em_xcalib_jet_eta_p->at(jI);
      if(jetEtaBinsDoAbs) jetEtaVal = TMath::Abs(jetEtaVal);

      if(jetEtaVal < jetEtaBinsLow) continue;
      if(jetEtaVal >= jetEtaBinsHigh) continue;
      
      goodJetPts.push_back(aktRhi_em_xcalib_jet_pt_p->at(jI));
      goodJetPhis.push_back(aktRhi_em_xcalib_jet_phi_p->at(jI));
      goodJetEtas.push_back(aktRhi_em_xcalib_jet_eta_p->at(jI));
    }

    for(unsigned int pI = 0; pI < goodPhotonPts.size(); ++pI){
      for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
	if(goodPhotonPts[pI] < gammaPtBins[gI]) break;

	++(nGammaPerMinPt[gI]);

	std::vector<std::vector<TLorentzVector> > jetsPerDPhi;
	for(unsigned dI = 0; dI < dphiCuts.size(); ++dI){
	  jetsPerDPhi.push_back({});
	}
	
	for(unsigned int jI = 0; jI < goodJetPts.size(); ++jI){
	  double deltaPhi = TMath::Abs(getDPHI(goodPhotonPhis[pI], goodJetPhis[jI]));

	  double deltaR = getDR(goodPhotonEtas[pI], goodPhotonPhis[pI], goodJetEtas[jI], goodJetPhis[jI]);
	  if(deltaR < 0.4) continue;
	  

	  for(unsigned dI = 0; dI < dphiCuts.size(); ++dI){
	    if(deltaPhi < dphiCuts[dI]) break;

	    TLorentzVector tL;
	    tL.SetPtEtaPhiM(goodJetPts[jI], goodJetEtas[jI], goodJetPhis[jI], 0.0);
	    
	    jetsPerDPhi[dI].push_back(tL);
	    jetPtPerGammaPtDPhi_h[gI][dI]->Fill(goodJetPts[jI]);
	  }
	}
     
	for(unsigned dI = 0; dI < dphiCuts.size(); ++dI){
	  for(Int_t jI = 0; jI < nTwoJetPtBins; ++jI){
	    Int_t counter = 0;
	    
	    for(unsigned int dI2 = 0; dI2 < jetsPerDPhi[dI].size(); ++dI2){
	      if(jetsPerDPhi[dI][dI2].Pt() >= twoJetPtBins[jI]) ++counter;
	    }

	    if(counter >= 2){
	      TLorentzVector tL;
	      tL.SetPtEtaPhiM(0,0,0,0);

	      for(unsigned int dI2 = 0; dI2 < jetsPerDPhi[dI].size(); ++dI2){
		tL += jetsPerDPhi[dI][dI2];
		jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Fill(jetsPerDPhi[dI][dI2].Pt());
	      }

	      jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Fill(tL.Pt());
	    }
	  }
	}	 
      }
    }

  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    for(unsigned int dI = 0; dI < dphiCuts.size(); ++dI){
      jetPtPerGammaPtDPhi_h[gI][dI]->Scale(1./(double)nGammaPerMinPt[gI]);
      
      jetPtPerGammaPtDPhi_h[gI][dI]->Write("", TObject::kOverwrite);
      delete jetPtPerGammaPtDPhi_h[gI][dI];

      for(Int_t jI = 0; jI < nJetPtBins; ++jI){
	jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Scale(1./(double)nGammaPerMinPt[gI]);

	jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Write("", TObject::kOverwrite);
	delete jetPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI];	

	jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Scale(1./(double)nGammaPerMinPt[gI]);

	jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI]->Write("", TObject::kOverwrite);
	delete jetVectPtPerGammaPtDPhi_TwoJets_h[gI][dI][jI];	
      }
    }
  }

  inConfig_p->SetValue("NDPHICUTS", std::to_string((int)dphiCuts.size()).c_str());
  inConfig_p->Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;  

  std::cout << "GDJNTUPLETOSIGNALHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNTupleToSignalHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjNTupleToSignalHist(argv[1]);
  return retVal;
}

