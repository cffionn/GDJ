//Author: Chris McGinn (2020.02.19)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/histDefUtility.h"
#include "include/stringUtil.h"

int gdjNTupleToHist(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;

  configParser config(inConfigFileName);
  std::string inROOTFileName = config.GetConfigVal("INFILENAME");
  std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
  std::string outFileName = config.GetConfigVal("OUTFILENAME");
  
  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  centralityFromInput centTable(inCentFileName);
  
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + ".root";

  const bool isPP = std::stoi(config.GetConfigVal("ISPP"));
  const Int_t nMaxCentBins = 10;
  Int_t nCentBins = 1;
  std::vector<int> centBins;
  std::vector<std::string> centBinsStr = {"PP"};
  if(!isPP){
    centBins = strToVectI(config.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;

    if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

    centBinsStr.clear();
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));
    }
  }

  const Int_t nMaxPtBins = 200;
  const Int_t nPtBins = std::stoi(config.GetConfigVal("NPTBINS"));
  if(!goodBinning(inConfigFileName, nMaxPtBins, nPtBins, "NPTBINS")) return 1;
  const Float_t ptBinsLow = std::stof(config.GetConfigVal("PTBINSLOW"));
  const Float_t ptBinsHigh = std::stof(config.GetConfigVal("PTBINSHIGH"));
  const Bool_t ptBinsDoLog = std::stof(config.GetConfigVal("PTBINSDOLOG"));
  Double_t ptBins[nMaxPtBins+1];
  if(ptBinsDoLog) getLogBins(ptBinsLow, ptBinsHigh, nPtBins, ptBins);
  else getLinBins(ptBinsLow, ptBinsHigh, nPtBins, ptBins);

  const Int_t nMaxEtaBins = 100;
  const Int_t nEtaBins = std::stoi(config.GetConfigVal("NETABINS"));
  if(!goodBinning(inConfigFileName, nMaxEtaBins, nEtaBins, "NETABINS")) return 1;
  const Float_t etaBinsLow = std::stof(config.GetConfigVal("ETABINSLOW"));
  const Float_t etaBinsHigh = std::stof(config.GetConfigVal("ETABINSHIGH"));
  Double_t etaBins[nMaxEtaBins+1];
  getLinBins(etaBinsLow, etaBinsHigh, nEtaBins, etaBins);

  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* centrality_p = nullptr;
  TH1F* photonPt_p[nMaxCentBins];
  TH1F* photonEta_p[nMaxCentBins];
  TH2F* photonEtaPt_p[nMaxCentBins];
  
  if(!isPP){
    centrality_p = new TH1F("centrality_h", ";Centrality (%);Counts", 100, -0.5, 99.5);
    centerTitles(centrality_p);
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    photonPt_p[cI] = new TH1F(("photonPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma p_{T} [GeV];Counts", nPtBins, ptBins);
    photonEta_p[cI] = new TH1F(("photonEta_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;Counts", nEtaBins, etaBins);

    photonEtaPt_p[cI] = new TH2F(("photonEtaPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;#gamma p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins);
    
    centerTitles({photonPt_p[cI], photonEta_p[cI]});
    centerTitles(photonEtaPt_p[cI]);
  }
  
  TFile* inFile_p = new TFile(inROOTFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  Float_t fcalA_et, fcalC_et;
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  
  inTree_p->SetBranchStatus("*", 0);
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
  }
  inTree_p->SetBranchStatus("photon_pt", 1);
  inTree_p->SetBranchStatus("photon_eta", 1);
  
  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
  }
  inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
  inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);

  const ULong64_t nEntries = inTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    Int_t centPos = -1;
    if(!isPP){
      Double_t cent = centTable.GetCent(fcalA_et + fcalC_et);
      centPos = ghostPos(centBins, cent, false);
      centrality_p->Fill(cent);      
    }
    else centPos = 0;

    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      photonPt_p[centPos]->Fill(photon_pt_p->at(pI));
      photonEta_p[centPos]->Fill(photon_eta_p->at(pI));
      photonEtaPt_p[centPos]->Fill(photon_eta_p->at(pI), photon_pt_p->at(pI));
    }
  }  
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  if(!isPP) centrality_p->Write("", TObject::kOverwrite);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    photonPt_p[cI]->Write("", TObject::kOverwrite);
    photonEta_p[cI]->Write("", TObject::kOverwrite);
    photonEtaPt_p[cI]->Write("", TObject::kOverwrite);
  }

  if(!isPP) delete centrality_p;
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    delete photonPt_p[cI];
    delete photonEta_p[cI];
    delete photonEtaPt_p[cI];
  }
  
  outFile_p->Close();
  delete outFile_p;
  
  std::cout << "GDJNTUPLETOHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNTupleToHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjNTupleToHist(argv[1]);
  return retVal;
}
