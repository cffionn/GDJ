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
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjNTupleToHist(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;

  globalDebugHandler gDebug;
  bool doGlobalDebug = gDebug.GetDoGlobalDebug();

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
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const bool isPP = std::stoi(config.GetConfigVal("ISPP"));
  const Int_t nMaxSubBins = 10;
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

  const Int_t nMaxEtaPhiBins = 100;
  const Int_t nEtaBins = std::stoi(config.GetConfigVal("NETABINS"));
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nEtaBins, "NETABINS")) return 1;
  const Float_t etaBinsLow = std::stof(config.GetConfigVal("ETABINSLOW"));
  const Float_t etaBinsHigh = std::stof(config.GetConfigVal("ETABINSHIGH"));
  const Bool_t etaBinsDoAbs = std::stoi(config.GetConfigVal("ETABINSDOABS"));
  if(etaBinsDoAbs && etaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsLow \'" << etaBinsLow << "\' less than 0 despite requested etaBinsDoAbs \'" << etaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  Double_t etaBins[nMaxEtaPhiBins+1];
  getLinBins(etaBinsLow, etaBinsHigh, nEtaBins, etaBins);

  const Int_t nPhiBins = std::stoi(config.GetConfigVal("NPHIBINS"));
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nPhiBins, "NPHIBINS")) return 1;
  Double_t phiBins[nMaxEtaPhiBins+1];
  getLinBins(-TMath::Pi()+0.01, TMath::Pi()+0.01, nPhiBins, phiBins);

  //Pt sub bins handling
  const Int_t nPtBinsSub = std::stoi(config.GetConfigVal("NPTBINSSUB"));
  if(!goodBinning(inConfigFileName, nMaxSubBins, nPtBinsSub, "NPTBINSSUB")) return 1;
  const Float_t ptBinsSubLow = std::stof(config.GetConfigVal("PTBINSLOW"));
  const Float_t ptBinsSubHigh = std::stof(config.GetConfigVal("PTBINSHIGH"));
  if(ptBinsSubLow < ptBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains ptBinsSubLow \'" << ptBinsSubLow << "\' less than ptBinsLow \'" << ptBinsLow << "\'. return 1" << std::endl;
  if(ptBinsSubHigh > ptBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains ptBinsSubHigh \'" << ptBinsSubHigh << "\' greater than ptBinsHigh \'" << ptBinsHigh << "\'. return 1" << std::endl;
  if(ptBinsSubLow < ptBinsLow || ptBinsSubHigh > ptBinsHigh) return 1;
  const Bool_t ptBinsSubDoLog = std::stof(config.GetConfigVal("PTBINSDOLOG"));
  Double_t ptBinsSub[nMaxSubBins+1];
  if(ptBinsSubDoLog) getLogBins(ptBinsSubLow, ptBinsSubHigh, nPtBinsSub, ptBinsSub);
  else getLinBins(ptBinsSubLow, ptBinsSubHigh, nPtBinsSub, ptBinsSub);
  std::vector<std::string> ptBinsSubStr;
  for(Int_t pI = 0; pI < nPtBinsSub; ++pI){
    ptBinsSubStr.push_back("Pt" + prettyString(ptBinsSub[pI], 1, true) + "to" + prettyString(ptBinsSub[pI+1], 1, true));
  }
  ptBinsSubStr.push_back("Pt" + prettyString(ptBinsSub[0], 1, true) + "to" + prettyString(ptBinsSub[nPtBinsSub], 1, true));

  //Eta sub bins handling
  const Int_t nEtaBinsSub = std::stoi(config.GetConfigVal("NETABINSSUB"));
  if(!goodBinning(inConfigFileName, nMaxSubBins, nEtaBinsSub, "NETABINSSUB")) return 1;
  const Float_t etaBinsSubLow = std::stof(config.GetConfigVal("ETABINSLOW"));
  const Float_t etaBinsSubHigh = std::stof(config.GetConfigVal("ETABINSHIGH"));
  if(etaBinsSubLow < etaBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than etaBinsLow \'" << etaBinsLow << "\'. return 1" << std::endl;
  if(etaBinsSubHigh > etaBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubHigh \'" << etaBinsSubHigh << "\' greater than etaBinsHigh \'" << etaBinsHigh << "\'. return 1" << std::endl;
  if(etaBinsSubLow < etaBinsLow || etaBinsSubHigh > etaBinsHigh) return 1;
  const	Bool_t etaBinsSubDoAbs = std::stoi(config.GetConfigVal("ETABINSSUBDOABS"));
  if(etaBinsSubDoAbs && etaBinsSubLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than 0 despite requested etaBinsSubDoAbs \'" << etaBinsSubDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t etaBinsSub[nMaxSubBins+1];
  getLinBins(etaBinsSubLow, etaBinsSubHigh, nEtaBinsSub, etaBinsSub);
  std::vector<std::string> etaBinsSubStr;
  std::string preStr = "";
  if(etaBinsSubDoAbs) preStr = "Abs";
  for(Int_t eI = 0; eI < nEtaBinsSub; ++eI){
    etaBinsSubStr.push_back(preStr + "Eta" + prettyString(etaBinsSub[eI], 2, true) + "to" + prettyString(etaBinsSub[eI+1], 2, true));
  }
  etaBinsSubStr.push_back(preStr + "Eta" + prettyString(etaBinsSub[0], 2, true) + "to" + prettyString(etaBinsSub[nEtaBinsSub], 2, true));
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* centrality_p = nullptr;
  TH1F* photonPtVCentEta_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonEtaVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH1F* photonPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  TH2F* photonEtaPt_p[nMaxCentBins];
  TH2F* photonEtaPhiVCentPt_p[nMaxCentBins][nMaxSubBins+1];
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!isPP){
    centrality_p = new TH1F("centrality_h", ";Centrality (%);Counts", 100, -0.5, 99.5);
    centerTitles(centrality_p);
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      photonPtVCentEta_p[cI][eI] = new TH1F(("photonPtVCentEta" + centBinsStr[cI] + "_" + etaBinsSubStr[eI]+ "_h").c_str(), ";#gamma p_{T} [GeV];Counts", nPtBins, ptBins);
      centerTitles(photonPtVCentEta_p[cI][eI]);
    }

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      photonEtaVCentPt_p[cI][pI] = new TH1F(("photonEtaVCentPt" + centBinsStr[cI] + "_" + ptBinsSubStr[pI] + "_h").c_str(), ";#gamma #eta;Counts", nEtaBins, etaBins);
      photonPhiVCentPt_p[cI][pI] = new TH1F(("photonPhi_" + centBinsStr[cI] + "_" + ptBinsSubStr[pI] + "_h").c_str(), ";#gamma #phi;Counts", nPhiBins, phiBins);

      centerTitles({photonEtaVCentPt_p[cI][pI], photonPhiVCentPt_p[cI][pI]});
    }
    
    photonEtaPt_p[cI] = new TH2F(("photonEtaPt_" + centBinsStr[cI] + "_h").c_str(), ";#gamma #eta;#gamma p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins);

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI] = new TH2F(("photonEtaPhi_" + centBinsStr[cI] + "_" + ptBinsSubStr[pI] + "_h").c_str(), ";#gamma #eta;#gamma #phi", nEtaBins, etaBins, nPhiBins, phiBins);
      centerTitles(photonEtaPhiVCentPt_p[cI][pI]);
    }
    
    centerTitles(photonEtaPt_p[cI]);
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TFile* inFile_p = new TFile(inROOTFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  Float_t fcalA_et, fcalC_et;
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;

  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  
  inTree_p->SetBranchStatus("*", 0);
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
  }
  inTree_p->SetBranchStatus("photon_pt", 1);
  inTree_p->SetBranchStatus("photon_eta", 1);
  inTree_p->SetBranchStatus("photon_phi", 1);
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_pt", 1);
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_eta", 1);
  inTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_phi", 1);
  
  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
  }
  inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
  inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
  inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);

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
      if(photon_pt_p->at(pI) < ptBins[0]) continue;
      if(photon_pt_p->at(pI) >= ptBins[nPtBins]) continue;
      if(photon_eta_p->at(pI) <= etaBins[0]) continue;
      if(photon_eta_p->at(pI) >= etaBins[nEtaBins]) continue;

      Int_t ptPos = ghostPos(nPtBinsSub, ptBinsSub, photon_pt_p->at(pI), false);
      Int_t etaPos = ghostPos(nEtaBinsSub, etaBinsSub, photon_eta_p->at(pI), false);

      if(etaPos >= 0){
	photonPtVCentEta_p[centPos][etaPos]->Fill(photon_pt_p->at(pI));
	photonPtVCentEta_p[centPos][nEtaBinsSub]->Fill(photon_pt_p->at(pI));
      }

      if(ptPos >= 0){
	photonEtaVCentPt_p[centPos][ptPos]->Fill(photon_eta_p->at(pI));
	photonPhiVCentPt_p[centPos][ptPos]->Fill(photon_phi_p->at(pI));
	photonEtaPhiVCentPt_p[centPos][ptPos]->Fill(photon_eta_p->at(pI), photon_phi_p->at(pI));

      	photonEtaVCentPt_p[centPos][nPtBinsSub]->Fill(photon_eta_p->at(pI));
	photonPhiVCentPt_p[centPos][nPtBinsSub]->Fill(photon_phi_p->at(pI));
	photonEtaPhiVCentPt_p[centPos][nPtBinsSub]->Fill(photon_eta_p->at(pI), photon_phi_p->at(pI));
      }
      
      photonEtaPt_p[centPos]->Fill(photon_eta_p->at(pI), photon_pt_p->at(pI));
    }
  }  
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  if(!isPP) centrality_p->Write("", TObject::kOverwrite);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      photonPtVCentEta_p[cI][eI]->Write("", TObject::kOverwrite);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      photonEtaVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
      photonPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
    }
    
    photonEtaPt_p[cI]->Write("", TObject::kOverwrite);

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      photonEtaPhiVCentPt_p[cI][pI]->Write("", TObject::kOverwrite);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }

  
  if(!isPP) delete centrality_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t eI = 0; eI < nEtaBinsSub+1; ++eI){
      delete photonPtVCentEta_p[cI][eI];
    }

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      delete photonEtaVCentPt_p[cI][pI];
      delete photonPhiVCentPt_p[cI][pI];
    }

    delete photonEtaPt_p[cI];

    for(Int_t pI = 0; pI < nPtBinsSub+1; ++pI){
      delete photonEtaPhiVCentPt_p[cI][pI];
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form          
  TEnv configEnv; // We will convert to tenv and store in file
  for(auto const & val : configMap){
    configEnv.SetValue(val.first.c_str(), val.second.c_str()); //Fill out the map
  }
  configEnv.Write("config", TObject::kOverwrite);  
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
