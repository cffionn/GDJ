//Author: Chris McGinn (2020.07.30)
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
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int gdjNTupleToMBHist(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"INDATAFILENAME",
					"OUTFILENAME",
					"NCENTBINS",
					"CENTBINSLOW",
					"CENTBINSHIGH",
					"NFCALETBINS",
					"FCALETBINSLOW",
					"FCALETBINSHIGH",
					"FCALETBINSDOLOG",
					"FCALETBINSDOLIN",
					"JETR",
					"NJETPTBINS",
					"JETPTBINSLOW",
					"JETPTBINSHIGH",
					"JETETABINSLOW",
					"JETETABINSHIGH",
					"JETETABINSDOABS"};
  
  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;
  std::string outFileName = inConfig_p->GetValue("OUTFILENAME", "");
  std::string dateStr = getDateStr();

  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  outFileName.replace(outFileName.rfind("."), outFileName.size(), "");
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE"); 

  //Grab the jet R
  const int jetR = inConfig_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);
  
  //Grab binnings from config
  const Int_t nMaxCentBins = 400; //will almost never exceed
  const Int_t nCentBins = inConfig_p->GetValue("NCENTBINS", 0);
  if(nCentBins <= 0 || nCentBins > nMaxCentBins){
    std::cout << "NCENTBINS \'" << nCentBins << "\' given in config \'" << inConfigFileName << "\' is not valid, please give a value > 0 and <= " << nMaxCentBins << ". return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Float_t centBinsLow = inConfig_p->GetValue("CENTBINSLOW", 0);
  const Float_t centBinsHigh = inConfig_p->GetValue("CENTBINSHIGH", 0);
  if(centBinsHigh <= centBinsLow){
    std::cout << "CENTBINSHIGH \'" << centBinsHigh << "\' is less than or equal to CENTBINSLOW \'" << centBinsLow << "\' as specified in \'" << inConfigFileName << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nFCalEtBins = inConfig_p->GetValue("NFCALETBINS", 0);
  if(nFCalEtBins <= 0 || nFCalEtBins > nMaxCentBins){
    std::cout << "NCENTBINS \'" << nFCalEtBins << "\' given in config \'" << inConfigFileName << "\' is not valid, please give a value > 0 and <= " << nMaxCentBins << ". return 1" << std::endl;
    return 1;
  }
  
  const Float_t fcalEtBinsLow = inConfig_p->GetValue("FCALETBINSLOW", 0);
  const Float_t fcalEtBinsHigh = inConfig_p->GetValue("FCALETBINSHIGH", 0);
  if(fcalEtBinsHigh <= fcalEtBinsLow){
    std::cout << "FCALETBINSHIGH \'" << fcalEtBinsHigh << "\' is less than or equal to FCALETBINSLOW \'" << fcalEtBinsLow << "\' as specified in \'" << inConfigFileName << "\'. return 1" << std::endl;
    return 1;
  }

  const Bool_t fcalEtBinsDoLin = inConfig_p->GetValue("FCALETBINSDOLIN", 0);
  const Bool_t fcalEtBinsDoLog = inConfig_p->GetValue("FCALETBINSDOLOG", 0);

  if((fcalEtBinsDoLin && fcalEtBinsDoLog) || (!fcalEtBinsDoLin && !fcalEtBinsDoLog)){
    std::cout << "Both FCALETBINSDOLIN, " << fcalEtBinsDoLin << ", and FCALETBINSDOLOG, " << fcalEtBinsDoLog << ", have the same value from input file \'" << inConfigFileName << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nJetPtBins = inConfig_p->GetValue("NJETPTBINS", 20);
  const Float_t jetPtBinsLow = inConfig_p->GetValue("JETPTBINSLOW", 30);
  const Float_t jetPtBinsHigh = inConfig_p->GetValue("JETPTBINSHIGH", 50);

  //  const Int_t nJetEtaBins = inConfig_p->GetValue("NJETETABINS", 28);
  const Float_t jetEtaBinsLow = inConfig_p->GetValue("JETETABINSLOW", -2.8);
  const Float_t jetEtaBinsHigh = inConfig_p->GetValue("JETETABINSHIGH", 2.8);
  const Bool_t jetEtaBinsDoAbs = inConfig_p->GetValue("JETETABINSDOABS", 0);

  if(jetEtaBinsDoAbs && jetEtaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains jetEtaBinsLow \'" << jetEtaBinsLow << "\' less than 0 despite requested jetEtaBinsDoAbs \'" << jetEtaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }

  
  Double_t centBins[nMaxCentBins+1];
  getLinBins(centBinsLow, centBinsHigh, nCentBins, centBins);

  Double_t fcalEtBins[nMaxCentBins+1];
  if(fcalEtBinsDoLin) getLinBins(fcalEtBinsLow, fcalEtBinsHigh, nFCalEtBins, fcalEtBins);
  else getLogBins(fcalEtBinsLow, fcalEtBinsHigh, nFCalEtBins, fcalEtBins);

  Double_t jetPtBins[nMaxCentBins+1];
  getLinBins(jetPtBinsLow, jetPtBinsHigh, nJetPtBins, jetPtBins);

  TH1F* cent_p = new TH1F("cent_h", ";Centrality (%); Counts", nCentBins, centBins);
  TH1F* centData_p = new TH1F("centData_h", ";Centrality (%); Counts", nCentBins, centBins);
  TH1F* fcalEt_p = new TH1F("fcalEt_h", ";#Sigma_{FCal} E_{T};Counts", nFCalEtBins, fcalEtBins);
  TH2F* centVFCalEt_p = new TH2F("centVFCalEt_h", ";Centrality (%);#Sigma_{FCal} E_{T}", nCentBins, centBins, nFCalEtBins, fcalEtBins);

  TH1F* jtMultPerCent_p[nMaxCentBins];
  std::map<unsigned long long, unsigned long long> centCounter;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "Cent" + std::to_string((Int_t)centBins[cI]) + "to" + std::to_string((Int_t)centBins[cI+1]);
    std::string centStr2 = std::to_string((Int_t)centBins[cI]) + "-" + std::to_string((Int_t)centBins[cI+1]) + "%";
    jtMultPerCent_p[cI] = new TH1F(("jtMultPerCent_" + centStr + "_h").c_str(), (";Reco. Jet p_{T} [GeV];N_{Jet,R=" + jetRStr + "}/N_{Event," +  centStr2 + "}").c_str(), nJetPtBins, jetPtBins);

    jtMultPerCent_p[cI]->Sumw2();
    centCounter[cI] = 0;
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string inDataFileName = inConfig_p->GetValue("INDATAFILENAME", "");
  if(!check.checkFileExt(inDataFileName, ".root")) return 1;
  
  TFile* inFile_p = new TFile(inDataFileName.c_str(), "READ");
  TEnv* inFileConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> inFileReqParams = {"ISPP",
					      "ISCOLLISIONS"};
  if(!checkEnvForParams(inFileConfig_p, inFileReqParams)) return 1; 

  bool isPP = inFileConfig_p->GetValue("ISPP", 0);
  bool isCollisions = inFileConfig_p->GetValue("ISCOLLISIONS", 0);
  
  if(isPP){
    std::cout << "Given input \'" << inDataFileName << "\' is not Pb+Pb. return 1" << std::endl;
    return 1;
  }
  if(!isCollisions){
    std::cout << "Given input \'" << inDataFileName << "\' is not data. return 1" << std::endl;
    return 1;
  }

  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
  Int_t cent_;
  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("cent", 1);

  inTree_p->SetBranchAddress("cent", &cent_);

  for(ULong64_t entry = 0; entry < (ULong64_t)inTree_p->GetEntries(); ++entry){
    inTree_p->GetEntry(entry);
    centData_p->Fill(cent_);
  }
  
  inFile_p->Close();
  delete inFile_p;
  
  inFile_p = new TFile(inFileName.c_str(), "READ");
  inFileConfig_p = (TEnv*)inFile_p->Get("config");
  if(!checkEnvForParams(inFileConfig_p, inFileReqParams)) return 1; 

  isPP = inFileConfig_p->GetValue("ISPP", 0);
  isCollisions = inFileConfig_p->GetValue("ISCOLLISIONS", 0);

  if(isPP){
    std::cout << "Given input \'" << inFileName << "\' is not Pb+Pb. return 1" << std::endl;
    return 1;
  }
  if(!isCollisions){
    std::cout << "Given input \'" << inFileName << "\' is not data. return 1" << std::endl;
    return 1;
  }

  inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  Float_t fcalA_et_, fcalC_et_;
  std::vector<float>* aktRhi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_eta_p=nullptr;
  //  std::vector<float>* aktRhi_em_xcalib_jet_phi_p=nullptr;

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("cent", 1);
  inTree_p->SetBranchStatus("fcalA_et", 1);
  inTree_p->SetBranchStatus("fcalC_et", 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
  //  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);

  inTree_p->SetBranchAddress("cent", &cent_);
  inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
  inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
  //  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);


  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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

    unsigned long long centPos = ghostPos(nCentBins, centBins, cent_); 
    if(centPos > (unsigned long long)nCentBins) continue;
    
    cent_p->Fill(cent_);
    fcalEt_p->Fill(fcalA_et_ + fcalC_et_);
    centVFCalEt_p->Fill(cent_, fcalA_et_ + fcalC_et_);

    ++(centCounter[centPos]);
    
    if(doPrint && doGlobalDebug) std::cout << "NJET: " << aktRhi_em_xcalib_jet_pt_p->size() << std::endl;
    for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
      if(doPrint && doGlobalDebug) std::cout << " PT " << jI << ": " << aktRhi_em_xcalib_jet_pt_p->at(jI) << std::endl;      
      
      if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jetPtBinsLow) continue;
      if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jetPtBinsHigh) continue;

      double jetEtaVal = aktRhi_em_xcalib_jet_eta_p->at(jI);
      if(jetEtaBinsDoAbs) jetEtaVal = TMath::Abs(jetEtaVal);

      if(jetEtaVal < jetEtaBinsLow) continue;
      if(jetEtaVal >= jetEtaBinsHigh) continue;
      
      if(doPrint && doGlobalDebug) std::cout << "  JET FILLS" << std::endl;
      jtMultPerCent_p[centPos]->Fill(aktRhi_em_xcalib_jet_pt_p->at(jI));
    }
  }
  std::cout << "Processing complete!" << std::endl;

  /*
  for(unsigned long long cI = 0; cI < (unsigned long long)nCentBins; ++cI){
    std::cout << "CENT " << centBins[cI] << "-" << centBins[cI+1] << ": " << centCounter[cI] << std::endl;
  }
  */

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
 
  cent_p->Write("", TObject::kOverwrite);
  delete cent_p;

  centData_p->Write("", TObject::kOverwrite);
  delete centData_p;

  fcalEt_p->Write("", TObject::kOverwrite);
  delete fcalEt_p;

  centVFCalEt_p->Write("", TObject::kOverwrite);
  delete centVFCalEt_p;

  for(unsigned long long cI = 0; cI < (unsigned long long)nCentBins; ++cI){
    jtMultPerCent_p[cI]->Scale(1./(double)centCounter[cI]);
    jtMultPerCent_p[cI]->Write("", TObject::kOverwrite);
    delete jtMultPerCent_p[cI];
  }
  
  inConfig_p->Write("config", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete inConfig_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "GDJNTUPLETOMBHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNTupleToMBHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjNTupleToMBHist(argv[1]);
  return retVal;
}
