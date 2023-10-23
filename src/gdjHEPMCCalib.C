//Author: Chris McGinn (2022.11.21)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/photonUtil.h"
#include "include/stringUtil.h"
#include "include/varUtil.h"

int gdjHEPMCCalib(std::string inConfigFileName)
{
  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Global debug handler - do 'export DOGLOBALDEBUGROOT=1' at command line to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INFILENAME",
					      "OUTFILENAME",
					      "INHEPDATAFILENAME"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  std::string inFileName = config_p->GetValue("INFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string inHEPDataFileName = config_p->GetValue("INHEPDATAFILENAME", "");
  
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  if(!check.checkFileExt(inHEPDataFileName, ".root")) return 1;

  TFile* hepFile_p = new TFile(inHEPDataFileName.c_str(), "READ");
  TH1D* hepRAAHist_p = (TH1D*)hepFile_p->Get("Table 19/Hist1D_y1");
  TH1D* hepRAAHist_SysUp_p = (TH1D*)hepFile_p->Get("Table 19/Hist1D_y1_e1plus");
  TH1D* hepRAAHist_SysDown_p = (TH1D*)hepFile_p->Get("Table 19/Hist1D_y1_e1minus");
  TH1D* hepRAAHist_Stat_p = (TH1D*)hepFile_p->Get("Table 19/Hist1D_y1_e2");

  const Int_t nMaxBins = 1000;
  Int_t nBins = hepRAAHist_p->GetNbinsX();
  if(nBins > nMaxBins){
    std::cout << "Needed nBins \'" << nBins << "\' exceeds max bins \'" << nMaxBins << "\'. please fix. return 1" << std::endl;
    return 1;
  }
  Double_t bins[nMaxBins+1];
  for(Int_t bIX = 0; bIX < hepRAAHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    bins[bIX] = hepRAAHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
  }   
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* jewelTree_p = (TTree*)inFile_p->Get("jewelTree");
  TEnv* inFileConfig_p = (TEnv*)inFile_p->Get("config");  

  const Int_t nMaxJetRs = 7;
  std::string jetRStr = inFileConfig_p->GetValue("JETRS", "");
  std::vector<int> jetRVect = strToVectI(jetRStr);
  if(jetRVect.size() <= 0 || jetRVect.size() > (unsigned int)nMaxJetRs){
    std::cout << "Jet R query returned \'" << jetRStr << "\', either 0 valid jet Rs or in excess of cap, " << nMaxJetRs << ". return 1" << std::endl;
    return 1;
  }
  
  //Add some parameters from infileconfig to config for writeout
  config_p->SetValue("CENT", inFileConfig_p->GetValue("CENT", ""));
  config_p->SetValue("JETRS", jetRStr.c_str());
  
  Double_t evtWeight_;
  
  const Int_t nMaxPart = 10000;
  Int_t nPart_;
  Float_t pt_[nMaxPart];
  Float_t eta_[nMaxPart];
  Float_t phi_[nMaxPart];
  Int_t pid_[nMaxPart];

  const Int_t nMaxJets = 100;
  Int_t nJt_[nMaxJetRs];
  Float_t jtpt_[nMaxJetRs][nMaxJets];
  Float_t jtphi_[nMaxJetRs][nMaxJets];
  Float_t jteta_[nMaxJetRs][nMaxJets];
  Float_t jtm_[nMaxJetRs][nMaxJets];

  jewelTree_p->SetBranchStatus("*", 0);  
  jewelTree_p->SetBranchStatus("evtWeight", 1);
  jewelTree_p->SetBranchStatus("nPart", 1);  
  jewelTree_p->SetBranchStatus("pt", 1);  
  jewelTree_p->SetBranchStatus("eta", 1);  
  jewelTree_p->SetBranchStatus("phi", 1);  
  jewelTree_p->SetBranchStatus("pid", 1);  

  jewelTree_p->SetBranchAddress("evtWeight", &evtWeight_);
  jewelTree_p->SetBranchAddress("nPart", &nPart_);
  jewelTree_p->SetBranchAddress("pt", pt_);
  jewelTree_p->SetBranchAddress("eta", eta_);
  jewelTree_p->SetBranchAddress("phi", phi_);
  jewelTree_p->SetBranchAddress("pid", pid_);

  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    jewelTree_p->SetBranchStatus(("nJtR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtptR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtetaR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtphiR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtmR" + std::to_string(jetRVect[rI])).c_str(), 1);

    jewelTree_p->SetBranchAddress(("nJtR" + std::to_string(jetRVect[rI])).c_str(), &nJt_[rI]);
    jewelTree_p->SetBranchAddress(("jtptR" + std::to_string(jetRVect[rI])).c_str(), jtpt_[rI]);
    jewelTree_p->SetBranchAddress(("jtetaR" + std::to_string(jetRVect[rI])).c_str(), jteta_[rI]);
    jewelTree_p->SetBranchAddress(("jtphiR" + std::to_string(jetRVect[rI])).c_str(), jtphi_[rI]);
    jewelTree_p->SetBranchAddress(("jtmR" + std::to_string(jetRVect[rI])).c_str(), jtm_[rI]);
  }
  
  //Check the outrootfile directory exists or add directory path to outfilename
  if(outFileName.find("/") == std::string::npos){
    outFileName = outputStr + "/" + outFileName;
  }
  else{
    std::string dirName = outFileName.substr(0, outFileName.rfind("/"));
    check.doCheckMakeDir(dirName);
  }
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* histPt_p[nMaxJetRs];
  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    histPt_p[rI] = new TH1D(("jetPt_R" + std::to_string(jetRVect[rI]) + "_h").c_str(), ";Jet p_{T};#frac{1}{N_{Evt}} #frac{dN_{Jet}}{dp_{T}}", nBins, bins);
    histPt_p[rI]->Sumw2();
  }
  
  const ULong64_t nEntries = jewelTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry % nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    jewelTree_p->GetEntry(entry);

    for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
      for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	if(jtpt_[rI][jI] < bins[0]) continue;
	if(jtpt_[rI][jI] > bins[nBins]) continue;
	if(TMath::Abs(jteta_[rI][jI]) > 2.8) continue;
	
	histPt_p[rI]->Fill(jtpt_[rI][jI], evtWeight_);
      }
    }
  }
  std::cout << "Event processing complete." << std::endl;    
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  
  config_p->Write("config", TObject::kOverwrite);

  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    histPt_p[rI]->Write("", TObject::kOverwrite);
    delete histPt_p[rI];
  }
  
  hepRAAHist_p->Write("", TObject::kOverwrite);  
  hepRAAHist_SysUp_p->Write("", TObject::kOverwrite);  
  hepRAAHist_SysDown_p->Write("", TObject::kOverwrite);  
  hepRAAHist_Stat_p->Write("", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;
  
  hepFile_p->Close();
  delete hepFile_p;
  
  //Clean all new's
  delete config_p;  

  std::cout << "GDJHEPMCCALIB complete. return 0" << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHEPMCCalib.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHEPMCCalib(argv[1]);
  return retVal;
}
