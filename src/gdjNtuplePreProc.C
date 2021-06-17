//Author: Chris McGinn (2020.02.26)
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
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/keyHandler.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/plotUtilities.h"
#include "include/returnFileList.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"


int gdjNtuplePreProc(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  //  configParser config(inConfig_p);
  std::vector<std::string> necessaryParams = {"MCPREPROCDIRNAME",
					      "OUTFILENAME",
					      "CENTFILENAME",
					      "ISPP",
					      "ISMC",
					      "KEEPTRUTH",
					      "DOMINGAMMAPT",
					      "MINGAMMAPT"};  
  if(!checkEnvForParams(inConfig_p, necessaryParams)) return 1;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const bool isPP = inConfig_p->GetValue("ISPP", 0);
  const bool isMC = inConfig_p->GetValue("ISMC", 0);
  const bool keepTruth = inConfig_p->GetValue("KEEPTRUTH", 0);

  const bool doMinGammaPt = inConfig_p->GetValue("DOMINGAMMAPT", 0);
  const double minGammaPt = inConfig_p->GetValue("MINGAMMAPT", -1.0);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string inDirStr = inConfig_p->GetValue("MCPREPROCDIRNAME", "");
  const std::string inCentFileName = inConfig_p->GetValue("CENTFILENAME", "");
  if(!check.checkDir(inDirStr)){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }
  if(!check.checkFileExt(inCentFileName, "txt")) return 1;

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();

  std::string weirdStr = "";

  //NColl weights
  std::vector<double> ncollWeights;
  std::vector<double> centCounts;
  for(Int_t cI = 0; cI < 100; ++cI){
    ncollWeights.push_back(findAvgNColl_Cent(cI, cI+1));
    centCounts.push_back(0.0);
  }
  for(Int_t cI = 0; cI < 100; ++cI){
    ncollWeights[99-cI] /= ncollWeights[0];
  }  

  std::vector<std::string> fileList = returnFileList(inDirStr, ".root");
  if(fileList.size() == 0){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' contains no root files. return 1" << std::endl;
    return 1;
  }
  TFile* inFile_p = new TFile(fileList[0].c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");
  
  //  configParser inConfigGlobal(fileConfig_p);

  inFile_p->Close();
  delete inFile_p;

  bool getTracks = false;
  if(checkEnvForParams(fileConfig_p, {"GETTRACKS"})) getTracks = fileConfig_p->GetValue("GETTRACKS", 0);

  std::string topOutDir = "output";
  if(checkEnvForParams(inConfig_p, {"OUTDIRNAME"})){
    topOutDir = inConfig_p->GetValue("OUTDIRNAME", "");  
    if(!check.checkDir(topOutDir)){
      std::cout << "Given output directory \'" << topOutDir << "\' does not exist. Please create it or remove param to use local default \'output\' directory. return 1" << std::endl;
      return 1;
    }
  }

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir(topOutDir);
  check.doCheckMakeDir(topOutDir + "/" + dateStr);

  std::string outFileName = inConfig_p->GetValue("OUTFILENAME", "");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = topOutDir + "/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("gammaJetTree_p", "");
  
  std::vector<std::string> outBranchesToAdd = {"cent",
					       "sampleTag",
					       "sampleWeight",
					       "ncollWeight",
					       "fullWeight",
					       "truthPhotonPt",
					       "truthPhotonEta",
					       "truthPhotonPhi",
					       "truthPhotonType",
					       "truthPhotonOrigin",
					       "truthPhotonIso2",
					       "truthPhotonIso3",
					       "truthPhotonIso4"};

  Int_t runNumber_, eventNumber_;
  UInt_t lumiBlock_;
  Bool_t passesToroid_;

  Float_t pthat_;
  Int_t cent_;
  Int_t sampleTag_;
  Float_t sampleWeight_;
  Float_t ncollWeight_;
  Float_t fullWeight_;

  Bool_t is_pileup_, is_oo_pileup_;
  
  const Int_t nTreeParton_ = 2;
  Float_t treePartonPt_[nTreeParton_];
  Float_t treePartonEta_[nTreeParton_];
  Float_t treePartonPhi_[nTreeParton_];
  Int_t treePartonId_[nTreeParton_];

  Float_t actualInteractionsPerCrossing_, averageInteractionsPerCrossing_;
  Int_t nvert_;
  std::vector<float>* vert_x_p=nullptr;
  std::vector<float>* vert_y_p=nullptr;
  std::vector<float>* vert_z_p=nullptr;
  std::vector<int>* vert_type_p=nullptr;
  std::vector<int>* vert_ntrk_p=nullptr;

  Float_t fcalA_et_, fcalC_et_, fcalA_et_Cos2_, fcalC_et_Cos2_, fcalA_et_Sin2_, fcalC_et_Sin2_, fcalA_et_Cos3_, fcalC_et_Cos3_, fcalA_et_Sin3_, fcalC_et_Sin3_, fcalA_et_Cos4_, fcalC_et_Cos4_, fcalA_et_Sin4_, fcalC_et_Sin4_, evtPlane2Phi_, evtPlane3Phi_, evtPlane4Phi_;

  Int_t ntrk_;
  std::vector<float> *trk_pt_p=nullptr;
  std::vector<float> *trk_eta_p=nullptr;
  std::vector<float> *trk_phi_p=nullptr;
  std::vector<float> *trk_charge_p=nullptr;
  std::vector<bool> *trk_tightPrimary_p=nullptr;
  std::vector<bool> *trk_HILoose_p=nullptr;
  std::vector<bool> *trk_HITight_p=nullptr;
  std::vector<float> *trk_d0_p=nullptr;
  std::vector<float> *trk_z0_p=nullptr;
  std::vector<float> *trk_vz_p=nullptr;
  std::vector<float> *trk_theta_p=nullptr;
  std::vector<int> *trk_nPixelHits_p=nullptr;
  std::vector<int> *trk_nSCTHits_p=nullptr;
  std::vector<int> *trk_nBlayerHits_p=nullptr;
  
  Int_t truth_n_;
  std::vector<float>* truth_charge_p=nullptr;
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_e_p=nullptr;
  std::vector<float>* truth_pdg_p=nullptr;
  std::vector<int>* truth_type_p=nullptr;
  std::vector<int>* truth_origin_p=nullptr;
  std::vector<int>* truth_status_p=nullptr;

  Int_t truthOut_n_;
  std::vector<float>* truthOut_charge_p=new std::vector<float>;
  std::vector<float>* truthOut_pt_p=new std::vector<float>;
  std::vector<float>* truthOut_eta_p=new std::vector<float>;
  std::vector<float>* truthOut_phi_p=new std::vector<float>;
  std::vector<float>* truthOut_e_p=new std::vector<float>;
  std::vector<float>* truthOut_pdg_p=new std::vector<float>;
  std::vector<int>* truthOut_type_p=new std::vector<int>;
  std::vector<int>* truthOut_origin_p=new std::vector<int>;
  std::vector<int>* truthOut_status_p=new std::vector<int>;

  Float_t truthPhotonPt_;
  Float_t truthPhotonEta_;
  Float_t truthPhotonPhi_;
  Int_t truthPhotonType_;
  Int_t truthPhotonOrigin_;
  Float_t truthPhotonIso2_;
  Float_t truthPhotonIso3_;
  Float_t truthPhotonIso4_;
  
  Int_t akt2hi_jet_n_;
  const Int_t nJES = 18;
  std::vector<std::vector<float>* > akt2hi_etajes_jet_pt_sys_JES_p;
  if(isMC){
    akt2hi_etajes_jet_pt_sys_JES_p.reserve(nJES);
    for(unsigned int eI = 0; eI < nJES; ++eI){
      akt2hi_etajes_jet_pt_sys_JES_p[eI] = nullptr;
    }
  }

  const Int_t nJER = 9;
  std::vector<std::vector<float>* > akt2hi_etajes_jet_pt_sys_JER_p;
  if(isMC){
    akt2hi_etajes_jet_pt_sys_JER_p.reserve(nJER);
    for(unsigned int eI = 0; eI < nJER; ++eI){
      akt2hi_etajes_jet_pt_sys_JER_p[eI] = nullptr;
    }
  }

  std::vector<float>* akt2hi_etajes_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_etajes_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_etajes_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_etajes_jet_e_p=nullptr;
  std::vector<float>* akt2hi_etajes_jet_m_p=nullptr;
  std::vector<float>* akt2hi_insitu_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_insitu_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_insitu_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_insitu_jet_e_p=nullptr;
  std::vector<float>* akt2hi_insitu_jet_m_p=nullptr;

  std::vector<bool>* akt2hi_jet_clean_p=nullptr;
  std::vector<int>* akt2hi_truthpos_p=nullptr;
  
  Int_t akt4hi_jet_n_;
  std::vector<float>* akt4hi_etajes_jet_pt_p=nullptr;
  std::vector<std::vector<float>* > akt4hi_etajes_jet_pt_sys_JES_p;
  if(isMC){
    akt4hi_etajes_jet_pt_sys_JES_p.reserve(nJES);
    for(unsigned int eI = 0; eI < nJES; ++eI){
      akt4hi_etajes_jet_pt_sys_JES_p[eI] = nullptr;
    }
  }
  std::vector<std::vector<float>* > akt4hi_etajes_jet_pt_sys_JER_p;
  if(isMC){
    akt4hi_etajes_jet_pt_sys_JER_p.reserve(nJER);
    for(unsigned int eI = 0; eI < nJER; ++eI){
      akt4hi_etajes_jet_pt_sys_JER_p[eI] = nullptr;
    }
  }

  std::vector<float>* akt4hi_etajes_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_etajes_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_etajes_jet_e_p=nullptr;
  std::vector<float>* akt4hi_etajes_jet_m_p=nullptr;
  std::vector<float>* akt4hi_insitu_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_insitu_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_insitu_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_insitu_jet_e_p=nullptr;
  std::vector<float>* akt4hi_insitu_jet_m_p=nullptr;

  std::vector<bool>* akt4hi_jet_clean_p=nullptr;
  std::vector<int>* akt4hi_truthpos_p=nullptr;

  
  Int_t akt10hi_jet_n_;
  std::vector<float>* akt10hi_etajes_jet_pt_p=nullptr;
  std::vector<float>* akt10hi_etajes_jet_eta_p=nullptr;
  std::vector<float>* akt10hi_etajes_jet_phi_p=nullptr;
  std::vector<float>* akt10hi_etajes_jet_e_p=nullptr;
  std::vector<float>* akt10hi_etajes_jet_m_p=nullptr;
  std::vector<float>* akt10hi_insitu_jet_pt_p=nullptr;
  std::vector<float>* akt10hi_insitu_jet_eta_p=nullptr;
  std::vector<float>* akt10hi_insitu_jet_phi_p=nullptr;
  std::vector<float>* akt10hi_insitu_jet_e_p=nullptr;
  std::vector<float>* akt10hi_insitu_jet_m_p=nullptr;

  std::vector<bool>* akt10hi_jet_clean_p=nullptr;
  std::vector<int>* akt10hi_truthpos_p=nullptr;

  
  Int_t akt2to10hi_jet_n_;
  std::vector<float>* akt2to10hi_etajes_jet_pt_p=nullptr;
  std::vector<float>* akt2to10hi_etajes_jet_eta_p=nullptr;
  std::vector<float>* akt2to10hi_etajes_jet_phi_p=nullptr;
  std::vector<float>* akt2to10hi_etajes_jet_e_p=nullptr;
  std::vector<float>* akt2to10hi_etajes_jet_m_p=nullptr;
  std::vector<int>* akt2to10hi_truthpos_p=nullptr;  

  Int_t photon_n_;
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_pt_precali_p=nullptr;
  std::vector<float>* photon_pt_sys1_p=nullptr;
  std::vector<float>* photon_pt_sys2_p=nullptr;
  std::vector<float>* photon_pt_sys3_p=nullptr;
  std::vector<float>* photon_pt_sys4_p=nullptr;  
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;
  std::vector<bool>* photon_tight_b4FudgeTool_p=nullptr;  
  std::vector<bool>* photon_loose_p=nullptr;
  std::vector<unsigned int>* photon_isem_p=nullptr;
  std::vector<int>* photon_convFlag_p=nullptr;
  std::vector<float>* photon_Rconv_p=nullptr;
  std::vector<float>* photon_etcone20_p=nullptr;
  std::vector<float>* photon_etcone30_p=nullptr;
  std::vector<float>* photon_etcone40_p=nullptr;
  std::vector<float>* photon_topoetcone20_p=nullptr;
  std::vector<float>* photon_topoetcone30_p=nullptr;
  std::vector<float>* photon_topoetcone40_p=nullptr;
  std::vector<float>* photon_Rhad1_p=nullptr;
  std::vector<float>* photon_Rhad_p=nullptr;
  std::vector<float>* photon_e277_p=nullptr;
  std::vector<float>* photon_Reta_p=nullptr;
  std::vector<float>* photon_Rphi_p=nullptr;
  std::vector<float>* photon_weta1_p=nullptr;
  std::vector<float>* photon_weta2_p=nullptr;
  std::vector<float>* photon_wtots1_p=nullptr;
  std::vector<float>* photon_f1_p=nullptr;
  std::vector<float>* photon_f3_p=nullptr;
  std::vector<float>* photon_fracs1_p=nullptr;
  std::vector<float>* photon_DeltaE_p=nullptr;
  std::vector<float>* photon_Eratio_p=nullptr;

  Int_t akt2_truth_jet_n_;
  std::vector<float>* akt2_truth_jet_pt_p=nullptr;
  std::vector<float>* akt2_truth_jet_eta_p=nullptr;
  std::vector<float>* akt2_truth_jet_phi_p=nullptr;
  std::vector<float>* akt2_truth_jet_e_p=nullptr;
  std::vector<float>* akt2_truth_jet_m_p=nullptr;
  std::vector<int>* akt2_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt2_truth_jet_recopos_p=nullptr;

  Int_t akt4_truth_jet_n_;
  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;
  std::vector<float>* akt4_truth_jet_e_p=nullptr;
  std::vector<float>* akt4_truth_jet_m_p=nullptr;
  std::vector<int>* akt4_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt4_truth_jet_recopos_p=nullptr;

  
  Int_t akt10_truth_jet_n_;
  std::vector<float>* akt10_truth_jet_pt_p=nullptr;
  std::vector<float>* akt10_truth_jet_eta_p=nullptr;
  std::vector<float>* akt10_truth_jet_phi_p=nullptr;
  std::vector<float>* akt10_truth_jet_e_p=nullptr;
  std::vector<float>* akt10_truth_jet_m_p=nullptr;
  std::vector<int>* akt10_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt10_truth_jet_recopos_p=nullptr;

  Int_t akt2to10_truth_jet_n_;
  std::vector<float>* akt2to10_truth_jet_pt_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_eta_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_phi_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_e_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_m_p=nullptr;
  std::vector<int>* akt2to10_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt2to10_truth_jet_recopos_p=nullptr;
  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outTree_p->Branch("runNumber", &runNumber_, "runNumber/I");
  outTree_p->Branch("eventNumber", &eventNumber_, "eventNumber/I");
  outTree_p->Branch("lumiBlock", &lumiBlock_, "lumiBlock/i");
  outTree_p->Branch("passesToroid", &passesToroid_, "passesToroid/i");

  if(isMC) outTree_p->Branch("pthat", &pthat_, "pthat/F");
  if(!isPP) outTree_p->Branch("cent", &cent_, "cent/I");
  outTree_p->Branch("sampleTag", &sampleTag_, "sampleTag/I");
  if(isMC) outTree_p->Branch("sampleWeight", &sampleWeight_, "sampleWeight/F");
  if(!isPP && isMC) outTree_p->Branch("ncollWeight", &ncollWeight_, "ncollWeight/F");
  if(isMC) outTree_p->Branch("fullWeight", &fullWeight_, "fullWeight/F");

  outTree_p->Branch("is_pileup", &is_pileup_, "is_pileup/O");
  outTree_p->Branch("is_oo_pileup", &is_oo_pileup_, "is_oo_pileup/O");
  
  if(isMC){
    outTree_p->Branch("treePartonPt", treePartonPt_, ("treePartonPt[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonEta", treePartonEta_, ("treePartonEta[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonPhi", treePartonPhi_, ("treePartonPhi[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonId", treePartonId_, ("treePartonId[" + std::to_string(nTreeParton_) + "]/I").c_str());
  }
  
  outTree_p->Branch("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_, "actualInteractionsPerCrossing/F");
  outTree_p->Branch("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_, "averageInteractionsPerCrossing/F");
  outTree_p->Branch("nvert", &nvert_, "nvert/I");
  
  outTree_p->Branch("vert_x", &vert_x_p);
  outTree_p->Branch("vert_y", &vert_y_p);
  outTree_p->Branch("vert_z", &vert_z_p);
  outTree_p->Branch("vert_type", &vert_type_p);
  outTree_p->Branch("vert_ntrk", &vert_ntrk_p);
  
  outTree_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  outTree_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F"); 
  outTree_p->Branch("fcalA_et_Cos2", &fcalA_et_Cos2_, "fcalA_et_Cos2/F");
  outTree_p->Branch("fcalC_et_Cos2", &fcalC_et_Cos2_, "fcalC_et_Cos2/F");
  outTree_p->Branch("fcalA_et_Sin2", &fcalA_et_Sin2_, "fcalA_et_Sin2/F");
  outTree_p->Branch("fcalC_et_Sin2", &fcalC_et_Sin2_, "fcalC_et_Sin2/F");

  outTree_p->Branch("fcalA_et_Cos3", &fcalA_et_Cos3_, "fcalA_et_Cos3/F");
  outTree_p->Branch("fcalC_et_Cos3", &fcalC_et_Cos3_, "fcalC_et_Cos3/F");
  outTree_p->Branch("fcalA_et_Sin3", &fcalA_et_Sin3_, "fcalA_et_Sin3/F");
  outTree_p->Branch("fcalC_et_Sin3", &fcalC_et_Sin3_, "fcalC_et_Sin3/F");
    
  outTree_p->Branch("fcalA_et_Cos4", &fcalA_et_Cos4_, "fcalA_et_Cos4/F");
  outTree_p->Branch("fcalC_et_Cos4", &fcalC_et_Cos4_, "fcalC_et_Cos4/F");
  outTree_p->Branch("fcalA_et_Sin4", &fcalA_et_Sin4_, "fcalA_et_Sin4/F");
  outTree_p->Branch("fcalC_et_Sin4", &fcalC_et_Sin4_, "fcalC_et_Sin4/F");
  
  outTree_p->Branch("evtPlane2Phi", &evtPlane2Phi_, "evtPlane2Phi/F");
  outTree_p->Branch("evtPlane3Phi", &evtPlane3Phi_, "evtPlane3Phi/F");
  outTree_p->Branch("evtPlane4Phi", &evtPlane4Phi_, "evtPlane4Phi/F");
  
  if(getTracks){
    outTree_p->Branch("ntrk", &ntrk_, "ntrk/I");
    outTree_p->Branch("trk_pt", &trk_pt_p);
    outTree_p->Branch("trk_eta", &trk_eta_p);
    outTree_p->Branch("trk_phi", &trk_phi_p);
    
    outTree_p->Branch("trk_charge", &trk_charge_p);
    outTree_p->Branch("trk_tightPrimary", &trk_tightPrimary_p);
    outTree_p->Branch("trk_HILoose", &trk_HILoose_p);
    outTree_p->Branch("trk_HITight", &trk_HITight_p);
    outTree_p->Branch("trk_d0", &trk_d0_p);
    outTree_p->Branch("trk_z0", &trk_z0_p);
    outTree_p->Branch("trk_vz", &trk_vz_p);
    outTree_p->Branch("trk_theta", &trk_theta_p);
    outTree_p->Branch("trk_nPixelHits", &trk_nPixelHits_p);
    outTree_p->Branch("trk_nSCTHits", &trk_nSCTHits_p);
    outTree_p->Branch("trk_nBlayerHits", &trk_nBlayerHits_p);
  }
  
  if(isMC){
    if(keepTruth){
      outTree_p->Branch("truth_n", &truthOut_n_, "truth_n/I");
      outTree_p->Branch("truth_charge", &truthOut_charge_p);
      outTree_p->Branch("truth_pt", &truthOut_pt_p);
      outTree_p->Branch("truth_eta", &truthOut_eta_p);
      outTree_p->Branch("truth_phi", &truthOut_phi_p);
      outTree_p->Branch("truth_pdg", &truthOut_pdg_p);
      outTree_p->Branch("truth_e", &truthOut_e_p);
      outTree_p->Branch("truth_type", &truthOut_type_p);
      outTree_p->Branch("truth_origin", &truthOut_origin_p);
      outTree_p->Branch("truth_status", &truthOut_status_p);
    }

    outTree_p->Branch("truthPhotonPt", &truthPhotonPt_, "truthPhotonPt/F");
    outTree_p->Branch("truthPhotonPhi", &truthPhotonPhi_, "truthPhotonPhi/F");
    outTree_p->Branch("truthPhotonEta", &truthPhotonEta_, "truthPhotonEta/F");
    outTree_p->Branch("truthPhotonType", &truthPhotonType_, "truthPhotonType/I");
    outTree_p->Branch("truthPhotonOrigin", &truthPhotonOrigin_, "truthPhotonOrigin/I");
    outTree_p->Branch("truthPhotonIso2", &truthPhotonIso2_, "truthPhotonIso2/F");
    outTree_p->Branch("truthPhotonIso3", &truthPhotonIso3_, "truthPhotonIso3/F");
    outTree_p->Branch("truthPhotonIso4", &truthPhotonIso4_, "truthPhotonIso4/F");
  }
  
  outTree_p->Branch("akt2hi_jet_n", &akt2hi_jet_n_, "akt2hi_jet_n/I");

  if(isMC){
    for(Int_t eI = 0; eI < nJES; ++eI){
      outTree_p->Branch(("akt2hi_etajes_jet_pt_sys_JES_" + std::to_string(eI)).c_str(), &(akt2hi_etajes_jet_pt_sys_JES_p[eI]));
    }
    
    for(Int_t eI = 0; eI < nJER; ++eI){
      outTree_p->Branch(("akt2hi_etajes_jet_pt_sys_JER_" + std::to_string(eI)).c_str(), &(akt2hi_etajes_jet_pt_sys_JER_p[eI]));
    }
  }

  outTree_p->Branch("akt2hi_etajes_jet_pt", &akt2hi_etajes_jet_pt_p);
  outTree_p->Branch("akt2hi_etajes_jet_eta", &akt2hi_etajes_jet_eta_p);
  outTree_p->Branch("akt2hi_etajes_jet_phi", &akt2hi_etajes_jet_phi_p);
  outTree_p->Branch("akt2hi_etajes_jet_e", &akt2hi_etajes_jet_e_p);
  outTree_p->Branch("akt2hi_etajes_jet_m", &akt2hi_etajes_jet_m_p);
  outTree_p->Branch("akt2hi_insitu_jet_pt", &akt2hi_insitu_jet_pt_p);
  outTree_p->Branch("akt2hi_insitu_jet_eta", &akt2hi_insitu_jet_eta_p);
  outTree_p->Branch("akt2hi_insitu_jet_phi", &akt2hi_insitu_jet_phi_p);
  outTree_p->Branch("akt2hi_insitu_jet_e", &akt2hi_insitu_jet_e_p);
  outTree_p->Branch("akt2hi_insitu_jet_m", &akt2hi_insitu_jet_m_p);

  outTree_p->Branch("akt2hi_jet_clean", &akt2hi_jet_clean_p);
  if(isMC) outTree_p->Branch("akt2hi_truthpos", &akt2hi_truthpos_p);
  
  outTree_p->Branch("akt4hi_jet_n", &akt4hi_jet_n_, "akt4hi_jet_n/I");

  if(isMC){
    for(Int_t eI = 0; eI < nJES; ++eI){
      outTree_p->Branch(("akt4hi_etajes_jet_pt_sys_JES_" + std::to_string(eI)).c_str(), &(akt4hi_etajes_jet_pt_sys_JES_p[eI]));
    }
    
    for(Int_t eI = 0; eI < nJER; ++eI){
      outTree_p->Branch(("akt4hi_etajes_jet_pt_sys_JER_" + std::to_string(eI)).c_str(), &(akt4hi_etajes_jet_pt_sys_JER_p[eI]));
    }
  }

  outTree_p->Branch("akt4hi_etajes_jet_pt", &akt4hi_etajes_jet_pt_p);
  outTree_p->Branch("akt4hi_etajes_jet_eta", &akt4hi_etajes_jet_eta_p);
  outTree_p->Branch("akt4hi_etajes_jet_phi", &akt4hi_etajes_jet_phi_p);
  outTree_p->Branch("akt4hi_etajes_jet_e", &akt4hi_etajes_jet_e_p);
  outTree_p->Branch("akt4hi_etajes_jet_m", &akt4hi_etajes_jet_m_p);
  outTree_p->Branch("akt4hi_insitu_jet_pt", &akt4hi_insitu_jet_pt_p);
  outTree_p->Branch("akt4hi_insitu_jet_eta", &akt4hi_insitu_jet_eta_p);
  outTree_p->Branch("akt4hi_insitu_jet_phi", &akt4hi_insitu_jet_phi_p);
  outTree_p->Branch("akt4hi_insitu_jet_e", &akt4hi_insitu_jet_e_p);
  outTree_p->Branch("akt4hi_insitu_jet_m", &akt4hi_insitu_jet_m_p);

  outTree_p->Branch("akt4hi_jet_clean", &akt4hi_jet_clean_p);
  if(isMC) outTree_p->Branch("akt4hi_truthpos", &akt4hi_truthpos_p);

  if(!isPP){   
    outTree_p->Branch("akt10hi_jet_n", &akt10hi_jet_n_, "akt10hi_jet_n/I");
    outTree_p->Branch("akt10hi_etajes_jet_pt", &akt10hi_etajes_jet_pt_p);
    outTree_p->Branch("akt10hi_etajes_jet_eta", &akt10hi_etajes_jet_eta_p);
    outTree_p->Branch("akt10hi_etajes_jet_phi", &akt10hi_etajes_jet_phi_p);
    outTree_p->Branch("akt10hi_etajes_jet_e", &akt10hi_etajes_jet_e_p);
    outTree_p->Branch("akt10hi_etajes_jet_m", &akt10hi_etajes_jet_m_p);
    outTree_p->Branch("akt10hi_insitu_jet_pt", &akt10hi_insitu_jet_pt_p);
    outTree_p->Branch("akt10hi_insitu_jet_eta", &akt10hi_insitu_jet_eta_p);
    outTree_p->Branch("akt10hi_insitu_jet_phi", &akt10hi_insitu_jet_phi_p);
    outTree_p->Branch("akt10hi_insitu_jet_e", &akt10hi_insitu_jet_e_p);
    outTree_p->Branch("akt10hi_insitu_jet_m", &akt10hi_insitu_jet_m_p);

    outTree_p->Branch("akt10hi_jet_clean", &akt10hi_jet_clean_p);
    if(isMC) outTree_p->Branch("akt10hi_truthpos", &akt10hi_truthpos_p);    
  }

  
  outTree_p->Branch("akt2to10hi_jet_n", &akt2to10hi_jet_n_, "akt2to10hi_jet_n/I");
  outTree_p->Branch("akt2to10hi_etajes_jet_pt", &akt2to10hi_etajes_jet_pt_p);
  outTree_p->Branch("akt2to10hi_etajes_jet_eta", &akt2to10hi_etajes_jet_eta_p);
  outTree_p->Branch("akt2to10hi_etajes_jet_phi", &akt2to10hi_etajes_jet_phi_p);
  outTree_p->Branch("akt2to10hi_etajes_jet_e", &akt2to10hi_etajes_jet_e_p);
  outTree_p->Branch("akt2to10hi_etajes_jet_m", &akt2to10hi_etajes_jet_m_p);
  if(isMC) outTree_p->Branch("akt2to10hi_truthpos", &akt2to10hi_truthpos_p);
  

  outTree_p->Branch("photon_n", &photon_n_, "photon_n/I");
  outTree_p->Branch("photon_pt", &photon_pt_p);
  outTree_p->Branch("photon_pt_precali", &photon_pt_precali_p);
  if(isMC){
    outTree_p->Branch("photon_pt_sys1", &photon_pt_sys1_p);
    outTree_p->Branch("photon_pt_sys2", &photon_pt_sys2_p);
    outTree_p->Branch("photon_pt_sys3", &photon_pt_sys3_p);
    outTree_p->Branch("photon_pt_sys4", &photon_pt_sys4_p);
  }
  outTree_p->Branch("photon_eta", &photon_eta_p);
  outTree_p->Branch("photon_phi", &photon_phi_p);
  outTree_p->Branch("photon_tight", &photon_tight_p);
  outTree_p->Branch("photon_tight_b4FudgeTool", &photon_tight_b4FudgeTool_p);
  outTree_p->Branch("photon_loose", &photon_loose_p);
  outTree_p->Branch("photon_isem", &photon_isem_p);
  outTree_p->Branch("photon_convFlag", &photon_convFlag_p);
  outTree_p->Branch("photon_Rconv", &photon_Rconv_p);
  outTree_p->Branch("photon_etcone20", &photon_etcone20_p);
  outTree_p->Branch("photon_etcone30", &photon_etcone30_p);
  outTree_p->Branch("photon_etcone40", &photon_etcone40_p);
  outTree_p->Branch("photon_topoetcone20", &photon_topoetcone20_p);
  outTree_p->Branch("photon_topoetcone30", &photon_topoetcone30_p);
  outTree_p->Branch("photon_topoetcone40", &photon_topoetcone40_p);
  outTree_p->Branch("photon_Rhad1", &photon_Rhad1_p);
  outTree_p->Branch("photon_Rhad", &photon_Rhad_p);
  outTree_p->Branch("photon_e277", &photon_e277_p);
  outTree_p->Branch("photon_Reta", &photon_Reta_p);
  outTree_p->Branch("photon_Rphi", &photon_Rphi_p);
  outTree_p->Branch("photon_weta1", &photon_weta1_p);
  outTree_p->Branch("photon_weta2", &photon_weta2_p);
  outTree_p->Branch("photon_wtots1", &photon_wtots1_p);
  outTree_p->Branch("photon_f1", &photon_f1_p);
  outTree_p->Branch("photon_f3", &photon_f3_p);
  outTree_p->Branch("photon_fracs1", &photon_fracs1_p);
  outTree_p->Branch("photon_DeltaE", &photon_DeltaE_p);
  outTree_p->Branch("photon_Eratio", &photon_Eratio_p);
  
  if(isMC){
    outTree_p->Branch("akt2_truth_jet_n", &akt2_truth_jet_n_, "akt2_truth_jet_n/I");
    outTree_p->Branch("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
    outTree_p->Branch("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
    outTree_p->Branch("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
    outTree_p->Branch("akt2_truth_jet_e", &akt2_truth_jet_e_p);
    outTree_p->Branch("akt2_truth_jet_m", &akt2_truth_jet_m_p);
    outTree_p->Branch("akt2_truth_jet_partonid", &akt2_truth_jet_partonid_p);
    outTree_p->Branch("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);
    
    outTree_p->Branch("akt4_truth_jet_n", &akt4_truth_jet_n_, "akt4_truth_jet_n/I");
    outTree_p->Branch("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
    outTree_p->Branch("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
    outTree_p->Branch("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
    outTree_p->Branch("akt4_truth_jet_e", &akt4_truth_jet_e_p);
    outTree_p->Branch("akt4_truth_jet_m", &akt4_truth_jet_m_p);
    outTree_p->Branch("akt4_truth_jet_partonid", &akt4_truth_jet_partonid_p);
    outTree_p->Branch("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);

    if(!isPP){    
      outTree_p->Branch("akt10_truth_jet_n", &akt10_truth_jet_n_, "akt10_truth_jet_n/I");
      outTree_p->Branch("akt10_truth_jet_pt", &akt10_truth_jet_pt_p);
      outTree_p->Branch("akt10_truth_jet_eta", &akt10_truth_jet_eta_p);
      outTree_p->Branch("akt10_truth_jet_phi", &akt10_truth_jet_phi_p);
      outTree_p->Branch("akt10_truth_jet_e", &akt10_truth_jet_e_p);
      outTree_p->Branch("akt10_truth_jet_m", &akt10_truth_jet_m_p);
      outTree_p->Branch("akt10_truth_jet_partonid", &akt10_truth_jet_partonid_p);
      outTree_p->Branch("akt10_truth_jet_recopos", &akt10_truth_jet_recopos_p);
    }
   
    outTree_p->Branch("akt2to10_truth_jet_n", &akt2to10_truth_jet_n_, "akt2to10_truth_jet_n/I");
    outTree_p->Branch("akt2to10_truth_jet_pt", &akt2to10_truth_jet_pt_p);
    outTree_p->Branch("akt2to10_truth_jet_eta", &akt2to10_truth_jet_eta_p);
    outTree_p->Branch("akt2to10_truth_jet_phi", &akt2to10_truth_jet_phi_p);
    outTree_p->Branch("akt2to10_truth_jet_e", &akt2to10_truth_jet_e_p);
    outTree_p->Branch("akt2to10_truth_jet_m", &akt2to10_truth_jet_m_p);
    outTree_p->Branch("akt2to10_truth_jet_partonid", &akt2to10_truth_jet_partonid_p);
    outTree_p->Branch("akt2to10_truth_jet_recopos", &akt2to10_truth_jet_recopos_p);

  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<int, unsigned long long> tagToCounts;
  std::map<int, double> tagToXSec;
  std::map<int, double> tagToFilterEff;
  std::vector<int> minPthats;
  std::map<int, int> minPthatToTag;
  std::map<int, double> tagToWeight;
  
  std::map<std::string, std::vector<std::string> > configMap;
  configMap["MCPREPROCDIRNAME"] = {inDirStr};

  
  //Basic pre-processing for output config
  std::vector<std::string> listOfBranchesIn, listOfBranchesHLT, listOfBranchesHLTPre;
  std::vector<std::string> listOfBranchesOut = getVectBranchList(outTree_p);
  ULong64_t totalNEntries = 0;
  for(auto const & file : fileList){
    inFile_p = new TFile(file.c_str(), "READ");
    inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    totalNEntries += inTree_p->GetEntries();
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");

    if(!isPP){
      inTree_p->SetBranchStatus("*", 0);
      inTree_p->SetBranchStatus("fcalA_et", 1);
      inTree_p->SetBranchStatus("fcalC_et", 1);
      
      inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
      inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);
      
      for(Long64_t entry = 0; entry < inTree_p->GetEntries(); ++entry){
	inTree_p->GetEntry(entry);

	cent_ = centTable.GetCent(fcalA_et_ + fcalC_et_);
	++(centCounts[cent_]);
      }
    }
    
    std::map<std::string, std::string> tempConfigMap = GetMapFromEnv(inConfig_p);

    if(isMC){
      sampleHandler sHandler;
      
      std::string inDataSetName = inConfig_p->GetValue("INDATASET", "");
      
      if(inDataSetName.size() == 0 || !sHandler.Init(inDataSetName)){
	std::cout << "GDJMCNTUPLEPREPROC ERROR - Given input \'" << file << "\' contains INDATASET \'" << inDataSetName << "\' that is not valid. return 1" << std::endl;
	
	inFile_p->Close();
	delete inFile_p;
	
	outFile_p->Close();
	delete outFile_p;
	
	return 1;
      }
      
      sampleTag_ = sHandler.GetTag();
      
      if(tagToCounts.count(sampleTag_) == 0){
	tagToCounts[sampleTag_] = inTree_p->GetEntries();
	tagToXSec[sampleTag_] = sHandler.GetXSection();
	tagToFilterEff[sampleTag_] = sHandler.GetFilterEff();
	int minPthat = sHandler.GetMinPthat();
	minPthats.push_back(minPthat);
	minPthatToTag[sHandler.GetMinPthat()] = sampleTag_;
      }
      else tagToCounts[sampleTag_] += inTree_p->GetEntries();
    }
    
    for(auto const & val : tempConfigMap){
      bool isFound = false;
      std::vector<std::string> tempVect = configMap[val.first];
      for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
	if(isStrSame(tempVect[vI], val.second)){
	  isFound = true;
	  break;
	}
      }

      if(!isFound){
	tempVect.push_back(val.second);
	configMap[val.first] = tempVect;
      }      
    }

    std::vector<std::string> tempBranches = getVectBranchList(inTree_p);
    for(auto const & branch : tempBranches){
      if(branch.size() >= 4){
	if(branch.substr(0,4).find("HLT_") != std::string::npos){
	  if(branch.find("presc") != std::string::npos){
	    if(!vectContainsStr(branch, &listOfBranchesHLTPre)) listOfBranchesHLTPre.push_back(branch);
	  }
	  else{
	    if(!vectContainsStr(branch, &listOfBranchesHLT)) listOfBranchesHLT.push_back(branch);
	  }	    
	  
	  continue;
	}
      }
	
      if(!vectContainsStr(branch, &listOfBranchesIn)){
	listOfBranchesIn.push_back(branch);
      }
    }
        
    inFile_p->Close();
    delete inFile_p;
  }

  if(isMC){
    std::sort(std::begin(minPthats), std::end(minPthats));
    std::cout << "Min pthats: " << std::endl;
    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      int tag = minPthatToTag[minPthats[pI]];
      double weight = tagToXSec[tag]*tagToFilterEff[tag]/(double)tagToCounts[tag];
      tagToWeight[tag] = weight;
    }
    
    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      if(pI == 0) continue;
      
      int tag = minPthatToTag[minPthats[pI]];
      int tag0 = minPthatToTag[minPthats[0]];
      
      tagToWeight[tag] /= tagToWeight[tag0];
    }
    
    int tag0 = minPthatToTag[minPthats[0]];
    tagToWeight[tag0] = 1.0;
  
    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      int tag = minPthatToTag[minPthats[pI]];

      std::cout << " " << pI << "/" << minPthats.size() << ": " << minPthats[pI] << ", " << tagToWeight[tag] << std::endl;
    }
  }

  std::vector<Bool_t*> hltVect;
  std::vector<Float_t*> hltPreVect;
  hltVect.reserve(listOfBranchesHLT.size());
  hltPreVect.reserve(listOfBranchesHLTPre.size());

  for(auto const & branchHLT : listOfBranchesHLT){
    hltVect.push_back(new Bool_t(false));

    outTree_p->Branch(branchHLT.c_str(), hltVect[hltVect.size()-1], (branchHLT + "/O").c_str());
  }
  for(auto const & branchHLTPre : listOfBranchesHLTPre){
    hltPreVect.push_back(new Float_t(0.0));

    outTree_p->Branch(branchHLTPre.c_str(), hltPreVect[hltPreVect.size()-1], (branchHLTPre + "/F").c_str());
  }
  
  std::cout << "BRANCHES: " << std::endl;
  for(auto const &branch : listOfBranchesIn){
    std::cout << " BRANCH '" << branch << "'" << std::endl;
  }

  bool allBranchesGood = true;
  for(auto const & branchIn : listOfBranchesIn){
    bool containsBranch = vectContainsStr(branchIn, &listOfBranchesOut);
    if(!containsBranch){
      if(!keepTruth){
	if(isStrSame(branchIn, "truth_n")) continue;
	if(isStrSame(branchIn, "truth_charge")) continue;
	if(isStrSame(branchIn, "truth_pt")) continue;
	if(isStrSame(branchIn, "truth_eta")) continue;
	if(isStrSame(branchIn, "truth_phi")) continue;
	if(isStrSame(branchIn, "truth_e")) continue;
	if(isStrSame(branchIn, "truth_pdg")) continue;
	if(isStrSame(branchIn, "truth_type")) continue;
	if(isStrSame(branchIn, "truth_origin")) continue;
	if(isStrSame(branchIn, "truth_status")) continue;
      }

      std::cout << "GDJMCNTUPLEPREPROC ERROR - \'" << branchIn << "\' branch is not part of output tree. please fix macro, return 1" << std::endl;
    }

    allBranchesGood = allBranchesGood && containsBranch;    
  }

  bool allBranchesGood2 = true;
  for(auto const & branchOut : listOfBranchesOut){
    if(vectContainsStr(branchOut, &outBranchesToAdd)) continue;//Dont want to throw error on branches we added
    bool containsBranch = vectContainsStr(branchOut, &listOfBranchesIn);

    if(!containsBranch){
      if(!keepTruth){
	if(isStrSame(branchOut, "truth_n")) continue;
	if(isStrSame(branchOut, "truth_charge")) continue;
	if(isStrSame(branchOut, "truth_pt")) continue;
	if(isStrSame(branchOut, "truth_eta")) continue;
	if(isStrSame(branchOut, "truth_phi")) continue;
	if(isStrSame(branchOut, "truth_e")) continue;
	if(isStrSame(branchOut, "truth_pdg")) continue;
	if(isStrSame(branchOut, "truth_type")) continue;
	if(isStrSame(branchOut, "truth_origin")) continue;
	if(isStrSame(branchOut, "truth_status")) continue;
      }

      std::cout << "GDJMCNTUPLEPREPROC ERROR - \'" << branchOut << "\' branch is not part of input tree. please fix macro, return 1" << std::endl;
    }

    allBranchesGood2 = allBranchesGood2 && containsBranch;
  }

  
  if(!allBranchesGood || !allBranchesGood2){
    outFile_p->Close();
    delete outFile_p;
    return 1;
  }
 

  if(!isPP && isMC){
    for(unsigned int cI = 0; cI < ncollWeights.size(); ++cI){      
      if(centCounts[cI] < 0.5) continue;
      ncollWeights[cI] /= centCounts[cI];
    }
  }
  
  sampleHandler sHandler;

  const Int_t nDivInt = inConfig_p->GetValue("NDIVPERMSG", 20);  
  ULong64_t nDiv = TMath::Max((ULong64_t)1, totalNEntries/nDivInt);
  ULong64_t currTotalEntries = 0;
  UInt_t nFile = 0;

  const Int_t nTerm = inConfig_p->GetValue("NTERM", -1);
  cppWatch timer;
  cppWatch subTimer1;
  cppWatch subTimer2;
  cppWatch subTimer3;
  cppWatch subTimer4;
  
  double prevTimeCPU = 0.0;
  double prevTimeWall = 0.0;

  int minNVert = 99999999;
  int maxNVert = 0;
  
  int minNTruth = 99999999;
  int maxNTruth = 0;
  
  int minNPhoton = 99999999;
  int maxNPhoton = 0;
  
  int minNRecoR2 = 99999999;
  int maxNRecoR2 = 0;
  
  int minNRecoR4 = 99999999;
  int maxNRecoR4 = 0;
  
  int minNRecoR10 = 99999999;
  int maxNRecoR10 = 0;
  
  int minNRecoR2to10 = 99999999;
  int maxNRecoR2to10 = 0;
  
  int minNTruthR2 = 99999999;
  int maxNTruthR2 = 0;
  
  int minNTruthR4 = 99999999;
  int maxNTruthR4 = 0;
  
  int minNTruthR10 = 99999999;
  int maxNTruthR10 = 0;
  
  int minNTruthR2to10 = 99999999;
  int maxNTruthR2to10 = 0;

  for(auto const & file : fileList){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TFile* inFile_p = new TFile(file.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
    std::string inDataSetName = inConfig_p->GetValue("INDATASET", "");

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(isMC){
      if(inDataSetName.size() == 0 || !sHandler.Init(inDataSetName)){
	std::cout << "GDJMCNTUPLEPREPROC ERROR - Given input \'" << file << "\' contains INDATASET \'" << inDataSetName << "\' that is not valid. return 1" << std::endl;
	
	inFile_p->Close();
	delete inFile_p;
	
	outFile_p->Close();
	delete outFile_p;
	
	return 1;
      }
    
      sampleTag_ = sHandler.GetTag();
      sampleWeight_ = tagToWeight[sampleTag_];
    }

    for(unsigned int bI = 0; bI < listOfBranchesHLT.size(); ++bI){
      inTree_p->SetBranchAddress(listOfBranchesHLT[bI].c_str(), hltVect[bI]);
    }

    for(unsigned int bI = 0; bI < listOfBranchesHLTPre.size(); ++bI){
      inTree_p->SetBranchAddress(listOfBranchesHLTPre[bI].c_str(), hltPreVect[bI]);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    
    inTree_p->SetBranchAddress("runNumber", &runNumber_);
    inTree_p->SetBranchAddress("eventNumber", &eventNumber_);
    inTree_p->SetBranchAddress("lumiBlock", &lumiBlock_);
    inTree_p->SetBranchAddress("passesToroid", &passesToroid_);

    if(isMC){
      inTree_p->SetBranchAddress("pthat", &pthat_);

      inTree_p->SetBranchAddress("treePartonPt", treePartonPt_);
      inTree_p->SetBranchAddress("treePartonEta", treePartonEta_);
      inTree_p->SetBranchAddress("treePartonPhi", treePartonPhi_);
      inTree_p->SetBranchAddress("treePartonId", treePartonId_);
    }

    inTree_p->SetBranchAddress("is_pileup", &is_pileup_);
    inTree_p->SetBranchAddress("is_oo_pileup", &is_oo_pileup_);
    
    inTree_p->SetBranchAddress("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_);
    inTree_p->SetBranchAddress("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_);

    inTree_p->SetBranchAddress("nvert", &nvert_);
    inTree_p->SetBranchAddress("vert_x", &vert_x_p);
    inTree_p->SetBranchAddress("vert_y", &vert_y_p);
    inTree_p->SetBranchAddress("vert_z", &vert_z_p);
    inTree_p->SetBranchAddress("vert_type", &vert_type_p);
    inTree_p->SetBranchAddress("vert_ntrk", &vert_ntrk_p);

    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);
    inTree_p->SetBranchAddress("fcalA_et_Cos2", &fcalA_et_Cos2_);
    inTree_p->SetBranchAddress("fcalC_et_Cos2", &fcalC_et_Cos2_);
    inTree_p->SetBranchAddress("fcalA_et_Sin2", &fcalA_et_Sin2_);
    inTree_p->SetBranchAddress("fcalC_et_Sin2", &fcalC_et_Sin2_);
    
    inTree_p->SetBranchAddress("fcalA_et_Cos3", &fcalA_et_Cos3_);
    inTree_p->SetBranchAddress("fcalC_et_Cos3", &fcalC_et_Cos3_);
    inTree_p->SetBranchAddress("fcalA_et_Sin3", &fcalA_et_Sin3_);
    inTree_p->SetBranchAddress("fcalC_et_Sin3", &fcalC_et_Sin3_);
    
    inTree_p->SetBranchAddress("fcalA_et_Cos4", &fcalA_et_Cos4_);
    inTree_p->SetBranchAddress("fcalC_et_Cos4", &fcalC_et_Cos4_);
    inTree_p->SetBranchAddress("fcalA_et_Sin4", &fcalA_et_Sin4_);
    inTree_p->SetBranchAddress("fcalC_et_Sin4", &fcalC_et_Sin4_);
    
    inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi_);
    inTree_p->SetBranchAddress("evtPlane3Phi", &evtPlane3Phi_);
    inTree_p->SetBranchAddress("evtPlane4Phi", &evtPlane4Phi_);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(getTracks){
      inTree_p->SetBranchAddress("ntrk", &ntrk_);
      inTree_p->SetBranchAddress("trk_pt", &trk_pt_p);
      inTree_p->SetBranchAddress("trk_eta", &trk_eta_p);
      inTree_p->SetBranchAddress("trk_phi", &trk_phi_p);
      
      inTree_p->SetBranchAddress("trk_charge", &trk_charge_p);
      inTree_p->SetBranchAddress("trk_tightPrimary", &trk_tightPrimary_p);
      inTree_p->SetBranchAddress("trk_HILoose", &trk_HILoose_p);
      inTree_p->SetBranchAddress("trk_HITight", &trk_HITight_p);
      inTree_p->SetBranchAddress("trk_d0", &trk_d0_p);
      inTree_p->SetBranchAddress("trk_z0", &trk_z0_p);
      inTree_p->SetBranchAddress("trk_vz", &trk_vz_p);
      inTree_p->SetBranchAddress("trk_theta", &trk_theta_p);
      inTree_p->SetBranchAddress("trk_nPixelHits", &trk_nPixelHits_p);
      inTree_p->SetBranchAddress("trk_nSCTHits", &trk_nSCTHits_p);
      inTree_p->SetBranchAddress("trk_nBlayerHits", &trk_nBlayerHits_p);
    }
    
    if(isMC){
      inTree_p->SetBranchAddress("truth_n", &truth_n_);
      inTree_p->SetBranchAddress("truth_charge", &truth_charge_p);
      inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
      inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
      inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
      inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);
      inTree_p->SetBranchAddress("truth_e", &truth_e_p);
      inTree_p->SetBranchAddress("truth_type", &truth_type_p);
      inTree_p->SetBranchAddress("truth_origin", &truth_origin_p);
      inTree_p->SetBranchAddress("truth_status", &truth_status_p);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    inTree_p->SetBranchAddress("akt2hi_jet_n", &akt2hi_jet_n_);
       
    inTree_p->SetBranchAddress("akt2hi_etajes_jet_pt", &akt2hi_etajes_jet_pt_p);
    if(isMC){
      for(Int_t eI = 0; eI < nJES; ++eI){
	inTree_p->SetBranchAddress(("akt2hi_etajes_jet_pt_sys_JES_" + std::to_string(eI)).c_str(), &(akt2hi_etajes_jet_pt_sys_JES_p[eI]));
      }

      for(Int_t eI = 0; eI < nJER; ++eI){
	inTree_p->SetBranchAddress(("akt2hi_etajes_jet_pt_sys_JER_" + std::to_string(eI)).c_str(), &(akt2hi_etajes_jet_pt_sys_JER_p[eI]));
      }
    }

    inTree_p->SetBranchAddress("akt2hi_etajes_jet_eta", &akt2hi_etajes_jet_eta_p);
    inTree_p->SetBranchAddress("akt2hi_etajes_jet_phi", &akt2hi_etajes_jet_phi_p);
    inTree_p->SetBranchAddress("akt2hi_etajes_jet_e", &akt2hi_etajes_jet_e_p);
    inTree_p->SetBranchAddress("akt2hi_etajes_jet_m", &akt2hi_etajes_jet_m_p);

    inTree_p->SetBranchAddress("akt2hi_insitu_jet_pt", &akt2hi_insitu_jet_pt_p);
    inTree_p->SetBranchAddress("akt2hi_insitu_jet_eta", &akt2hi_insitu_jet_eta_p);
    inTree_p->SetBranchAddress("akt2hi_insitu_jet_phi", &akt2hi_insitu_jet_phi_p);
    inTree_p->SetBranchAddress("akt2hi_insitu_jet_e", &akt2hi_insitu_jet_e_p);
    inTree_p->SetBranchAddress("akt2hi_insitu_jet_m", &akt2hi_insitu_jet_m_p);

    inTree_p->SetBranchAddress("akt2hi_jet_clean", &akt2hi_jet_clean_p);
    if(isMC) inTree_p->SetBranchAddress("akt2hi_truthpos", &akt2hi_truthpos_p);

    inTree_p->SetBranchAddress("akt4hi_jet_n", &akt4hi_jet_n_);
    inTree_p->SetBranchAddress("akt4hi_etajes_jet_pt", &akt4hi_etajes_jet_pt_p);
    if(isMC){
      for(Int_t eI = 0; eI < nJES; ++eI){
	inTree_p->SetBranchAddress(("akt4hi_etajes_jet_pt_sys_JES_" + std::to_string(eI)).c_str(), &(akt4hi_etajes_jet_pt_sys_JES_p[eI]));
      }

      for(Int_t eI = 0; eI < nJER; ++eI){
	inTree_p->SetBranchAddress(("akt4hi_etajes_jet_pt_sys_JER_" + std::to_string(eI)).c_str(), &(akt4hi_etajes_jet_pt_sys_JER_p[eI]));
      }
    }
    inTree_p->SetBranchAddress("akt4hi_etajes_jet_eta", &akt4hi_etajes_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_etajes_jet_phi", &akt4hi_etajes_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_etajes_jet_e", &akt4hi_etajes_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_etajes_jet_m", &akt4hi_etajes_jet_m_p);

    inTree_p->SetBranchAddress("akt4hi_insitu_jet_pt", &akt4hi_insitu_jet_pt_p);
    inTree_p->SetBranchAddress("akt4hi_insitu_jet_eta", &akt4hi_insitu_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_insitu_jet_phi", &akt4hi_insitu_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_insitu_jet_e", &akt4hi_insitu_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_insitu_jet_m", &akt4hi_insitu_jet_m_p);

    inTree_p->SetBranchAddress("akt4hi_jet_clean", &akt4hi_jet_clean_p);
    if(isMC) inTree_p->SetBranchAddress("akt4hi_truthpos", &akt4hi_truthpos_p);

    if(!isPP){      
      inTree_p->SetBranchAddress("akt10hi_jet_n", &akt10hi_jet_n_);
      inTree_p->SetBranchAddress("akt10hi_etajes_jet_pt", &akt10hi_etajes_jet_pt_p);
      inTree_p->SetBranchAddress("akt10hi_etajes_jet_eta", &akt10hi_etajes_jet_eta_p);
      inTree_p->SetBranchAddress("akt10hi_etajes_jet_phi", &akt10hi_etajes_jet_phi_p);
      inTree_p->SetBranchAddress("akt10hi_etajes_jet_e", &akt10hi_etajes_jet_e_p);
      inTree_p->SetBranchAddress("akt10hi_etajes_jet_m", &akt10hi_etajes_jet_m_p);

      inTree_p->SetBranchAddress("akt10hi_insitu_jet_pt", &akt10hi_insitu_jet_pt_p);
      inTree_p->SetBranchAddress("akt10hi_insitu_jet_eta", &akt10hi_insitu_jet_eta_p);
      inTree_p->SetBranchAddress("akt10hi_insitu_jet_phi", &akt10hi_insitu_jet_phi_p);
      inTree_p->SetBranchAddress("akt10hi_insitu_jet_e", &akt10hi_insitu_jet_e_p);
      inTree_p->SetBranchAddress("akt10hi_insitu_jet_m", &akt10hi_insitu_jet_m_p);

      inTree_p->SetBranchAddress("akt10hi_jet_clean", &akt10hi_jet_clean_p);
      if(isMC) inTree_p->SetBranchAddress("akt10hi_truthpos", &akt10hi_truthpos_p);      
    }
    
    inTree_p->SetBranchAddress("akt2to10hi_jet_n", &akt2to10hi_jet_n_);
    inTree_p->SetBranchAddress("akt2to10hi_etajes_jet_pt", &akt2to10hi_etajes_jet_pt_p);
    inTree_p->SetBranchAddress("akt2to10hi_etajes_jet_eta", &akt2to10hi_etajes_jet_eta_p);
    inTree_p->SetBranchAddress("akt2to10hi_etajes_jet_phi", &akt2to10hi_etajes_jet_phi_p);
    inTree_p->SetBranchAddress("akt2to10hi_etajes_jet_e", &akt2to10hi_etajes_jet_e_p);
    inTree_p->SetBranchAddress("akt2to10hi_etajes_jet_m", &akt2to10hi_etajes_jet_m_p);
    if(isMC) inTree_p->SetBranchAddress("akt2to10hi_truthpos", &akt2to10hi_truthpos_p);    

    inTree_p->SetBranchAddress("photon_n", &photon_n_);
    inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
    inTree_p->SetBranchAddress("photon_pt_precali", &photon_pt_precali_p);
    if(isMC){
      inTree_p->SetBranchAddress("photon_pt_sys1", &photon_pt_sys1_p);
      inTree_p->SetBranchAddress("photon_pt_sys2", &photon_pt_sys2_p);
      inTree_p->SetBranchAddress("photon_pt_sys3", &photon_pt_sys3_p);
      inTree_p->SetBranchAddress("photon_pt_sys4", &photon_pt_sys4_p);
    }
    inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
    inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
    inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
    inTree_p->SetBranchAddress("photon_tight_b4FudgeTool", &photon_tight_b4FudgeTool_p);
    inTree_p->SetBranchAddress("photon_loose", &photon_loose_p);
    inTree_p->SetBranchAddress("photon_isem", &photon_isem_p);
    inTree_p->SetBranchAddress("photon_convFlag", &photon_convFlag_p);
    inTree_p->SetBranchAddress("photon_Rconv", &photon_Rconv_p);
    inTree_p->SetBranchAddress("photon_etcone20", &photon_etcone20_p);
    inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);
    inTree_p->SetBranchAddress("photon_etcone40", &photon_etcone40_p);
    inTree_p->SetBranchAddress("photon_topoetcone20", &photon_topoetcone20_p);
    inTree_p->SetBranchAddress("photon_topoetcone30", &photon_topoetcone30_p);
    inTree_p->SetBranchAddress("photon_topoetcone40", &photon_topoetcone40_p);
    inTree_p->SetBranchAddress("photon_Rhad1", &photon_Rhad1_p);
    inTree_p->SetBranchAddress("photon_Rhad", &photon_Rhad_p);
    inTree_p->SetBranchAddress("photon_e277", &photon_e277_p);
    inTree_p->SetBranchAddress("photon_Reta", &photon_Reta_p);
    inTree_p->SetBranchAddress("photon_Rphi", &photon_Rphi_p);
    inTree_p->SetBranchAddress("photon_weta1", &photon_weta1_p);
    inTree_p->SetBranchAddress("photon_weta2", &photon_weta2_p);
    inTree_p->SetBranchAddress("photon_wtots1", &photon_wtots1_p);
    inTree_p->SetBranchAddress("photon_f1", &photon_f1_p);
    inTree_p->SetBranchAddress("photon_f3", &photon_f3_p);
    inTree_p->SetBranchAddress("photon_fracs1", &photon_fracs1_p);
    inTree_p->SetBranchAddress("photon_DeltaE", &photon_DeltaE_p);
    inTree_p->SetBranchAddress("photon_Eratio", &photon_Eratio_p);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(isMC){
      inTree_p->SetBranchAddress("akt2_truth_jet_n", &akt2_truth_jet_n_);
      inTree_p->SetBranchAddress("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_e", &akt2_truth_jet_e_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_m", &akt2_truth_jet_m_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_partonid", &akt2_truth_jet_partonid_p);
      inTree_p->SetBranchAddress("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);

      inTree_p->SetBranchAddress("akt4_truth_jet_n", &akt4_truth_jet_n_);
      inTree_p->SetBranchAddress("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_e", &akt4_truth_jet_e_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_m", &akt4_truth_jet_m_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_partonid", &akt4_truth_jet_partonid_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);


      if(!isPP){	
	inTree_p->SetBranchAddress("akt10_truth_jet_n", &akt10_truth_jet_n_);
	inTree_p->SetBranchAddress("akt10_truth_jet_pt", &akt10_truth_jet_pt_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_eta", &akt10_truth_jet_eta_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_phi", &akt10_truth_jet_phi_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_e", &akt10_truth_jet_e_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_m", &akt10_truth_jet_m_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_partonid", &akt10_truth_jet_partonid_p);
	inTree_p->SetBranchAddress("akt10_truth_jet_recopos", &akt10_truth_jet_recopos_p);       
      }
      
      inTree_p->SetBranchAddress("akt2to10_truth_jet_n", &akt2to10_truth_jet_n_);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_pt", &akt2to10_truth_jet_pt_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_eta", &akt2to10_truth_jet_eta_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_phi", &akt2to10_truth_jet_phi_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_e", &akt2to10_truth_jet_e_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_m", &akt2to10_truth_jet_m_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_partonid", &akt2to10_truth_jet_partonid_p);
      inTree_p->SetBranchAddress("akt2to10_truth_jet_recopos", &akt2to10_truth_jet_recopos_p);            
    }
  
    const ULong64_t nEntries = inTree_p->GetEntries();
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
      subTimer1.start();
      if(currTotalEntries%nDiv == 0){
	std::cout << " Entry " << currTotalEntries << "/" << totalNEntries << "... (File " << nFile << "/" << fileList.size() << ", disp on every " << nDiv << " events)"  << std::endl;

	if(currTotalEntries != 0){
	  timer.stop();

	  subTimer1.stop();
	  double subTime1 = subTimer1.totalCPU();
	  subTimer1.start();
	  double subTime2 = subTimer2.totalCPU();
	  double subTime3 = subTimer3.totalCPU();
	  double subTime4 = subTimer4.totalCPU();

	  double currTotalCPU = timer.totalCPU();
	  double currTotalWall = timer.totalWall();
	  double deltaTotalCPU = currTotalCPU - prevTimeCPU;
	  double deltaTotalWall = currTotalWall - prevTimeWall;

	  std::cout << "Current cpu time: " << prettyString(currTotalCPU, 1, false) << " (Delta: " << prettyString(deltaTotalCPU, 1, false) << ")" << std::endl;
	  std::cout << "Current wall time: " << prettyString(currTotalWall, 1, false) << " (Delta: " << prettyString(deltaTotalWall, 1, false) << ")" << std::endl;
	  
	  std::cout << "TIMING FRACTIONS: " << std::endl;
	  std::cout << " 1: " << subTime1/currTotalCPU << std::endl;
	  std::cout << " 2: " << subTime2/currTotalCPU << std::endl;
	  std::cout << " 3: " << subTime3/currTotalCPU << std::endl;
	  std::cout << " 4: " << subTime4/currTotalCPU << std::endl;

	  double totalTimeGuess = timer.totalWall()*totalNEntries/currTotalEntries;
	    
	  std::cout << " TOOK " << timer.totalWall() << " TO PROCESS " << currTotalEntries << " events..." << std::endl;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " SECONDS..." << std::endl;
	  totalTimeGuess /= 60;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " MINUTES..." << std::endl;
	  totalTimeGuess /= 60;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " HOURS..." << std::endl;
	  std::cout << std::endl;

	  prevTimeCPU = currTotalCPU;
	  prevTimeWall = currTotalWall;
	}	
	timer.start();
      }

      if(nTerm > 0){
	if(currTotalEntries > (ULong64_t)nTerm){
	  std::cout << "TERMINATING" << std::endl;
	  timer.stop();
	  
	  double totalTimeGuess = timer.totalWall()*totalNEntries/currTotalEntries;
	  
	  std::cout << " TOOK " << timer.totalWall() << " TO PROCESS " << currTotalEntries << " events..." << std::endl;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " SECONDS..." << std::endl;
	  totalTimeGuess /= 60;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " MINUTES..." << std::endl;
	  totalTimeGuess /= 60;
	  std::cout << " ESTIMATE TOTAL TIME TO BE " << totalTimeGuess << " HOURS..." << std::endl;
      	  
	  break;
	}
      }

      inTree_p->GetEntry(entry);

      subTimer1.stop();
      subTimer2.start();

      if(doMinGammaPt){
	bool isGoodGamma = false;
	for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
	  if(photon_pt_p->at(pI) > minGammaPt){
	    isGoodGamma = true;
	    break;
	  }
	}
	if(!isGoodGamma){
	  ++currTotalEntries;
	  continue;
	}
      }

      if(!isPP){
	cent_ = centTable.GetCent(fcalA_et_ + fcalC_et_);
	ncollWeight_ = ncollWeights[cent_];
      }
      else ncollWeight_ = 1.0;

      fullWeight_ = sampleWeight_*ncollWeight_;

      subTimer2.stop();
      subTimer3.start();
    
      

      if(minNVert > nvert_) minNVert = nvert_;
      if(maxNVert < nvert_) maxNVert = nvert_;
      
      if(isMC){
	if(minNTruth > truth_n_) minNTruth = truth_n_;
	if(maxNTruth < truth_n_) maxNTruth = truth_n_;
      }

      if(minNPhoton > photon_n_) minNPhoton = photon_n_;
      if(maxNPhoton < photon_n_) maxNPhoton = photon_n_;

      if(minNRecoR2 > akt2hi_jet_n_) minNRecoR2 = akt2hi_jet_n_;
      if(maxNRecoR2 < akt2hi_jet_n_) maxNRecoR2 = akt2hi_jet_n_;

      if(minNRecoR4 > akt4hi_jet_n_) minNRecoR4 = akt4hi_jet_n_;
      if(maxNRecoR4 < akt4hi_jet_n_) maxNRecoR4 = akt4hi_jet_n_;

      if(!isPP){
	if(minNRecoR10 > akt10hi_jet_n_) minNRecoR10 = akt10hi_jet_n_;
	if(maxNRecoR10 < akt10hi_jet_n_) maxNRecoR10 = akt10hi_jet_n_;
      }

      if(minNRecoR2to10 > akt2to10hi_jet_n_) minNRecoR2to10 = akt2to10hi_jet_n_;
      if(maxNRecoR2to10 < akt2to10hi_jet_n_) maxNRecoR2to10 = akt2to10hi_jet_n_;

      if(isMC){
	if(minNTruthR2 > akt2_truth_jet_n_) minNTruthR2 = akt2_truth_jet_n_;
	if(maxNTruthR2 < akt2_truth_jet_n_) maxNTruthR2 = akt2_truth_jet_n_;
	
	if(minNTruthR4 > akt4_truth_jet_n_) minNTruthR4 = akt4_truth_jet_n_;
	if(maxNTruthR4 < akt4_truth_jet_n_) maxNTruthR4 = akt4_truth_jet_n_;

	if(!isPP){
	  if(minNTruthR10 > akt10_truth_jet_n_) minNTruthR10 = akt10_truth_jet_n_;
	  if(maxNTruthR10 < akt10_truth_jet_n_) maxNTruthR10 = akt10_truth_jet_n_;
	}
	
	if(minNTruthR2to10 > akt2to10_truth_jet_n_) minNTruthR2to10 = akt2to10_truth_jet_n_;
	if(maxNTruthR2to10 < akt2to10_truth_jet_n_) maxNTruthR2to10 = akt2to10_truth_jet_n_;
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      //Vertex size checks
      if(vert_x_p->size() != (unsigned int)nvert_) std::cout << "VERT VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(vert_y_p->size() != (unsigned int)nvert_) std::cout << "VERT VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(vert_z_p->size() != (unsigned int)nvert_) std::cout << "VERT VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(vert_type_p->size() != (unsigned int)nvert_) std::cout << "VERT VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(vert_ntrk_p->size() != (unsigned int)nvert_) std::cout << "VERT VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      //truth size checks
      if(isMC){	
	if(truth_charge_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_pt_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_e_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_eta_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_phi_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_status_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_pdg_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_type_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(truth_origin_p->size() != (unsigned int)truth_n_) std::cout << "TRUTH VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      //akt2 jet collection size checks
      if(akt2hi_etajes_jet_pt_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(isMC){
	for(Int_t eI = 0; eI < nJES; ++eI){
	  if(akt2hi_etajes_jet_pt_sys_JES_p[eI]->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}

	for(Int_t eI = 0; eI < nJER; ++eI){
	  if(akt2hi_etajes_jet_pt_sys_JER_p[eI]->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}
      }

      if(akt2hi_etajes_jet_eta_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_etajes_jet_phi_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_etajes_jet_e_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_etajes_jet_m_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(akt2hi_insitu_jet_pt_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_insitu_jet_eta_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_insitu_jet_phi_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_insitu_jet_e_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2hi_insitu_jet_m_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(akt2hi_jet_clean_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(isMC){
	if(akt2hi_truthpos_p->size() != (unsigned int)akt2hi_jet_n_) std::cout << "AKT2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }

      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      //akt4 jet collection size checks
      if(akt4hi_etajes_jet_pt_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(isMC){
	for(Int_t eI = 0; eI < nJES; ++eI){
	  if(akt4hi_etajes_jet_pt_sys_JES_p[eI]->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}

	for(Int_t eI = 0; eI < nJER; ++eI){
	  if(akt4hi_etajes_jet_pt_sys_JER_p[eI]->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}
      }

      if(akt4hi_etajes_jet_eta_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_etajes_jet_phi_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_etajes_jet_e_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_etajes_jet_m_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(akt4hi_insitu_jet_pt_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_insitu_jet_eta_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_insitu_jet_phi_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_insitu_jet_e_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt4hi_insitu_jet_m_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(akt4hi_jet_clean_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(isMC){
	if(akt4hi_truthpos_p->size() != (unsigned int)akt4hi_jet_n_) std::cout << "AKT4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }

      //akt10 jet collection size checks
      if(!isPP){
	if(akt10hi_etajes_jet_pt_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_etajes_jet_eta_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_etajes_jet_phi_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_etajes_jet_e_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_etajes_jet_m_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

	if(akt10hi_insitu_jet_pt_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_insitu_jet_eta_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_insitu_jet_phi_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_insitu_jet_e_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt10hi_insitu_jet_m_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

	if(akt10hi_jet_clean_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(isMC){
	  if(akt10hi_truthpos_p->size() != (unsigned int)akt10hi_jet_n_) std::cout << "AKT10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}
      }

      //akt2to10 jet collection size checks
      if(akt2to10hi_etajes_jet_pt_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2to10hi_etajes_jet_eta_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2to10hi_etajes_jet_phi_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2to10hi_etajes_jet_e_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(akt2to10hi_etajes_jet_m_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(isMC){
	if(akt2to10hi_truthpos_p->size() != (unsigned int)akt2to10hi_jet_n_) std::cout << "AKT2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }

      //Burned a few times on photon vector mismatches so adding a check here
      if(photon_pt_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_pt_precali_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(isMC){
	if(photon_pt_sys1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(photon_pt_sys2_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(photon_pt_sys3_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(photon_pt_sys4_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }
      if(photon_eta_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_phi_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_tight_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_tight_b4FudgeTool_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_loose_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_isem_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_convFlag_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Rconv_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_etcone20_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_etcone30_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_etcone40_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_topoetcone20_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_topoetcone30_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_topoetcone40_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Rhad1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Rhad_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_e277_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Reta_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Rphi_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_weta1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_weta2_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_wtots1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_f1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_f3_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_fracs1_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_DeltaE_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      if(photon_Eratio_p->size() != (unsigned int)photon_n_) std::cout << "PHOTON VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

      if(isMC){
	if(akt2_truth_jet_pt_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_eta_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_phi_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_e_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_m_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_recopos_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2_truth_jet_partonid_p->size() != (unsigned int)akt2_truth_jet_n_) std::cout << "TRUTH JET R2 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      

	if(akt4_truth_jet_pt_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_eta_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_phi_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_e_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_m_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_recopos_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt4_truth_jet_partonid_p->size() != (unsigned int)akt4_truth_jet_n_) std::cout << "TRUTH JET R4 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;

	if(!isPP){
	  if(akt10_truth_jet_pt_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_eta_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_phi_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_e_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_m_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_recopos_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	  if(akt10_truth_jet_partonid_p->size() != (unsigned int)akt10_truth_jet_n_) std::cout << "TRUTH JET R10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	}
      

	if(akt2to10_truth_jet_pt_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_eta_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_phi_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_e_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_m_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_recopos_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
	if(akt2to10_truth_jet_partonid_p->size() != (unsigned int)akt2to10_truth_jet_n_) std::cout << "TRUTH JET R2to10 VECTOR WARNING: VECTOR SIZE MISMATCH" << std::endl;
      }

      if(isMC){
	truthPhotonPt_ = -999.;     
	truthPhotonPhi_ = -999.;      
	truthPhotonEta_ = -999.;
	truthPhotonType_ = -999;
	truthPhotonOrigin_ = -999;
	truthPhotonIso2_ = -999.;
	truthPhotonIso3_ = -999.;
	truthPhotonIso4_ = -999.;
	
	if(keepTruth){
	  truthOut_charge_p->clear();
	  truthOut_pt_p->clear();
	  truthOut_eta_p->clear();
	  truthOut_phi_p->clear();
	  truthOut_e_p->clear();
	  truthOut_pdg_p->clear();
	  truthOut_type_p->clear();
	  truthOut_origin_p->clear();
	  truthOut_status_p->clear();
	  truthOut_n_ = 0;
	}

	for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	  if(truth_status_p->at(tI) != 1) continue;	

	  if(keepTruth){
	    truthOut_charge_p->push_back(truth_charge_p->at(tI));
	    truthOut_pt_p->push_back(truth_pt_p->at(tI));
	    truthOut_eta_p->push_back(truth_eta_p->at(tI));
	    truthOut_phi_p->push_back(truth_phi_p->at(tI));
	    truthOut_e_p->push_back(truth_e_p->at(tI));
	    truthOut_pdg_p->push_back(truth_pdg_p->at(tI));
	    truthOut_type_p->push_back(truth_type_p->at(tI));
	    truthOut_origin_p->push_back(truth_origin_p->at(tI));
	    truthOut_status_p->push_back(truth_status_p->at(tI));	  
	    ++truthOut_n_;
	  }

	  //The following needs to be replaced according to what YJ found in MC truth photon definition
	  //Need to take more photons than just 14 type origin 37	 
	  //	  if(truth_type_p->at(tI) != 14) continue;
	  //	  if(truth_origin_p->at(tI) != 37) continue;

          int truthType = truth_type_p->at(tI);
          if(!( truthType >= 14 && truthType <= 16)) continue;
	  //  IsoPhoton         =  14,
	  //  NonIsoPhoton      =  15,
	  //  BkgPhoton         =  16,

	  int truthOrigin = truth_origin_p->at(tI);
          //https://gitlab.cern.ch/atlas/athena/-/blob/21.2/PhysicsAnalysis/MCTruthClassifier/MCTruthClassifier/MCTruthClassifierDefs.h
          if((truthOrigin<=35) || (truthOrigin>=41)) continue;
	  // BremPhot      = 36,
	  // PromptPhot    = 37,
	  // UndrPhot      = 38,
	  // ISRPhot       = 39,
	  // FSRPhot       = 40,

	  if(truth_pdg_p->at(tI) != 22) continue;

	  //Truth isolation via Yeonju Go (username YeonjuGo on github)
	  float genEtSum2 = 0;
	  float genEtSum3 = 0;
	  float genEtSum4 = 0;
	  for(unsigned int tI2 = 0; tI2 <truth_pt_p->size() ; ++tI2){
	    if( tI2==tI) continue;
	    if( truth_status_p->at(tI2) != 1) continue;

	    //	    if( (truth_type_p->at(tI2) == 14 && truth_origin_p->at(tI2)==37 && truth_pdg_p->at(tI2) == 22)) continue;

	    TLorentzVector temp;
	    temp.SetPtEtaPhiE(truth_pt_p->at(tI2), truth_eta_p->at(tI2), truth_phi_p->at(tI2), truth_e_p->at(tI2));
	    Float_t dR = getDR(truth_eta_p->at(tI), truth_phi_p->at(tI), truth_eta_p->at(tI2), truth_phi_p->at(tI2));
	    if(dR < 0.4) genEtSum4 += temp.Et();
	    if(dR < 0.3) genEtSum3 += temp.Et();
	    if(dR < 0.2) genEtSum2 += temp.Et();
	  }
	  //End truth isolation calculation

	  if(truthPhotonPt_ > 0){
	    if(doGlobalDebug) std::cout << "WARNING: MORE THAN ONE GEN LEVEL PHOTON FOUND IN FILE, ENTRY: " << file << ", " << entry << std::endl;
	    if(truth_pt_p->at(tI) > truthPhotonPt_){
	      truthPhotonPt_ = truth_pt_p->at(tI);
	      truthPhotonPhi_ = truth_phi_p->at(tI);
	      truthPhotonEta_ = truth_eta_p->at(tI);
	      truthPhotonType_ = truth_type_p->at(tI);
	      truthPhotonOrigin_ = truth_origin_p->at(tI);
	      truthPhotonIso2_ = genEtSum2; 
	      truthPhotonIso3_ = genEtSum3; 
	      truthPhotonIso4_ = genEtSum4; 
	    }
	  }
	  else{
	    truthPhotonPt_ = truth_pt_p->at(tI);
	    truthPhotonPhi_ = truth_phi_p->at(tI);
	    truthPhotonEta_ = truth_eta_p->at(tI);
	    truthPhotonType_ = truth_type_p->at(tI);
	    truthPhotonOrigin_ = truth_origin_p->at(tI);
	    truthPhotonIso2_ = genEtSum2; 
	    truthPhotonIso3_ = genEtSum3; 
	    truthPhotonIso4_ = genEtSum4; 
	  }
	}
      }

      subTimer3.stop();
      subTimer4.start();

      outTree_p->Fill();
      ++currTotalEntries;

      subTimer4.stop();
    }

    inFile_p->Close();
    delete inFile_p;

    if(nTerm > 0){
      if(currTotalEntries > (ULong64_t)nTerm) break;
    }

    ++nFile;
  }

  std::cout << "RANGES: " << std::endl;
  std::cout << " RANGE NVERT: " << minNVert << "-" << maxNVert << std::endl;
  std::cout << " RANGE NTRUTH: " << minNTruth << "-" << maxNTruth << std::endl;
  std::cout << " RANGE NPHOTON: " << minNPhoton << "-" << maxNPhoton << std::endl;
  std::cout << " RANGE NRECOR2: " << minNRecoR2 << "-" << maxNRecoR2 << std::endl;
  std::cout << " RANGE NRECOR4: " << minNRecoR4 << "-" << maxNRecoR4 << std::endl;
  std::cout << " RANGE NRECOR10: " << minNRecoR10 << "-" << maxNRecoR10 << std::endl;
  std::cout << " RANGE NRECOR2TO10: " << minNRecoR2to10 << "-" << maxNRecoR2to10 << std::endl;
  std::cout << " RANGE NTRUTHR2: " << minNTruthR2 << "-" << maxNTruthR2 << std::endl;
  std::cout << " RANGE NTRUTHR4: " << minNTruthR4 << "-" << maxNTruthR4 << std::endl;
  std::cout << " RANGE NTRUTHR10: " << minNTruthR10 << "-" << maxNTruthR10 << std::endl;
  std::cout << " RANGE NTRUTHR2TO10: " << minNTruthR2to10 << "-" << maxNTruthR2to10 << std::endl;
  
  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  
  TEnv outConfig;
  for(auto const & val : configMap){  
    for(unsigned int vI = 0; vI < val.second.size(); ++vI){
      std::string configVal = val.first;
      if(val.second.size() != 1) configVal = configVal + "_" + std::to_string(vI);
      outConfig.SetValue(configVal.c_str(), val.second[vI].c_str());
    }    
  }

  outConfig.Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;

  delete inConfig_p;
  
  std::cout << "GDJMCNTUPLEPREPROC COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNtuplePreProc.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjNtuplePreProc(argv[1]);
  return retVal;
}
