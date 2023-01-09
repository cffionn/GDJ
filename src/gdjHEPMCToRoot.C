//Author: Chris McGinn (2022.11.18)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//HepMC
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

//FastJet
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/stringUtil.h"

/// Utility function copied from example_UsingIterators.cc in HepMC
/// this predicate returns true if the input has no decay vertex
class IsStateFinal {
public:
    /// returns true if the GenParticle does not decay
    bool operator()( const HepMC::GenParticle* p ) {
        if ( !p->end_vertex() && p->status()==1 ) return 1;
        return 0;
    }
};

int gdjHEPMCToRoot(std::string inConfigFileName)
{
  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INHEPFILENAME",
					      "OUTROOTFILENAME",
					      "JETRS",
					      "JETPTMIN"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  std::string inHEPFileName = config_p->GetValue("INHEPFILENAME", "");
  std::string outROOTFileName = config_p->GetValue("OUTROOTFILENAME", "");
  std::string jetRs = config_p->GetValue("JETRS", "");

  const Double_t jtPtMin = config_p->GetValue("JETPTMIN", 15.0);
  
  std::vector<int> jetRVect = strToVectI(jetRs);
  const Int_t nMaxJetRs = 7;
  std::vector<int> validJetRs = {2, 3, 4, 5, 6, 8, 10};
  for(unsigned int i = 0; i < jetRVect.size(); ++i){
    if(!vectContainsInt(jetRVect[i], &validJetRs)){
      std::cout << "Error: JetR given \'" << jetRVect[i] << "\', " << std::endl;
    }
  }

  std::cout << "JET Rs, #: " << jetRVect.size() << std::endl;
  for(unsigned int jI = 0; jI < jetRVect.size(); ++jI){
    std::cout << " " << jI << ": " << jetRVect[jI] << std::endl;
  }

  
  //Check the outrootfile directory exists or add directory path to outfilename
  if(outROOTFileName.find("/") == std::string::npos){
    outROOTFileName = outputStr + "/" + outROOTFileName;
  }
  else{
    std::string dirName = outROOTFileName.substr(0, outROOTFileName.rfind("/"));
    check.doCheckMakeDir(dirName);
  }

  Double_t evtWeight_;
  
  const Int_t nMaxPart = 10000;
  Int_t nPart_;
  Float_t pt_[nMaxPart];
  Float_t eta_[nMaxPart];
  Float_t phi_[nMaxPart];
  Float_t m_[nMaxPart];
  Int_t pid_[nMaxPart];

  const Int_t nMaxJets = 100;
  Int_t nJt_[nMaxJetRs];
  Float_t jtpt_[nMaxJetRs][nMaxJets];
  Float_t jtphi_[nMaxJetRs][nMaxJets];
  Float_t jteta_[nMaxJetRs][nMaxJets];
  Float_t jtm_[nMaxJetRs][nMaxJets];
      
  TFile* outFile_p = new TFile(outROOTFileName.c_str(), "RECREATE");
  TTree* jewelTree_p = new TTree("jewelTree", "");
  
  jewelTree_p->Branch("evtWeight", &evtWeight_, "evtWeight/D");
  jewelTree_p->Branch("nPart", &nPart_, "nPart/I");
  jewelTree_p->Branch("pt", pt_, "pt[nPart]/F");
  jewelTree_p->Branch("eta", eta_, "eta[nPart]/F");
  jewelTree_p->Branch("phi", phi_, "phi[nPart]/F");
  jewelTree_p->Branch("m", m_, "m[nPart]/F");
  jewelTree_p->Branch("pid", pid_, "pid[nPart]/I");  				  

  for(unsigned int jI = 0; jI < jetRVect.size(); ++jI){
    std::string jetR = std::to_string(jetRVect[jI]);
    
    jewelTree_p->Branch(("nJtR" + jetR).c_str(), &nJt_[jI], ("nJtR" + jetR + "/I").c_str());
    jewelTree_p->Branch(("jtptR" + jetR).c_str(), jtpt_[jI], ("jtptR" + jetR + "[nJtR" + jetR + "]/F").c_str());
    jewelTree_p->Branch(("jtetaR" + jetR).c_str(), jteta_[jI], ("jtetaR" + jetR + "[nJtR" + jetR + "]/F").c_str());
    jewelTree_p->Branch(("jtphiR" + jetR).c_str(), jtphi_[jI], ("jtphiR" + jetR + "[nJtR" + jetR + "]/F").c_str());
    jewelTree_p->Branch(("jtmR" + jetR).c_str(), jtm_[jI], ("jtmR" + jetR + "[nJtR" + jetR + "]/F").c_str());
  }

  //Quickly analyze hep file name for centrality
  std::string centStr = "PP";
  if(inHEPFileName.find("Cent") != std::string::npos){
    centStr = inHEPFileName.substr(inHEPFileName.find("Cent"), inHEPFileName.size());
    centStr = centStr.substr(0, centStr.rfind("."));
  }
  
  //Following hepmc example_EventSelection.cc
  HepMC::IO_GenEvent ascii_in(inHEPFileName.c_str(),std::ios::in);
  HepMC::GenEvent* evt = ascii_in.read_next_event();
  ULong64_t nEvents = 0;
  while(evt){
    evtWeight_ = evt->weights().weights()[0];
    
    //grab all final state particles
    IsStateFinal isFinal;
    std::vector<HepMC::GenParticle*> finalStateParticles;
    for(HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p){
      if(isFinal(*p)) finalStateParticles.push_back(*p);
    }
    
    nPart_ = 0;
    std::vector<fastjet::PseudoJet> particlesForCluster;
    for(auto const & particle : finalStateParticles){
      pt_[nPart_] = particle->momentum().perp();
      eta_[nPart_] = particle->momentum().eta();
      phi_[nPart_] = particle->momentum().phi();
      m_[nPart_] = particle->momentum().m();
      pid_[nPart_] = particle->pdg_id();

      TLorentzVector tL;
      tL.SetPtEtaPhiM(pt_[nPart_], eta_[nPart_], phi_[nPart_], m_[nPart_]);
      particlesForCluster.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));
      
      ++nPart_;
    }

    if(nPart_ > nMaxPart){
      std::cout << "ERROR: nPart=" << nPart_ << " exceeds maximum allowed value nMaxPart=" << nMaxPart << ". Please extend array size. return 1" << std::endl;
      return 1;
    }

    //Cluster jets
    for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
      double jetR = jetRVect[rI];
      jetR /= 10.0;

      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
      fastjet::ClusterSequence cs(particlesForCluster, jetDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(jtPtMin));

      nJt_[rI] = 0;

      for(unsigned int jI = 0; jI < jets.size(); ++jI){
	jtpt_[rI][jI] = jets[jI].pt();
	jtphi_[rI][jI] = jets[jI].phi_std();
	jteta_[rI][jI] = jets[jI].eta();
	jtm_[rI][jI] = jets[jI].m();

	++nJt_[rI];
      }
    }
    
    jewelTree_p->Fill();
    
    //iterate event counter
    ++nEvents;

    //Clean and next event
    delete evt;
    ascii_in >> evt;
  }
  
  std::cout << "NEVENTS: " << nEvents << std::endl;

  outFile_p->cd();

  jewelTree_p->Write("", TObject::kOverwrite);
  delete jewelTree_p;

  config_p->SetValue("CENT", centStr.c_str());
  config_p->Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;
  
  //Clean all news
  delete config_p;  

  std::cout << "GDJHEPMCTOROOT complete. return 0" << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHEPMCToRoot.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHEPMCToRoot(argv[1]);
  return retVal;
}
