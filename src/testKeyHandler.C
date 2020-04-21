//Author: Chris McGinn (2020.02.27)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>

//ROOT
#include "TMath.h"

//Local
#include "include/getLinBins.h"
#include "include/ghostUtil.h"
#include "include/keyHandler.h"

int testKeyHandler(float cent, float vz, float evtPlane)
{
  const Int_t nCentBins = 100;
  Float_t centBinsLow = 0;
  Float_t centBinsHigh = 100;
  Double_t centBins[nCentBins+1];
  getLinBins(centBinsLow, centBinsHigh, nCentBins, centBins);

  const Int_t nVzBins = 30;
  Float_t vzBinsLow = -15.0;
  Float_t vzBinsHigh = 15.0;
  Double_t vzBins[nVzBins+1];
  getLinBins(vzBinsLow, vzBinsHigh, nVzBins, vzBins);
  
  const Int_t nEvtPlaneBins = 16;
  Float_t evtPlaneBinsLow = -TMath::Pi()/2.;
  Float_t evtPlaneBinsHigh = TMath::Pi()/2.;
  Double_t evtPlaneBins[nEvtPlaneBins+1];
  getLinBins(evtPlaneBinsLow, evtPlaneBinsHigh, nEvtPlaneBins, evtPlaneBins);

  keyHandler test;
  test.Init({nCentBins-1, nVzBins-1, nEvtPlaneBins-1});

  std::cout << "Key for cent, vz, evtPlane: " << cent << ", " << vz << ", " << evtPlane << std::endl;
  unsigned long long centPos = ghostPos(nCentBins, centBins, cent);
  unsigned long long vzPos = ghostPos(nVzBins, vzBins, vz);
  unsigned long long evtPlanePos = ghostPos(nEvtPlaneBins, evtPlaneBins, evtPlane);

  unsigned long long key = test.GetKey({centPos, vzPos, evtPlanePos});
  std::cout << " " << key << std::endl;
  std::cout << " BINS: " << centPos << ", " << vzPos << ", " << evtPlanePos << std::endl;
  std::cout << "Inverting key: " << std::endl;
  std::vector<unsigned long long> invertKey = test.InvertKey(key);

  for(unsigned int eI = 0; eI < invertKey.size(); ++eI){
    std::cout << " " << invertKey[eI] << std::endl;
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/testKeyHandler.exe <cent 0 to 100> <vz -15 to 15> <evt plane -pi/2 to pi/2>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += testKeyHandler(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]));
  return retVal;
}
