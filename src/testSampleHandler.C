//Author: Chris McGinn (2020.02.27)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//Local
#include "include/sampleHandler.h"

int testSampleHandler(bool isPP, bool isMC, int year, int minPthat)
{
  sampleHandler test;
  test.Init(isPP, isMC, year, minPthat);
  test.PrintTags();

  int tag = test.GetTag();
  double xsec = test.GetXSection();
  double filterEff = test.GetFilterEff();

  std::cout << std::endl;
  std::cout << "Your tag, xsec, eff: " << tag << ", " << xsec << ", " << filterEff << std::endl;
  
  std::cout << "TESTSAMPLEHANDLER COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/testSampleHandler.exe <isPP> <isMC> <year> <minPthat>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += testSampleHandler(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]));
  return retVal;
}
