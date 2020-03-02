//Author: Chris McGinn (2020.02.27)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>

//Local
#include "include/sampleHandler.h"
#include "include/stringUtil.h"

sampleHandler::sampleHandler(bool in_isPP, bool in_isMC, int in_year, int in_minPthat)
{
  Init(in_isPP, in_isMC, in_year, in_minPthat);
  return;
}
 
sampleHandler::~sampleHandler()
{
  Clean();
  return;
}

void sampleHandler::PreInit()
{
  validMinPthatsByYear[2015] = {};
  validMinPthatsByYear[2017] = {35, 50, 70, 140, 280};
  validMinPthatsByYear[2018] = {50, 70, 140};

  //PYTHIA8 + Overlay
  dataSetNameToIsPP["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 50;

  dataSetNameToIsPP["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 70;
  
  dataSetNameToIsPP["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 140;

  //PYTHIA8
  dataSetNameToIsPP["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 35;
 
  dataSetNameToIsPP["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 50;

  dataSetNameToIsPP["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 70;

  dataSetNameToIsPP["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 140;

  dataSetNameToIsPP["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 280;

  for(int ppI = 0; ppI < 2; ++ppI){
    for(int mcI = 0; mcI < 2; ++mcI){
      for(auto const & year : validYears){
	if(ppI == 0 && year == 2017) continue;
	if(ppI == 1 && (year == 2015 || year == 2018)) continue;

	for(auto const & min : validMinPthatsByYear[year]){
	  int tag = CreateTag(ppI, mcI, year, min);

	  //Values as extracted in AMI
	  //PYT.8 search: mc16_5TeV.%.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP%.merge.AOD.e5094_s3238_r10441_r10210
	  //PYT.8+Overlay search: mc16_5TeV.%.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP%.merge.AOD.e5094_d1516_r11439_r11217
	  
	  if(tag == 35201711){
	    tagToXSec[tag] = 351620;//in nanobarn
	    tagToFilterEff[tag] = 0.000029108;
	  }
	  else if(tag == 50201711){
	    tagToXSec[tag] = 85898;//in nanobarn
	    tagToFilterEff[tag] = 0.00003339;
	  }
	  else if(tag == 70201711){
	    tagToXSec[tag] = 21551;//in nanobarn
	    tagToFilterEff[tag] = 0.000045787;
	  }
	  else if(tag == 140201711){
	    tagToXSec[tag] = 1044;//in nanobarn
	    tagToFilterEff[tag] = 0.000050981;
	  }
	  else if(tag == 280201711){
	    tagToXSec[tag] = 37.592;//in nanobarn
	    tagToFilterEff[tag] = 0.000043848;
	  }
	  else if(tag == 50201810){
	    tagToXSec[tag] = 85898;//in nanobar
	    tagToFilterEff[tag] = 0.00003339;
	  }
	  else if(tag == 70201810){
	    tagToXSec[tag] = 21551;//in nanobarn
	    tagToFilterEff[tag] = 0.000045787;
	  }
	  else if(tag == 140201810){
	    tagToXSec[tag] = 1044;//in nanobarn
	    tagToFilterEff[tag] = 0.000050981;
	  }
	  else{
	    tagToXSec[tag] = 0.0;//in nanobarn
	    tagToFilterEff[tag] = 0.0;
	  }
	}
      }      
    }
  }

  return;
}

bool sampleHandler::Init(std::string sampleString)
{
  //Pre-Initialize the maps
  PreInit();
     
  m_isInit = false;

  if(dataSetNameToIsPP.count(sampleString) == 0){
    std::cout << "sampleHandler error - Given sample string \'" << sampleString << "\' is not valid. Initialization failed." << std::endl;
    std::cout << "To fix, pick from: " << std::endl;
    for(auto const & sample : dataSetNameToIsPP){
      std::cout << " " << sample.first << std::endl;
    }
    std::cout << "Or modify sampleHandler class. return false" << std::endl;
    Clean();    
    return m_isInit;
  }

  m_isInit = true;

  m_isPP = dataSetNameToIsPP[sampleString];
  m_isMC = dataSetNameToIsMC[sampleString];
  m_year = dataSetNameToYear[sampleString];
  m_minPthat = dataSetNameToMinPthat[sampleString];

  m_tagVal = CreateTag();
  
  return m_isInit;
}

bool sampleHandler::Init(bool in_isPP, bool in_isMC, int in_year, int in_minPthat)
{
  //Pre-Initialize the maps
  PreInit();
     
  m_isInit = false;

  if(!vectContainsInt(in_year, &validYears)){
    std::cout << "sampleHandler error - Given year val \'" << in_year << "\' is not valid. Initialization failed." << std::endl;
    std::cout << " To fix, pick from: ";
    for(auto const & year : validYears){
      std::cout << year << ", ";
    }
    std::cout << " or modify sampleHandler class. return false" << std::endl;
    Clean();
    return m_isInit;
  }
  
  if(in_isMC){
    if(!vectContainsInt(in_minPthat, &(validMinPthatsByYear[in_year]))){
      std::cout << "sampleHandler error - Given year/minPthat vals \'" << in_year << "/" << in_minPthat << "\' is not valid. Initialization failed." << std::endl;
      std::cout << " To fix, pick from: ";
      for(auto const & iter : validMinPthatsByYear){
	std::cout << iter.first << "{";
	std::string goodStuff = "";
	for(auto const & iter2 : iter.second){
	  goodStuff = goodStuff + std::to_string(iter2) + ", ";
	}
	if(goodStuff.find(",") != std::string::npos) goodStuff.replace(goodStuff.rfind(","), goodStuff.size(), "");

	std::cout << goodStuff << "}, ";
      }
      std::cout << " or modify sampleHandler class. return false" << std::endl;
      
      Clean();
      return m_isInit;
    }    
  }

  m_isInit = true;

  m_isPP = in_isPP;
  m_isMC = in_isMC;
  m_year = in_year;
  m_minPthat = in_minPthat;

  m_tagVal = CreateTag();
  
  return m_isInit;
}

int sampleHandler::GetTag()
{
  if(!m_isInit){
    std::cout << "sampleHandler GetTag Error - Called w/o proper initialization. return tag -1" << std::endl;
  }
  
  return m_tagVal;
}

double sampleHandler::GetXSection()
{
  double xsec = 0.0;
  if(m_isInit) xsec = tagToXSec[m_tagVal];
  else std::cout << "SAMPLEHANDLER ERROR - GetXSection() called despite no init. returning 0.0" << std::endl;
  return xsec;
}

double sampleHandler::GetFilterEff()
{
  double eff = 0.0;
  if(m_isInit) eff = tagToFilterEff[m_tagVal];
  else std::cout << "SAMPLEHANDLER ERROR - GetFilterEff() called despite no init. returning 0.0" << std::endl;
  return eff;
}

void sampleHandler::Clean()
{
  m_isInit = false;

  m_isPP = false;
  m_isMC = false;
  m_year = -1;
  m_minPthat = -1;

  m_tagVal = -1;
  
  return;
}

void sampleHandler::PrintTags()
{
  std::cout << "SAMPLEHANDLER - PrintTags: " << std::endl;
  for(auto const & tag : tagToXSec){
    std::cout << " " << tag.first << ": " << tag.second << ", " << tagToFilterEff[tag.first] << std::endl;
  }
  
  return;
}

int sampleHandler::CreateTag()
{
  int retVal = -1;
  if(m_isInit) retVal = CreateTag(m_isPP, m_isMC, m_year, m_minPthat);
  else std::cout << "sampleHandler CreateTag Error - Called w/o proper initialization. return tag -1" << std::endl;
  
  return retVal;
}

int sampleHandler::CreateTag(bool in_isPP, bool in_isMC, int in_year, int in_minPthat)
{
  int isPP = (int)in_isPP;
  int isMC = (int)in_isMC;
  
  int retVal = isPP + isMC*10 + in_year*100 + in_minPthat*1000000;
  
  return retVal;
}
