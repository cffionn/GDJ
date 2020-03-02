//Author: Chris McGinn (2020.01.24)

//cpp
#include <iostream>
#include <vector>

//ROOT
#include "THashList.h"

//Local
#include "include/configParser.h"

configParser::configParser(std::string inConfigFileName)
{
  Init(inConfigFileName);
  return;
}

configParser::configParser(TEnv* inConfigEnv_p)
{
  Init(inConfigEnv_p);
  return;
}

configParser::configParser(TEnv* inConfigEnv_p, TEnv* inDefinitionEnv_p)
{
  Init(inConfigEnv_p, inDefinitionEnv_p);
  return;
}

bool configParser::Init(std::string inConfigFileName)
{
  Clean(); //Make sure all is reset in-case we forget a Clean call before invoking Init

  m_configFileName = inConfigFileName;
  if(!check.checkFile(m_configFileName)){
    std::cout << "configParser::Init - Given inConfigFileName \'" << inConfigFileName << "\' is invalid. Init failed, return false" << std::endl;
    Clean();
    return false;
  }
  
  m_configFile = new std::ifstream(m_configFileName.c_str());//Process our config file
  std::string tempStr;
  std::vector<std::string> tempVect;
  while(std::getline(*m_configFile, tempStr)){
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "),1,"");}
    if(tempStr.size()==0) continue;
    if(tempStr.substr(0,1).find("#") != std::string::npos) continue;//# is comment

    tempVect = commaSepStringToVect(tempStr);//We need to parse string list

    if(tempVect.size() < 2){
      std::cout << "configParser::Init - Given inConfigFileName \'" << inConfigFileName << "\' contains invalid line \'" << tempStr << "\'. Please fix (At least 2 items, comma sep). Clean and return false." << std::endl;
      Clean();
      return false;
    }
    std::string valStr = tempVect[1];
    for(unsigned int tI = 2; tI < tempVect.size(); ++tI){
      valStr = valStr + "," + tempVect[tI];
    }

    m_configVals[tempVect[0]] = valStr;
  }
  
  return true;
}

bool configParser::Init(TEnv* inConfigEnv_p)
{
  Clean();
  
  THashList* hash_p = (THashList*)inConfigEnv_p->GetTable();
  for(Int_t entry = 0; entry < hash_p->GetEntries(); ++entry){
    m_configVals[hash_p->At(entry)->GetName()] = inConfigEnv_p->GetValue(hash_p->At(entry)->GetName(), "");
  }

  for(auto const & vals : m_configVals){
    std::string sec = vals.second;
    while(sec.find(" ") != std::string::npos){
      sec.replace(sec.find(" "), 1, "");
    }
    m_configVals[vals.first] = sec;
  }

  return true;
}

bool configParser::Init(TEnv* inConfigEnv_p, TEnv* inDefinitionEnv_p)
{
  Init(inConfigEnv_p);

  THashList* hash_p = (THashList*)inDefinitionEnv_p->GetTable();
  for(Int_t entry = 0; entry < hash_p->GetEntries(); ++entry){
    m_definitionVals[hash_p->At(entry)->GetName()] = inDefinitionEnv_p->GetValue(hash_p->At(entry)->GetName(), "");
  }
  return true;
}

bool configParser::ContainsParam(std::string inStr)
{
  if(m_configVals.count(inStr) == 0){
    std::cout << "configParser::ContainsParam - Given inStr \'" << inStr << "\' not found. Please define in given config \'" << m_configFileName << "\'. return false" << std::endl;
    return false;
  }
  return true;
}

std::string configParser::GetConfigVal(std::string inStr)
{
  if(m_configVals.count(inStr) == 0){
    std::cout << "configParser::GetConfigVal - Given inStr \'" << inStr << "\' is not found. Please define in given config \'" << m_configFileName << "\'. return empty string." << std::endl;
    return "";
  }
  return m_configVals[inStr];
}

std::string configParser::GetDefinitionVal(std::string inStr)
{
  if(m_definitionVals.count(inStr) == 0){
    std::cout << "configParser::GetDefinitionVal - Given inStr \'" << inStr << "\' is not found. Please define in given config \'" << m_configFileName << "\'. return empty string." << std::endl;
    return "";
  }
  return m_definitionVals[inStr];
}

std::map<std::string, std::string> configParser::GetConfigMap(){return m_configVals;}

void configParser::Clean()
{
  m_configFileName = "";
  if(m_configFile != nullptr){
    m_configFile->close();
    delete m_configFile;
    m_configFile = nullptr;
  }
  m_configVals.clear();
  m_definitionVals.clear();
  
  return;
}

void configParser::SetConfigVal(std::string inStrParam, std::string inStrVal)
{
  if(m_configVals.count(inStrParam) != 0){
    std::cout << "CONFIGPARSER WARNING - SetConfigVal called to override, \'" << inStrParam << "\' with new value \'" << inStrVal << "\'" << std::endl;
  }
  m_configVals[inStrParam] = inStrVal;
  return;
}
