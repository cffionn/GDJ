//Author: Chris McGinn (2020.01.24)
#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

//cpp
#include <fstream>
#include <map>
#include <string>

//ROOT
#include "TEnv.h"

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"

class configParser
{
 public:
  configParser(){};
  configParser(std::string inConfigFileName);
  configParser(TEnv* inConfigEnv_p);
  ~configParser(){};

  bool Init(std::string inConfigFileName);
  bool Init(TEnv* inConfigEnv_p);
  std::string GetConfigVal(std::string inStr);
  std::map<std::string, std::string> GetConfigMap();
  void Clean();
  
 private:
  checkMakeDir check;
  std::string m_configFileName;
  std::ifstream* m_configFile = nullptr;
  std::map<std::string, std::string> m_configVals;  
};

#endif
