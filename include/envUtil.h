#ifndef ENVUTIL_H
#define ENVUTIL_H

//cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "THashList.h"

//Local
#include "include/stringUtil.h"

inline bool checkEnvForParams(TEnv* inEnv_p, std::vector<std::string> inParams)
{
  THashList* hash_p = (THashList*)inEnv_p->GetTable();
  bool retVal = true;

  for(auto const & param : inParams){
    bool isFound = false;

    for(Int_t entry = 0; entry < hash_p->GetEntries(); ++entry){
      std::string name = hash_p->At(entry)->GetName();
      if(isStrSame(name, param)){
        isFound = true;
        break;
      }
    }

    if(!isFound){
      std::cout << "PLOTMAKING ERROR - Missing required parameter \'" << param << "\'. return false" << std::endl;
      retVal = false;
    }
  }

  return retVal;
}

#endif
