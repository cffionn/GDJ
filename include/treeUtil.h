#ifndef TREEUTIL_H
#define TREEUTIL_H

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TObjArray.h"
#include "TTree.h"

inline std::vector<std::string> getVectBranchList(TTree* inTree_p)
{
  std::vector<std::string> branchList;
  TObjArray* temp_p = (TObjArray*)inTree_p->GetListOfBranches();
  for(Int_t bI = 0; bI < temp_p->GetEntries(); ++bI){
    std::string branchStr = temp_p->At(bI)->GetName();
    branchList.push_back(branchStr);
  }  
  return branchList;
}

#endif
