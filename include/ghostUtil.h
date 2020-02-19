#ifndef GHOSTUTIL_H
#define GHOSTUTIL_H

//cpp
#include <iostream>
#include <vector>

inline int ghostPos(std::vector<float> bins_, double ghostVal, bool printWarnings)
{
  if(bins_.size() == 0){
    std::cout << "MAKECLUSTERTREE GHOSTPOS ERROR: Given bins have size \'0\'. returning int32 max for absurd result" << std::endl;
    return 2147483647;//lol                                                                  
  }

  int ghostPos = -1;
  if(ghostVal<bins_.at(0)){
    if(printWarnings) std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' falls below binning (low edge \'" << bins_[0] << "\'). return pos 0" << std::endl;
    ghostPos = 0;
  }
  else if(ghostVal>=bins_.at(bins_.size()-1)){
    ghostPos = bins_.size()-2;//-2 because etabins are 1 greater than rhobins                
    if(printWarnings) std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' is above binning (high edge \'" << bins_[bins_.size()-1] << "\'). return pos " << bins_.size()-2 << ", " << bins_[bins_.size()-2] << "-" << bins_[bins_.size()-1] << std::endl;
  }
  else{
    for(unsigned int ie = 0; ie < bins_.size()-1; ++ie){
      if(ghostVal>=bins_.at(ie) && ghostVal<bins_.at(ie+1)){
        ghostPos = ie;
        break;
      }
    }
  }
  return ghostPos;
}

inline int ghostPos(Int_t nBins, Double_t bins[], double ghostVal, bool printWarnings)
{
  std::vector<float> binsVect;
  for(Int_t bI = 0; bI < nBins+1; ++bI){
    binsVect.push_back(bins[bI]);
  }
  return ghostPos(binsVect, ghostVal, printWarnings);
}


inline int ghostPos(std::vector<int> bins_, double ghostVal, bool printWarnings)
{
  if(bins_.size() == 0){
    std::cout << "MAKECLUSTERTREE GHOSTPOS ERROR: Given bins have size \'0\'. returning int32 max for absurd result" << std::endl;
    return 2147483647;//lol                                                                  
  }

  int ghostPos = -1;
  if(ghostVal<bins_.at(0)){
    if(printWarnings) std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' falls below binning (low edge \'" << bins_[0] << "\'). return pos 0" << std::endl;
    ghostPos = 0;
  }
  else if(ghostVal>=bins_.at(bins_.size()-1)){
    ghostPos = bins_.size()-2;//-2 because etabins are 1 greater than rhobins                
    if(printWarnings) std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' is above binning (high edge \'" << bins_[bins_.size()-1] << "\'). return pos " << bins_.size()-2 << ", " << bins_[bins_.size()-2] << "-" << bins_[bins_.size()-1] << std::endl;
  }
  else{
    for(unsigned int ie = 0; ie < bins_.size()-1; ++ie){
      if(ghostVal>=bins_.at(ie) && ghostVal<bins_.at(ie+1)){
        ghostPos = ie;
        break;
      }
    }
  }
  return ghostPos;
}

#endif
