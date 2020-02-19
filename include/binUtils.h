#ifndef BINUTILS_H
#define BINUTILS_H

bool goodBinning(std::string config, Int_t nMaxBins, Int_t nBins, std::string binType)
{
  if(nBins > nMaxBins){
    std::cout << "HISTMAKING ERROR IN GOODBINNING - config \'" << config << "\' contains " << binType << "  of size \'" << nBins << "', greater than max val \'" << nMaxBins << "\'. Please fix. return false" << std::endl;
    return false;
  }
  return true;
}

#endif
