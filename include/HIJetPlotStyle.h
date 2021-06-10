#ifndef HISTMAKINGTOOL_HIJETPLOTSTYLE_H
#define HISTMAKINGTOOL_HIJETPLOTSTYLE_H

#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "TGraph.h"
#include <TF1.h>

namespace HIJet
{
  namespace Style
  {
    const int NSTYLES=15;
    const int color_scheme[NSTYLES]={1,kRed-4,kAzure-3,kGreen+2,kMagenta+2,
				     kOrange+2,kRed-4,kAzure-3,kGreen+2,kMagenta+2,
				     kOrange+2,kCyan+3,28,41,kGray};

    const int style_scheme[NSTYLES]={20,21,33,34,29,
    				     24,25,27,28,30,
    				     23,20,21,33,34};
  
    const float size_scheme[NSTYLES]={1,1,1.6,1.2,1.6,1,
				      1,1,1.6,1,1,
				      1,1,1.6,1.2};

    const float opacity_scheme[NSTYLES]={0.3,0.3,0.3,0.3,0.3,
					 0.3,0.3,0.3,0.3,0.3,
					 0.3,0.3,0.3,0.3,0.3};
    
    const Width_t line_width=3;
    const int Npx=1000;

    inline void EquipHistogram(TH1* h1, unsigned int index)
    {
      unsigned int index_mod=index%NSTYLES;
      h1->SetMarkerColor(color_scheme[index_mod]);
      h1->SetLineColor(color_scheme[index_mod]);
      h1->SetLineWidth(line_width);
      h1->SetMarkerStyle(style_scheme[index_mod]);
      h1->SetMarkerSize(size_scheme[index_mod]);
    }
    inline void EquipGraph(TGraph* g1, unsigned int index)
    {
      unsigned int index_mod=index%NSTYLES;
      g1->SetMarkerColor(color_scheme[index_mod]);
      g1->SetLineColor(color_scheme[index_mod]);
      g1->SetLineWidth(line_width);
      g1->SetMarkerStyle(style_scheme[index_mod]);
      g1->SetMarkerSize(size_scheme[index_mod]);
    }
    inline void EquipTF1(TF1* f1, unsigned int index)
    {
      unsigned int index_mod=index%NSTYLES;
      f1->SetLineColor(color_scheme[index_mod]);
      f1->SetLineStyle( (index_mod+1)%2 + 1);
      f1->SetNpx(Npx);
    }//Next 2 functions added by CFM for error box opacities tied to marker/color scheme    
    inline int GetColor(unsigned int index){return color_scheme[index%NSTYLES];}
    inline float GetOpacity(unsigned int index){return opacity_scheme[index%NSTYLES];}    
    
    inline void SetAtlasStyle(std::string inStyleFileName)
    {
      std::cout << "INSTYLEFILENAME: " << inStyleFileName << std::endl;
      ///
      TFile* inFile_p = TFile::Open(inStyleFileName.c_str());
      if(!inFile_p)
      {
	std::cerr << "HIJet::SetAtlasStyle: Could not open file " << inStyleFileName << ". Skipping style settings." << std::endl;
	return;
      }
      TStyle* atlas_style=(TStyle*)inFile_p->Get("ATLAS");
      if(!atlas_style)
      {
	std::cerr << "HIJet::SetAtlasStyle: Could find TStyle object named ATLAS in " << inStyleFileName << ". Skipping style settings." << std::endl;
	return;
      }
      std::cout << "Applying ATLAS style settings...\n" << std::endl;
      gROOT->SetStyle("ATLAS");
      gROOT->ForceStyle();
    }

    inline void SetAtlasStyle()
    {
      std::string calibEnvVar(gSystem->Getenv("CALIBPATH"));
      std::string calibPath=calibEnvVar.substr(0,calibEnvVar.find_first_of(":"));
      calibPath.append("/HistMakingTool/");

      std::string inStyleFileName("AtlasStyle.root");
      if(calibPath.size()>1) inStyleFileName.insert(0,calibPath);
      SetAtlasStyle(inStyleFileName);
    }

  }
}
#endif
