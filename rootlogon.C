{
  gInterpreter->AddIncludePath(gSystem->Getenv("GDJDIR"));  
  gInterpreter->AddIncludePath(gSystem->Getenv("ROOUNFOLDDIR"));  
  gStyle->SetOptStat(0);
}
