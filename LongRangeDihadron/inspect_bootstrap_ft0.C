#include <TFile.h>
#include <TH1D.h>
#include <cstdio>

void stat(const char* tag,const char* path,const char* hname){
  TFile* f=TFile::Open(path,"READ");
  if(!f||!f->IsOpen()){printf("%s: MISSING file %s\n",tag,path);return;}
  TH1D* h=(TH1D*)f->Get(hname);
  if(!h){printf("%s: missing hist %s\n",tag,hname); f->Close(); delete f; return;}
  printf("%s: integral=%.8g, min=%.6g, max=%.6g, mean=%.6g\n",tag,h->Integral(),h->GetMinimum(),h->GetMaximum(),h->GetMean());
  f->Close(); delete f;
}

void inspect_bootstrap_ft0(){
  stat("648799 data 0-20","./ProcessOutput/BootstrapSample_LHC25ae_pass2_648799_Cent_0_20_FT0A_FT0C.root","bsSample_hPhiSameOverMixed_0_20");
  stat("648799 templ 80-100","./ProcessOutput/BootstrapSample_LHC25ae_pass2_648799_Cent_80_100_FT0A_FT0C.root","bsSample_hPhiSameOverMixed_80_100");
  stat("648800 data 0-20","./ProcessOutput/BootstrapSample_LHC25ae_pass2_648800_Cent_0_20_FT0A_FT0C.root","bsSample_hPhiSameOverMixed_0_20");
  stat("648800 templ 80-100","./ProcessOutput/BootstrapSample_LHC25ae_pass2_648800_Cent_80_100_FT0A_FT0C.root","bsSample_hPhiSameOverMixed_80_100");
}
