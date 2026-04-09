#include <TFile.h>
#include <TKey.h>
#include <cstdio>

void dumpKeys(const char* path){
  TFile* f=TFile::Open(path,"READ");
  if(!f||!f->IsOpen()){printf("MISSING %s\n",path);return;}
  printf("\n== %s ==\n", path);
  TIter next(f->GetListOfKeys());
  while(TKey* k=(TKey*)next()) printf("%s\n", k->GetName());
  f->Close(); delete f;
}

void list_keys(){
  dumpKeys("./ProcessOutput/BootstrapSample_LHC25ae_pass2_648799_Cent_0_20_FT0A_FT0C.root");
  dumpKeys("./ProcessOutput/BootstrapSample_LHC25ae_pass2_648799_Cent_80_100_FT0A_FT0C.root");
  dumpKeys("./ProcessOutput/BootstrapSample_LHC25ae_pass2_648800_Cent_0_20_FT0A_FT0C.root");
  dumpKeys("./ProcessOutput/BootstrapSample_LHC25ae_pass2_648800_Cent_80_100_FT0A_FT0C.root");
}
