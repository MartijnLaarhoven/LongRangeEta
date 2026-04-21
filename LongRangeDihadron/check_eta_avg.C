#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include <cstdio>

void check(const char* path){
  TFile* f=TFile::Open(path);
  if(!f||f->IsZombie()){printf("open fail %s\n",path); return;}
  TH1D* h=(TH1D*)f->Get("hV2");
  if(!h){printf("no hV2 %s\n",path); f->Close(); return;}
  double sum=0, sw=0, swv=0, bwSum=0, bwv=0; int n=0;
  for(int i=1;i<=h->GetNbinsX();++i){
    double v=h->GetBinContent(i), e=h->GetBinError(i), bw=h->GetXaxis()->GetBinWidth(i);
    if(std::isfinite(v) && std::isfinite(e) && e>0 && e<9.9 && v>=0){
      n++; sum+=v; double w=1.0/(e*e); sw+=w; swv+=w*v; bwSum+=bw; bwv+=bw*v;
    }
  }
  printf("%s\n bins=%d mean=%g ivw=%g bw-mean=%g\n", path,n,sum/n,swv/sw,bwv/bwSum);
  f->Close();
}

void check_eta_avg(){
  check("TemplateFit/EtaDiff/VnDelta_LHC25af_pass2_650317_nch10_50_Mult_10_50_TPC_FT0C.root");
  check("TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_653254_Mult_10_50_TPC_FT0C.root");
  check("TemplateFit/EtaDiff/VnDelta_LHC25af_pass2_650316_nch10_50_Mult_10_50_TPC_FT0A.root");
  check("TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_653254_Mult_10_50_TPC_FT0A.root");
}

check_eta_avg();
