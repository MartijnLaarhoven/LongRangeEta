#include <TFile.h>
#include <TH1.h>
#include <iostream>

void printOne(const char* tag, const char* pathFull, const char* pathInner, const char* pathOuter, int n){
  TFile* fF=TFile::Open(pathFull,"READ");
  TFile* fI=TFile::Open(pathInner,"READ");
  TFile* fO=TFile::Open(pathOuter,"READ");
  if(!fF||!fI||!fO){ std::cout<<"MISSING FILES for "<<tag<<" v"<<n<<"\n"; return; }
  TH1* hF=(TH1*)fF->Get(Form("hV%d_Combined",n));
  TH1* hI=(TH1*)fI->Get(Form("hV%d_Sides",n));
  TH1* hO=(TH1*)fO->Get(Form("hV%d_Sides",n));
  TH1* hIcmb=(TH1*)fI->Get(Form("hV%d_Combined",n));
  TH1* hOcmb=(TH1*)fO->Get(Form("hV%d_Combined",n));
  std::cout<<"\n["<<tag<<" v"<<n<<"]\n";
  std::cout<<"full has hVn_Combined: "<<(hF?"yes":"no")<<"\n";
  std::cout<<"inner has hVn_Sides:   "<<(hI?"yes":"no")<<" , hVn_Combined: "<<(hIcmb?"yes":"no")<<"\n";
  std::cout<<"outer has hVn_Sides:   "<<(hO?"yes":"no")<<" , hVn_Combined: "<<(hOcmb?"yes":"no")<<"\n";
  if(hF){
    std::cout<<"full FT0C(bin1)="<<hF->GetBinContent(1)<<" +/- "<<hF->GetBinError(1)
             <<", FT0A(last)="<<hF->GetBinContent(hF->GetNbinsX())<<" +/- "<<hF->GetBinError(hF->GetNbinsX())<<"\n";
  }
  auto p=[&](const char* name, TH1* h){
    if(!h) return;
    std::cout<<name<<" sides bin1="<<h->GetBinContent(1)<<" +/- "<<h->GetBinError(1)
             <<", bin2="<<h->GetBinContent(2)<<" +/- "<<h->GetBinError(2)<<"\n";
  };
  p("inner",hI); p("outer",hO);
  fF->Close(); fI->Close(); fO->Close();
}

void check_ring_full(){
  const char* ooF="./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
  const char* ooI="./3times2PC/Vn_LHC25ae_pass2_innerRing_Cent_0_20.root";
  const char* ooO="./3times2PC/Vn_LHC25ae_pass2_outerRing_Cent_0_20.root";
  const char* neF="./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
  const char* neI="./3times2PC/Vn_LHC25af_pass2_innerRing_Cent_0_20.root";
  const char* neO="./3times2PC/Vn_LHC25af_pass2_outerRing_Cent_0_20.root";
  for(int n=2;n<=4;++n){
    printOne("O-O",ooF,ooI,ooO,n);
    printOne("Ne-Ne",neF,neI,neO,n);
  }
}
