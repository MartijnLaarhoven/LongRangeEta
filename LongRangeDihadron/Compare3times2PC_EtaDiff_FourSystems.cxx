#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <string>

void Compare3times2PC_EtaDiff_FourSystems() {
    const char* fileNeNe = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileOO   = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
    const char* filePO   = "./3times2PC/Vn_LHC25ad_pass2_644389_Cent_0_20.root";
    const char* filePP   = "./3times2PC/Vn_LHC24af_pass1_644663_Cent_0_20.root";

    TFile* fNeNe = TFile::Open(fileNeNe, "READ");
    TFile* fOO   = TFile::Open(fileOO,   "READ");
    TFile* fPO   = TFile::Open(filePO,   "READ");
    TFile* fPP   = TFile::Open(filePP,   "READ");

    if (!fNeNe || !fNeNe->IsOpen()) { std::cerr << "Cannot open Ne-Ne file!" << std::endl; return; }
    if (!fOO   || !fOO->IsOpen())   { std::cerr << "Cannot open O-O file!" << std::endl; return; }
    if (!fPO   || !fPO->IsOpen())   { std::cerr << "Cannot open p-O file!" << std::endl; return; }
    if (!fPP   || !fPP->IsOpen())   { std::cerr << "Cannot open pp file!" << std::endl; return; }

    TH1D* hNeNe = (TH1D*)fNeNe->Get("hV2");
    TH1D* hOO   = (TH1D*)fOO->Get("hV2");
    TH1D* hPO   = (TH1D*)fPO->Get("hV2");
    TH1D* hPP   = (TH1D*)fPP->Get("hV2_Combined");
    if (!hPP) hPP = (TH1D*)fPP->Get("hV2");

    if (!hNeNe || !hOO || !hPO || !hPP) {
        std::cerr << "Missing v2 histogram in one or more files." << std::endl;
        return;
    }

    hNeNe->SetMarkerStyle(20); hNeNe->SetMarkerColor(kMagenta+2); hNeNe->SetLineColor(kMagenta+2);
    hOO->SetMarkerStyle(21);   hOO->SetMarkerColor(kBlue+1);     hOO->SetLineColor(kBlue+1);
    hPO->SetMarkerStyle(24);   hPO->SetMarkerColor(kGreen+2);    hPO->SetLineColor(kGreen+2);
    hPP->SetMarkerStyle(25);   hPP->SetMarkerColor(kRed+1);      hPP->SetLineColor(kRed+1);

    double maxVal = 0;
    for (int i = 1; i <= hNeNe->GetNbinsX(); ++i) maxVal = std::max(maxVal, hNeNe->GetBinContent(i) + hNeNe->GetBinError(i));
    for (int i = 1; i <= hOO->GetNbinsX();   ++i) maxVal = std::max(maxVal, hOO->GetBinContent(i) + hOO->GetBinError(i));
    for (int i = 1; i <= hPO->GetNbinsX();   ++i) maxVal = std::max(maxVal, hPO->GetBinContent(i) + hPO->GetBinError(i));
    for (int i = 1; i <= hPP->GetNbinsX();   ++i) maxVal = std::max(maxVal, hPP->GetBinContent(i) + hPP->GetBinError(i));
    if (maxVal <= 0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff4Sys", "3x2PC Eta-diff v2: Ne-Ne vs O-O vs p-O vs pp", 900, 650);
    hNeNe->SetTitle("3x2PC #eta-differential v_{2}: Ne-Ne vs O-O vs p-O vs pp; #eta; v_{2}");
    hNeNe->GetYaxis()->SetRangeUser(0, maxVal * 1.25);
    hNeNe->Draw("E");
    hOO->Draw("E same");
    hPO->Draw("E same");
    hPP->Draw("E same");

    TLegend* leg = new TLegend(0.46, 0.66, 0.90, 0.90);
    leg->AddEntry(hNeNe, "Ne-Ne (LHC25af_pass2_645746)", "lep");
    leg->AddEntry(hOO,   "O-O (LHC25ae_pass2_645657)",   "lep");
    leg->AddEntry(hPO,   "p-O (LHC25ad_pass2_644389)",   "lep");
    leg->AddEntry(hPP,   "pp (LHC24af_pass1_644663)",    "lep");
    leg->Draw();

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp.root");
    std::cout << "Saved eta-diff v2 comparison plot for Ne-Ne, O-O, p-O, pp." << std::endl;

    fNeNe->Close();
    fOO->Close();
    fPO->Close();
    fPP->Close();
    delete fNeNe;
    delete fOO;
    delete fPO;
    delete fPP;
    delete c;
}
