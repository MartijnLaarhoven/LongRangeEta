#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <iostream>

void Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe() {
    const char* fileNeNe = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileOO   = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";

    TFile* fNeNe = TFile::Open(fileNeNe, "READ");
    TFile* fOO   = TFile::Open(fileOO,   "READ");

    if (!fNeNe || !fNeNe->IsOpen()) {
        std::cerr << "Cannot open Ne-Ne file: " << fileNeNe << std::endl;
        return;
    }
    if (!fOO || !fOO->IsOpen()) {
        std::cerr << "Cannot open O-O file: " << fileOO << std::endl;
        fNeNe->Close();
        delete fNeNe;
        return;
    }

    TH1D* hNeNe = (TH1D*)fNeNe->Get("hV2_Combined");
    TH1D* hOO   = (TH1D*)fOO->Get("hV2_Combined");
    if (!hNeNe) hNeNe = (TH1D*)fNeNe->Get("hV2");
    if (!hOO)   hOO   = (TH1D*)fOO->Get("hV2");

    if (!hNeNe || !hOO) {
        std::cerr << "Missing hV2_Combined/hV2 in one or both files." << std::endl;
        fNeNe->Close();
        fOO->Close();
        delete fNeNe;
        delete fOO;
        return;
    }

    hNeNe->SetMarkerStyle(20);
    hNeNe->SetMarkerColor(kMagenta + 2);
    hNeNe->SetLineColor(kMagenta + 2);

    hOO->SetMarkerStyle(21);
    hOO->SetMarkerColor(kBlue + 1);
    hOO->SetLineColor(kBlue + 1);

    double maxVal = 0.0;
    for (int i = 1; i <= hNeNe->GetNbinsX(); ++i) {
        maxVal = std::max(maxVal, hNeNe->GetBinContent(i) + hNeNe->GetBinError(i));
    }
    for (int i = 1; i <= hOO->GetNbinsX(); ++i) {
        maxVal = std::max(maxVal, hOO->GetBinContent(i) + hOO->GetBinError(i));
    }
    if (maxVal <= 0.0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff2Sys", "3x2PC Eta-diff v2: O-O vs Ne-Ne", 900, 650);
    hNeNe->SetTitle("3x2PC #eta-differential v_{2}: Ne-Ne vs O-O (full-range combination); #eta; v_{2}");
    hNeNe->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
    hNeNe->Draw("E");
    hOO->Draw("E same");

    auto drawEdgePoints = [](TH1D* h, int color, int style) {
        if (!h || h->GetNbinsX() < 2) return;
        TGraphErrors* g = new TGraphErrors(2);
        g->SetPoint(0, 1.0, h->GetBinContent(1));
        g->SetPointError(0, 0.0, h->GetBinError(1));
        g->SetPoint(1, h->GetNbinsX(), h->GetBinContent(h->GetNbinsX()));
        g->SetPointError(1, 0.0, h->GetBinError(h->GetNbinsX()));
        g->SetMarkerStyle(style);
        g->SetMarkerSize(1.6);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->Draw("P same");
    };

    drawEdgePoints(hNeNe, kMagenta + 2, 29);
    drawEdgePoints(hOO,   kBlue + 1,    29);

    TLegend* leg = new TLegend(0.56, 0.76, 0.90, 0.90);
    leg->AddEntry(hNeNe, "Ne-Ne full range combination", "lep");
    leg->AddEntry(hOO,   "O-O full range combination",   "lep");
    leg->Draw();

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_OO_vs_NeNe.root");
    std::cout << "Saved eta-diff v2 comparison plot for O-O vs Ne-Ne." << std::endl;

    fNeNe->Close();
    fOO->Close();
    delete fNeNe;
    delete fOO;
    delete c;
}
