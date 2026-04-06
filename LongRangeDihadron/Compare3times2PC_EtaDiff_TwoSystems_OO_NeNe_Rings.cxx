#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace {
TH1D* LoadFullV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

TH1D* LoadRingV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Sides");
    if (!h) h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

void Style(TH1D* h, int color, int marker) {
    if (!h) return;
    h->SetMarkerStyle(marker);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
}

void Style(TGraphErrors* g, int color, int marker) {
    if (!g) return;
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
}

TH1D* BuildEtaScaleHistogram(const char* name, const char* title) {
    TH1D* h = new TH1D(name, title, 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("v_{2}");
    h->SetStats(0);
    h->GetXaxis()->SetNdivisions(510);
    return h;
}

TGraphErrors* BuildEtaGraph(const char* name, TH1D* hFull, TH1D* hInnerSides, TH1D* hOuterSides, bool useRings) {
    if (!hFull) return nullptr;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ex;
    std::vector<double> ey;

    auto pushPoint = [&](double xval, double yval, double yerr) {
        x.push_back(xval);
        y.push_back(yval);
        ex.push_back(0.0);
        ey.push_back(yerr);
    };

    if (useRings && hInnerSides && hOuterSides) {
        const double ft0cInnerY = hInnerSides->GetBinContent(1);
        const double ft0cInnerErr = hInnerSides->GetBinError(1);
        const double ft0cOuterY = hOuterSides->GetBinContent(1);
        const double ft0cOuterErr = hOuterSides->GetBinError(1);
        const double ft0aOuterY = hOuterSides->GetBinContent(2);
        const double ft0aOuterErr = hOuterSides->GetBinError(2);
        const double ft0aInnerY = hInnerSides->GetBinContent(2);
        const double ft0aInnerErr = hInnerSides->GetBinError(2);

        pushPoint(-3.0, ft0cInnerY, ft0cInnerErr);
        pushPoint(-2.4, ft0cOuterY, ft0cOuterErr);
        pushPoint(2.4, ft0cOuterY, ft0cOuterErr);
        pushPoint(3.0, ft0cInnerY, ft0cInnerErr);
        pushPoint(-4.55, ft0aInnerY, ft0aInnerErr);
        pushPoint(-3.85, ft0aOuterY, ft0aOuterErr);
        pushPoint(3.85, ft0aOuterY, ft0aOuterErr);
        pushPoint(4.55, ft0aInnerY, ft0aInnerErr);
    } else {
        const double ft0cY = hFull->GetBinContent(1);
        const double ft0cErr = hFull->GetBinError(1);
        const double ft0aY = hFull->GetBinContent(hFull->GetNbinsX());
        const double ft0aErr = hFull->GetBinError(hFull->GetNbinsX());
        pushPoint(-2.4, ft0cY, ft0cErr);
        pushPoint(2.4, ft0cY, ft0cErr);
        pushPoint(-4.55, ft0aY, ft0aErr);
        pushPoint(4.55, ft0aY, ft0aErr);
    }

    for (int itpc = 0; itpc < 16; ++itpc) {
        const int sourceBin = itpc + 2;
        pushPoint(-0.75 + 0.1 * itpc, hFull->GetBinContent(sourceBin), hFull->GetBinError(sourceBin));
    }

    TGraphErrors* graph = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    graph->SetName(name);
    return graph;
}
}

void Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe_Rings() {
    const char* fileOOFull      = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
    const char* fileNeNeFull    = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileNeNeInner   = "./3times2PC/Vn_LHC25af_pass2_innerRing_Cent_0_20.root";
    const char* fileNeNeOuter   = "./3times2PC/Vn_LHC25af_pass2_outerRing_Cent_0_20.root";

    TFile* fOOFull    = TFile::Open(fileOOFull, "READ");
    TFile* fNeNeFull  = TFile::Open(fileNeNeFull, "READ");
    TFile* fNeNeInner = TFile::Open(fileNeNeInner, "READ");
    TFile* fNeNeOuter = TFile::Open(fileNeNeOuter, "READ");

    if (!fOOFull || !fOOFull->IsOpen()) {
        std::cerr << "Cannot open O-O full-range file: " << fileOOFull << std::endl;
        return;
    }
    if (!fNeNeFull || !fNeNeFull->IsOpen()) {
        std::cerr << "Cannot open Ne-Ne full-range file: " << fileNeNeFull << std::endl;
        fOOFull->Close();
        delete fOOFull;
        return;
    }
    if (!fNeNeInner || !fNeNeInner->IsOpen()) {
        std::cerr << "Warning: Cannot open Ne-Ne inner-ring file: " << fileNeNeInner << std::endl;
    }
    if (!fNeNeOuter || !fNeNeOuter->IsOpen()) {
        std::cerr << "Warning: Cannot open Ne-Ne outer-ring file: " << fileNeNeOuter << std::endl;
    }

    TH1D* hOOFull    = LoadFullV2(fOOFull);
    TH1D* hNeNeFull  = LoadFullV2(fNeNeFull);
    TH1D* hNeNeInner = LoadRingV2(fNeNeInner);
    TH1D* hNeNeOuter = LoadRingV2(fNeNeOuter);

    if (!hOOFull || !hNeNeFull) {
        std::cerr << "Missing hV2_Combined/hV2 in O-O or Ne-Ne full-range file." << std::endl;
        fOOFull->Close();
        fNeNeFull->Close();
        if (fNeNeInner) fNeNeInner->Close();
        if (fNeNeOuter) fNeNeOuter->Close();
        delete fOOFull;
        delete fNeNeFull;
        if (fNeNeInner) delete fNeNeInner;
        if (fNeNeOuter) delete fNeNeOuter;
        return;
    }

    TH1D* hFrame = BuildEtaScaleHistogram("hFrame", "3x2PC v_{2}: FT0 ring ranges on a physical #eta scale;#eta;v_{2}");
    TGraphErrors* gNeNe = BuildEtaGraph("gNeNe", hNeNeFull, hNeNeInner, hNeNeOuter, true);
    TGraphErrors* gOO = BuildEtaGraph("gOO", hOOFull, nullptr, nullptr, false);
    if (!hFrame || !gNeNe || !gOO) {
        std::cerr << "Failed to build eta-scale frame or graphs." << std::endl;
        fOOFull->Close();
        fNeNeFull->Close();
        if (fNeNeInner) fNeNeInner->Close();
        if (fNeNeOuter) fNeNeOuter->Close();
        delete fOOFull;
        delete fNeNeFull;
        if (fNeNeInner) delete fNeNeInner;
        if (fNeNeOuter) delete fNeNeOuter;
        return;
    }

    Style(gOO,    kBlue + 1, 21);
    Style(gNeNe,  kMagenta + 2, 20);

    double maxVal = 0.0;
    TGraphErrors* gs[] = {gOO, gNeNe};
    for (TGraphErrors* g : gs) {
        if (!g) continue;
        for (int i = 0; i < g->GetN(); ++i) {
            double xval = 0.0;
            double yval = 0.0;
            g->GetPoint(i, xval, yval);
            maxVal = std::max(maxVal, yval + g->GetErrorY(i));
        }
    }
    if (maxVal <= 0.0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff2SysRings", "3x2PC with explicit FT0 eta scale", 1150, 700);
    hFrame->SetTitle("3x2PC v_{2}: FT0 ring ranges on a physical #eta scale;#eta;v_{2}");
    hFrame->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
    hFrame->Draw("AXIS");
    gNeNe->Draw("P same");
    gOO->Draw("P same");

    TLine* sepLeft = new TLine(-2.1, 0.0, -2.1, maxVal * 1.22);
    TLine* sepTPCLeft = new TLine(-0.8, 0.0, -0.8, maxVal * 1.22);
    TLine* sepTPCRight = new TLine(0.8, 0.0, 0.8, maxVal * 1.22);
    TLine* sepRight = new TLine(2.1, 0.0, 2.1, maxVal * 1.22);
    sepLeft->SetLineStyle(2);
    sepTPCLeft->SetLineStyle(2);
    sepTPCRight->SetLineStyle(2);
    sepRight->SetLineStyle(2);
    sepLeft->SetLineColor(kGray + 2);
    sepTPCLeft->SetLineColor(kGray + 2);
    sepTPCRight->SetLineColor(kGray + 2);
    sepRight->SetLineColor(kGray + 2);
    sepLeft->Draw("same");
    sepTPCLeft->Draw("same");
    sepTPCRight->Draw("same");
    sepRight->Draw("same");

    TLegend* leg = new TLegend(0.46, 0.64, 0.90, 0.90);
    leg->AddEntry(gOO,    "O-O (full FT0 baseline)",    "lep");
    leg->AddEntry(gNeNe,  "Ne-Ne (FT0 ring ranges)",  "lep");
    leg->Draw();

    c->SaveAs("./3times2PC/Compare_v2_FT0Side_OO_vs_NeNe_inner_outer_full.root");
    std::cout << "Saved comparison plot with full-range TPC bins and FT0 ring points for O-O vs Ne-Ne." << std::endl;

    fOOFull->Close();
    fNeNeFull->Close();
    if (fNeNeInner) fNeNeInner->Close();
    if (fNeNeOuter) fNeNeOuter->Close();
    delete fOOFull;
    delete fNeNeFull;
    if (fNeNeInner) delete fNeNeInner;
    if (fNeNeOuter) delete fNeNeOuter;
}
