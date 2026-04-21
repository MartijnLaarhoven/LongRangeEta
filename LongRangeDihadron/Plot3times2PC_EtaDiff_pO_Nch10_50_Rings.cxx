#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct EtaPoint {
    double x;
    double y;
    double exLow;
    double exHigh;
    double eyLow;
    double eyHigh;
};

TH1D* LoadV2Combined(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

TFile* OpenFirstExistingFile(const std::vector<std::string>& candidates, std::string& selectedPath) {
    for (const auto& path : candidates) {
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (f && f->IsOpen() && !f->IsZombie()) {
            selectedPath = path;
            return f;
        }
        if (f) {
            f->Close();
            delete f;
        }
    }
    selectedPath.clear();
    return nullptr;
}

TH1D* LoadV2Side(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Side");
    if (!h) h = (TH1D*)f->Get("hV2_Sides");
    if (!h) h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

TH1D* BuildSidesFromFiles(TFile* fFT0C, TFile* fFT0A, const char* name) {
    TH1D* hFT0C = LoadV2Side(fFT0C);
    TH1D* hFT0A = LoadV2Side(fFT0A);
    if (!hFT0C || !hFT0A) return nullptr;

    TH1D* h = new TH1D(name, "", 2, 0.5, 2.5);
    h->SetDirectory(nullptr);
    h->SetBinContent(1, hFT0C->GetBinContent(1));
    h->SetBinError(1, hFT0C->GetBinError(1));
    h->SetBinContent(2, hFT0A->GetBinContent(1));
    h->SetBinError(2, hFT0A->GetBinError(1));
    return h;
}

TH1D* BuildEtaFrame(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
    h->SetDirectory(nullptr);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("v_{2}");
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.042);
    h->GetXaxis()->SetTitleOffset(1.00);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetTitleOffset(1.08);
    h->GetXaxis()->SetNdivisions(510);
    h->GetYaxis()->SetNdivisions(508);
    return h;
}

std::vector<EtaPoint> BuildEtaPoints(TFile* fFull, TFile* fInnerC, TFile* fInnerA, TFile* fOuterC, TFile* fOuterA) {
    std::vector<EtaPoint> points;
    TH1D* hFull = LoadV2Combined(fFull);
    TH1D* hInnerCSide = LoadV2Side(fInnerC);
    TH1D* hInnerASide = LoadV2Side(fInnerA);
    TH1D* hOuterCSide = LoadV2Side(fOuterC);
    TH1D* hOuterASide = LoadV2Side(fOuterA);
    if (!hFull || !hInnerCSide || !hInnerASide || !hOuterCSide || !hOuterASide) return points;

    auto pushPoint = [&](double xval,
                         double yval,
                         double xerrLow,
                         double xerrHigh,
                         double yerrLow,
                         double yerrHigh) {
        points.push_back({xval, yval, xerrLow, xerrHigh, yerrLow, yerrHigh});
    };

    const double ft0cInnerY = hInnerCSide->GetBinContent(1);
    const double ft0cInnerErr = hInnerCSide->GetBinError(1);
    const double ft0cOuterY = hOuterCSide->GetBinContent(1);
    const double ft0cOuterErr = hOuterCSide->GetBinError(1);
    const double ft0aOuterY = hOuterASide->GetBinContent(1);
    const double ft0aOuterErr = hOuterASide->GetBinError(1);
    const double ft0aInnerY = hInnerASide->GetBinContent(1);
    const double ft0aInnerErr = hInnerASide->GetBinError(1);

    pushPoint(-3.0, ft0cInnerY, 0.3, 0.3, ft0cInnerErr, ft0cInnerErr);
    pushPoint(-2.4, ft0cOuterY, 0.3, 0.3, ft0cOuterErr, ft0cOuterErr);
    pushPoint(3.85, ft0aInnerY, 0.35, 0.35, ft0aInnerErr, ft0aInnerErr);
    pushPoint(4.55, ft0aOuterY, 0.35, 0.35, ft0aOuterErr, ft0aOuterErr);

    for (int itpc = 0; itpc < 16; ++itpc) {
        const int sourceBin = itpc + 2;
        pushPoint(-0.75 + 0.1 * itpc,
                  hFull->GetBinContent(sourceBin),
                  0.05,
                  0.05,
                  hFull->GetBinError(sourceBin),
                  hFull->GetBinError(sourceBin));
    }

    return points;
}

TGraphAsymmErrors* BuildGraph(const char* name, const std::vector<EtaPoint>& points) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> exLow;
    std::vector<double> exHigh;
    std::vector<double> eyLow;
    std::vector<double> eyHigh;

    for (const auto& p : points) {
        x.push_back(p.x);
        y.push_back(p.y);
        exLow.push_back(p.exLow);
        exHigh.push_back(p.exHigh);
        eyLow.push_back(p.eyLow);
        eyHigh.push_back(p.eyHigh);
    }

    if (x.empty()) return nullptr;

    TGraphAsymmErrors* g = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), y.data(), exLow.data(), exHigh.data(), eyLow.data(), eyHigh.data());
    g->SetName(name);
    return g;
}

void StyleGraph(TGraphAsymmErrors* g, int color, int markerStyle) {
    if (!g) return;
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(1.35);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

double GetGraphMax(TGraphAsymmErrors* g) {
    if (!g) return 0.0;
    double maxVal = 0.0;
    for (int i = 0; i < g->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        g->GetPoint(i, x, y);
        maxVal = std::max(maxVal, y + g->GetErrorYhigh(i));
    }
    return maxVal;
}

double GetGraphMin(TGraphAsymmErrors* g) {
    if (!g) return 0.0;
    double minVal = 1e9;
    for (int i = 0; i < g->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        g->GetPoint(i, x, y);
        minVal = std::min(minVal, y - g->GetErrorYlow(i));
    }
    if (!std::isfinite(minVal) || minVal > 1e8) minVal = 0.0;
    return minVal;
}
}

void Plot3times2PC_EtaDiff_pO_Nch10_50_Rings() {
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    std::string pathFull, pathInnerC, pathInnerA, pathOuterC, pathOuterA;
    TFile* fFull = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_650299_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_650299_Mult_10_50.root"
    }, pathFull);
    TFile* fInnerC = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0C_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0C_Mult_10_50.root"
    }, pathInnerC);
    TFile* fInnerA = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0A_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0A_Mult_10_50.root"
    }, pathInnerA);
    TFile* fOuterC = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0C_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0C_Mult_10_50.root"
    }, pathOuterC);
    TFile* fOuterA = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0A_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0A_Mult_10_50.root"
    }, pathOuterA);

    if (!fFull) {
        std::cerr << "Cannot open p-O Nch full file (tried Nch/Mult candidates)." << std::endl;
        return;
    }

    std::cout << "[Plot p-O Nch rings] Full file: " << pathFull << std::endl;
    std::cout << "[Plot p-O Nch rings] Inner FT0C file: " << (pathInnerC.empty() ? "<missing>" : pathInnerC) << std::endl;
    std::cout << "[Plot p-O Nch rings] Inner FT0A file: " << (pathInnerA.empty() ? "<missing>" : pathInnerA) << std::endl;
    std::cout << "[Plot p-O Nch rings] Outer FT0C file: " << (pathOuterC.empty() ? "<missing>" : pathOuterC) << std::endl;
    std::cout << "[Plot p-O Nch rings] Outer FT0A file: " << (pathOuterA.empty() ? "<missing>" : pathOuterA) << std::endl;

    if (!fInnerC || !fInnerA || !fOuterC || !fOuterA) {
        std::cerr << "Missing p-O ring files for ring-combination plot." << std::endl;
        if (fFull) { fFull->Close(); delete fFull; }
        if (fInnerC) { fInnerC->Close(); delete fInnerC; }
        if (fInnerA) { fInnerA->Close(); delete fInnerA; }
        if (fOuterC) { fOuterC->Close(); delete fOuterC; }
        if (fOuterA) { fOuterA->Close(); delete fOuterA; }
        return;
    }

    TH1D* hFull = LoadV2Combined(fFull);
    if (!hFull) {
        std::cerr << "Missing full v2 histogram in p-O Nch file." << std::endl;
        if (fFull) { fFull->Close(); delete fFull; }
        if (fInnerC) { fInnerC->Close(); delete fInnerC; }
        if (fInnerA) { fInnerA->Close(); delete fInnerA; }
        if (fOuterC) { fOuterC->Close(); delete fOuterC; }
        if (fOuterA) { fOuterA->Close(); delete fOuterA; }
        return;
    }

    std::vector<EtaPoint> points = BuildEtaPoints(fFull, fInnerC, fInnerA, fOuterC, fOuterA);
    TGraphAsymmErrors* g = BuildGraph("gpONchRings", points);
    if (!g) {
        std::cerr << "Failed to build p-O Nch ring-combination graph." << std::endl;
        if (fFull) { fFull->Close(); delete fFull; }
        if (fInnerC) { fInnerC->Close(); delete fInnerC; }
        if (fInnerA) { fInnerA->Close(); delete fInnerA; }
        if (fOuterC) { fOuterC->Close(); delete fOuterC; }
        if (fOuterA) { fOuterA->Close(); delete fOuterA; }
        return;
    }

    StyleGraph(g, kBlue + 1, 24);

    TH1D* hFrame = BuildEtaFrame("hFramepONchRings");
    double maxVal = GetGraphMax(g);
    if (maxVal <= 0.0) maxVal = 0.08;
    double minVal = GetGraphMin(g);
    if (minVal > 0.0) minVal = 0.0;
    hFrame->GetYaxis()->SetRangeUser(minVal - 0.02 * maxVal, 1.30 * maxVal);

    TCanvas* c = new TCanvas("cpONchRings", "p-O 3times2PC Nch 10-50 with rings", 1150, 700);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.12);
    c->SetTicks(1, 1);

    hFrame->Draw("AXIS");
    g->Draw("PZ same");

    TLegend* leg = new TLegend(0.74, 0.80, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(g, "p-O (Nch 10-50)", "lep");
    leg->Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextColor(kBlack);
    lat.SetTextSize(0.042);
    lat.DrawLatex(0.13, 0.88, "ALICE Preliminary");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.13, 0.825, "3x2PC, N_{ch} 10-50");
    lat.DrawLatex(0.13, 0.775, "FT0 rings + TPC acceptance");

    c->SaveAs("./3times2PC/Plot_v2_FT0Side_pO_Nch10_50_inner_outer_full.root");
    c->SaveAs("./3times2PC/Plot_v2_FT0Side_pO_Nch10_50_inner_outer_full.png");

    if (fFull) { fFull->Close(); delete fFull; }
    if (fInnerC) { fInnerC->Close(); delete fInnerC; }
    if (fInnerA) { fInnerA->Close(); delete fInnerA; }
    if (fOuterC) { fOuterC->Close(); delete fOuterC; }
    if (fOuterA) { fOuterA->Close(); delete fOuterA; }
    delete c;
}