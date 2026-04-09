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
#include <vector>

namespace {
struct EtaPoint {
    double x;
    double y;
    double exLow;
    double exHigh;
    double eyLow;
    double eyHigh;
    bool mirrored;
};

TH1D* LoadV2Combined(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

TH1D* BuildEtaFrame(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
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

std::vector<EtaPoint> BuildEtaPoints(TH1D* h) {
    std::vector<EtaPoint> points;
    if (!h) return points;

    auto pushPoint = [&](double xval,
                         double yval,
                         double xerrLow,
                         double xerrHigh,
                         double yerrLow,
                         double yerrHigh,
                         bool mirrored) {
        points.push_back({xval, yval, xerrLow, xerrHigh, yerrLow, yerrHigh, mirrored});
    };

    // FT0 detector-side points; mirrored copies are displayed with open markers.
    const double ft0cY = h->GetBinContent(1);
    const double ft0cErr = h->GetBinError(1);
    const double ft0aY = h->GetBinContent(h->GetNbinsX());
    const double ft0aErr = h->GetBinError(h->GetNbinsX());

    pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, false);
    pushPoint(2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, true);
    pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, false);
    pushPoint(-4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, true);

    for (int itpc = 0; itpc < 16; ++itpc) {
        const int sourceBin = itpc + 2;
        pushPoint(-0.75 + 0.1 * itpc,
                  h->GetBinContent(sourceBin),
                  0.05,
                  0.05,
                  h->GetBinError(sourceBin),
                  h->GetBinError(sourceBin),
                  false);
    }

    return points;
}

TGraphAsymmErrors* BuildGraph(const char* name, const std::vector<EtaPoint>& points, bool mirroredOnly) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> exLow;
    std::vector<double> exHigh;
    std::vector<double> eyLow;
    std::vector<double> eyHigh;

    for (const auto& p : points) {
        if (mirroredOnly && !p.mirrored) continue;
        if (!mirroredOnly && p.mirrored) continue;
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

void Plot3times2PC_EtaDiff_NeNe_Nch10_50() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const char* inputFile = "./3times2PC/Vn_LHC25af_pass2_650641_nch10_50_Mult_10_50.root";

    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || !f->IsOpen()) {
        std::cerr << "Cannot open input file: " << inputFile << std::endl;
        return;
    }

    TH1D* hV2 = LoadV2Combined(f);
    if (!hV2) {
        std::cerr << "Missing hV2_Combined/hV2 in: " << inputFile << std::endl;
        f->Close();
        delete f;
        return;
    }

    std::vector<EtaPoint> points = BuildEtaPoints(hV2);
    TGraphAsymmErrors* gFilled = BuildGraph("gNeNeNchFilled", points, false);
    TGraphAsymmErrors* gMirror = BuildGraph("gNeNeNchMirror", points, true);

    if (!gFilled) {
        std::cerr << "Failed to build Ne-Ne Nch point graph." << std::endl;
        f->Close();
        delete f;
        return;
    }

    StyleGraph(gFilled, kGreen + 2, 20);
    StyleGraph(gMirror, kGreen + 2, 24);

    TH1D* hFrame = BuildEtaFrame("hFrameNeNeNch");

    double maxVal = std::max(GetGraphMax(gFilled), GetGraphMax(gMirror));
    if (maxVal <= 0.0) maxVal = 0.08;

    double minVal = std::min(GetGraphMin(gFilled), GetGraphMin(gMirror));
    if (minVal > 0.0) minVal = 0.0;

    hFrame->GetYaxis()->SetRangeUser(minVal - 0.02 * maxVal, 1.30 * maxVal);

    TCanvas* c = new TCanvas("cNeNeNch10_50", "Ne-Ne 3times2PC Nch 10-50", 980, 760);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.07);

    hFrame->Draw();
    gFilled->Draw("PZ same");
    if (gMirror) gMirror->Draw("PZ same");

    TLegend* leg = new TLegend(0.64, 0.80, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);
    leg->AddEntry(gFilled, "Ne-Ne (Nch 10-50)", "lep");
    leg->Draw();

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.14, 0.93, "ALICE Preliminary");
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.14, 0.875, "Ne-Ne #sqrt{s_{NN}} = 5.36 TeV, N_{ch}: 10-50");
    lat.DrawLatex(0.14, 0.83, "3times2PC, full-range combination");

    c->SaveAs("./3times2PC/Plot_v2_FT0Side_NeNe_Nch10_50.root");
    c->SaveAs("./3times2PC/Plot_v2_FT0Side_NeNe_Nch10_50.png");

    f->Close();
    delete f;
    delete c;
}
