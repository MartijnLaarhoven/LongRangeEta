#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TLine.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>

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

TH1D* LoadFullV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

void Style(TGraphAsymmErrors* g, int color, int marker) {
    if (!g) return;
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.35);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
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

TH1D* BuildRatioFrame(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("Ne-Ne / O-O");
    h->GetXaxis()->SetTitleSize(0.11);
    h->GetXaxis()->SetLabelSize(0.10);
    h->GetXaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetTitleSize(0.10);
    h->GetYaxis()->SetLabelSize(0.09);
    h->GetYaxis()->SetTitleOffset(0.48);
    h->GetXaxis()->SetNdivisions(510);
    h->GetYaxis()->SetNdivisions(505);
    return h;
}

std::vector<EtaPoint> BuildEtaPoints(TH1D* hFull, bool mirrorDetectors) {
    std::vector<EtaPoint> points;
    if (!hFull) return points;

    auto pushPoint = [&](double xval,
                         double yval,
                         double exLow,
                         double exHigh,
                         double eyLow,
                         double eyHigh,
                         bool mirrored) {
        points.push_back({xval, yval, exLow, exHigh, eyLow, eyHigh, mirrored});
    };

    const int nBins = hFull->GetNbinsX();
    if (nBins < 4) return points;

    const double ft0cY = hFull->GetBinContent(1);
    const double ft0cErr = hFull->GetBinError(1);
    const double ft0aY = hFull->GetBinContent(nBins);
    const double ft0aErr = hFull->GetBinError(nBins);

    if (mirrorDetectors) {
        // Physical-side points
        pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, false);
        pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, false);
        // Mirrored points drawn with open symbols
        pushPoint(2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, true);
        pushPoint(-4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, true);
    } else {
        // Asymmetric p-O: keep full-range detector points without mirroring.
        pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, false);
        pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, false);
    }

    for (int itpc = 0; itpc < 16; ++itpc) {
        const int sourceBin = itpc + 2;
        pushPoint(-0.75 + 0.1 * itpc,
                  hFull->GetBinContent(sourceBin),
                  0.05,
                  0.05,
                  hFull->GetBinError(sourceBin),
                  hFull->GetBinError(sourceBin),
                  false);
    }

    return points;
}

TGraphAsymmErrors* BuildGraphFromPoints(const char* name, const std::vector<EtaPoint>& points, bool includeMirrored, bool onlyMirrored) {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> exLow;
    std::vector<double> exHigh;
    std::vector<double> eyLow;
    std::vector<double> eyHigh;

    for (const auto& p : points) {
        if (!includeMirrored && p.mirrored) continue;
        if (onlyMirrored && !p.mirrored) continue;
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

std::vector<EtaPoint> BuildRatioPoints(const std::vector<EtaPoint>& numPoints, const std::vector<EtaPoint>& denPoints) {
    std::map<int, EtaPoint> numByX;
    std::map<int, EtaPoint> denByX;
    auto keyX = [](double x) { return static_cast<int>(std::lround(1000.0 * x)); };

    for (const auto& p : numPoints) numByX[keyX(p.x)] = p;
    for (const auto& p : denPoints) denByX[keyX(p.x)] = p;

    std::vector<EtaPoint> ratio;
    for (const auto& it : numByX) {
        if (denByX.find(it.first) == denByX.end()) continue;
        const EtaPoint& a = it.second;
        const EtaPoint& b = denByX[it.first];
        if (!std::isfinite(a.y) || !std::isfinite(b.y) || a.y <= 0.0 || b.y <= 0.0) continue;

        const double r = a.y / b.y;
        const double ea = std::max(a.eyLow, a.eyHigh);
        const double eb = std::max(b.eyLow, b.eyHigh);
        const double er = r * std::sqrt((ea / a.y) * (ea / a.y) + (eb / b.y) * (eb / b.y));
        ratio.push_back({a.x,
                         r,
                         std::max(a.exLow, b.exLow),
                         std::max(a.exHigh, b.exHigh),
                         er,
                         er,
                         a.mirrored || b.mirrored});
    }
    return ratio;
}

double GetGraphYMax(TGraphAsymmErrors* g) {
    if (!g) return 0.0;
    double maxVal = 0.0;
    for (int i = 0; i < g->GetN(); ++i) {
        double xval = 0.0;
        double yval = 0.0;
        g->GetPoint(i, xval, yval);
        maxVal = std::max(maxVal, yval + g->GetErrorYhigh(i));
    }
    return maxVal;
}

double GetGraphYMin(TGraphAsymmErrors* g) {
    if (!g) return 1.0;
    double minVal = 1e9;
    for (int i = 0; i < g->GetN(); ++i) {
        double xval = 0.0;
        double yval = 0.0;
        g->GetPoint(i, xval, yval);
        minVal = std::min(minVal, yval - g->GetErrorYlow(i));
    }
    if (!std::isfinite(minVal) || minVal > 1e8) minVal = 1.0;
    return minVal;
}
}

void Compare3times2PC_EtaDiff_ThreeSystems() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const char* fileNeNe = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileOO   = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
    const char* filePO   = "./3times2PC/Vn_LHC25ad_pass2_644389_Cent_0_20.root";

    TFile* fNeNe = TFile::Open(fileNeNe, "READ");
    TFile* fOO   = TFile::Open(fileOO,   "READ");
    TFile* fPO   = TFile::Open(filePO,   "READ");

    if (!fNeNe || !fNeNe->IsOpen()) { std::cerr << "Cannot open Ne-Ne file!" << std::endl; return; }
    if (!fOO   || !fOO->IsOpen())   { std::cerr << "Cannot open O-O file!" << std::endl; return; }
    if (!fPO   || !fPO->IsOpen())   { std::cerr << "Warning: cannot open p-O file, skip p-O." << std::endl; }

    TH1D* hNeNe = LoadFullV2(fNeNe);
    TH1D* hOO   = LoadFullV2(fOO);
    TH1D* hPO   = LoadFullV2(fPO);
    if (!hNeNe || !hOO) {
        std::cerr << "Missing mandatory Ne-Ne or O-O v2 histogram." << std::endl;
        return;
    }

    std::vector<EtaPoint> nePoints = BuildEtaPoints(hNeNe, true);
    std::vector<EtaPoint> ooPoints = BuildEtaPoints(hOO, true);
    std::vector<EtaPoint> poPoints = BuildEtaPoints(hPO, false);

    TH1D* hFrame = BuildEtaFrame("hFrame3Sys");

    TGraphAsymmErrors* gNeFilled = BuildGraphFromPoints("gNeFilled", nePoints, false, false);
    TGraphAsymmErrors* gNeMirror = BuildGraphFromPoints("gNeMirror", nePoints, true, true);
    TGraphAsymmErrors* gOOFilled = BuildGraphFromPoints("gOOFilled", ooPoints, false, false);
    TGraphAsymmErrors* gOOMirror = BuildGraphFromPoints("gOOMirror", ooPoints, true, true);
    TGraphAsymmErrors* gPOFilled = BuildGraphFromPoints("gPOFilled", poPoints, true, false);

    if (!hFrame || !gNeFilled || !gOOFilled) {
        std::cerr << "Failed to build eta-scale frame or graphs." << std::endl;
        return;
    }

    Style(gOOFilled, kRed + 1, 21);
    Style(gOOMirror, kRed + 1, 25);
    Style(gNeFilled, kGreen + 2, 20);
    Style(gNeMirror, kGreen + 2, 24);
    Style(gPOFilled, kBlue + 1, 24);

    double maxVal = 0.0;
    maxVal = std::max(maxVal, GetGraphYMax(gNeFilled));
    maxVal = std::max(maxVal, GetGraphYMax(gNeMirror));
    maxVal = std::max(maxVal, GetGraphYMax(gOOFilled));
    maxVal = std::max(maxVal, GetGraphYMax(gOOMirror));
    maxVal = std::max(maxVal, GetGraphYMax(gPOFilled));
    if (maxVal <= 0.0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff3Sys", "3x2PC three systems", 1150, 700);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.07);
    c->SetBottomMargin(0.12);
    c->SetTicks(1, 1);

    hFrame->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
    hFrame->Draw("AXIS");
    gNeFilled->Draw("PZ same");
    if (gNeMirror) gNeMirror->Draw("PZ same");
    gOOFilled->Draw("PZ same");
    if (gOOMirror) gOOMirror->Draw("PZ same");
    if (gPOFilled) gPOFilled->Draw("PZ same");

    TLegend* leg = new TLegend(0.64, 0.70, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gOOFilled, "O-O", "lep");
    leg->AddEntry(gNeFilled, "Ne-Ne", "lep");
    if (gPOFilled) leg->AddEntry(gPOFilled, "p-O", "lep");
    if (gOOMirror || gNeMirror) leg->AddEntry(gOOMirror ? gOOMirror : gNeMirror, "Mirrored points", "p");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.13, 0.91, "ALICE Preliminary");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.13, 0.855, "3x2PC, 0-20%");
    latex.DrawLatex(0.13, 0.805, "Full FT0 ranges + TPC acceptance");

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO.root");
    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO.png");
    std::cout << "Saved eta-diff v2 comparison plot for available systems." << std::endl;

    if (fNeNe) fNeNe->Close();
    if (fOO) fOO->Close();
    if (fPO) fPO->Close();
    delete fNeNe;
    delete fOO;
    delete fPO;
    delete c;
}
