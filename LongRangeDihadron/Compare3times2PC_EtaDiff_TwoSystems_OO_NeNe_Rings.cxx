#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TLine.h>
#include <TStyle.h>
#include <cmath>
#include <algorithm>
#include <iostream>
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

TH1D* LoadRingV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Sides");
    if (!h) h = (TH1D*)f->Get("hV2_Combined");
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

TH1D* BuildEtaScaleHistogram(const char* name, const char* title) {
    TH1D* h = new TH1D(name, title, 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("v_{2}");
    h->SetStats(0);
    h->SetTitle("");
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

TH1D* BuildRatioHistogram(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("Ne-Ne / O-O");
    h->SetStats(0);
    h->SetTitle("");
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

std::vector<EtaPoint> BuildEtaPoints(TH1D* hFull, TH1D* hInnerSides, TH1D* hOuterSides, bool useRings) {
    std::vector<EtaPoint> points;
    if (!hFull) return points;

    auto pushPoint = [&](double xval,
                         double yval,
                         double xerrLow,
                         double xerrHigh,
                         double yerrLow,
                         double yerrHigh,
                         bool mirrored) {
        points.push_back({xval, yval, xerrLow, xerrHigh, yerrLow, yerrHigh, mirrored});
    };

    if (useRings && hInnerSides && hOuterSides) {
        const double fullFt0cY = hFull->GetBinContent(1);
        const double fullFt0cErr = hFull->GetBinError(1);
        const double fullFt0aY = hFull->GetBinContent(hFull->GetNbinsX());
        const double fullFt0aErr = hFull->GetBinError(hFull->GetNbinsX());

        auto sanitizeRingPoint = [&](double y, double err, double fallbackY, double fallbackErr, const char* label) {
            const bool invalid = (!std::isfinite(y) || !std::isfinite(err) || y < 0.0 || err >= 9.9);
            if (invalid) {
                std::cerr << "[Compare rings] Warning: invalid " << label
                          << " (y=" << y << ", err=" << err << "), fallback to full-range value" << std::endl;
                return std::make_pair(fallbackY, fallbackErr);
            }
            return std::make_pair(y, err);
        };

        auto ft0cInner = sanitizeRingPoint(hInnerSides->GetBinContent(1), hInnerSides->GetBinError(1), fullFt0cY, fullFt0cErr, "FT0C inner");
        auto ft0cOuter = sanitizeRingPoint(hOuterSides->GetBinContent(1), hOuterSides->GetBinError(1), fullFt0cY, fullFt0cErr, "FT0C outer");
        auto ft0aOuter = sanitizeRingPoint(hOuterSides->GetBinContent(2), hOuterSides->GetBinError(2), fullFt0aY, fullFt0aErr, "FT0A outer");
        auto ft0aInner = sanitizeRingPoint(hInnerSides->GetBinContent(2), hInnerSides->GetBinError(2), fullFt0aY, fullFt0aErr, "FT0A inner");

        // Physical-side points are filled; mirrored copies are open-symbol points.
        pushPoint(-3.0, ft0cInner.first, 0.3, 0.3, ft0cInner.second, ft0cInner.second, false);
        pushPoint(-2.4, ft0cOuter.first, 0.3, 0.3, ft0cOuter.second, ft0cOuter.second, false);
        pushPoint(2.4, ft0cOuter.first, 0.3, 0.3, ft0cOuter.second, ft0cOuter.second, true);
        pushPoint(3.0, ft0cInner.first, 0.3, 0.3, ft0cInner.second, ft0cInner.second, true);
        pushPoint(3.85, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, false);
        pushPoint(4.55, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, false);
        pushPoint(-4.55, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, true);
        pushPoint(-3.85, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, true);
    } else {
        const double ft0cY = hFull->GetBinContent(1);
        const double ft0cErr = hFull->GetBinError(1);
        const double ft0aY = hFull->GetBinContent(hFull->GetNbinsX());
        const double ft0aErr = hFull->GetBinError(hFull->GetNbinsX());
        pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, false);
        pushPoint(2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr, true);
        pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, false);
        pushPoint(-4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr, true);
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

    TGraphAsymmErrors* graph = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), y.data(), exLow.data(), exHigh.data(), eyLow.data(), eyHigh.data());
    graph->SetName(name);
    return graph;
}

std::vector<EtaPoint> BuildRatioPoints(const std::vector<EtaPoint>& numPoints, const std::vector<EtaPoint>& denPoints) {
    std::map<int, EtaPoint> numByX;
    std::map<int, EtaPoint> denByX;

    auto keyX = [](double x) { return static_cast<int>(std::lround(1000.0 * x)); };

    for (const auto& p : numPoints) {
        numByX[keyX(p.x)] = p;
    }
    for (const auto& p : denPoints) {
        denByX[keyX(p.x)] = p;
    }

    std::vector<EtaPoint> ratioPoints;
    for (const auto& it : numByX) {
        const int key = it.first;
        if (denByX.find(key) == denByX.end()) continue;

        const EtaPoint& a = it.second;
        const EtaPoint& b = denByX[key];
        if (!std::isfinite(a.y) || !std::isfinite(b.y) || a.y <= 0.0 || b.y <= 0.0) continue;

        const double r = a.y / b.y;
        const double errA = std::max(a.eyLow, a.eyHigh);
        const double errB = std::max(b.eyLow, b.eyHigh);
        const double relA = (a.y > 0.0) ? (errA / a.y) : 0.0;
        const double relB = (b.y > 0.0) ? (errB / b.y) : 0.0;
        const double errR = r * std::sqrt(relA * relA + relB * relB);

        ratioPoints.push_back({a.x,
                               r,
                               std::max(a.exLow, b.exLow),
                               std::max(a.exHigh, b.exHigh),
                               errR,
                               errR,
                               a.mirrored || b.mirrored});
    }

    return ratioPoints;
}

double GetGraphMax(TGraphAsymmErrors* g) {
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

double GetGraphMin(TGraphAsymmErrors* g) {
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

void Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe_Rings() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const char* fileOOFull      = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
    const char* fileOOInner     = "./3times2PC/Vn_LHC25ae_pass2_innerRing_Cent_0_20.root";
    const char* fileOOOuter     = "./3times2PC/Vn_LHC25ae_pass2_outerRing_Cent_0_20.root";
    const char* fileNeNeFull    = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileNeNeInner   = "./3times2PC/Vn_LHC25af_pass2_innerRing_Cent_0_20.root";
    const char* fileNeNeOuter   = "./3times2PC/Vn_LHC25af_pass2_outerRing_Cent_0_20.root";

    TFile* fOOFull    = TFile::Open(fileOOFull, "READ");
    TFile* fOOInner   = TFile::Open(fileOOInner, "READ");
    TFile* fOOOuter   = TFile::Open(fileOOOuter, "READ");
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
    if (!fOOInner || !fOOInner->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O inner-ring file: " << fileOOInner << std::endl;
    }
    if (!fOOOuter || !fOOOuter->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O outer-ring file: " << fileOOOuter << std::endl;
    }
    if (!fNeNeInner || !fNeNeInner->IsOpen()) {
        std::cerr << "Warning: Cannot open Ne-Ne inner-ring file: " << fileNeNeInner << std::endl;
    }
    if (!fNeNeOuter || !fNeNeOuter->IsOpen()) {
        std::cerr << "Warning: Cannot open Ne-Ne outer-ring file: " << fileNeNeOuter << std::endl;
    }

    TH1D* hOOFull    = LoadFullV2(fOOFull);
    TH1D* hOOInner   = LoadRingV2(fOOInner);
    TH1D* hOOOuter   = LoadRingV2(fOOOuter);
    TH1D* hNeNeFull  = LoadFullV2(fNeNeFull);
    TH1D* hNeNeInner = LoadRingV2(fNeNeInner);
    TH1D* hNeNeOuter = LoadRingV2(fNeNeOuter);

    if (!hOOFull || !hNeNeFull) {
        std::cerr << "Missing hV2_Combined/hV2 in O-O or Ne-Ne full-range file." << std::endl;
        fOOFull->Close();
        fNeNeFull->Close();
        if (fOOInner) fOOInner->Close();
        if (fOOOuter) fOOOuter->Close();
        if (fNeNeInner) fNeNeInner->Close();
        if (fNeNeOuter) fNeNeOuter->Close();
        delete fOOFull;
        delete fNeNeFull;
        if (fOOInner) delete fOOInner;
        if (fOOOuter) delete fOOOuter;
        if (fNeNeInner) delete fNeNeInner;
        if (fNeNeOuter) delete fNeNeOuter;
        return;
    }

    const bool ooHasRings = (hOOInner && hOOOuter);
    std::vector<EtaPoint> ooPoints = BuildEtaPoints(hOOFull, hOOInner, hOOOuter, ooHasRings);
    std::vector<EtaPoint> nePoints = BuildEtaPoints(hNeNeFull, hNeNeInner, hNeNeOuter, true);
    std::vector<EtaPoint> ratioPoints = BuildRatioPoints(nePoints, ooPoints);

    TH1D* hFrame = BuildEtaScaleHistogram("hFrame", "");
    TH1D* hRatioFrame = BuildRatioHistogram("hRatioFrame");

    TGraphAsymmErrors* gOOFilled = BuildGraphFromPoints("gOOFilled", ooPoints, false, false);
    TGraphAsymmErrors* gOOMirror = BuildGraphFromPoints("gOOMirror", ooPoints, true, true);
    TGraphAsymmErrors* gNeNeFilled = BuildGraphFromPoints("gNeNeFilled", nePoints, false, false);
    TGraphAsymmErrors* gNeNeMirror = BuildGraphFromPoints("gNeNeMirror", nePoints, true, true);
    TGraphAsymmErrors* gRatioFilled = BuildGraphFromPoints("gRatioFilled", ratioPoints, false, false);
    TGraphAsymmErrors* gRatioMirror = BuildGraphFromPoints("gRatioMirror", ratioPoints, true, true);

    if (!hFrame || !hRatioFrame || !gOOFilled || !gNeNeFilled || !gRatioFilled) {
        std::cerr << "Failed to build eta-scale frame/graphs or ratio graph." << std::endl;
        fOOFull->Close();
        fNeNeFull->Close();
        if (fOOInner) fOOInner->Close();
        if (fOOOuter) fOOOuter->Close();
        if (fNeNeInner) fNeNeInner->Close();
        if (fNeNeOuter) fNeNeOuter->Close();
        delete fOOFull;
        delete fNeNeFull;
        if (fOOInner) delete fOOInner;
        if (fOOOuter) delete fOOOuter;
        if (fNeNeInner) delete fNeNeInner;
        if (fNeNeOuter) delete fNeNeOuter;
        return;
    }

    Style(gOOFilled,    kRed + 1, 21);
    Style(gOOMirror,    kRed + 1, 25);
    Style(gNeNeFilled,  kGreen + 2, 20);
    Style(gNeNeMirror,  kGreen + 2, 24);
    Style(gRatioFilled, kBlack, 20);
    Style(gRatioMirror, kBlack, 24);

    double maxVal = 0.0;
    maxVal = std::max(maxVal, GetGraphMax(gOOFilled));
    maxVal = std::max(maxVal, GetGraphMax(gOOMirror));
    maxVal = std::max(maxVal, GetGraphMax(gNeNeFilled));
    maxVal = std::max(maxVal, GetGraphMax(gNeNeMirror));
    if (maxVal <= 0.0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff2SysRings", "3x2PC with explicit FT0 eta scale and ratio", 1150, 850);

    TPad* padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    padTop->SetLeftMargin(0.11);
    padTop->SetRightMargin(0.03);
    padTop->SetTopMargin(0.07);
    padTop->SetBottomMargin(0.02);
    padTop->SetTicks(1, 1);
    padTop->Draw();

    TPad* padBottom = new TPad("padBottom", "padBottom", 0.0, 0.00, 1.0, 0.30);
    padBottom->SetLeftMargin(0.11);
    padBottom->SetRightMargin(0.03);
    padBottom->SetTopMargin(0.02);
    padBottom->SetBottomMargin(0.34);
    padBottom->SetTicks(1, 1);
    padBottom->Draw();

    padTop->cd();
    hFrame->GetXaxis()->SetLabelSize(0.0);
    hFrame->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
    hFrame->Draw("AXIS");
    gNeNeFilled->Draw("PZ same");
    if (gNeNeMirror) gNeNeMirror->Draw("PZ same");
    gOOFilled->Draw("PZ same");
    if (gOOMirror) gOOMirror->Draw("PZ same");

    TLegend* leg = new TLegend(0.74, 0.78, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gOOFilled,    ooHasRings ? "O-O" : "O-O (full FT0 baseline)", "lep");
    leg->AddEntry(gNeNeFilled,  "Ne-Ne", "lep");
    if (gOOMirror) leg->AddEntry(gOOMirror, "Mirrored points", "p");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.13, 0.91, "ALICE Preliminary");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.13, 0.855, "3x2PC, 0-20%");
    latex.DrawLatex(0.13, 0.805, "FT0 rings + TPC acceptance (+ ratio)");

    padBottom->cd();
    double ratioMin = std::min(GetGraphMin(gRatioFilled), GetGraphMin(gRatioMirror));
    double ratioMax = std::max(GetGraphMax(gRatioFilled), GetGraphMax(gRatioMirror));
    if (!std::isfinite(ratioMin) || !std::isfinite(ratioMax) || ratioMax <= ratioMin) {
        ratioMin = 0.6;
        ratioMax = 1.4;
    }
    const double pad = 0.15 * (ratioMax - ratioMin);
    hRatioFrame->GetYaxis()->SetRangeUser(std::max(0.0, ratioMin - pad), ratioMax + pad);
    hRatioFrame->Draw("AXIS");
    TLine* unity = new TLine(-5.0, 1.0, 5.0, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kGray + 2);
    unity->Draw("same");
    gRatioFilled->Draw("PZ same");
    if (gRatioMirror) gRatioMirror->Draw("PZ same");

    c->SaveAs("./3times2PC/Compare_v2_FT0Side_OO_vs_NeNe_inner_outer_full.root");
    c->SaveAs("./3times2PC/Compare_v2_FT0Side_OO_vs_NeNe_inner_outer_full.png");

    auto DrawSingleSystem = [&](const char* canvasName,
                                const char* canvasTitle,
                                const char* outStem,
                                TGraphAsymmErrors* graphFilled,
                                TGraphAsymmErrors* graphMirror,
                                const char* legendLabel,
                                const char* systemLabel) {
        TCanvas* cSingle = new TCanvas(canvasName, canvasTitle, 1150, 700);
        cSingle->SetLeftMargin(0.11);
        cSingle->SetRightMargin(0.03);
        cSingle->SetTopMargin(0.06);
        cSingle->SetBottomMargin(0.12);
        cSingle->SetTicks(1, 1);

        TH1D* hFrameSingle = BuildEtaScaleHistogram(Form("hFrame_%s", canvasName), "");
        hFrameSingle->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
        hFrameSingle->Draw("AXIS");
        graphFilled->Draw("PZ same");
        if (graphMirror) graphMirror->Draw("PZ same");

        TLegend* legSingle = new TLegend(0.71, 0.79, 0.91, 0.90);
        legSingle->SetBorderSize(0);
        legSingle->SetFillStyle(0);
        legSingle->SetTextFont(42);
        legSingle->SetTextSize(0.030);
        legSingle->AddEntry(graphFilled, legendLabel, "lep");
        if (graphMirror) legSingle->AddEntry(graphMirror, "Mirrored points", "p");
        legSingle->Draw();

        TLatex latexSingle;
        latexSingle.SetNDC();
        latexSingle.SetTextFont(42);
        latexSingle.SetTextColor(kBlack);
        latexSingle.SetTextSize(0.042);
        latexSingle.DrawLatex(0.13, 0.91, "ALICE Preliminary");
        latexSingle.SetTextSize(0.036);
        latexSingle.DrawLatex(0.13, 0.855, "3x2PC, 0-20%");
        latexSingle.DrawLatex(0.13, 0.805, systemLabel);

        cSingle->SaveAs(Form("./3times2PC/%s.root", outStem));
        cSingle->SaveAs(Form("./3times2PC/%s.png", outStem));
    };

    DrawSingleSystem("cEtaDiff_OOOnly",
                     "3x2PC O-O only",
                     "Compare_v2_FT0Side_OO_only_inner_outer_full",
                     gOOFilled,
                     gOOMirror,
                     ooHasRings ? "O-O" : "O-O (full FT0 baseline)",
                     "O-O, FT0 rings + TPC acceptance");

    DrawSingleSystem("cEtaDiff_NeNeOnly",
                     "3x2PC Ne-Ne only",
                     "Compare_v2_FT0Side_NeNe_only_inner_outer_full",
                     gNeNeFilled,
                     gNeNeMirror,
                     "Ne-Ne",
                     "Ne-Ne, FT0 rings + TPC acceptance");

    if (!ooHasRings) {
        std::cout << "Saved comparison plot. Note: O-O ring files missing, so O-O uses full-range FT0 points." << std::endl;
    } else {
        std::cout << "Saved comparison plot with FT0 ring points for both O-O and Ne-Ne." << std::endl;
    }
    std::cout << "Saved individual O-O and Ne-Ne plots in the same style." << std::endl;

    fOOFull->Close();
    fNeNeFull->Close();
    if (fOOInner) fOOInner->Close();
    if (fOOOuter) fOOOuter->Close();
    if (fNeNeInner) fNeNeInner->Close();
    if (fNeNeOuter) fNeNeOuter->Close();
    delete fOOFull;
    delete fNeNeFull;
    if (fOOInner) delete fOOInner;
    if (fOOOuter) delete fOOOuter;
    if (fNeNeInner) delete fNeNeInner;
    if (fNeNeOuter) delete fNeNeOuter;
}
