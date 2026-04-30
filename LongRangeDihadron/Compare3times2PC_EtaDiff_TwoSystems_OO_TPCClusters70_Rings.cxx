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

TH1D* LoadFullVn(TFile* f, int harmonic) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get(Form("hV%d_Combined", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d", harmonic));
    return h;
}

TH1D* LoadRingVn(TFile* f, int harmonic) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get(Form("hV%d_Side", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d_Sides", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d_Combined", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d", harmonic));
    return h;
}

TH1D* LoadSideValue(TFile* f, int harmonic) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get(Form("hV%d_Side", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d_Sides", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d_Combined", harmonic));
    if (!h) h = (TH1D*)f->Get(Form("hV%d", harmonic));
    return h;
}

TH1D* BuildRingSidesFromSideFiles(TFile* fFT0C, TFile* fFT0A, int harmonic, const char* name, const char* title) {
    TH1D* hFT0C = LoadSideValue(fFT0C, harmonic);
    TH1D* hFT0A = LoadSideValue(fFT0A, harmonic);
    if (!hFT0C || !hFT0A) return nullptr;

    TH1D* h = new TH1D(name, title, 2, 0.5, 2.5);
    h->GetXaxis()->SetBinLabel(1, "FT0C [-3.3,-2.1]");
    h->GetXaxis()->SetBinLabel(2, "FT0A [3.5,4.9]");
    h->SetBinContent(1, hFT0C->GetBinContent(1));
    h->SetBinError(1, hFT0C->GetBinError(1));
    h->SetBinContent(2, hFT0A->GetBinContent(1));
    h->SetBinError(2, hFT0A->GetBinError(1));
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

TH1D* BuildEtaScaleHistogram(const char* name, const char* title, int harmonic) {
    TH1D* h = new TH1D(name, title, 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle(Form("v_{%d}", harmonic));
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
    h->GetYaxis()->SetTitle("O-O TPCClusters70 / O-O");
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

TH1D* BuildRingToFullHistogram(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("ring / full FT0");
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

TH1D* BuildRingToFullDeltaPercentHistogram(const char* name) {
    TH1D* h = new TH1D(name, "", 20, -5.0, 5.0);
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("100#times(ring/full - 1) [%]");
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

std::vector<EtaPoint> BuildEtaPoints(TH1D* hFull,
                                     TH1D* hInnerSides,
                                     TH1D* hOuterSides,
                                     bool useRings,
                                     bool ft0aInnerAtLargeAbsEta) {
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
        if (ft0aInnerAtLargeAbsEta) {
            pushPoint(4.55, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, false);
            pushPoint(3.85, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, false);
            pushPoint(-3.85, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, true);
            pushPoint(-4.55, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, true);
        } else {
            pushPoint(3.85, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, false);
            pushPoint(4.55, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, false);
            pushPoint(-4.55, ft0aOuter.first, 0.35, 0.35, ft0aOuter.second, ft0aOuter.second, true);
            pushPoint(-3.85, ft0aInner.first, 0.35, 0.35, ft0aInner.second, ft0aInner.second, true);
        }
    } else {
        const double ft0cY = hFull->GetBinContent(1);
        const double ft0cErr = hFull->GetBinError(1);
        const double ft0aY = hFull->GetBinContent(hFull->GetNbinsX());
        const double ft0aErr = hFull->GetBinError(hFull->GetNbinsX());
        pushPoint(-2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr, false);
        pushPoint(2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr, true);
        pushPoint(4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr, false);
        pushPoint(-4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr, true);
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

std::vector<EtaPoint> BuildRingToFullPoints(TH1D* hFull,
                                            TH1D* hInnerSides,
                                            TH1D* hOuterSides,
                                            bool ft0aInnerAtLargeAbsEta) {
    std::vector<EtaPoint> points;
    if (!hFull || !hInnerSides || !hOuterSides) return points;

    const double fullFt0cY = hFull->GetBinContent(1);
    const double fullFt0cErr = hFull->GetBinError(1);
    const double fullFt0aY = hFull->GetBinContent(hFull->GetNbinsX());
    const double fullFt0aErr = hFull->GetBinError(hFull->GetNbinsX());

    auto pushRatioPoint = [&](double x,
                              double ringY,
                              double ringErr,
                              double fullY,
                              double fullErr,
                              bool mirrored,
                              double ex) {
        if (!std::isfinite(ringY) || !std::isfinite(ringErr) || !std::isfinite(fullY) || !std::isfinite(fullErr)) return;
        if (ringY <= 0.0 || fullY <= 0.0) return;
        const double r = ringY / fullY;
        const double relRing = ringErr / ringY;
        const double relFull = fullErr / fullY;
        const double errR = r * std::sqrt(relRing * relRing + relFull * relFull);
        points.push_back({x, r, ex, ex, errR, errR, mirrored});
    };

    const double ft0cInnerY = hInnerSides->GetBinContent(1);
    const double ft0cInnerErr = hInnerSides->GetBinError(1);
    const double ft0cOuterY = hOuterSides->GetBinContent(1);
    const double ft0cOuterErr = hOuterSides->GetBinError(1);
    const double ft0aOuterY = hOuterSides->GetBinContent(2);
    const double ft0aOuterErr = hOuterSides->GetBinError(2);
    const double ft0aInnerY = hInnerSides->GetBinContent(2);
    const double ft0aInnerErr = hInnerSides->GetBinError(2);

    pushRatioPoint(-3.0, ft0cInnerY, ft0cInnerErr, fullFt0cY, fullFt0cErr, false, 0.3);
    pushRatioPoint(-2.4, ft0cOuterY, ft0cOuterErr, fullFt0cY, fullFt0cErr, false, 0.3);
    pushRatioPoint(2.4, ft0cOuterY, ft0cOuterErr, fullFt0cY, fullFt0cErr, true, 0.3);
    pushRatioPoint(3.0, ft0cInnerY, ft0cInnerErr, fullFt0cY, fullFt0cErr, true, 0.3);

    if (ft0aInnerAtLargeAbsEta) {
        pushRatioPoint(4.55, ft0aInnerY, ft0aInnerErr, fullFt0aY, fullFt0aErr, false, 0.35);
        pushRatioPoint(3.85, ft0aOuterY, ft0aOuterErr, fullFt0aY, fullFt0aErr, false, 0.35);
        pushRatioPoint(-3.85, ft0aOuterY, ft0aOuterErr, fullFt0aY, fullFt0aErr, true, 0.35);
        pushRatioPoint(-4.55, ft0aInnerY, ft0aInnerErr, fullFt0aY, fullFt0aErr, true, 0.35);
    } else {
        pushRatioPoint(3.85, ft0aInnerY, ft0aInnerErr, fullFt0aY, fullFt0aErr, false, 0.35);
        pushRatioPoint(4.55, ft0aOuterY, ft0aOuterErr, fullFt0aY, fullFt0aErr, false, 0.35);
        pushRatioPoint(-4.55, ft0aOuterY, ft0aOuterErr, fullFt0aY, fullFt0aErr, true, 0.35);
        pushRatioPoint(-3.85, ft0aInnerY, ft0aInnerErr, fullFt0aY, fullFt0aErr, true, 0.35);
    }

    return points;
}

std::vector<EtaPoint> BuildFullOnlyPointsWithRingAverage(TH1D* hFull, TH1D* hInnerSides, TH1D* hOuterSides) {
    // Keep TPC bins from full histogram; replace FT0 side values by (inner+outer)/2 when available.
    if (!hFull) return {};
    std::vector<EtaPoint> points;

    auto pushPoint = [&](double xval,
                         double yval,
                         double xerrLow,
                         double xerrHigh,
                         double yerrLow,
                         double yerrHigh,
                         bool mirrored) {
        points.push_back({xval, yval, xerrLow, xerrHigh, yerrLow, yerrHigh, mirrored});
    };

    double ft0cY = hFull->GetBinContent(1);
    double ft0cErr = hFull->GetBinError(1);
    double ft0aY = hFull->GetBinContent(hFull->GetNbinsX());
    double ft0aErr = hFull->GetBinError(hFull->GetNbinsX());

    if (hInnerSides && hOuterSides) {
        const double innerC = hInnerSides->GetBinContent(1);
        const double outerC = hOuterSides->GetBinContent(1);
        const double innerCErr = hInnerSides->GetBinError(1);
        const double outerCErr = hOuterSides->GetBinError(1);
        const double innerA = hInnerSides->GetBinContent(2);
        const double outerA = hOuterSides->GetBinContent(2);
        const double innerAErr = hInnerSides->GetBinError(2);
        const double outerAErr = hOuterSides->GetBinError(2);

        const bool validC = std::isfinite(innerC) && std::isfinite(outerC) && innerC > 0.0 && outerC > 0.0;
        const bool validA = std::isfinite(innerA) && std::isfinite(outerA) && innerA > 0.0 && outerA > 0.0;

        if (validC) {
            ft0cY = 0.5 * (innerC + outerC);
            ft0cErr = 0.5 * std::sqrt(innerCErr * innerCErr + outerCErr * outerCErr);
        }
        if (validA) {
            ft0aY = 0.5 * (innerA + outerA);
            ft0aErr = 0.5 * std::sqrt(innerAErr * innerAErr + outerAErr * outerAErr);
        }
    }

    pushPoint(-2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr, false);
    pushPoint(2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr, true);
    pushPoint(4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr, false);
    pushPoint(-4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr, true);

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

void Compare3times2PC_EtaDiff_TwoSystems_OO_TPCClusters70_Rings() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const char* fileOOFull      = "./3times2PC/Vn_LHC25ae_pass2_648798_Cent_0_20.root";
    const char* fileOOInnerC    = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_Cent_0_20.root";
    const char* fileOOInnerA    = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_Cent_0_20.root";
    const char* fileOOOuterC    = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_Cent_0_20.root";
    const char* fileOOOuterA    = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_Cent_0_20.root";
    const char* fileNeNeFull    = "./3times2PC/Vn_LHC25ae_pass2_659588_TPCClusters70_Cent_0_20.root";
    const char* fileNeNeInnerC  = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_TPCClusters70_Cent_0_20.root";
    const char* fileNeNeInnerA  = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_TPCClusters70_Cent_0_20.root";
    const char* fileNeNeOuterC  = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_TPCClusters70_Cent_0_20.root";
    const char* fileNeNeOuterA  = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_TPCClusters70_Cent_0_20.root";

    TFile* fOOFull    = TFile::Open(fileOOFull, "READ");
    TFile* fOOInnerC  = TFile::Open(fileOOInnerC, "READ");
    TFile* fOOInnerA  = TFile::Open(fileOOInnerA, "READ");
    TFile* fOOOuterC  = TFile::Open(fileOOOuterC, "READ");
    TFile* fOOOuterA  = TFile::Open(fileOOOuterA, "READ");
    TFile* fNeNeFull  = TFile::Open(fileNeNeFull, "READ");
    TFile* fNeNeInnerC = TFile::Open(fileNeNeInnerC, "READ");
    TFile* fNeNeInnerA = TFile::Open(fileNeNeInnerA, "READ");
    TFile* fNeNeOuterC = TFile::Open(fileNeNeOuterC, "READ");
    TFile* fNeNeOuterA = TFile::Open(fileNeNeOuterA, "READ");

    std::cout << "[Compare rings] O-O inputs: full, inner FT0C, inner FT0A, outer FT0C, outer FT0A" << std::endl;
    std::cout << "[Compare rings]   full  : " << fileOOFull << std::endl;
    std::cout << "[Compare rings]   innerC: " << fileOOInnerC << std::endl;
    std::cout << "[Compare rings]   innerA: " << fileOOInnerA << std::endl;
    std::cout << "[Compare rings]   outerC: " << fileOOOuterC << std::endl;
    std::cout << "[Compare rings]   outerA: " << fileOOOuterA << std::endl;
    std::cout << "[Compare rings] O-O TPCClusters70 inputs: full, inner FT0C, inner FT0A, outer FT0C, outer FT0A" << std::endl;
    std::cout << "[Compare rings]   full  : " << fileNeNeFull << std::endl;
    std::cout << "[Compare rings]   innerC: " << fileNeNeInnerC << std::endl;
    std::cout << "[Compare rings]   innerA: " << fileNeNeInnerA << std::endl;
    std::cout << "[Compare rings]   outerC: " << fileNeNeOuterC << std::endl;
    std::cout << "[Compare rings]   outerA: " << fileNeNeOuterA << std::endl;

    if (!fOOFull || !fOOFull->IsOpen()) {
        std::cerr << "Cannot open O-O full-range file: " << fileOOFull << std::endl;
        return;
    }
    if (!fNeNeFull || !fNeNeFull->IsOpen()) {
        std::cerr << "Cannot open O-O TPCClusters70 full-range file: " << fileNeNeFull << std::endl;
        fOOFull->Close();
        delete fOOFull;
        return;
    }
    if (!fOOInnerC || !fOOInnerC->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O inner-ring FT0C file: " << fileOOInnerC << std::endl;
    }
    if (!fOOInnerA || !fOOInnerA->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O inner-ring FT0A file: " << fileOOInnerA << std::endl;
    }
    if (!fOOOuterC || !fOOOuterC->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O outer-ring FT0C file: " << fileOOOuterC << std::endl;
    }
    if (!fOOOuterA || !fOOOuterA->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O outer-ring FT0A file: " << fileOOOuterA << std::endl;
    }
    if (!fNeNeInnerC || !fNeNeInnerC->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O TPCClusters70 inner-ring FT0C file: " << fileNeNeInnerC << std::endl;
    }
    if (!fNeNeInnerA || !fNeNeInnerA->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O TPCClusters70 inner-ring FT0A file: " << fileNeNeInnerA << std::endl;
    }
    if (!fNeNeOuterC || !fNeNeOuterC->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O TPCClusters70 outer-ring FT0C file: " << fileNeNeOuterC << std::endl;
    }
    if (!fNeNeOuterA || !fNeNeOuterA->IsOpen()) {
        std::cerr << "Warning: Cannot open O-O TPCClusters70 outer-ring FT0A file: " << fileNeNeOuterA << std::endl;
    }

    auto DrawSingleSystem = [&](int harmonic,
                                const char* canvasName,
                                const char* canvasTitle,
                                const char* outStem,
                                TGraphAsymmErrors* graphFilled,
                                TGraphAsymmErrors* graphMirror,
                                const char* legendLabel,
                                const char* systemLabel,
                                double yMax) {
        TCanvas* cSingle = new TCanvas(canvasName, canvasTitle, 1150, 700);
        cSingle->SetLeftMargin(0.11);
        cSingle->SetRightMargin(0.03);
        cSingle->SetTopMargin(0.06);
        cSingle->SetBottomMargin(0.12);
        cSingle->SetTicks(1, 1);

        TH1D* hFrameSingle = BuildEtaScaleHistogram(Form("hFrame_%s", canvasName), "", harmonic);
        hFrameSingle->GetYaxis()->SetRangeUser(0.0, yMax * 1.25);
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
        latexSingle.DrawLatex(0.13, 0.88, "ALICE Preliminary");
        latexSingle.SetTextSize(0.036);
        latexSingle.DrawLatex(0.13, 0.825, "3x2PC, 0-20%");
        latexSingle.DrawLatex(0.13, 0.775, systemLabel);

        cSingle->SaveAs(Form("./3times2PC/%s.root", outStem));
        cSingle->SaveAs(Form("./3times2PC/%s.png", outStem));
    };

    auto DrawRingToFullSystem = [&](int harmonic,
                                    const char* canvasName,
                                    const char* canvasTitle,
                                    const char* outStem,
                                    const char* systemLabel,
                                    int color,
                                    TH1D* hFull,
                                    TH1D* hInner,
                                    TH1D* hOuter,
                                    bool ft0aInnerAtLargeAbsEta) {
        if (!hFull || !hInner || !hOuter) {
            std::cout << "[Compare rings] Skip ring/full for " << systemLabel
                      << " v" << harmonic << ": missing full or ring histogram." << std::endl;
            return;
        }

        std::vector<EtaPoint> ringToFullPoints = BuildRingToFullPoints(hFull, hInner, hOuter, ft0aInnerAtLargeAbsEta);
        if (ringToFullPoints.empty()) {
            std::cout << "[Compare rings] Skip ring/full for " << systemLabel
                      << " v" << harmonic << ": no valid ratio points." << std::endl;
            return;
        }

        std::cout << "[Compare rings] " << systemLabel << " v" << harmonic << " ring/full values:" << std::endl;
        for (const auto& p : ringToFullPoints) {
            std::cout << "  eta=" << p.x << " ratio=" << p.y << " +/- " << p.eyHigh
                      << (p.mirrored ? " (mirror)" : "") << std::endl;
        }

        TGraphAsymmErrors* gFilled = BuildGraphFromPoints(Form("gRingToFullFilled_%s_v%d", systemLabel, harmonic),
                                                           ringToFullPoints,
                                                           false,
                                                           false);
        TGraphAsymmErrors* gMirror = BuildGraphFromPoints(Form("gRingToFullMirror_%s_v%d", systemLabel, harmonic),
                                                           ringToFullPoints,
                                                           true,
                                                           true);
        if (!gFilled) {
            std::cout << "[Compare rings] Skip ring/full for " << systemLabel
                      << " v" << harmonic << ": graph build failed." << std::endl;
            return;
        }

        std::vector<EtaPoint> ringToFullDeltaPercent;
        ringToFullDeltaPercent.reserve(ringToFullPoints.size());
        for (const auto& p : ringToFullPoints) {
            ringToFullDeltaPercent.push_back({p.x,
                                              100.0 * (p.y - 1.0),
                                              p.exLow,
                                              p.exHigh,
                                              100.0 * p.eyLow,
                                              100.0 * p.eyHigh,
                                              p.mirrored});
        }

        TGraphAsymmErrors* gDeltaFilled = BuildGraphFromPoints(Form("gRingToFullDeltaFilled_%s_v%d", systemLabel, harmonic),
                                                                ringToFullDeltaPercent,
                                                                false,
                                                                false);
        TGraphAsymmErrors* gDeltaMirror = BuildGraphFromPoints(Form("gRingToFullDeltaMirror_%s_v%d", systemLabel, harmonic),
                                                                ringToFullDeltaPercent,
                                                                true,
                                                                true);

        Style(gFilled, color, 20);
        Style(gMirror, color, 24);
        Style(gDeltaFilled, color, 20);
        Style(gDeltaMirror, color, 24);

        TH1D* hFrameRing = BuildRingToFullHistogram(Form("hFrameRingToFull_%s_v%d", systemLabel, harmonic));
        if (!hFrameRing) return;

        const double yMinRaw = std::min(GetGraphMin(gFilled), GetGraphMin(gMirror));
        const double yMaxRaw = std::max(GetGraphMax(gFilled), GetGraphMax(gMirror));
        double yMin = yMinRaw;
        double yMax = yMaxRaw;
        if (!std::isfinite(yMin) || !std::isfinite(yMax) || yMax <= yMin) {
            yMin = 0.7;
            yMax = 1.3;
        }
        const double pad = 0.20 * (yMax - yMin);
        hFrameRing->GetYaxis()->SetRangeUser(std::max(0.0, yMin - pad), yMax + pad);

        TCanvas* cRing = new TCanvas(canvasName, canvasTitle, 1100, 700);
        cRing->SetLeftMargin(0.11);
        cRing->SetRightMargin(0.03);
        cRing->SetTopMargin(0.06);
        cRing->SetBottomMargin(0.12);
        cRing->SetTicks(1, 1);

        hFrameRing->Draw("AXIS");
        TLine* unity = new TLine(-5.0, 1.0, 5.0, 1.0);
        unity->SetLineStyle(2);
        unity->SetLineColor(kGray + 2);
        unity->Draw("same");
        gFilled->Draw("PZ same");
        if (gMirror) gMirror->Draw("PZ same");

        TLegend* leg = new TLegend(0.67, 0.76, 0.91, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.030);
        leg->AddEntry(gFilled, Form("%s rings / %s full", systemLabel, systemLabel), "lep");
        if (gMirror) leg->AddEntry(gMirror, "Mirrored points", "p");
        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextColor(kBlack);
        latex.SetTextSize(0.042);
        latex.DrawLatex(0.13, 0.88, "ALICE Preliminary");
        latex.SetTextSize(0.036);
        latex.DrawLatex(0.13, 0.825, "3x2PC, 0-20%");
        latex.DrawLatex(0.13, 0.775, Form("%s: v_{%d} ring/full FT0 check", systemLabel, harmonic));

        cRing->SaveAs(Form("./3times2PC/%s.root", outStem));
        cRing->SaveAs(Form("./3times2PC/%s.png", outStem));

        if (gDeltaFilled) {
            TH1D* hFrameDelta = BuildRingToFullDeltaPercentHistogram(Form("hFrameRingToFullDelta_%s_v%d", systemLabel, harmonic));
            if (hFrameDelta) {
                const double dMinRaw = std::min(GetGraphMin(gDeltaFilled), GetGraphMin(gDeltaMirror));
                const double dMaxRaw = std::max(GetGraphMax(gDeltaFilled), GetGraphMax(gDeltaMirror));
                double dMin = dMinRaw;
                double dMax = dMaxRaw;
                if (!std::isfinite(dMin) || !std::isfinite(dMax) || dMax <= dMin) {
                    dMin = -10.0;
                    dMax = 10.0;
                }
                const double dPad = 0.20 * (dMax - dMin);
                hFrameDelta->GetYaxis()->SetRangeUser(dMin - dPad, dMax + dPad);

                TCanvas* cDelta = new TCanvas(Form("%s_delta", canvasName),
                                              Form("%s delta", canvasTitle),
                                              1100,
                                              700);
                cDelta->SetLeftMargin(0.11);
                cDelta->SetRightMargin(0.03);
                cDelta->SetTopMargin(0.06);
                cDelta->SetBottomMargin(0.12);
                cDelta->SetTicks(1, 1);

                hFrameDelta->Draw("AXIS");
                TLine* zero = new TLine(-5.0, 0.0, 5.0, 0.0);
                zero->SetLineStyle(2);
                zero->SetLineColor(kGray + 2);
                zero->Draw("same");
                gDeltaFilled->Draw("PZ same");
                if (gDeltaMirror) gDeltaMirror->Draw("PZ same");

                TLatex latexDelta;
                latexDelta.SetNDC();
                latexDelta.SetTextFont(42);
                latexDelta.SetTextColor(kBlack);
                latexDelta.SetTextSize(0.042);
                latexDelta.DrawLatex(0.13, 0.88, "ALICE Preliminary");
                latexDelta.SetTextSize(0.036);
                latexDelta.DrawLatex(0.13, 0.825, "3x2PC, 0-20%");
                latexDelta.DrawLatex(0.13, 0.775, Form("%s: v_{%d} ring/full deviation", systemLabel, harmonic));

                cDelta->SaveAs(Form("./3times2PC/%s_DeltaPercent.root", outStem));
                cDelta->SaveAs(Form("./3times2PC/%s_DeltaPercent.png", outStem));
            }
        }
    };

    const int harmonics[] = {2, 3, 4};
    bool anySaved = false;
    for (const int harmonic : harmonics) {
        TH1D* hOOFull    = LoadFullVn(fOOFull, harmonic);
        TH1D* hOOInner   = BuildRingSidesFromSideFiles(fOOInnerC, fOOInnerA, harmonic, Form("hOOInnerSides_v%d", harmonic), Form("O-O inner sides;side;v_{%d}", harmonic));
        TH1D* hOOOuter   = BuildRingSidesFromSideFiles(fOOOuterC, fOOOuterA, harmonic, Form("hOOOuterSides_v%d", harmonic), Form("O-O outer sides;side;v_{%d}", harmonic));
        TH1D* hNeNeFull  = LoadFullVn(fNeNeFull, harmonic);
        TH1D* hNeNeInner = BuildRingSidesFromSideFiles(fNeNeInnerC, fNeNeInnerA, harmonic, Form("hNeNeInnerSides_v%d", harmonic), Form("O-O TPCClusters70 inner sides;side;v_{%d}", harmonic));
        TH1D* hNeNeOuter = BuildRingSidesFromSideFiles(fNeNeOuterC, fNeNeOuterA, harmonic, Form("hNeNeOuterSides_v%d", harmonic), Form("O-O TPCClusters70 outer sides;side;v_{%d}", harmonic));

        if (!hOOInner || !hOOOuter || !hNeNeInner || !hNeNeOuter) {
            std::cout << "[Compare rings] Falling back to legacy combined ring files for v" << harmonic << std::endl;
            delete hOOInner;
            delete hOOOuter;
            delete hNeNeInner;
            delete hNeNeOuter;
            hOOInner = LoadRingVn(fOOInnerC, harmonic);
            hOOOuter = LoadRingVn(fOOOuterC, harmonic);
            hNeNeInner = LoadRingVn(fNeNeInnerC, harmonic);
            hNeNeOuter = LoadRingVn(fNeNeOuterC, harmonic);
        }

        if (!hOOFull || !hNeNeFull) {
            std::cout << "[Compare rings] Skip v" << harmonic
                      << ": missing full-range histogram in O-O or O-O TPCClusters70 file." << std::endl;
            continue;
        }

        const bool ooHasRings = (hOOInner && hOOOuter);
        std::vector<EtaPoint> ooPoints = BuildEtaPoints(hOOFull, hOOInner, hOOOuter, ooHasRings, false);
        std::vector<EtaPoint> nePoints = BuildEtaPoints(hNeNeFull, hNeNeInner, hNeNeOuter, true, false);
        std::vector<EtaPoint> ratioPoints = BuildRatioPoints(nePoints, ooPoints);

        TH1D* hFrame = BuildEtaScaleHistogram(Form("hFrame_v%d", harmonic), "", harmonic);
        TH1D* hRatioFrame = BuildRatioHistogram(Form("hRatioFrame_v%d", harmonic));

        TGraphAsymmErrors* gOOFilled = BuildGraphFromPoints(Form("gOOFilled_v%d", harmonic), ooPoints, false, false);
        TGraphAsymmErrors* gOOMirror = BuildGraphFromPoints(Form("gOOMirror_v%d", harmonic), ooPoints, true, true);
        TGraphAsymmErrors* gNeNeFilled = BuildGraphFromPoints(Form("gNeNeFilled_v%d", harmonic), nePoints, false, false);
        TGraphAsymmErrors* gNeNeMirror = BuildGraphFromPoints(Form("gNeNeMirror_v%d", harmonic), nePoints, true, true);
        TGraphAsymmErrors* gRatioFilled = BuildGraphFromPoints(Form("gRatioFilled_v%d", harmonic), ratioPoints, false, false);
        TGraphAsymmErrors* gRatioMirror = BuildGraphFromPoints(Form("gRatioMirror_v%d", harmonic), ratioPoints, true, true);

        if (!hFrame || !hRatioFrame || !gOOFilled || !gNeNeFilled || !gRatioFilled) {
            std::cout << "[Compare rings] Skip v" << harmonic
                      << ": failed to build frame/graphs for this harmonic." << std::endl;
            continue;
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

        TCanvas* c = new TCanvas(Form("cEtaDiff2SysRings_v%d", harmonic),
                                 Form("3x2PC v%d with explicit FT0 eta scale and ratio", harmonic),
                                 1150,
                                 850);

        TPad* padTop = new TPad(Form("padTop_v%d", harmonic), "padTop", 0.0, 0.30, 1.0, 1.0);
        padTop->SetLeftMargin(0.11);
        padTop->SetRightMargin(0.03);
        padTop->SetTopMargin(0.07);
        padTop->SetBottomMargin(0.02);
        padTop->SetTicks(1, 1);
        padTop->Draw();

        TPad* padBottom = new TPad(Form("padBottom_v%d", harmonic), "padBottom", 0.0, 0.00, 1.0, 0.30);
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
        leg->AddEntry(gNeNeFilled,  "O-O TPCClusters70", "lep");
        if (gOOMirror) leg->AddEntry(gOOMirror, "Mirrored points", "p");
        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextColor(kBlack);
        latex.SetTextSize(0.042);
        latex.DrawLatex(0.13, 0.88, "ALICE Preliminary");
        latex.SetTextSize(0.036);
        latex.DrawLatex(0.13, 0.825, "3x2PC, 0-20%");
        latex.DrawLatex(0.13, 0.775, "FT0 rings + TPC acceptance (+ ratio)");

        padBottom->cd();
        hRatioFrame->GetYaxis()->SetRangeUser(0.95, 1.05);
        hRatioFrame->Draw("AXIS");
        TLine* unity = new TLine(-5.0, 1.0, 5.0, 1.0);
        unity->SetLineStyle(2);
        unity->SetLineColor(kGray + 2);
        unity->Draw("same");
        gRatioFilled->Draw("PZ same");
        if (gRatioMirror) gRatioMirror->Draw("PZ same");

        c->SaveAs(Form("./3times2PC/Compare_v%d_FT0Side_OO_vs_OO_TPCClusters70_inner_outer_full.root", harmonic));
        c->SaveAs(Form("./3times2PC/Compare_v%d_FT0Side_OO_vs_OO_TPCClusters70_inner_outer_full.png", harmonic));

        std::vector<EtaPoint> ooFullOnlyPoints = BuildFullOnlyPointsWithRingAverage(hOOFull, hOOInner, hOOOuter);
        std::vector<EtaPoint> neFullOnlyPoints = BuildFullOnlyPointsWithRingAverage(hNeNeFull, hNeNeInner, hNeNeOuter);
        std::vector<EtaPoint> fullOnlyRatioPoints = BuildRatioPoints(neFullOnlyPoints, ooFullOnlyPoints);

        TGraphAsymmErrors* gOOFullOnlyFilled = BuildGraphFromPoints(Form("gOOFullOnlyFilled_v%d", harmonic), ooFullOnlyPoints, false, false);
        TGraphAsymmErrors* gOOFullOnlyMirror = BuildGraphFromPoints(Form("gOOFullOnlyMirror_v%d", harmonic), ooFullOnlyPoints, true, true);
        TGraphAsymmErrors* gNeNeFullOnlyFilled = BuildGraphFromPoints(Form("gNeNeFullOnlyFilled_v%d", harmonic), neFullOnlyPoints, false, false);
        TGraphAsymmErrors* gNeNeFullOnlyMirror = BuildGraphFromPoints(Form("gNeNeFullOnlyMirror_v%d", harmonic), neFullOnlyPoints, true, true);
        TGraphAsymmErrors* gFullOnlyRatioFilled = BuildGraphFromPoints(Form("gFullOnlyRatioFilled_v%d", harmonic), fullOnlyRatioPoints, false, false);
        TGraphAsymmErrors* gFullOnlyRatioMirror = BuildGraphFromPoints(Form("gFullOnlyRatioMirror_v%d", harmonic), fullOnlyRatioPoints, true, true);

        if (gOOFullOnlyFilled && gNeNeFullOnlyFilled && gFullOnlyRatioFilled) {
            Style(gOOFullOnlyFilled, kRed + 1, 21);
            Style(gOOFullOnlyMirror, kRed + 1, 25);
            Style(gNeNeFullOnlyFilled, kGreen + 2, 20);
            Style(gNeNeFullOnlyMirror, kGreen + 2, 24);
            Style(gFullOnlyRatioFilled, kBlack, 20);
            Style(gFullOnlyRatioMirror, kBlack, 24);

            TH1D* hFullOnlyFrame = BuildEtaScaleHistogram(Form("hFullOnlyFrame_v%d", harmonic), "", harmonic);
            TH1D* hFullOnlyRatioFrame = BuildRatioHistogram(Form("hFullOnlyRatioFrame_v%d", harmonic));
            if (hFullOnlyFrame && hFullOnlyRatioFrame) {
                double fullOnlyMax = 0.0;
                fullOnlyMax = std::max(fullOnlyMax, GetGraphMax(gOOFullOnlyFilled));
                fullOnlyMax = std::max(fullOnlyMax, GetGraphMax(gOOFullOnlyMirror));
                fullOnlyMax = std::max(fullOnlyMax, GetGraphMax(gNeNeFullOnlyFilled));
                fullOnlyMax = std::max(fullOnlyMax, GetGraphMax(gNeNeFullOnlyMirror));
                if (fullOnlyMax <= 0.0) fullOnlyMax = 0.1;

                TCanvas* cFullOnly = new TCanvas(Form("cEtaDiff2SysFullOnly_v%d", harmonic),
                                                 Form("3x2PC v%d full FT0 comparison without rings", harmonic),
                                                 1150,
                                                 850);

                TPad* padTopFullOnly = new TPad(Form("padTopFullOnly_v%d", harmonic), "padTopFullOnly", 0.0, 0.30, 1.0, 1.0);
                padTopFullOnly->SetLeftMargin(0.11);
                padTopFullOnly->SetRightMargin(0.03);
                padTopFullOnly->SetTopMargin(0.07);
                padTopFullOnly->SetBottomMargin(0.02);
                padTopFullOnly->SetTicks(1, 1);
                padTopFullOnly->Draw();

                TPad* padBottomFullOnly = new TPad(Form("padBottomFullOnly_v%d", harmonic), "padBottomFullOnly", 0.0, 0.00, 1.0, 0.30);
                padBottomFullOnly->SetLeftMargin(0.11);
                padBottomFullOnly->SetRightMargin(0.03);
                padBottomFullOnly->SetTopMargin(0.02);
                padBottomFullOnly->SetBottomMargin(0.34);
                padBottomFullOnly->SetTicks(1, 1);
                padBottomFullOnly->Draw();

                padTopFullOnly->cd();
                hFullOnlyFrame->GetXaxis()->SetLabelSize(0.0);
                hFullOnlyFrame->GetYaxis()->SetRangeUser(0.0, fullOnlyMax * 1.25);
                hFullOnlyFrame->Draw("AXIS");
                gNeNeFullOnlyFilled->Draw("PZ same");
                if (gNeNeFullOnlyMirror) gNeNeFullOnlyMirror->Draw("PZ same");
                gOOFullOnlyFilled->Draw("PZ same");
                if (gOOFullOnlyMirror) gOOFullOnlyMirror->Draw("PZ same");

                TLegend* legFullOnly = new TLegend(0.74, 0.78, 0.91, 0.90);
                legFullOnly->SetBorderSize(0);
                legFullOnly->SetFillStyle(0);
                legFullOnly->SetTextFont(42);
                legFullOnly->SetTextSize(0.030);
                legFullOnly->AddEntry(gOOFullOnlyFilled, "O-O (full FT0)", "lep");
                legFullOnly->AddEntry(gNeNeFullOnlyFilled, "O-O TPCClusters70 (full FT0)", "lep");
                if (gOOFullOnlyMirror) legFullOnly->AddEntry(gOOFullOnlyMirror, "Mirrored points", "p");
                legFullOnly->Draw();

                TLatex latexFullOnly;
                latexFullOnly.SetNDC();
                latexFullOnly.SetTextFont(42);
                latexFullOnly.SetTextColor(kBlack);
                latexFullOnly.SetTextSize(0.042);
                latexFullOnly.DrawLatex(0.13, 0.88, "ALICE Preliminary");
                latexFullOnly.SetTextSize(0.036);
                latexFullOnly.DrawLatex(0.13, 0.825, "3x2PC, 0-20%");
                latexFullOnly.DrawLatex(0.13, 0.775, "Full FT0 ranges + TPC acceptance (+ratio)");

                padBottomFullOnly->cd();
                hFullOnlyRatioFrame->GetYaxis()->SetRangeUser(0.95, 1.05);
                hFullOnlyRatioFrame->Draw("AXIS");
                TLine* fullOnlyUnity = new TLine(-5.0, 1.0, 5.0, 1.0);
                fullOnlyUnity->SetLineStyle(2);
                fullOnlyUnity->SetLineColor(kGray + 2);
                fullOnlyUnity->Draw("same");
                gFullOnlyRatioFilled->Draw("PZ same");
                if (gFullOnlyRatioMirror) gFullOnlyRatioMirror->Draw("PZ same");

                cFullOnly->SaveAs(Form("./3times2PC/Compare_v%d_FT0Full_OO_vs_OO_TPCClusters70_ratio.root", harmonic));
                cFullOnly->SaveAs(Form("./3times2PC/Compare_v%d_FT0Full_OO_vs_OO_TPCClusters70_ratio.png", harmonic));
                std::cout << "[Compare rings] Saved full-FT0-only O-O vs O-O TPCClusters70 comparison for v"
                          << harmonic << " as Compare_v" << harmonic
                          << "_FT0Full_OO_vs_OO_TPCClusters70_ratio.png" << std::endl;
            }
        }

        DrawSingleSystem(harmonic,
                         Form("cEtaDiff_OOOnly_v%d", harmonic),
                         Form("3x2PC O-O only v%d", harmonic),
                         Form("Compare_v%d_FT0Side_OO_only_inner_outer_full", harmonic),
                         gOOFilled,
                         gOOMirror,
                         ooHasRings ? "O-O" : "O-O (full FT0 baseline)",
                         "O-O, FT0 rings + TPC acceptance",
                         maxVal);

        DrawSingleSystem(harmonic,
                         Form("cEtaDiff_NeNeOnly_v%d", harmonic),
                         Form("3x2PC O-O TPCClusters70 only v%d", harmonic),
                         Form("Compare_v%d_FT0Side_OO_TPCClusters70_only_inner_outer_full", harmonic),
                         gNeNeFilled,
                         gNeNeMirror,
                         "O-O TPCClusters70",
                         "O-O TPCClusters70, FT0 rings + TPC acceptance",
                         maxVal);

        DrawRingToFullSystem(harmonic,
                             Form("cRingToFull_OO_v%d", harmonic),
                             Form("FT0 ring/full check O-O v%d", harmonic),
                             Form("Compare_v%d_FT0RingToFull_OO", harmonic),
                             "O-O",
                             kRed + 1,
                             hOOFull,
                             hOOInner,
                             hOOOuter,
                             true);

        DrawRingToFullSystem(harmonic,
                             Form("cRingToFull_NeNe_v%d", harmonic),
                             Form("FT0 ring/full check O-O TPCClusters70 v%d", harmonic),
                             Form("Compare_v%d_FT0RingToFull_OO_TPCClusters70", harmonic),
                             "O-O TPCClusters70",
                             kGreen + 2,
                             hNeNeFull,
                             hNeNeInner,
                             hNeNeOuter,
                             false);

        anySaved = true;
        if (!ooHasRings) {
            std::cout << "Saved v" << harmonic
                      << " comparison plot. Note: O-O ring files missing, so O-O uses full-range FT0 points." << std::endl;
        } else {
            std::cout << "Saved v" << harmonic
                      << " systematic-check comparison plot with FT0 ring points for O-O and O-O TPCClusters70." << std::endl;
        }
            std::cout << "Saved individual O-O and O-O TPCClusters70 v" << harmonic << " plots in the same style." << std::endl;
    }

    if (!anySaved) {
        std::cerr << "No harmonic plots were produced (v2/v3/v4 missing in inputs)." << std::endl;
    }

    fOOFull->Close();
    fNeNeFull->Close();
    if (fOOInnerC) fOOInnerC->Close();
    if (fOOInnerA) fOOInnerA->Close();
    if (fOOOuterC) fOOOuterC->Close();
    if (fOOOuterA) fOOOuterA->Close();
    if (fNeNeInnerC) fNeNeInnerC->Close();
    if (fNeNeInnerA) fNeNeInnerA->Close();
    if (fNeNeOuterC) fNeNeOuterC->Close();
    if (fNeNeOuterA) fNeNeOuterA->Close();
    delete fOOFull;
    delete fNeNeFull;
    if (fOOInnerC) delete fOOInnerC;
    if (fOOInnerA) delete fOOInnerA;
    if (fOOOuterC) delete fOOOuterC;
    if (fOOOuterA) delete fOOOuterA;
    if (fNeNeInnerC) delete fNeNeInnerC;
    if (fNeNeInnerA) delete fNeNeInnerA;
    if (fNeNeOuterC) delete fNeNeOuterC;
    if (fNeNeOuterA) delete fNeNeOuterA;
}
