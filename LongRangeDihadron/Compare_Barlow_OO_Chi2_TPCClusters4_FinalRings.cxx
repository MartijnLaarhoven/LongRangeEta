#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace TPCC90Barlow {
constexpr double kBarlowAcceptThreshold = 1.0;

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
    if (!f || !f->IsOpen() || f->IsZombie()) return nullptr;
    TH1D* h = static_cast<TH1D*>(f->Get(Form("hV%d_Combined", harmonic)));
    if (!h) h = static_cast<TH1D*>(f->Get(Form("hV%d", harmonic)));
    return h;
}

TH1D* LoadSideValue(TFile* f, int harmonic) {
    if (!f || !f->IsOpen() || f->IsZombie()) return nullptr;
    TH1D* h = static_cast<TH1D*>(f->Get(Form("hV%d_Side", harmonic)));
    if (!h) h = static_cast<TH1D*>(f->Get(Form("hV%d_Sides", harmonic)));
    if (!h) h = static_cast<TH1D*>(f->Get(Form("hV%d_Combined", harmonic)));
    if (!h) h = static_cast<TH1D*>(f->Get(Form("hV%d", harmonic)));
    return h;
}

TH1D* BuildRingSidesFromSideFiles(TFile* fFT0C, TFile* fFT0A, int harmonic, const char* name, const char* title) {
    TH1D* hFT0C = LoadSideValue(fFT0C, harmonic);
    TH1D* hFT0A = LoadSideValue(fFT0A, harmonic);
    if (!hFT0C || !hFT0A) return nullptr;

    TH1D* h = new TH1D(name, title, 2, 0.5, 2.5);
    h->GetXaxis()->SetBinLabel(1, "FT0C");
    h->GetXaxis()->SetBinLabel(2, "FT0A");
    h->SetBinContent(1, hFT0C->GetBinContent(1));
    h->SetBinError(1, hFT0C->GetBinError(1));
    h->SetBinContent(2, hFT0A->GetBinContent(1));
    h->SetBinError(2, hFT0A->GetBinError(1));
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

        auto sanitizeRingPoint = [&](double y, double err, double fallbackY, double fallbackErr) {
            const bool invalid = (!std::isfinite(y) || !std::isfinite(err) || y <= 0.0 || err >= 9.9);
            if (invalid) {
                return std::make_pair(fallbackY, fallbackErr);
            }
            return std::make_pair(y, err);
        };

        auto ft0cInner = sanitizeRingPoint(hInnerSides->GetBinContent(1), hInnerSides->GetBinError(1), fullFt0cY, fullFt0cErr);
        auto ft0cOuter = sanitizeRingPoint(hOuterSides->GetBinContent(1), hOuterSides->GetBinError(1), fullFt0cY, fullFt0cErr);
        auto ft0aOuter = sanitizeRingPoint(hOuterSides->GetBinContent(2), hOuterSides->GetBinError(2), fullFt0aY, fullFt0aErr);
        auto ft0aInner = sanitizeRingPoint(hInnerSides->GetBinContent(2), hInnerSides->GetBinError(2), fullFt0aY, fullFt0aErr);

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

    std::sort(points.begin(), points.end(), [](const EtaPoint& a, const EtaPoint& b) {
        if (a.x == b.x) return a.mirrored < b.mirrored;
        return a.x < b.x;
    });

    return points;
}

double ComputeBarlow(double nominal, double nominalErr, double varied, double variedErr) {
    const double diff = std::fabs(varied - nominal);
    const double sigmaVaried2 = variedErr * variedErr;
    const double sigmaNominal2 = nominalErr * nominalErr;
    const double sigmaBarlow2 = std::fabs(sigmaVaried2 - sigmaNominal2);
    const double sigmaIndependent2 = sigmaVaried2 + sigmaNominal2;

    if (sigmaIndependent2 <= 1e-16) return 0.0;

    const double relSeparation = sigmaBarlow2 / sigmaIndependent2;
    const bool useIndependent = (relSeparation < 0.20);
    const double den2 = useIndependent ? sigmaIndependent2 : sigmaBarlow2;
    const double den = std::sqrt(den2);
    if (den <= 0.0) return 0.0;
    return diff / den;
}

bool AcceptByBarlow(double barlow) {
    return std::isfinite(barlow) && barlow < kBarlowAcceptThreshold;
}

std::string EtaLabel(double x) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << x;
    return ss.str();
}

}

void Compare_Barlow_OO_Chi2_TPCClusters4_FinalRings() {
    using namespace TPCC90Barlow;
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir("./3times2PC", kTRUE);
    gSystem->mkdir("./3times2PC/PDFs", kTRUE);

    const char* nominalFullPath   = "./3times2PC/Vn_LHC25ae_pass2_648798_Cent_0_20.root";
    const char* nominalInnerCPath = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_Cent_0_20.root";
    const char* nominalInnerAPath = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_Cent_0_20.root";
    const char* nominalOuterCPath = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_Cent_0_20.root";
    const char* nominalOuterAPath = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_Cent_0_20.root";

    const char* variedFullPath   = "./3times2PC/Vn_LHC25ae_pass2_659596_Chi2_TPCClusters4_Cent_0_20.root";
    const char* variedInnerCPath = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_Chi2_TPCClusters4_Cent_0_20.root";
    const char* variedInnerAPath = "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_Chi2_TPCClusters4_Cent_0_20.root";
    const char* variedOuterCPath = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_Chi2_TPCClusters4_Cent_0_20.root";
    const char* variedOuterAPath = "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_Chi2_TPCClusters4_Cent_0_20.root";

    TFile* fNominalFull = TFile::Open(nominalFullPath, "READ");
    TFile* fNominalInnerC = TFile::Open(nominalInnerCPath, "READ");
    TFile* fNominalInnerA = TFile::Open(nominalInnerAPath, "READ");
    TFile* fNominalOuterC = TFile::Open(nominalOuterCPath, "READ");
    TFile* fNominalOuterA = TFile::Open(nominalOuterAPath, "READ");

    TFile* fVariedFull = TFile::Open(variedFullPath, "READ");
    TFile* fVariedInnerC = TFile::Open(variedInnerCPath, "READ");
    TFile* fVariedInnerA = TFile::Open(variedInnerAPath, "READ");
    TFile* fVariedOuterC = TFile::Open(variedOuterCPath, "READ");
    TFile* fVariedOuterA = TFile::Open(variedOuterAPath, "READ");

    if (!fNominalFull || !fVariedFull || !fNominalFull->IsOpen() || !fVariedFull->IsOpen()) {
        std::cerr << "[Barlow final-rings] Missing full-range nominal or varied file." << std::endl;
        return;
    }
    if (!fNominalInnerC || !fNominalInnerA || !fNominalOuterC || !fNominalOuterA ||
        !fVariedInnerC || !fVariedInnerA || !fVariedOuterC || !fVariedOuterA ||
        !fNominalInnerC->IsOpen() || !fNominalInnerA->IsOpen() || !fNominalOuterC->IsOpen() || !fNominalOuterA->IsOpen() ||
        !fVariedInnerC->IsOpen() || !fVariedInnerA->IsOpen() || !fVariedOuterC->IsOpen() || !fVariedOuterA->IsOpen()) {
        std::cerr << "[Barlow final-rings] Missing one or more ring files; final ring-included Barlow cannot run." << std::endl;
        return;
    }

    TFile* out = TFile::Open("./3times2PC/Barlow_OO_vs_OO_Chi2_TPCClusters4_FinalRings_Cent_0_20.root", "RECREATE");
    if (!out || !out->IsOpen() || out->IsZombie()) {
        std::cerr << "[Barlow final-rings] Cannot create output ROOT file." << std::endl;
        return;
    }

    for (int harmonic = 2; harmonic <= 4; ++harmonic) {
        TH1D* hNominalFull = LoadFullVn(fNominalFull, harmonic);
        TH1D* hVariedFull = LoadFullVn(fVariedFull, harmonic);
        TH1D* hNominalInner = BuildRingSidesFromSideFiles(fNominalInnerC, fNominalInnerA, harmonic,
                                                          Form("hNominalInnerSides_v%d", harmonic), "");
        TH1D* hNominalOuter = BuildRingSidesFromSideFiles(fNominalOuterC, fNominalOuterA, harmonic,
                                                          Form("hNominalOuterSides_v%d", harmonic), "");
        TH1D* hVariedInner = BuildRingSidesFromSideFiles(fVariedInnerC, fVariedInnerA, harmonic,
                                                         Form("hVariedInnerSides_v%d", harmonic), "");
        TH1D* hVariedOuter = BuildRingSidesFromSideFiles(fVariedOuterC, fVariedOuterA, harmonic,
                                                         Form("hVariedOuterSides_v%d", harmonic), "");

        if (!hNominalFull || !hVariedFull || !hNominalInner || !hNominalOuter || !hVariedInner || !hVariedOuter) {
            std::cerr << "[Barlow final-rings] Missing histograms for v" << harmonic << ", skipping." << std::endl;
            delete hNominalInner;
            delete hNominalOuter;
            delete hVariedInner;
            delete hVariedOuter;
            continue;
        }

        const auto nominalPoints = BuildEtaPoints(hNominalFull, hNominalInner, hNominalOuter, true, false);
        const auto variedPoints = BuildEtaPoints(hVariedFull, hVariedInner, hVariedOuter, true, false);

        std::map<int, EtaPoint> variedByX;
        auto keyX = [](double x) { return static_cast<int>(std::lround(1000.0 * x)); };
        for (const auto& p : variedPoints) {
            variedByX[keyX(p.x)] = p;
        }

        std::vector<EtaPoint> matchedNominal;
        std::vector<EtaPoint> matchedVaried;
        for (const auto& p : nominalPoints) {
            const int k = keyX(p.x);
            auto it = variedByX.find(k);
            if (it == variedByX.end()) continue;
            matchedNominal.push_back(p);
            matchedVaried.push_back(it->second);
        }

        if (matchedNominal.empty()) {
            std::cerr << "[Barlow final-rings] No matched eta points for v" << harmonic << ", skipping." << std::endl;
            delete hNominalInner;
            delete hNominalOuter;
            delete hVariedInner;
            delete hVariedOuter;
            continue;
        }

        const int nBins = static_cast<int>(matchedNominal.size());
        TH1D* hNominal = new TH1D(Form("hNominalFinalRings_v%d", harmonic),
                                  Form("Nominal final ring-included v_{%d};eta point;v_{%d}", harmonic, harmonic),
                                  nBins, 0.5, nBins + 0.5);
        TH1D* hVaried = new TH1D(Form("hVariedFinalRings_v%d", harmonic),
                                 Form("Varied final ring-included v_{%d};eta point;v_{%d}", harmonic, harmonic),
                                 nBins, 0.5, nBins + 0.5);
        TH1D* hBarlow = new TH1D(Form("hBarlowFinalRings_v%d", harmonic),
                                 Form("Barlow final ring-included v_{%d};eta point;|#Delta|/#sigma_{Barlow}", harmonic),
                                 nBins, 0.5, nBins + 0.5);
        TH1D* hSystAbs = new TH1D(Form("hSystematicAbsFinalRings_v%d", harmonic),
                                  Form("Accepted systematic uncertainty v_{%d};eta point;|#Delta v_{%d}|", harmonic, harmonic),
                                  nBins, 0.5, nBins + 0.5);
        TH1D* hSystRel = new TH1D(Form("hSystematicRelFinalRings_v%d", harmonic),
                                  Form("Accepted systematic uncertainty v_{%d};eta point;100#times|#Delta|/v_{%d}", harmonic, harmonic),
                                  nBins, 0.5, nBins + 0.5);
        TH1D* hTotalAbs = new TH1D(Form("hTotalAbsFinalRings_v%d", harmonic),
                                   Form("Total uncertainty v_{%d};eta point;#sigma_{tot}(v_{%d})", harmonic, harmonic),
                                   nBins, 0.5, nBins + 0.5);
        TH1D* hTotalRel = new TH1D(Form("hTotalRelFinalRings_v%d", harmonic),
                                   Form("Total uncertainty v_{%d};eta point;100#times#sigma_{tot}/v_{%d}", harmonic, harmonic),
                                   nBins, 0.5, nBins + 0.5);
        TH1D* hDecision = new TH1D(Form("hBarlowDecisionFinalRings_v%d", harmonic),
                                   Form("Barlow decision v_{%d};eta point;1=accept, 0=overrule", harmonic),
                                   nBins, 0.5, nBins + 0.5);

        for (int i = 0; i < nBins; ++i) {
            const int bin = i + 1;
            const double nominal = matchedNominal[i].y;
            const double nominalErr = std::max(matchedNominal[i].eyLow, matchedNominal[i].eyHigh);
            const double varied = matchedVaried[i].y;
            const double variedErr = std::max(matchedVaried[i].eyLow, matchedVaried[i].eyHigh);
            const double diff = varied - nominal;
            const double barlow = ComputeBarlow(nominal, nominalErr, varied, variedErr);
            const bool accepted = AcceptByBarlow(barlow);
            const double systAbs = accepted ? std::fabs(diff) : 0.0;
            const double systRel = (accepted && std::fabs(nominal) > 1e-16) ? 100.0 * systAbs / std::fabs(nominal) : 0.0;
            const double totalAbs = std::sqrt(nominalErr * nominalErr + systAbs * systAbs);
            const double totalRel = (std::fabs(nominal) > 1e-16) ? 100.0 * totalAbs / std::fabs(nominal) : 0.0;

            const std::string label = EtaLabel(matchedNominal[i].x);

            hNominal->SetBinContent(bin, nominal);
            hNominal->SetBinError(bin, nominalErr);
            hVaried->SetBinContent(bin, varied);
            hVaried->SetBinError(bin, variedErr);
            hBarlow->SetBinContent(bin, barlow);
            hSystAbs->SetBinContent(bin, systAbs);
            hSystRel->SetBinContent(bin, systRel);
            hTotalAbs->SetBinContent(bin, totalAbs);
            hTotalRel->SetBinContent(bin, totalRel);
            hDecision->SetBinContent(bin, accepted ? 1.0 : 0.0);

            hNominal->GetXaxis()->SetBinLabel(bin, label.c_str());
            hVaried->GetXaxis()->SetBinLabel(bin, label.c_str());
            hBarlow->GetXaxis()->SetBinLabel(bin, label.c_str());
            hSystAbs->GetXaxis()->SetBinLabel(bin, label.c_str());
            hSystRel->GetXaxis()->SetBinLabel(bin, label.c_str());
            hTotalAbs->GetXaxis()->SetBinLabel(bin, label.c_str());
            hTotalRel->GetXaxis()->SetBinLabel(bin, label.c_str());
            hDecision->GetXaxis()->SetBinLabel(bin, label.c_str());

            std::cout << "[Barlow final-rings] v" << harmonic
                      << " eta=" << matchedNominal[i].x
                      << " nominal=" << nominal << " +/- " << nominalErr
                      << " varied=" << varied << " +/- " << variedErr
                      << " barlow=" << barlow
                      << " decision=" << (accepted ? "ACCEPT" : "OVERRULE")
                      << " systAbs=" << systAbs
                      << " totalAbs=" << totalAbs
                      << std::endl;
        }

        out->cd();
        hNominal->Write();
        hVaried->Write();
        hBarlow->Write();
        hSystAbs->Write();
        hSystRel->Write();
        hTotalAbs->Write();
        hTotalRel->Write();
        hDecision->Write();

        TH1D* hNominalBand = static_cast<TH1D*>(hNominal->Clone(Form("hNominalBandFinalRings_v%d", harmonic)));
        hNominalBand->SetDirectory(nullptr);
        for (int bin = 1; bin <= nBins; ++bin) {
            hNominalBand->SetBinError(bin, hTotalAbs->GetBinContent(bin));
        }

        TCanvas* c = new TCanvas(Form("cBarlowFinalRings_v%d", harmonic), "", 1050, 860);
        TPad* padTop = new TPad(Form("padTop_v%d", harmonic), "", 0.0, 0.35, 1.0, 1.0);
        TPad* padBottom = new TPad(Form("padBottom_v%d", harmonic), "", 0.0, 0.0, 1.0, 0.35);
        padTop->SetBottomMargin(0.02);
        padTop->SetLeftMargin(0.09);
        padTop->SetRightMargin(0.03);
        padBottom->SetTopMargin(0.04);
        padBottom->SetBottomMargin(0.28);
        padBottom->SetLeftMargin(0.09);
        padBottom->SetRightMargin(0.03);
        padTop->Draw();
        padBottom->Draw();

        padTop->cd();
        hNominalBand->SetStats(0);
        hNominalBand->GetXaxis()->SetLabelSize(0.0);
        hNominalBand->GetYaxis()->SetTitleSize(0.06);
        hNominalBand->GetYaxis()->SetLabelSize(0.05);
        hNominalBand->SetFillColor(kGray + 1);
        hNominalBand->SetFillStyle(3004);
        hNominalBand->SetLineColor(kGray + 2);

        double yMax = 0.0;
        double yMin = 1e9;
        for (int bin = 1; bin <= nBins; ++bin) {
            yMax = std::max(yMax, hVaried->GetBinContent(bin) + hVaried->GetBinError(bin));
            yMax = std::max(yMax, hNominal->GetBinContent(bin) + hNominalBand->GetBinError(bin));
            yMin = std::min(yMin, hVaried->GetBinContent(bin) - hVaried->GetBinError(bin));
            yMin = std::min(yMin, hNominal->GetBinContent(bin) - hNominalBand->GetBinError(bin));
        }
        if (!std::isfinite(yMin)) yMin = 0.0;
        if (!std::isfinite(yMax) || yMax <= yMin) yMax = yMin + 0.05;
        hNominalBand->GetYaxis()->SetRangeUser(yMin - 0.20 * std::fabs(yMin), yMax * 1.25);
        hNominalBand->Draw("E2");

        hNominal->SetLineColor(kBlue + 1);
        hNominal->SetMarkerColor(kBlue + 1);
        hNominal->SetMarkerStyle(20);
        hNominal->Draw("E1 SAME");

        hVaried->SetLineColor(kRed + 1);
        hVaried->SetMarkerColor(kRed + 1);
        hVaried->SetMarkerStyle(24);
        hVaried->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.11, 0.74, 0.43, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.032);
        leg->AddEntry(hNominal, "O-O nominal", "lep");
        leg->AddEntry(hVaried, "O-O Chi2_TPCClusters4", "lep");
        leg->AddEntry(hNominalBand, Form("Total uncertainty (Barlow < %.1f)", kBarlowAcceptThreshold), "f");
        leg->Draw();

        TLatex latex;
        latex.SetNDC(true);
        latex.SetTextSize(0.045);
        latex.DrawLatex(0.10, 0.92, Form("Final ring-included Barlow, v_{%d}", harmonic));

        padBottom->cd();
        hBarlow->SetStats(0);
        hBarlow->GetXaxis()->LabelsOption("v");
        hBarlow->GetXaxis()->SetTitle("eta point");
        hBarlow->GetYaxis()->SetTitle("|#Delta|/#sigma_{Barlow}");
        hBarlow->GetXaxis()->SetTitleSize(0.10);
        hBarlow->GetXaxis()->SetLabelSize(0.08);
        hBarlow->GetYaxis()->SetTitleSize(0.09);
        hBarlow->GetYaxis()->SetLabelSize(0.08);
        hBarlow->GetYaxis()->SetTitleOffset(0.58);
        hBarlow->SetMinimum(0.0);
        hBarlow->SetMaximum(std::max(1.5, hBarlow->GetMaximum() * 1.25));
        hBarlow->SetLineColor(kBlack);
        hBarlow->SetMarkerStyle(20);
        hBarlow->Draw("E1");

        TLine* lineAccept = new TLine(0.5, kBarlowAcceptThreshold, nBins + 0.5, kBarlowAcceptThreshold);
        lineAccept->SetLineStyle(7);
        lineAccept->SetLineWidth(2);
        lineAccept->SetLineColor(kGreen + 2);
        lineAccept->Draw();

        c->SaveAs(Form("./3times2PC/PDFs/Barlow_OO_Chi2_TPCClusters4_FinalRings_v%d.pdf", harmonic));

        delete lineAccept;
        delete leg;
        delete padTop;
        delete padBottom;
        delete c;
        delete hNominalBand;
        delete hNominalInner;
        delete hNominalOuter;
        delete hVariedInner;
        delete hVariedOuter;
    }

    out->Close();
    delete out;

    fNominalFull->Close();
    fNominalInnerC->Close();
    fNominalInnerA->Close();
    fNominalOuterC->Close();
    fNominalOuterA->Close();
    fVariedFull->Close();
    fVariedInnerC->Close();
    fVariedInnerA->Close();
    fVariedOuterC->Close();
    fVariedOuterA->Close();

    delete fNominalFull;
    delete fNominalInnerC;
    delete fNominalInnerA;
    delete fNominalOuterC;
    delete fNominalOuterA;
    delete fVariedFull;
    delete fVariedInnerC;
    delete fVariedInnerA;
    delete fVariedOuterC;
    delete fVariedOuterA;

    std::cout << "[Barlow final-rings] Wrote final ring-included Barlow outputs for O-O vs O-O Chi2_TPCClusters4." << std::endl;
}
