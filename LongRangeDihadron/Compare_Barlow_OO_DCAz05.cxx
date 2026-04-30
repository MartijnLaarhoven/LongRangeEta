#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
constexpr double kBarlowAcceptThreshold = 1.0;

struct BarlowTask {
    std::string label;
    std::string nominalPath;
    std::string variedPath;
    std::vector<std::string> histNames;
};

TH1D* LoadFirstExistingHistogram(TFile* file, const std::vector<std::string>& names) {
    if (!file || !file->IsOpen() || file->IsZombie()) return nullptr;
    for (const auto& name : names) {
        TObject* obj = file->Get(name.c_str());
        if (!obj) continue;
        TH1D* h = dynamic_cast<TH1D*>(obj);
        if (!h) continue;
        TH1D* clone = static_cast<TH1D*>(h->Clone(Form("%s_clone", h->GetName())));
        clone->SetDirectory(nullptr);
        return clone;
    }
    return nullptr;
}

double ComputeBarlow(double nominal, double nominalErr, double varied, double variedErr) {
    const double diff = std::fabs(varied - nominal);
    const double sigmaVaried2 = variedErr * variedErr;
    const double sigmaNominal2 = nominalErr * nominalErr;
    const double sigmaBarlow2 = std::fabs(sigmaVaried2 - sigmaNominal2);
    const double sigmaIndependent2 = sigmaVaried2 + sigmaNominal2;

    if (sigmaIndependent2 <= 1e-16) return 0.0;

    // Standard Barlow uses |sigma_var^2 - sigma_nom^2|. For nearly equal errors this
    // can become artificially tiny and inflate significance, so fall back to the
    // independent-error denominator in that regime.
    const double relSeparation = sigmaBarlow2 / sigmaIndependent2;
    const bool useIndependent = (relSeparation < 0.20);
    const double den2 = useIndependent ? sigmaIndependent2 : sigmaBarlow2;
    const double den = std::sqrt(den2);
    if (den <= 0.0) return 0.0;
    return diff / den;
}

void StyleHistogram(TH1D* h, int color, int marker) {
    if (!h) return;
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.0);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
}

bool AcceptByBarlow(double barlow) {
    return std::isfinite(barlow) && barlow < kBarlowAcceptThreshold;
}
}

void Compare_Barlow_OO_DCAz05() {
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir("./3times2PC", kTRUE);
    gSystem->mkdir("./3times2PC/PDFs", kTRUE);

    const std::vector<BarlowTask> tasks = {
        {"FULL_3x2", "./3times2PC/Vn_LHC25ae_pass2_648798_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_659116_DCAz05_Cent_0_20.root", {"hV2", "hV2_Combined", "hV2_Sides", "hV2_Side"}},
        {"INNER_RING_FT0C", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_DCAz05_Cent_0_20.root", {"hV2", "hV2_Side", "hV2_Sides", "hV2_Combined"}},
        {"INNER_RING_FT0A", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_DCAz05_Cent_0_20.root", {"hV2", "hV2_Side", "hV2_Sides", "hV2_Combined"}},
        {"OUTER_RING_FT0C", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_DCAz05_Cent_0_20.root", {"hV2", "hV2_Side", "hV2_Sides", "hV2_Combined"}},
        {"OUTER_RING_FT0A", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_DCAz05_Cent_0_20.root", {"hV2", "hV2_Side", "hV2_Sides", "hV2_Combined"}}
    };

    const char* outputPath = "./3times2PC/Barlow_OO_vs_OO_DCAz05_Cent_0_20.root";
    TFile* outputFile = nullptr;
    bool wroteAny = false;
    int processedTasks = 0;
    std::vector<std::string> missingTasks;

    for (const auto& task : tasks) {
        TFile* fNominal = TFile::Open(task.nominalPath.c_str(), "READ");
        TFile* fVaried = TFile::Open(task.variedPath.c_str(), "READ");

        if (!fNominal || !fNominal->IsOpen() || fNominal->IsZombie()) {
            std::cerr << "[Barlow] Missing nominal file for " << task.label << ": " << task.nominalPath << std::endl;
            missingTasks.push_back(task.label);
            if (fNominal) { fNominal->Close(); delete fNominal; }
            if (fVaried) { fVaried->Close(); delete fVaried; }
            continue;
        }
        if (!fVaried || !fVaried->IsOpen() || fVaried->IsZombie()) {
            std::cerr << "[Barlow] Missing varied file for " << task.label << ": " << task.variedPath << std::endl;
            missingTasks.push_back(task.label);
            fNominal->Close();
            delete fNominal;
            if (fVaried) { fVaried->Close(); delete fVaried; }
            continue;
        }

        TH1D* hNominal = LoadFirstExistingHistogram(fNominal, task.histNames);
        TH1D* hVaried = LoadFirstExistingHistogram(fVaried, task.histNames);

        if (!hNominal || !hVaried) {
            std::cerr << "[Barlow] Missing histogram for " << task.label << std::endl;
            missingTasks.push_back(task.label);
            delete hNominal;
            delete hVaried;
            fNominal->Close();
            fVaried->Close();
            delete fNominal;
            delete fVaried;
            continue;
        }

        const int nBins = std::min(hNominal->GetNbinsX(), hVaried->GetNbinsX());
        TH1D* hDiff = new TH1D(Form("hDiff_%s", task.label.c_str()), Form("%s: varied-nominal;bin;#Delta v_{2}", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hBarlow = new TH1D(Form("hBarlow_%s", task.label.c_str()), Form("%s: Barlow significance;bin;|#Delta| / #sigma_{Barlow}", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hSystematicAbs = new TH1D(Form("hSystematicAbs_%s", task.label.c_str()), Form("%s: accepted systematic uncertainty;bin;|#Delta v_{2}|", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hSystematicRel = new TH1D(Form("hSystematicRel_%s", task.label.c_str()), Form("%s: accepted systematic uncertainty;bin;100#times|#Delta v_{2}|/v_{2}", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hTotalAbs = new TH1D(Form("hTotalAbs_%s", task.label.c_str()), Form("%s: total uncertainty;bin;#sigma_{tot}(v_{2})", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hTotalRel = new TH1D(Form("hTotalRel_%s", task.label.c_str()), Form("%s: total uncertainty;bin;100#times#sigma_{tot}(v_{2})/v_{2}", task.label.c_str()), nBins, 0.5, nBins + 0.5);
        TH1D* hBarlowDecision = new TH1D(Form("hBarlowDecision_%s", task.label.c_str()), Form("%s: Barlow decision;bin;1=accepted, 0=overruled", task.label.c_str()), nBins, 0.5, nBins + 0.5);

        for (int bin = 1; bin <= nBins; ++bin) {
            const double nominal = hNominal->GetBinContent(bin);
            const double nominalErr = hNominal->GetBinError(bin);
            const double varied = hVaried->GetBinContent(bin);
            const double variedErr = hVaried->GetBinError(bin);
            const double diff = varied - nominal;
            const double diffErr = std::sqrt(nominalErr * nominalErr + variedErr * variedErr);
            const double barlow = ComputeBarlow(nominal, nominalErr, varied, variedErr);
            const bool accepted = AcceptByBarlow(barlow);
            const double systAbs = accepted ? std::fabs(diff) : 0.0;
            const double systRel = (accepted && std::fabs(nominal) > 1e-16) ? 100.0 * systAbs / std::fabs(nominal) : 0.0;
            const double totalAbs = std::sqrt(nominalErr * nominalErr + systAbs * systAbs);
            const double totalRel = (std::fabs(nominal) > 1e-16) ? 100.0 * totalAbs / std::fabs(nominal) : 0.0;

            hDiff->SetBinContent(bin, diff);
            hDiff->SetBinError(bin, diffErr);
            hBarlow->SetBinContent(bin, barlow);
            hBarlow->SetBinError(bin, 0.0);
            hSystematicAbs->SetBinContent(bin, systAbs);
            hSystematicAbs->SetBinError(bin, 0.0);
            hSystematicRel->SetBinContent(bin, systRel);
            hSystematicRel->SetBinError(bin, 0.0);
            hTotalAbs->SetBinContent(bin, totalAbs);
            hTotalAbs->SetBinError(bin, 0.0);
            hTotalRel->SetBinContent(bin, totalRel);
            hTotalRel->SetBinError(bin, 0.0);
            hBarlowDecision->SetBinContent(bin, accepted ? 1.0 : 0.0);
            hBarlowDecision->SetBinError(bin, 0.0);

            std::cout << "[Barlow] " << task.label
                      << " bin=" << bin
                      << " nominal=" << nominal << " +/- " << nominalErr
                      << " varied=" << varied << " +/- " << variedErr
                      << " diff=" << diff
                      << " barlow=" << barlow
                      << " decision=" << (accepted ? "ACCEPT" : "OVERRULE")
                      << " systAbs=" << systAbs
                      << " systRel=" << systRel << "%"
                      << " totalAbs=" << totalAbs
                      << " totalRel=" << totalRel << "%"
                      << std::endl;
        }

        StyleHistogram(hDiff, kBlue + 1, 20);
        StyleHistogram(hBarlow, kRed + 1, 24);
        StyleHistogram(hSystematicAbs, kGreen + 2, 20);
        StyleHistogram(hSystematicRel, kGreen + 2, 20);
        StyleHistogram(hTotalAbs, kOrange + 7, 20);
        StyleHistogram(hTotalRel, kOrange + 7, 20);
        StyleHistogram(hBarlowDecision, kBlack, 20);

        if (!outputFile) {
            outputFile = TFile::Open(outputPath, "RECREATE");
            if (!outputFile || !outputFile->IsOpen() || outputFile->IsZombie()) {
                std::cerr << "[Barlow] Failed to create output file." << std::endl;
                delete hDiff;
                delete hBarlow;
                delete hNominal;
                delete hVaried;
                fNominal->Close();
                fVaried->Close();
                delete fNominal;
                delete fVaried;
                if (outputFile) {
                    outputFile->Close();
                    delete outputFile;
                }
                return;
            }
        }

        outputFile->cd();
        hDiff->Write();
        hBarlow->Write();
        hSystematicAbs->Write();
        hSystematicRel->Write();
        hTotalAbs->Write();
        hTotalRel->Write();
        hBarlowDecision->Write();
        wroteAny = true;
        ++processedTasks;

        TH1D* hNominalBand = static_cast<TH1D*>(hNominal->Clone(Form("hNominalBand_%s", task.label.c_str())));
        hNominalBand->SetDirectory(nullptr);
        for (int bin = 1; bin <= nBins; ++bin) {
            hNominalBand->SetBinError(bin, hTotalAbs->GetBinContent(bin));
        }

        TCanvas* cOverlay = new TCanvas(Form("cBarlowOverlay_%s", task.label.c_str()), task.label.c_str(), 950, 850);
        TPad* padOverlayTop = new TPad(Form("padOverlayTop_%s", task.label.c_str()), "", 0.0, 0.35, 1.0, 1.0);
        TPad* padOverlayBottom = new TPad(Form("padOverlayBottom_%s", task.label.c_str()), "", 0.0, 0.0, 1.0, 0.35);
        padOverlayTop->SetBottomMargin(0.02);
        padOverlayTop->SetLeftMargin(0.12);
        padOverlayTop->SetRightMargin(0.04);
        padOverlayBottom->SetTopMargin(0.04);
        padOverlayBottom->SetBottomMargin(0.25);
        padOverlayBottom->SetLeftMargin(0.12);
        padOverlayBottom->SetRightMargin(0.04);
        padOverlayTop->Draw();
        padOverlayBottom->Draw();

        double overlayMax = 0.0;
        for (int bin = 1; bin <= nBins; ++bin) {
            overlayMax = std::max(overlayMax, hNominal->GetBinContent(bin) + hNominalBand->GetBinError(bin));
            overlayMax = std::max(overlayMax, hVaried->GetBinContent(bin) + hVaried->GetBinError(bin));
        }
        if (overlayMax <= 0.0) overlayMax = 0.1;

        padOverlayTop->cd();
        hNominalBand->SetStats(0);
        hNominalBand->GetXaxis()->SetLabelSize(0.0);
        hNominalBand->GetYaxis()->SetTitleSize(0.06);
        hNominalBand->GetYaxis()->SetLabelSize(0.05);
        hNominalBand->GetYaxis()->SetRangeUser(0.0, overlayMax * 1.25);
        hNominalBand->SetFillColor(kGray + 1);
        hNominalBand->SetFillStyle(3004);
        hNominalBand->SetLineColor(kGray + 2);
        hNominalBand->SetMarkerSize(0);
        hNominalBand->Draw("E2");
        hNominal->SetLineColor(kBlue + 1);
        hNominal->SetMarkerColor(kBlue + 1);
        hNominal->Draw("HIST SAME");
        hVaried->SetLineColor(kRed + 1);
        hVaried->SetMarkerColor(kRed + 1);
        hVaried->Draw("HIST SAME");

        TLegend* overlayLegend = new TLegend(0.14, 0.74, 0.47, 0.90);
        overlayLegend->SetBorderSize(0);
        overlayLegend->SetFillStyle(0);
        overlayLegend->SetTextFont(42);
        overlayLegend->SetTextSize(0.032);
        overlayLegend->AddEntry(hNominal, "Nominal data", "l");
        overlayLegend->AddEntry(hVaried, "Systematic variation", "l");
        overlayLegend->AddEntry(hNominalBand, Form("Total uncertainty band (stat #oplus syst, Barlow < %.1f)", kBarlowAcceptThreshold), "f");
        overlayLegend->Draw();

        TLatex overlayLatex;
        overlayLatex.SetNDC(true);
        overlayLatex.SetTextSize(0.05);
        overlayLatex.DrawLatex(0.13, 0.92, Form("O-O 0-20%%: nominal vs DCA_{z}<0.5 cm (%s)", task.label.c_str()));

        padOverlayBottom->cd();
        hBarlow->SetMinimum(0.0);
        hBarlow->SetMaximum(std::max(1.5, hBarlow->GetMaximum() * 1.20));
        hBarlow->GetXaxis()->SetTitleSize(0.10);
        hBarlow->GetXaxis()->SetLabelSize(0.09);
        hBarlow->GetYaxis()->SetTitleSize(0.09);
        hBarlow->GetYaxis()->SetLabelSize(0.08);
        hBarlow->GetYaxis()->SetTitleOffset(0.60);
        hBarlow->Draw("E1");

        TLine* acceptLine = new TLine(0.5, kBarlowAcceptThreshold, nBins + 0.5, kBarlowAcceptThreshold);
        acceptLine->SetLineColor(kGreen + 2);
        acceptLine->SetLineStyle(7);
        acceptLine->SetLineWidth(2);
        acceptLine->Draw();

        TLine* zeroDecision = new TLine(0.5, 0.0, nBins + 0.5, 0.0);
        zeroDecision->SetLineStyle(2);
        zeroDecision->SetLineColor(kGray + 2);
        zeroDecision->Draw();

        TLine* oneDecision = new TLine(0.5, 1.0, nBins + 0.5, 1.0);
        oneDecision->SetLineStyle(2);
        oneDecision->SetLineColor(kBlue + 1);
        oneDecision->Draw();

        overlayLatex.SetTextSize(0.04);
        overlayLatex.DrawLatex(0.13, 0.86, "Accepted bins keep |#Delta| as the systematic uncertainty");
        overlayLatex.DrawLatex(0.13, 0.81, "Total uncertainty = stat #oplus accepted syst");

        cOverlay->SaveAs(Form("./3times2PC/PDFs/SystematicOverlay_OO_DCAz05_%s.pdf", task.label.c_str()));

        TCanvas* c = new TCanvas(Form("cBarlow_%s", task.label.c_str()), task.label.c_str(), 950, 850);
        TPad* padTop = new TPad(Form("padTop_%s", task.label.c_str()), "", 0.0, 0.35, 1.0, 1.0);
        TPad* padBottom = new TPad(Form("padBottom_%s", task.label.c_str()), "", 0.0, 0.0, 1.0, 0.35);
        padTop->SetBottomMargin(0.02);
        padTop->SetLeftMargin(0.12);
        padTop->SetRightMargin(0.04);
        padBottom->SetTopMargin(0.04);
        padBottom->SetBottomMargin(0.25);
        padBottom->SetLeftMargin(0.12);
        padBottom->SetRightMargin(0.04);
        padTop->Draw();
        padBottom->Draw();

        padTop->cd();
        hDiff->GetXaxis()->SetLabelSize(0.0);
        hDiff->GetYaxis()->SetTitleSize(0.06);
        hDiff->GetYaxis()->SetLabelSize(0.05);
        hDiff->Draw("E1");
        TLine* zeroLine = new TLine(0.5, 0.0, nBins + 0.5, 0.0);
        zeroLine->SetLineStyle(2);
        zeroLine->SetLineColor(kGray + 2);
        zeroLine->Draw();

        TLatex latex;
        latex.SetNDC(true);
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.13, 0.90, Form("O-O 0-20%%: nominal vs DCA_{z}<0.5 cm (%s)", task.label.c_str()));

        padBottom->cd();
        hBarlow->SetMinimum(0.0);
        hBarlow->SetMaximum(std::max(3.5, hBarlow->GetMaximum() * 1.25));
        hBarlow->GetXaxis()->SetTitleSize(0.10);
        hBarlow->GetXaxis()->SetLabelSize(0.09);
        hBarlow->GetYaxis()->SetTitleSize(0.09);
        hBarlow->GetYaxis()->SetLabelSize(0.08);
        hBarlow->GetYaxis()->SetTitleOffset(0.60);
        hBarlow->Draw("E1");

        TLine* line2sigma = new TLine(0.5, 2.0, nBins + 0.5, 2.0);
        line2sigma->SetLineColor(kOrange + 7);
        line2sigma->SetLineStyle(7);
        line2sigma->SetLineWidth(2);
        line2sigma->Draw();

        TLine* line3sigma = new TLine(0.5, 3.0, nBins + 0.5, 3.0);
        line3sigma->SetLineColor(kRed + 1);
        line3sigma->SetLineStyle(7);
        line3sigma->SetLineWidth(2);
        line3sigma->Draw();

        c->SaveAs(Form("./3times2PC/PDFs/Barlow_OO_DCAz05_%s.pdf", task.label.c_str()));

        delete line3sigma;
        delete line2sigma;
        delete zeroLine;
        delete padTop;
        delete padBottom;
        delete c;
        delete acceptLine;
        delete zeroDecision;
        delete oneDecision;
        delete overlayLegend;
        delete padOverlayTop;
        delete padOverlayBottom;
        delete cOverlay;
        delete hNominalBand;
        delete hDiff;
        delete hBarlow;
        delete hSystematicAbs;
        delete hSystematicRel;
        delete hTotalAbs;
        delete hTotalRel;
        delete hBarlowDecision;
        delete hNominal;
        delete hVaried;

        fNominal->Close();
        fVaried->Close();
        delete fNominal;
        delete fVaried;
    }

    if (outputFile) {
        outputFile->Close();
        delete outputFile;
    }

    if (!wroteAny) {
        gSystem->Unlink(outputPath);
        std::cerr << "[Barlow] No valid nominal/varied input pairs found. "
                  << "No output file written." << std::endl;
        return;
    }

    std::cout << "[Barlow] Processed " << processedTasks << " / " << tasks.size() << " tasks." << std::endl;
    if (!missingTasks.empty()) {
        std::cerr << "[Barlow] WARNING: skipped tasks:";
        for (const auto& label : missingTasks) {
            std::cerr << " " << label;
        }
        std::cerr << std::endl;
    }

    std::cout << "[Barlow] Wrote: " << outputPath << std::endl;
}
