#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
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
    const double diff = varied - nominal;
    const double den2 = std::fabs(variedErr * variedErr - nominalErr * nominalErr);
    const double denAlt2 = variedErr * variedErr + nominalErr * nominalErr;
    const double den = (den2 > 1e-16) ? std::sqrt(den2) : ((denAlt2 > 1e-16) ? std::sqrt(denAlt2) : 0.0);
    if (den <= 0.0) return 0.0;
    return std::fabs(diff) / den;
}

void StyleHistogram(TH1D* h, int color, int marker) {
    if (!h) return;
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.0);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
}
}

void Compare_Barlow_OO_DCAz05() {
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir("./3times2PC", kTRUE);
    gSystem->mkdir("./3times2PC/PDFs", kTRUE);

    const std::vector<BarlowTask> tasks = {
        {"FULL_3x2", "./3times2PC/Vn_LHC25ae_pass2_648798_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_659116_DCAz05_Cent_0_20.root", {"hV2_Combined", "hV2"}},
        {"INNER_RING_FT0C", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0C_DCAz05_Cent_0_20.root", {"hV2_Side", "hV2_Sides", "hV2_Combined", "hV2"}},
        {"INNER_RING_FT0A", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_innerRing_FT0A_DCAz05_Cent_0_20.root", {"hV2_Side", "hV2_Sides", "hV2_Combined", "hV2"}},
        {"OUTER_RING_FT0C", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0C_DCAz05_Cent_0_20.root", {"hV2_Side", "hV2_Sides", "hV2_Combined", "hV2"}},
        {"OUTER_RING_FT0A", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_Cent_0_20.root", "./3times2PC/Vn_LHC25ae_pass2_outerRing_FT0A_DCAz05_Cent_0_20.root", {"hV2_Side", "hV2_Sides", "hV2_Combined", "hV2"}}
    };

    const char* outputPath = "./3times2PC/Barlow_OO_vs_OO_DCAz05_Cent_0_20.root";
    TFile* outputFile = nullptr;
    bool wroteAny = false;

    for (const auto& task : tasks) {
        TFile* fNominal = TFile::Open(task.nominalPath.c_str(), "READ");
        TFile* fVaried = TFile::Open(task.variedPath.c_str(), "READ");

        if (!fNominal || !fNominal->IsOpen() || fNominal->IsZombie()) {
            std::cerr << "[Barlow] Missing nominal file for " << task.label << ": " << task.nominalPath << std::endl;
            if (fNominal) { fNominal->Close(); delete fNominal; }
            if (fVaried) { fVaried->Close(); delete fVaried; }
            continue;
        }
        if (!fVaried || !fVaried->IsOpen() || fVaried->IsZombie()) {
            std::cerr << "[Barlow] Missing varied file for " << task.label << ": " << task.variedPath << std::endl;
            fNominal->Close();
            delete fNominal;
            if (fVaried) { fVaried->Close(); delete fVaried; }
            continue;
        }

        TH1D* hNominal = LoadFirstExistingHistogram(fNominal, task.histNames);
        TH1D* hVaried = LoadFirstExistingHistogram(fVaried, task.histNames);

        if (!hNominal || !hVaried) {
            std::cerr << "[Barlow] Missing histogram for " << task.label << std::endl;
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

        for (int bin = 1; bin <= nBins; ++bin) {
            const double nominal = hNominal->GetBinContent(bin);
            const double nominalErr = hNominal->GetBinError(bin);
            const double varied = hVaried->GetBinContent(bin);
            const double variedErr = hVaried->GetBinError(bin);
            const double diff = varied - nominal;
            const double diffErr = std::sqrt(nominalErr * nominalErr + variedErr * variedErr);
            const double barlow = ComputeBarlow(nominal, nominalErr, varied, variedErr);

            hDiff->SetBinContent(bin, diff);
            hDiff->SetBinError(bin, diffErr);
            hBarlow->SetBinContent(bin, barlow);
            hBarlow->SetBinError(bin, 0.0);

            std::cout << "[Barlow] " << task.label
                      << " bin=" << bin
                      << " nominal=" << nominal << " +/- " << nominalErr
                      << " varied=" << varied << " +/- " << variedErr
                      << " diff=" << diff
                      << " barlow=" << barlow
                      << std::endl;
        }

        StyleHistogram(hDiff, kBlue + 1, 20);
        StyleHistogram(hBarlow, kRed + 1, 24);

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
        wroteAny = true;

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
        delete hDiff;
        delete hBarlow;
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

    std::cout << "[Barlow] Wrote: " << outputPath << std::endl;
}
