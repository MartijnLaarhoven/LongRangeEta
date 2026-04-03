#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

struct SystemInput {
    std::string label;
    std::string dataset;
    std::string filePath;
    int color;
    int marker;
};

void CheckTemplateFit_TPCFT0C_Trend() {
    gStyle->SetOptStat(0);

    std::vector<SystemInput> systems = {
        {"Ne-Ne", "LHC25af_pass2_637596", "./TemplateFit/EtaDiff/VnDelta_LHC25af_pass2_637596_Cent_0_20_TPC_FT0C.root", kRed + 1, 20},
        {"O-O",   "LHC25ae_pass2_644429", "./TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_644429_Cent_0_20_TPC_FT0C.root", kBlue + 1, 21},
        {"p-O",   "LHC25ad_pass2_644389", "./TemplateFit/EtaDiff/VnDelta_LHC25ad_pass2_644389_Cent_0_20_TPC_FT0C.root", kGreen + 2, 22},
        {"pp",    "LHC24af_pass1_644663", "./TemplateFit/EtaDiff/VnDelta_LHC24af_pass1_644663_Cent_0_20_TPC_FT0C.root", kBlack, 24}
    };

    std::vector<TH1D*> hists;
    hists.reserve(systems.size());

    for (const auto& s : systems) {
        TFile* f = TFile::Open(s.filePath.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[ERROR] Cannot open file: " << s.filePath << std::endl;
            return;
        }

        TH1D* h = dynamic_cast<TH1D*>(f->Get("hV2"));
        if (!h) {
            std::cerr << "[ERROR] Missing hV2 in file: " << s.filePath << std::endl;
            f->Close();
            delete f;
            return;
        }

        TH1D* cloned = dynamic_cast<TH1D*>(h->Clone(("hV2_" + s.dataset).c_str()));
        cloned->SetDirectory(nullptr);
        f->Close();
        delete f;

        cloned->SetLineColor(s.color);
        cloned->SetMarkerColor(s.color);
        cloned->SetMarkerStyle(s.marker);
        cloned->SetMarkerSize(1.0);
        cloned->SetLineWidth(2);

        hists.push_back(cloned);
    }

    const int nbins = hists[0]->GetNbinsX();

    std::ofstream report("./TemplateFit/EtaDiff/TPC_FT0C_TrendReport.txt");
    if (!report.is_open()) {
        std::cerr << "[ERROR] Cannot write report file." << std::endl;
        return;
    }

    report << "TemplateFit TPC-FT0C v2 trend report (EtaDiff, Cent 0-20)\n\n";
    report << std::left << std::setw(16) << "EtaRange";
    for (const auto& s : systems) {
        report << std::setw(24) << (s.label + " v2 +/- err");
    }
    report << "\n";

    std::cout << "\n=== TemplateFit TPC-FT0C v2 (per eta bin) ===\n";
    std::cout << std::left << std::setw(16) << "EtaRange";
    for (const auto& s : systems) {
        std::cout << std::setw(24) << (s.label + " v2 +/- err");
    }
    std::cout << "\n";

    for (int ibin = 1; ibin <= nbins; ++ibin) {
        const double etaLow = hists[0]->GetXaxis()->GetBinLowEdge(ibin);
        const double etaHigh = hists[0]->GetXaxis()->GetBinUpEdge(ibin);

        char etaLabel[64];
        snprintf(etaLabel, sizeof(etaLabel), "[%.1f, %.1f]", etaLow, etaHigh);

        report << std::left << std::setw(16) << etaLabel;
        std::cout << std::left << std::setw(16) << etaLabel;

        for (size_t isys = 0; isys < hists.size(); ++isys) {
            const double v = hists[isys]->GetBinContent(ibin);
            const double e = hists[isys]->GetBinError(ibin);

            char valueText[64];
            snprintf(valueText, sizeof(valueText), "%.6f +/- %.6f", v, e);
            report << std::setw(24) << valueText;
            std::cout << std::setw(24) << valueText;
        }

        report << "\n";
        std::cout << "\n";
    }

    report << "\nLinear slope check: fit v2(eta) with pol1 over full eta range.\n";
    std::cout << "\n=== Linear slope check (pol1 fit) ===\n";

    for (size_t isys = 0; isys < hists.size(); ++isys) {
        TF1 fitFunc("fitFunc", "pol1", -0.8, 0.8);
        hists[isys]->Fit(&fitFunc, "Q0");
        const double slope = fitFunc.GetParameter(1);
        const double slopeErr = fitFunc.GetParError(1);

        report << systems[isys].label << " (" << systems[isys].dataset << "): slope = "
               << slope << " +/- " << slopeErr << "\n";
        std::cout << systems[isys].label << " (" << systems[isys].dataset << "): slope = "
                  << slope << " +/- " << slopeErr << "\n";
    }

    report.close();

    TCanvas* c = new TCanvas("cTPCFT0CTrend", "TemplateFit TPC-FT0C Trend", 1000, 700);
    c->SetMargin(0.12, 0.04, 0.12, 0.04);

    double yMin = 1e9;
    double yMax = -1e9;
    for (auto* h : hists) {
        for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
            const double v = h->GetBinContent(ibin);
            const double e = h->GetBinError(ibin);
            yMin = std::min(yMin, v - e);
            yMax = std::max(yMax, v + e);
        }
    }
    yMin *= 0.95;
    yMax *= 1.05;

    hists[0]->SetTitle("TemplateFit EtaDiff: TPC-FT0C v_{2} trend;#eta^{trig} bin;v_{2}");
    hists[0]->GetYaxis()->SetRangeUser(yMin, yMax);
    hists[0]->Draw("E1");
    for (size_t i = 1; i < hists.size(); ++i) {
        hists[i]->Draw("E1 SAME");
    }

    TLegend* leg = new TLegend(0.14, 0.68, 0.45, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (size_t i = 0; i < systems.size(); ++i) {
        leg->AddEntry(hists[i], (systems[i].label + " (" + systems[i].dataset + ")").c_str(), "lep");
    }
    leg->Draw();

    TLatex txt;
    txt.SetNDC();
    txt.SetTextSize(0.032);
    txt.DrawLatex(0.56, 0.88, "Cent 0-20, TemplateFit EtaDiff");
    txt.DrawLatex(0.56, 0.83, "Input: hV2 from TemplateFit/EtaDiff");

    c->SaveAs("./TemplateFit/EtaDiff/TPC_FT0C_TrendCompare.root");

    std::cout << "\nSaved: ./TemplateFit/EtaDiff/TPC_FT0C_TrendReport.txt\n";
    std::cout << "Saved: ./TemplateFit/EtaDiff/TPC_FT0C_TrendCompare.root\n";

    for (auto* h : hists) {
        delete h;
    }
    delete leg;
    delete c;
}
