#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {
TFile* OpenFirstExistingFile(const std::vector<const char*>& candidates) {
    for (const auto* path : candidates) {
        if (gSystem && gSystem->AccessPathName(path)) {
            continue;
        }
        TFile* file = TFile::Open(path, "READ");
        if (file && file->IsOpen()) {
            std::cout << "Using file: " << path << std::endl;
            return file;
        }
        delete file;
    }
    return nullptr;
}

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

TGraphAsymmErrors* BuildEtaGraph(const char* name, TH1D* hFull, bool mirrorDetectors) {
    if (!hFull) return nullptr;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> exLow;
    std::vector<double> exHigh;
    std::vector<double> eyLow;
    std::vector<double> eyHigh;

    auto pushPoint = [&](double xval, double yval, double xerrLow, double xerrHigh, double yerrLow, double yerrHigh) {
        x.push_back(xval);
        y.push_back(yval);
        exLow.push_back(xerrLow);
        exHigh.push_back(xerrHigh);
        eyLow.push_back(yerrLow);
        eyHigh.push_back(yerrHigh);
    };

    const int nBins = hFull->GetNbinsX();
    if (nBins < 4) return nullptr;

    const double ft0cY = hFull->GetBinContent(1);
    const double ft0cErr = hFull->GetBinError(1);
    const double ft0aY = hFull->GetBinContent(nBins);
    const double ft0aErr = hFull->GetBinError(nBins);

    if (mirrorDetectors) {
        pushPoint(-2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr);
        pushPoint(2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr);
        pushPoint(-4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr);
        pushPoint(4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr);
    } else {
        pushPoint(-2.7, ft0cY, 0.6, 0.6, ft0cErr, ft0cErr);
        pushPoint(4.2, ft0aY, 0.7, 0.7, ft0aErr, ft0aErr);
    }

    const int nTpcBins = nBins - 2;
    for (int itpc = 0; itpc < nTpcBins; ++itpc) {
        const int sourceBin = itpc + 2;
        pushPoint(-0.75 + 0.1 * itpc,
                  hFull->GetBinContent(sourceBin),
                  0.05,
                  0.05,
                  hFull->GetBinError(sourceBin),
                  hFull->GetBinError(sourceBin));
    }

    TGraphAsymmErrors* graph = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), y.data(), exLow.data(), exHigh.data(), eyLow.data(), eyHigh.data());
    graph->SetName(name);
    return graph;
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
}

void Compare3times2PC_EtaDiff_FourSystems_Nch() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    TFile* fNeNe = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25af_pass2_650315_Nch_10_50.root",
        "./3times2PC/Vn_LHC25af_pass2_650315_nch10_50_Mult_10_50.root",
        "./3times2PC/Vn_LHC25af_pass2_650315_Mult_10_50.root"
    });
    TFile* fOO = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ae_pass2_653254_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ae_pass2_653254_Mult_10_50.root"
    });
    TFile* fPO = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC25ad_pass2_650299_Nch_10_50.root",
        "./3times2PC/Vn_LHC25ad_pass2_650299_Mult_10_50.root"
    });
    TFile* fPP = OpenFirstExistingFile({
        "./3times2PC/Vn_LHC24af_pass1_650588_Nch_10_50.root",
        "./3times2PC/Vn_LHC24af_pass1_650588_Mult_10_50.root"
    });

    if (!fNeNe || !fNeNe->IsOpen()) { std::cerr << "Cannot open Ne-Ne file!" << std::endl; return; }
    if (!fOO   || !fOO->IsOpen())   { std::cerr << "Cannot open O-O file!" << std::endl; return; }
    if (!fPO   || !fPO->IsOpen())   { std::cerr << "Warning: cannot open p-O file, skip p-O." << std::endl; }
    if (!fPP   || !fPP->IsOpen())   { std::cerr << "Warning: cannot open pp file, skip pp." << std::endl; }

    TH1D* hNeNe = LoadFullV2(fNeNe);
    TH1D* hOO   = LoadFullV2(fOO);
    TH1D* hPO   = LoadFullV2(fPO);
    TH1D* hPP   = LoadFullV2(fPP);

    if (!hNeNe || !hOO) {
        std::cerr << "Missing mandatory Ne-Ne or O-O v2 histogram." << std::endl;
        return;
    }
    if (fPO && fPO->IsOpen() && !hPO) {
        std::cerr << "Warning: p-O file opened but v2 histogram missing, skip p-O." << std::endl;
    }
    if (fPP && fPP->IsOpen() && !hPP) {
        std::cerr << "Warning: pp file opened but v2 histogram missing, skip pp." << std::endl;
    }

    TH1D* hFrame = BuildEtaScaleHistogram("hFrame4SysNch", "");
    TGraphAsymmErrors* gNeNe = BuildEtaGraph("gNeNe", hNeNe, true);
    TGraphAsymmErrors* gOO   = BuildEtaGraph("gOO",   hOO,   true);
    TGraphAsymmErrors* gPO   = hPO ? BuildEtaGraph("gPO", hPO, false) : nullptr;
    TGraphAsymmErrors* gPP   = hPP ? BuildEtaGraph("gPP", hPP, true) : nullptr;
    if (!hFrame || !gNeNe || !gOO) {
        std::cerr << "Failed to build eta-scale frame or graphs." << std::endl;
        return;
    }

    // Keep the same color assignment as the current 4-system comparison.
    Style(gOO,   kRed + 1,   21);
    Style(gNeNe, kGreen + 2, 20);
    Style(gPO,   kBlue + 1,  24);
    Style(gPP,   kBlack,     25);

    double maxVal = 0.0;
    maxVal = std::max(maxVal, GetGraphYMax(gNeNe));
    maxVal = std::max(maxVal, GetGraphYMax(gOO));
    maxVal = std::max(maxVal, GetGraphYMax(gPO));
    maxVal = std::max(maxVal, GetGraphYMax(gPP));
    if (maxVal <= 0.0) maxVal = 0.1;

    TCanvas* c = new TCanvas("cEtaDiff4SysNch", "3x2PC Eta-diff v2: four Nch systems", 1150, 700);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.12);
    c->SetTicks(1, 1);

    hFrame->GetYaxis()->SetRangeUser(0.0, maxVal * 1.25);
    hFrame->Draw("AXIS");
    if (gPP) gPP->Draw("PZ same");
    if (gPO) gPO->Draw("PZ same");
    gNeNe->Draw("PZ same");
    gOO->Draw("PZ same");

    TLegend* leg = new TLegend(0.74, 0.74, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gOO,   "O-O",   "lep");
    leg->AddEntry(gNeNe, "Ne-Ne", "lep");
    if (gPO) leg->AddEntry(gPO, "p-O", "lep");
    if (gPP) leg->AddEntry(gPP, "pp", "lep");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.13, 0.88, "ALICE Preliminary");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.13, 0.825, "3x2PC, 10-50 Nch");
    latex.DrawLatex(0.13, 0.775, "Full FT0 ranges + TPC acceptance (+ratio)");

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp_Nch.root");
    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp_Nch.png");
    std::cout << "Saved eta-diff v2 Nch comparison plot (available systems drawn)." << std::endl;

    if (fNeNe) fNeNe->Close();
    if (fOO) fOO->Close();
    if (fPO) fPO->Close();
    if (fPP) fPP->Close();
    delete fNeNe;
    delete fOO;
    delete fPO;
    delete fPP;
    delete c;
}
