#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {
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
        // Symmetric systems: mirror FT0C and FT0A points around eta=0.
        pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr);
        pushPoint(2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr);
        pushPoint(-4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr);
        pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr);
    } else {
        // Asymmetric p-O: keep only the physical detector-side points.
        pushPoint(-2.4, ft0cY, 0.3, 0.3, ft0cErr, ft0cErr);
        pushPoint(4.55, ft0aY, 0.35, 0.35, ft0aErr, ft0aErr);
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

void Compare3times2PC_EtaDiff_FourSystems() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const char* fileNeNe = "./3times2PC/Vn_LHC25af_pass2_645746_Cent_0_20.root";
    const char* fileOO   = "./3times2PC/Vn_LHC25ae_pass2_645657_Cent_0_20.root";
    const char* filePO   = "./3times2PC/Vn_LHC25ad_pass2_644389_Cent_0_20.root";
    const char* filePP   = "./3times2PC/Vn_LHC24af_pass1_644663_Cent_0_20.root";

    TFile* fNeNe = TFile::Open(fileNeNe, "READ");
    TFile* fOO   = TFile::Open(fileOO,   "READ");
    TFile* fPO   = TFile::Open(filePO,   "READ");
    TFile* fPP   = TFile::Open(filePP,   "READ");

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

    TH1D* hFrame = BuildEtaScaleHistogram("hFrame4Sys", "");
    TGraphAsymmErrors* gNeNe = BuildEtaGraph("gNeNe", hNeNe, true);
    TGraphAsymmErrors* gOO   = BuildEtaGraph("gOO",   hOO,   true);
    TGraphAsymmErrors* gPO   = hPO ? BuildEtaGraph("gPO", hPO, false) : nullptr;
    TGraphAsymmErrors* gPP   = hPP ? BuildEtaGraph("gPP", hPP, true) : nullptr;
    if (!hFrame || !gNeNe || !gOO) {
        std::cerr << "Failed to build eta-scale frame or graphs." << std::endl;
        return;
    }

    // Requested color mapping: O-O red, Ne-Ne green, p-O blue, pp black.
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

    TCanvas* c = new TCanvas("cEtaDiff4Sys", "3x2PC Eta-diff v2: four systems", 1150, 700);
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
    latex.DrawLatex(0.13, 0.91, "ALICE Preliminary");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.13, 0.855, "3x2PC, 0-20%");
    latex.DrawLatex(0.13, 0.805, "Full FT0 ranges + TPC acceptance");

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp.root");
    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp.png");
    std::cout << "Saved eta-diff v2 comparison plot (available systems drawn)." << std::endl;

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
