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
TFile* OpenFirstExistingFile(const std::vector<std::string>& candidates, std::string& selectedPath) {
    for (const auto& path : candidates) {
        TFile* file = TFile::Open(path.c_str(), "READ");
        if (file && file->IsOpen() && !file->IsZombie()) {
            selectedPath = path;
            return file;
        }
        if (file) {
            file->Close();
            delete file;
        }
    }
    selectedPath.clear();
    return nullptr;
}

TH1D* LoadCombinedV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Combined");
    return h;
}

TH1D* LoadSideV2(TFile* f) {
    if (!f || !f->IsOpen()) return nullptr;
    TH1D* h = (TH1D*)f->Get("hV2_Side");
    if (!h) h = (TH1D*)f->Get("hV2_Sides");
    if (!h) h = (TH1D*)f->Get("hV2_Combined");
    if (!h) h = (TH1D*)f->Get("hV2");
    return h;
}

TH1D* BuildSidesFromFiles(TFile* fFT0C, TFile* fFT0A, const char* name) {
    TH1D* hFT0C = LoadSideV2(fFT0C);
    TH1D* hFT0A = LoadSideV2(fFT0A);
    if (!hFT0C || !hFT0A) return nullptr;

    TH1D* h = new TH1D(name, "", 2, 0.5, 2.5);
    h->SetDirectory(nullptr);
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

TGraphAsymmErrors* BuildEtaGraphWithRings(const char* name, TH1D* hFull, TH1D* hInnerSides, TH1D* hOuterSides, bool mirrorDetectors) {
    if (!hFull || !hInnerSides || !hOuterSides) return nullptr;

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

    const double ft0cInnerY = hInnerSides->GetBinContent(1);
    const double ft0cInnerErr = hInnerSides->GetBinError(1);
    const double ft0cOuterY = hOuterSides->GetBinContent(1);
    const double ft0cOuterErr = hOuterSides->GetBinError(1);
    const double ft0aInnerY = hInnerSides->GetBinContent(2);
    const double ft0aInnerErr = hInnerSides->GetBinError(2);
    const double ft0aOuterY = hOuterSides->GetBinContent(2);
    const double ft0aOuterErr = hOuterSides->GetBinError(2);

    pushPoint(-3.0, ft0cInnerY, 0.3, 0.3, ft0cInnerErr, ft0cInnerErr);
    pushPoint(-2.4, ft0cOuterY, 0.3, 0.3, ft0cOuterErr, ft0cOuterErr);
    if (mirrorDetectors) {
        pushPoint(2.4, ft0cOuterY, 0.3, 0.3, ft0cOuterErr, ft0cOuterErr);
        pushPoint(3.0, ft0cInnerY, 0.3, 0.3, ft0cInnerErr, ft0cInnerErr);
    }

    pushPoint(3.85, ft0aInnerY, 0.35, 0.35, ft0aInnerErr, ft0aInnerErr);
    pushPoint(4.55, ft0aOuterY, 0.35, 0.35, ft0aOuterErr, ft0aOuterErr);
    if (mirrorDetectors) {
        pushPoint(-4.55, ft0aOuterY, 0.35, 0.35, ft0aOuterErr, ft0aOuterErr);
        pushPoint(-3.85, ft0aInnerY, 0.35, 0.35, ft0aInnerErr, ft0aInnerErr);
    }

    for (int itpc = 0; itpc < 16; ++itpc) {
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

TGraphAsymmErrors* BuildSystemGraph(const char* graphName,
                                    const char* systemLabel,
                                    const std::vector<std::string>& fullCandidates,
                                    const std::vector<std::string>& innerFT0CCandidates,
                                    const std::vector<std::string>& innerFT0ACandidates,
                                    const std::vector<std::string>& outerFT0CCandidates,
                                    const std::vector<std::string>& outerFT0ACandidates,
                                    bool mirrorDetectors) {
    std::string pFull;
    std::string pInnerC;
    std::string pInnerA;
    std::string pOuterC;
    std::string pOuterA;

    TFile* fFull = OpenFirstExistingFile(fullCandidates, pFull);
    TFile* fInnerC = OpenFirstExistingFile(innerFT0CCandidates, pInnerC);
    TFile* fInnerA = OpenFirstExistingFile(innerFT0ACandidates, pInnerA);
    TFile* fOuterC = OpenFirstExistingFile(outerFT0CCandidates, pOuterC);
    TFile* fOuterA = OpenFirstExistingFile(outerFT0ACandidates, pOuterA);

    if (!fFull || !fInnerC || !fInnerA || !fOuterC || !fOuterA) {
        std::cerr << "Missing required files for " << systemLabel << " ring comparison." << std::endl;
        if (fFull) { fFull->Close(); delete fFull; }
        if (fInnerC) { fInnerC->Close(); delete fInnerC; }
        if (fInnerA) { fInnerA->Close(); delete fInnerA; }
        if (fOuterC) { fOuterC->Close(); delete fOuterC; }
        if (fOuterA) { fOuterA->Close(); delete fOuterA; }
        return nullptr;
    }

    std::cout << "Using ring inputs for " << systemLabel << ":\n"
              << "  full   : " << pFull << "\n"
              << "  innerC : " << pInnerC << "\n"
              << "  innerA : " << pInnerA << "\n"
              << "  outerC : " << pOuterC << "\n"
              << "  outerA : " << pOuterA << std::endl;

    TH1D* hFull = LoadCombinedV2(fFull);
    TH1D* hInnerSides = BuildSidesFromFiles(fInnerC, fInnerA, (std::string("hInnerSides_") + systemLabel).c_str());
    TH1D* hOuterSides = BuildSidesFromFiles(fOuterC, fOuterA, (std::string("hOuterSides_") + systemLabel).c_str());
    TGraphAsymmErrors* g = BuildEtaGraphWithRings(graphName, hFull, hInnerSides, hOuterSides, mirrorDetectors);

    delete hInnerSides;
    delete hOuterSides;
    if (fFull) { fFull->Close(); delete fFull; }
    if (fInnerC) { fInnerC->Close(); delete fInnerC; }
    if (fInnerA) { fInnerA->Close(); delete fInnerA; }
    if (fOuterC) { fOuterC->Close(); delete fOuterC; }
    if (fOuterA) { fOuterA->Close(); delete fOuterA; }
    return g;
}
}

void Compare3times2PC_EtaDiff_FourSystems_Nch() {
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    TGraphAsymmErrors* gNeNe = BuildSystemGraph(
        "gNeNe",
        "Ne-Ne",
        {
        "./3times2PC/Vn_LHC25af_pass2_650315_nch10_50_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25af_pass2_innerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25af_pass2_innerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25af_pass2_outerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25af_pass2_outerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        true);

    TGraphAsymmErrors* gOO = BuildSystemGraph(
        "gOO",
        "O-O",
        {
        "./3times2PC/Vn_LHC25ae_pass2_650311_nch10_50_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ae_pass2_innerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ae_pass2_innerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ae_pass2_outerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ae_pass2_outerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        true);

    TGraphAsymmErrors* gPO = BuildSystemGraph(
        "gPO",
        "p-O",
        {
        "./3times2PC/Vn_LHC25ad_pass2_650299_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ad_pass2_innerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC25ad_pass2_outerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        false);

    TGraphAsymmErrors* gPP = BuildSystemGraph(
        "gPP",
        "pp",
        {
        "./3times2PC/Vn_LHC24af_pass1_650588_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC24af_pass1_innerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC24af_pass1_innerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC24af_pass1_outerRing_nch10_50_FT0C_Nch_10_50.root"
        },
        {
        "./3times2PC/Vn_LHC24af_pass1_outerRing_nch10_50_FT0A_Nch_10_50.root"
        },
        true);

    if (!gNeNe || !gOO || !gPO || !gPP) {
        std::cerr << "Cannot build four-system ring comparison: at least one system is missing required ring inputs." << std::endl;
        delete gNeNe;
        delete gOO;
        delete gPO;
        delete gPP;
        return;
    }

    TH1D* hFrame = BuildEtaScaleHistogram("hFrame4SysNch", "");
    if (!hFrame) {
        std::cerr << "Failed to build eta-scale frame or graphs." << std::endl;
        delete gNeNe;
        delete gOO;
        delete gPO;
        delete gPP;
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
    gPP->Draw("PZ same");
    gPO->Draw("PZ same");
    gNeNe->Draw("PZ same");
    gOO->Draw("PZ same");

    TLegend* leg = new TLegend(0.74, 0.74, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gOO,   "O-O",   "lep");
    leg->AddEntry(gNeNe, "Ne-Ne", "lep");
    leg->AddEntry(gPO, "p-O", "lep");
    leg->AddEntry(gPP, "pp", "lep");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.042);
    latex.DrawLatex(0.13, 0.88, "ALICE Preliminary");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.13, 0.825, "3x2PC, 10-50 Nch");
    latex.DrawLatex(0.13, 0.775, "FT0 rings + TPC acceptance (+ratio)");

    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp_Nch.root");
    c->SaveAs("./3times2PC/Compare_v2_EtaDiff_NeNe_vs_OO_vs_pO_vs_pp_Nch.png");
    std::cout << "Saved eta-diff v2 Nch comparison plot (available systems drawn)." << std::endl;

    delete gNeNe;
    delete gOO;
    delete gPO;
    delete gPP;
    delete hFrame;
    delete c;
}
