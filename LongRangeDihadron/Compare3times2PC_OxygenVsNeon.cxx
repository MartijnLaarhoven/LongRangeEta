/*
 * Compare3times2PC_OxygenVsNeon.cxx
 * Script to compare 3x2PC results between Oxygen-Oxygen and Neon-Neon systems.
 *
 * Author: [Your Name]
 * Date: 2026-03-31
 */

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>

// Helper to load a histogram from a file
TH1* LoadHist(const std::string& filename, const std::string& histname) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return nullptr;
    }
    TH1* h = dynamic_cast<TH1*>(file->Get(histname.c_str()));
    if (!h) {
        std::cerr << "Error: Cannot find histogram " << histname << " in " << filename << std::endl;
        file->Close();
        return nullptr;
    }
    h->SetDirectory(0); // Detach from file
    file->Close();
    return h;
}

void Compare3times2PC_OxygenVsNeon() {
    // Allow ROOT macro execution
    main();
}

int main() {
    gStyle->SetOptStat(0);
    // Define file and histogram names (update as needed)
    std::string neonFile = "./ProcessOutput/3times2PC_NeNe.root";
    std::string oxygenFile = "./ProcessOutput/3times2PC_OO.root";
    std::vector<std::string> histNames = {
        "h3times2PC_FT0AFT0C", // Example: ultra long range
        "h3times2PC_TPCFT0A",
        "h3times2PC_TPCFT0C"
    };
    std::vector<std::string> labels = {"FT0A-FT0C", "TPC-FT0A", "TPC-FT0C"};

    for (size_t i = 0; i < histNames.size(); ++i) {
        TH1* hNe = LoadHist(neonFile, histNames[i]);
        TH1* hO = LoadHist(oxygenFile, histNames[i]);
        if (!hNe || !hO) continue;

        TCanvas* c = new TCanvas(Form("c_%zu", i), labels[i].c_str(), 800, 600);
        hNe->SetLineColor(kBlue+1);
        hNe->SetLineWidth(2);
        hNe->SetMarkerColor(kBlue+1);
        hNe->SetMarkerStyle(20);
        hNe->SetTitle((labels[i] + ": Neon vs Oxygen").c_str());
        hNe->GetYaxis()->SetTitle("3x2PC value");
        hNe->GetXaxis()->SetTitle("Centrality or #eta bin");
        hNe->Draw("E1");

        hO->SetLineColor(kRed+1);
        hO->SetLineWidth(2);
        hO->SetMarkerColor(kRed+1);
        hO->SetMarkerStyle(24);
        hO->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(hNe, "Ne-Ne", "lp");
        leg->AddEntry(hO, "O-O", "lp");
        leg->Draw();

        c->SaveAs(Form("Compare3times2PC_%s.png", histNames[i].c_str()));
        delete c;
        delete hNe;
        delete hO;
    }
    std::cout << "Comparison plots saved as Compare3times2PC_*.png" << std::endl;
    return 0;
}
