#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

struct FT0AFT0CInput {
  std::string fileNameSuffix;
  int minRange;
  int maxRange;
};

void ProcessOneFT0AFT0C(const FT0AFT0CInput& input) {
  TFile* file = TFile::Open(Form("../../../AnalysisResultsROOTFiles/LongRangeEta/AnalysisResults_%s.root", input.fileNameSuffix.c_str()), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input: " << input.fileNameSuffix << std::endl;
    return;
  }

  TString flowDir = Form("flow-decorrelation_Cent_%d_%d", input.minRange, input.maxRange);
  TH2D* hSame = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_same_FT0A_FT0C", flowDir.Data()));
  TH2D* hMixed = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_mixed_FT0A_FT0C", flowDir.Data()));
  if (!hSame || !hMixed) {
    std::cerr << "[ERROR] Missing FT0A-FT0C 2D hist in " << input.fileNameSuffix
              << " cent [" << input.minRange << ", " << input.maxRange << "]" << std::endl;
    file->Close();
    delete file;
    return;
  }

  TH2D* hSameLocal = (TH2D*)hSame->Clone("hSameLocal");
  TH2D* hMixedLocal = (TH2D*)hMixed->Clone("hMixedLocal");

  const int binPhi1 = hMixedLocal->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
  const int binPhi2 = hMixedLocal->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);

  TH1D* hEtaMixNorm = hMixedLocal->ProjectionY("hEtaMixNorm");
  const int peakEtaBin = hEtaMixNorm ? hEtaMixNorm->GetMaximumBin() : hMixedLocal->GetYaxis()->FindBin(6.9);
  const int binEta1 = std::max(1, peakEtaBin - 1);
  const int binEta2 = std::min(hMixedLocal->GetNbinsY(), peakEtaBin + 1);
  const int nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);

  double norm = hMixedLocal->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
  if (!std::isfinite(norm) || norm <= 0.0) {
    const int nAllBins = hMixedLocal->GetNbinsX() * hMixedLocal->GetNbinsY();
    const double fullIntegral = hMixedLocal->Integral(1, hMixedLocal->GetNbinsX(), 1, hMixedLocal->GetNbinsY());
    norm = (nAllBins > 0) ? (fullIntegral / nAllBins) : -1.0;
  }
  delete hEtaMixNorm;

  if (!std::isfinite(norm) || norm <= 0.0) {
    std::cerr << "[ERROR] Invalid mixed-event norm for " << input.fileNameSuffix
              << " cent [" << input.minRange << ", " << input.maxRange << "]" << std::endl;
    delete hSameLocal;
    delete hMixedLocal;
    file->Close();
    delete file;
    return;
  }

  hMixedLocal->Scale(1.0 / norm);
  TH2D* hSM = (TH2D*)hSameLocal->Clone("dphideta_SM_full");
  hSM->Divide(hMixedLocal);

  TH1D* hPhiFull = hSM->ProjectionX("hPhiSameOverMixed_full", 1, hSM->GetNbinsY());
  hPhiFull->Scale(1.0 / hSM->GetNbinsY());

  const int oldBinLow = hSM->GetYaxis()->FindBin(6.5 + 0.001);
  const int oldBinHigh = hSM->GetYaxis()->FindBin(8.0 - 0.001);
  TH1D* hPhiOldWindow = hSM->ProjectionX("hPhiSameOverMixed_window6p5_8p0", oldBinLow, oldBinHigh);
  hPhiOldWindow->Scale(1.0 / (oldBinHigh - oldBinLow + 1));

  int nNan = 0;
  int nInf = 0;
  for (int ix = 1; ix <= hSM->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hSM->GetNbinsY(); ++iy) {
      double value = hSM->GetBinContent(ix, iy);
      if (std::isnan(value)) nNan++;
      if (!std::isfinite(value) && !std::isnan(value)) nInf++;
    }
  }

  gSystem->mkdir("./ProcessOutput/EtaDiffFullRange", kTRUE);
  TFile* fout = TFile::Open(Form("./ProcessOutput/EtaDiffFullRange/Mixed_%s_Cent_%d_%d_Eta_-0.8_0.8_FT0A_FT0C_FullDEta.root",
                                 input.fileNameSuffix.c_str(), input.minRange, input.maxRange), "RECREATE");
  if (!fout || !fout->IsOpen()) {
    std::cerr << "[ERROR] Cannot create output file for " << input.fileNameSuffix << std::endl;
    delete hPhiOldWindow;
    delete hPhiFull;
    delete hSM;
    delete hSameLocal;
    delete hMixedLocal;
    file->Close();
    delete file;
    return;
  }

  hSM->Write("dphideta_SM_0_100");
  hSameLocal->Write("dphideta_SE_0_100");
  hMixedLocal->Write("dphideta_ME_0_100");
  hPhiFull->Write();
  hPhiOldWindow->Write();

  fout->Close();
  delete fout;

  std::cout << "[OK] " << input.fileNameSuffix << " cent [" << input.minRange << ", " << input.maxRange << "]"
            << " norm=" << norm << " fullIntegral=" << hPhiFull->Integral()
            << " windowIntegral=" << hPhiOldWindow->Integral()
            << " nan=" << nNan << " inf=" << nInf << std::endl;

  delete hPhiOldWindow;
  delete hPhiFull;
  delete hSM;
  delete hSameLocal;
  delete hMixedLocal;
  file->Close();
  delete file;
}

void Process_dPhidEta_FT0AFT0C_FullRange() {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  std::vector<FT0AFT0CInput> inputs = {
    {"LHC25af_pass2_640018", 0, 20},
    {"LHC25af_pass2_640018", 80, 100}
  };

  for (const auto& input : inputs) {
    ProcessOneFT0AFT0C(input);
  }
}
