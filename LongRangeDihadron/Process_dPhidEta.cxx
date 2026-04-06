/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-15 21:14:52 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 14:49:02
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
R__LOAD_LIBRARY(libO2PhysicsPWGCFCore)

struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;
    Bool_t isMc;

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange, Bool_t _isMc=false) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange), isMc(_isMc) {}
};

void printAxesInfo(THnSparseF* sparseHist);
void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc);
void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Bool_t isMc);
std::string GetInputFileNameSuffix(const std::string &fileNameSuffix);
std::string GetDatasetDirectoryTag(const std::string &fileNameSuffix, Int_t corrType);
std::string ResolveInputDirectory(TFile *file, const std::string &prefix, const std::string &fileNameSuffix, Int_t corrType, const std::string &splitName, Int_t minRange, Int_t maxRange);

// global variables
std::string collisionSystemName = "peripheral PbPb";
std::string additionalSuffix = "";

std::string GetInputFileNameSuffix(const std::string &fileNameSuffix) {
    const std::string taggedPrefix = "LHC25af_pass2_646139_id";
    if (fileNameSuffix.rfind(taggedPrefix, 0) == 0) {
        return "LHC25af_pass2_646139";
    }
    return fileNameSuffix;
}

std::string GetDatasetDirectoryTag(const std::string &fileNameSuffix, Int_t corrType) {
    const std::string taggedPrefix = "LHC25af_pass2_646139_id";
    if (fileNameSuffix.rfind(taggedPrefix, 0) == 0) {
        return fileNameSuffix.substr(taggedPrefix.size() - 2); // keep "idXXXXX"
    }

    if (fileNameSuffix == "LHC25ae_pass2_644429") {
        if (corrType == kTPCFT0A) return "id50663";
        if (corrType == kTPCFT0C) return "id50664";
        if (corrType == kFT0AFT0C) return "id50586";
    }
    if (fileNameSuffix == "LHC25ad_pass2_644389") {
        if (corrType == kTPCFT0A) return "id50674";
        if (corrType == kTPCFT0C) return "id50675";
        if (corrType == kFT0AFT0C) return "id50587";
    }
    if (fileNameSuffix == "LHC24af_pass1_644663") {
        if (corrType == kTPCFT0A) return "id50684";
        if (corrType == kTPCFT0C) return "id50690";
        if (corrType == kFT0AFT0C) return "id50579";
    }
    if (fileNameSuffix == "LHC25af_pass2_631290") {
        if (corrType == kTPCFT0C) return "id47799";
    }
    return "";
}

std::string ResolveInputDirectory(TFile *file, const std::string &prefix, const std::string &fileNameSuffix, Int_t corrType, const std::string &splitName, Int_t minRange, Int_t maxRange) {
    std::vector<std::string> candidates;
    candidates.push_back(Form("%s_%s%s_%d_%d", prefix.c_str(), additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));

    const std::string datasetTag = GetDatasetDirectoryTag(fileNameSuffix, corrType);
    if (!datasetTag.empty()) {
        candidates.push_back(Form("%s_%s_%s_%d_%d", prefix.c_str(), datasetTag.c_str(), splitName.c_str(), minRange, maxRange));
        candidates.push_back(Form("%s_%s%s_%d_%d", prefix.c_str(), datasetTag.c_str(), splitName.c_str(), minRange, maxRange));
        candidates.push_back(Form("%s_%s%s%s_%d_%d", prefix.c_str(), datasetTag.c_str(), additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
        candidates.push_back(Form("%s_%s_%d_%d_%s", prefix.c_str(), splitName.c_str(), minRange, maxRange, datasetTag.c_str()));
        candidates.push_back(Form("%s_%s%s_%d_%d_%s", prefix.c_str(), additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, datasetTag.c_str()));

        // LHC24af FT0A-FT0C has two observed id tags in input files; try both.
        if (fileNameSuffix == "LHC24af_pass1_644663" && corrType == kFT0AFT0C) {
            candidates.push_back(Form("%s_%s_%s_%d_%d", prefix.c_str(), "id50588", splitName.c_str(), minRange, maxRange));
            candidates.push_back(Form("%s_%s%s_%d_%d", prefix.c_str(), "id50588", splitName.c_str(), minRange, maxRange));
            candidates.push_back(Form("%s_%s%s%s_%d_%d", prefix.c_str(), "id50588", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
            candidates.push_back(Form("%s_%s_%d_%d_%s", prefix.c_str(), splitName.c_str(), minRange, maxRange, "id50588"));
            candidates.push_back(Form("%s_%s%s_%d_%d_%s", prefix.c_str(), additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, "id50588"));
        }
    }

    for (const auto &candidate : candidates) {
        if (file->Get(candidate.c_str()) ||
            file->Get(Form("%s/sameEvent", candidate.c_str())) ||
            file->Get(Form("%s/sameEvent_%s", candidate.c_str(), DihadronCorrTypeName[corrType].c_str())) ||
            file->Get(Form("%s/deltaEta_deltaPhi_same_%s", candidate.c_str(), DihadronCorrTypeName[corrType].c_str())) ||
            file->Get(Form("%s/Trig_hist_%s", candidate.c_str(), DihadronCorrTypeName[corrType].c_str()))) {
            return candidate;
        }
    }

    return candidates.front();
}

void Process_dPhidEta() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);

    std::vector<InputUnit> inputList;
    additionalSuffix = "";
    
    collisionSystemName = "Ne-Ne";
    // Ne-Ne EtaDiff centrality-based
    // inputList.push_back(InputUnit("LHC25af_pass2_642734", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_642734", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // Dataset 1: LHC25af_pass2_632504 with TPC_FT0A
    // inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // Dataset 2: Ne-Ne full-acceptance TPC-FT0C (kept for baseline/full comparison)
    // inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // Ring-specific Ne-Ne TPC-FT0 inputs
    // inputList.push_back(InputUnit("LHC25af_pass2_637597", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_637597", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_631290", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_631290", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_637594", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_637594", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // Dataset 3: LHC25af_pass2_642734 with FT0A_FT0C (single full-range, NOT eta-differential)
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50585", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50585", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50559", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50559", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50560", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50560", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50561", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50561", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50562", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_646139_id50562", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // Dataset 4: O-O full-range baseline
    // inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));

    // Dataset 5: O-O ring-specific TPC-FT0 and ultra-long-range FT0A-FT0C side inputs
    // inputList.push_back(InputUnit("LHC25ae_pass2_638221", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_638221", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_634099", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_634099", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_634103", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_634103", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_637591", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_637591", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648799", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648799", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648800", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648800", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_644433", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_644433", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648788", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_648788", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    
    // p-O datasets (template: 80-100, signal: 0-20)
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ad_pass2_644389", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));
    
    // p-p datasets (LHC24af)
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC24af_pass1_644663", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));

    // Ne-Ne Nch-dependent datasets (template: 0-10, data: 10-50), full-range only
    inputList.push_back(InputUnit("LHC25af_pass2_650316", kTPCFT0A, kNch, kEtaDiffOn, 0, 10));
    inputList.push_back(InputUnit("LHC25af_pass2_650316", kTPCFT0A, kNch, kEtaDiffOn, 10, 50));
    inputList.push_back(InputUnit("LHC25af_pass2_650317", kTPCFT0C, kNch, kEtaDiffOn, 0, 10));
    inputList.push_back(InputUnit("LHC25af_pass2_650317", kTPCFT0C, kNch, kEtaDiffOn, 10, 50));
    inputList.push_back(InputUnit("LHC25af_pass2_650315", kFT0AFT0C, kNch, kEtaDiffOff, 0, 10));
    inputList.push_back(InputUnit("LHC25af_pass2_650315", kFT0AFT0C, kNch, kEtaDiffOff, 10, 50));
    

    for (auto input : inputList) {
        collisionSystemName = GetCollisionSystemNameFromDataset(input.fileNameSuffix);
        if (input.isEtadiff) {
            for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
                double etaMin = etaBins[iEta];
                double etaMax = etaBins[iEta + 1];
                Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.isMc);
            }
        }
        else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, input.isMc);
        }
    }
}

void printAxesInfo(THnSparseF* sparseHist) {
    // print the information of the axes of THnSparseF
    if (!sparseHist) {
        std::cerr << "Error: Null histogram pointer!" << std::endl;
        return;
    }

    const Int_t nDims = sparseHist->GetNdimensions();
    std::cout << "Sparse histogram has " << nDims << " dimensions\n";

    for (Int_t iDim = 0; iDim < nDims; ++iDim) {
        TAxis* axis = sparseHist->GetAxis(iDim);
        if (!axis) {
        std::cerr << "Error: Axis " << iDim << " not found!" << std::endl;
        continue;
        }
        if (iDim == 0) {
            axis->SetRangeUser(0,0.001);
              // 获取当前显示的区间参数
            const int firstBin = axis->GetFirst();
            const int lastBin = axis->GetLast();
            const double firstEdge = axis->GetBinLowEdge(firstBin);
            const double lastEdge = axis->GetBinUpEdge(lastBin);

            std::cout << "=== Axis 0" << " (" << axis->GetTitle() << ") ===\n"
                        << "Current displayed range: [" << firstEdge << ", " << lastEdge << "]\n"
                        << "Using bins: " << firstBin << " (" << firstEdge << ") -> "
                        << lastBin << " (" << lastEdge << ")\n\n";
        }

        std::cout << "Axis " << iDim << ": \n"
                << "  Title: " << axis->GetTitle() << "\n"
                << "  Bins:  " << axis->GetNbins() << "\n"
                << "  Range: [" << axis->GetXmin() 
                << ", " << axis->GetXmax() << "]\n"
                << "---------------------------------\n";
    }
}

void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc=false) {
    const std::string inputFileNameSuffix = GetInputFileNameSuffix(fileNameSuffix);
    TFile *file = TFile::Open(Form("../../../AnalysisResultsROOTFiles/LongRangeEta/AnalysisResults_%s.root", inputFileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TString longRangeDir = ResolveInputDirectory(file, "long-range-dihadron-cor", fileNameSuffix, corrType, splitName, minRange, maxRange).c_str();

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("%s/MCTrue", longRangeDir.Data()))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("%s/MCTrue", longRangeDir.Data()))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }

    TString flowDirProbe = ResolveInputDirectory(file, "flow-decorrelation", fileNameSuffix, corrType, splitName, minRange, maxRange).c_str();
    bool hasDirectFT0Input =
        file->Get(Form("%s/deltaEta_deltaPhi_same_%s", flowDirProbe.Data(), DihadronCorrTypeName[corrType].c_str())) &&
        file->Get(Form("%s/deltaEta_deltaPhi_mixed_%s", flowDirProbe.Data(), DihadronCorrTypeName[corrType].c_str())) &&
        file->Get(Form("%s/Trig_hist_%s", flowDirProbe.Data(), DihadronCorrTypeName[corrType].c_str()));

    if (!isMc && corrType == kFT0AFT0C) {
        std::cout << "[FT0 routing] " << fileNameSuffix << " [" << minRange << ", " << maxRange << "] hasDirectFT0Input=" << hasDirectFT0Input << " (flowDir=" << flowDirProbe.Data() << ")" << std::endl;
    }

    // Special handling for FT0A_FT0C in LongRangeEta input: direct 2D histograms.
    // Prefer direct objects for ultra-long-range FT0 EtaDiffOff and keep
    // CorrelationContainer path as fallback when direct objects are unavailable.
    
    if (!isMc && corrType == kFT0AFT0C && hasDirectFT0Input) {
        std::cout << "[FT0 routing] Using direct FT0 histograms for " << fileNameSuffix << " [" << minRange << ", " << maxRange << "]" << std::endl;
        TString flowDir = flowDirProbe;
        TH2D* hPhiEtaS = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_same_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        TH2D* hPhiEtaM = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_mixed_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        THnSparseD *trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));

        if (!hPhiEtaS || !hPhiEtaM || !trig) {
            std::cerr << "Error getting FT0A_FT0C histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
            file->Close();
            delete file;
            return;
        }

        TH2D* hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d", minRange, maxRange));
        TH2D* hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d", minRange, maxRange));

        TH2D* hPhiEtaMSafe = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_safe_%d_%d", minRange, maxRange));
        Double_t minPositive = 0.0;
        // for (Int_t ix = 1; ix <= hPhiEtaMSafe->GetNbinsX(); ++ix) {
        //     for (Int_t iy = 1; iy <= hPhiEtaMSafe->GetNbinsY(); ++iy) {
        //         Double_t value = hPhiEtaMSafe->GetBinContent(ix, iy);
        //         if (value > 0.0 && (minPositive == 0.0 || value < minPositive)) {
        //             minPositive = value;
        //         }
        //     }
        // }
        // if (minPositive <= 0.0) minPositive = 1e-12;
        // for (Int_t ix = 1; ix <= hPhiEtaMSafe->GetNbinsX(); ++ix) {
        //     for (Int_t iy = 1; iy <= hPhiEtaMSafe->GetNbinsY(); ++iy) {
        //         if (hPhiEtaMSafe->GetBinContent(ix, iy) <= 0.0) {
        //             hPhiEtaMSafe->SetBinContent(ix, iy, minPositive);
        //         }
        //     }
        // }

        Int_t Nbins = hPhiEtaMSafe->GetYaxis()->GetNbins();

        TH2D* hPhiEtaSMdiv = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SMdiv2_%d_%d", minRange, maxRange));
        TH2D* hPhiEtaMdiv = (TH2D*)hPhiEtaMSafe->Clone(Form("dphideta_MEdiv2_%d_%d", minRange, maxRange));
        hPhiEtaSMdiv->Rebin2D(3, Nbins);
        hPhiEtaMdiv->Rebin2D(3, Nbins);
        // hPhiEtaMdiv->GetYaxis()->SetRangeUser(6.2, 7);
        // hPhiEtaSMdiv->GetYaxis()->SetRangeUser(6.2, 7);

        TH2D* hPhiEtaSMsum = (TH2D*)hPhiEtaSMdiv->Clone(Form("dphideta_SM_%d_%d", minRange, maxRange));
        hPhiEtaSMsum->Divide(hPhiEtaSMdiv, hPhiEtaMdiv, 1.0, 1.0);

        Double_t nTriggersS = trig->Integral(true);
        if (nTriggersS > 0) {
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
            hPhiEtaMsum->Scale(1.0 / nTriggersS);
        }

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));

        hPhiEtaSMsum->Rebin2D(1, 1);
        hPhiEtaSsum->Rebin2D(1, 1);
        hPhiEtaMsum->Rebin2D(1, 1);

        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d", minRange, maxRange));
        hEta->SetTitle("#Delta#eta");

        TH1D* hPhiSameOverMixed = hPhiEtaSMsum->ProjectionX(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
        hPhiSameOverMixed->SetTitle(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
        hPhiSameOverMixed->GetXaxis()->SetTitle("#Delta#varphi");

        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        if (isNch) {
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        } else {
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));
        }

        TFile* fout = TFile::Open(Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()), "RECREATE");
        if (!fout || !fout->IsOpen()) {
            std::cerr << "Error: Cannot create output file for FT0A_FT0C case" << std::endl;
            delete hPhiEtaSMsum;
            delete hPhiEtaSsum;
            delete hPhiEtaMsum;
            delete hEta;
            delete hPhiSameOverMixed;
            file->Close();
            delete file;
            return;
        }

        TCanvas* c1 = new TCanvas("dPhidEta ", "dPhidEta ", 1200, 800);
        c1->Divide(2, 2);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);

        c1->cd(1);
        TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
        hPhiEtaSMsum_draw->Draw("surf1");
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        c1->cd(3);
        hPhiEtaSsum->Draw("surf1");
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        c1->cd(4);
        hPhiEtaMsum->Draw("surf1");
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        c1->Write();
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();
        hPhiSameOverMixed->Write();

        UInt_t seedBase = 12345 + 1000 * minRange + 10 * maxRange;
        TRandom3 randGen(seedBase);
        for (Int_t sample = 0; sample < maxSample; ++sample) {
            TH1D* hSample = (TH1D*)hPhiSameOverMixed->Clone(Form("hPhiSameOverMixed_%i_%i_%d", minRange, maxRange, sample));
            for (Int_t ibin = 1; ibin <= hSample->GetNbinsX(); ++ibin) {
                Double_t mean = hPhiSameOverMixed->GetBinContent(ibin);
                Double_t sigma = hPhiSameOverMixed->GetBinError(ibin);
                if (!std::isfinite(sigma) || sigma <= 0.0) {
                    sigma = std::max(1e-6, 1e-3 * std::fabs(mean));
                }
                Double_t fluctuated = randGen.Gaus(mean, sigma);
                if (!std::isfinite(fluctuated)) fluctuated = mean;
                hSample->SetBinContent(ibin, fluctuated);
                hSample->SetBinError(ibin, sigma);
            }
            hSample->Write();
            delete hSample;
        }

        fout->Close();
        delete fout;
        delete hPhiEtaSMsum;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiEtaMSafe;
        delete hPhiEtaSMsum_draw;
        delete c1;
        delete hEta;
        delete hPhiSameOverMixed;
        file->Close();
        delete file;
        std::cout << "Output file: " << Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()) << std::endl;
        std::cout << "Processing completed for all samples." << std::endl;
        return;
    }

        

    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("%s/sameEvent_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    if (!same) {
        same = (CorrelationContainer*)file->Get(Form("%s/sameEvent", longRangeDir.Data()));
    }
    if (!mixed) {
        mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent", longRangeDir.Data()));
    }
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("%s/MCTrue/MCTrig_hist_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    }

    if (!same || !mixed || !trig) {
        TString flowDir = ResolveInputDirectory(file, "flow-decorrelation", fileNameSuffix, corrType, splitName, minRange, maxRange).c_str();
        same = (CorrelationContainer*)file->Get(Form("%s/sameEvent_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        if (!same) {
            same = (CorrelationContainer*)file->Get(Form("%s/sameEvent", flowDir.Data()));
        }
        if (!mixed) {
            mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent", flowDir.Data()));
        }
        if (!isMc) {
            trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        } else {
            trig = (THnSparseD*)file->Get(Form("%s/MCTrue/MCTrig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        }
        if (same && mixed && trig) {
            longRangeDir = flowDir;
        }
    }

    if (!isMc && corrType == kFT0AFT0C) {
        std::cout << "[FT0 routing] Using CorrelationContainer fallback for " << fileNameSuffix << " [" << minRange << ", " << maxRange << "] from dir " << longRangeDir.Data() << std::endl;
    }
    

    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    THnSparseF *sparSig = nullptr;
    THnSparseF *sparMix = nullptr;
    Int_t selectedStep = isMc ? CorrelationContainer::kCFStepAll : 6;
    sparSig = (THnSparseF*)same->getPairHist()->getTHn(selectedStep);
    sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(selectedStep);
    if (!sparSig || !sparMix || sparSig->GetNbins() <= 0 || sparMix->GetNbins() <= 0) {
        std::cerr << "Error: CF step " << selectedStep << " unavailable/empty (zero filled bins) for " << fileNameSuffix << " in [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }
    std::cout << "Using CF step " << selectedStep << " for " << fileNameSuffix << " " << DihadronCorrTypeName[corrType] << " [" << minRange << ", " << maxRange << "]" << std::endl;

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Common axis settings for all samples
    sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    trig->GetAxis(trigAxis_pT)->SetRangeUser(minPt+0.001, maxPt-0.001);

    // Create output file
    TFile* fout = TFile::Open(Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()), "RECREATE");

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) {
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);

            // Normalize mixed event using peak-based method
            Double_t norm = 1.;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            
            // Find peak in eta distribution
            TH1D* hEtaNorm = hPhiEtaM->ProjectionY("hEtaNorm_temp");
            Int_t peakEtaBin = hEtaNorm->GetMaximumBin();
            Int_t binEta1 = std::max(1, peakEtaBin - 1);
            Int_t binEta2 = std::min(hPhiEtaM->GetNbinsY(), peakEtaBin + 1);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            delete hEtaNorm;
            
            // Fallback to full histogram if peak method fails
            if (!std::isfinite(norm) || norm <= 0.0) {
                const Int_t nAllBins = hPhiEtaM->GetNbinsX() * hPhiEtaM->GetNbinsY();
                const Double_t fullIntegral = hPhiEtaM->Integral(1, hPhiEtaM->GetNbinsX(), 1, hPhiEtaM->GetNbinsY());
                norm = (nAllBins > 0) ? (fullIntegral / nAllBins) : 1.0;
            }

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            Int_t dEtaNims = hPhiEtaM->GetYaxis()->GetNbins();
            TH2D* hPhiEtaMdiv = (TH2D*)hPhiEtaM->Clone(Form("dphideta_Mdiv_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            hPhiEtaMdiv->Rebin2D(1, dEtaNims);
            hPhiEtaSM->Rebin2D(1, dEtaNims);

            hPhiEtaSM->Divide(hPhiEtaMdiv);
            for (Int_t ix = 1; ix <= hPhiEtaSM->GetNbinsX(); ++ix) {
                for (Int_t iy = 1; iy <= hPhiEtaSM->GetNbinsY(); ++iy) {
                    const Double_t value = hPhiEtaSM->GetBinContent(ix, iy);
                    if (!std::isfinite(value)) {
                        hPhiEtaSM->SetBinContent(ix, iy, 0.0);
                        hPhiEtaSM->SetBinError(ix, iy, 0.0);
                    }
                }
            }

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");
        hPhiEtaSMsum->Rebin2D(1, 1);

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);
        if (corrType == kFT0AFT0C) {
            hPhiEtaMsum->Rebin2D(1, 1);
        }

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);
         if (corrType == kFT0AFT0C) {
            hPhiEtaSsum->Rebin2D(1, 1);
        }

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // set the ranges
        // hPhiEtaSMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaSsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Draw histograms
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

        c1->cd(1);
        TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
        // hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
        hPhiEtaSMsum_draw->Draw("surf1");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


        c1->cd(3);
        hPhiEtaSsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


        c1->cd(4);
        hPhiEtaMsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        // Saving to file
        fout->cd();
        c1->Write();
        // write canvas to file
        // c1->SaveAs(Form("./ProcessOutput/Mixed_%s_%s_%i_%i.png", fileNameSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)));
        // delete c1;

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        for (Int_t ibin = 1; ibin <= hPhiSameOverMixed_pos->GetNbinsX(); ++ibin) {
            const Double_t value = hPhiSameOverMixed_pos->GetBinContent(ibin);
            if (!std::isfinite(value)) {
                hPhiSameOverMixed_pos->SetBinContent(ibin, 0.0);
                hPhiSameOverMixed_pos->SetBinError(ibin, 0.0);
            }
        }
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete c1;
        
    }

    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    std::cout << "Output file: " << Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()) << std::endl; 
    std::cout << "Processing completed for all samples." << std::endl;
}

void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Bool_t isMc=false) {
    const std::string inputFileNameSuffix = GetInputFileNameSuffix(fileNameSuffix);
    TFile *file = TFile::Open(Form("../../../AnalysisResultsROOTFiles/LongRangeEta/AnalysisResults_%s.root", inputFileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TString longRangeDir = ResolveInputDirectory(file, "long-range-dihadron-cor", fileNameSuffix, corrType, splitName, minRange, maxRange).c_str();
    TString flowDir = ResolveInputDirectory(file, "flow-decorrelation", fileNameSuffix, corrType, splitName, minRange, maxRange).c_str();

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("%s/MCTrue", flowDir.Data()))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("%s/MCTrue", flowDir.Data()))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }

    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("%s/MCTrue/MCTrig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    }

    if (!isMc && corrType == kFT0AFT0C) {
        TH2D* hPhiEtaS = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_same_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        TH2D* hPhiEtaM = (TH2D*)file->Get(Form("%s/deltaEta_deltaPhi_mixed_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));

        if (!hPhiEtaS || !hPhiEtaM || !trig) {
            std::cerr << "Error getting FT0A_FT0C histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
            file->Close();
            delete file;
            return;
        }

        TH2D* hPhiEtaSLocal = (TH2D*)hPhiEtaS->Clone(Form("sameEvent_%i_%i", minRange, maxRange));
        TH2D* hPhiEtaMLocal = (TH2D*)hPhiEtaM->Clone(Form("mixedEvent_%i_%i", minRange, maxRange));

        Int_t binPhi1 = hPhiEtaMLocal->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
        Int_t binPhi2 = hPhiEtaMLocal->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
        Int_t binEta1 = hPhiEtaMLocal->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]-0.001);
        Int_t binEta2 = hPhiEtaMLocal->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]+0.001);
        Double_t norm = hPhiEtaMLocal->Integral(binPhi1, binPhi2, binEta1, binEta2);
        if (norm <= 0) {
            std::cerr << "Error: invalid mixed-event normalization for FT0A_FT0C in " << fileNameSuffix << std::endl;
            delete hPhiEtaSLocal;
            delete hPhiEtaMLocal;
            file->Close();
            delete file;
            return;
        }
        hPhiEtaMLocal->Scale(1./norm);

        TH2D* hPhiEtaSM = (TH2D*)hPhiEtaSLocal->Clone(Form("hPhiEtaSameOverMixed_%i_%i", minRange, maxRange));
        hPhiEtaSM->Divide(hPhiEtaMLocal);

        Int_t binEtaLow = hPhiEtaSM->GetYaxis()->FindBin(DihadrondEtaRange[corrType][0] + 0.001);
        Int_t binEtaHigh = hPhiEtaSM->GetYaxis()->FindBin(DihadrondEtaRange[corrType][1] - 0.001);
        TH1D* hPhiSMsum = (TH1D*)hPhiEtaSM->ProjectionX(Form("hPhiSameOverMixed_%i_%i", minRange, maxRange), binEtaLow, binEtaHigh);
        hPhiSMsum->Scale(1./(binEtaHigh - binEtaLow + 1));

        gSystem->mkdir("./ProcessOutput/EtaDiff", kTRUE);
        TFile* fout = TFile::Open(Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()), "RECREATE");
        if (!fout || !fout->IsOpen()) {
            std::cerr << "Error: Cannot create output file for FT0A_FT0C case" << std::endl;
            file->Close();
            delete file;
            return;
        }
        hPhiEtaSM->Write();
        hPhiSMsum->Write();

        for (Int_t sample = 0; sample < maxSample; ++sample) {
            TH1D* hSample = (TH1D*)hPhiSMsum->Clone(Form("hPhiSameOverMixed_%i_%i_%d", minRange, maxRange, sample));
            hSample->Write();
            delete hSample;
        }

        fout->Close();
        delete fout;
        delete hPhiEtaSM;
        delete hPhiSMsum;
        delete hPhiEtaSLocal;
        delete hPhiEtaMLocal;
        file->Close();
        delete file;
        std::cout << "Saved FT0A_FT0C direct-hist result to ./ProcessOutput/EtaDiff/Mixed_" << fileNameSuffix << additionalSuffix << "_" << splitName << "_" << int(minRange) << "_" << int(maxRange) << "_Eta_" << etaMin << "_" << etaMax << "_" << DihadronCorrTypeName[corrType] << ".root" << std::endl;
        return;
    }

    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("%s/sameEvent_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    if (!same) {
        same = (CorrelationContainer*)file->Get(Form("%s/sameEvent", longRangeDir.Data()));
    }
    if (!mixed) {
        mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent", longRangeDir.Data()));
    }

    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("%s/MCTrue/MCTrig_hist_%s", longRangeDir.Data(), DihadronCorrTypeName[corrType].c_str()));
    }

    if (!same || !mixed || !trig) {
        same = (CorrelationContainer*)file->Get(Form("%s/sameEvent_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        if (!same) {
            same = (CorrelationContainer*)file->Get(Form("%s/sameEvent", flowDir.Data()));
        }
        if (!mixed) {
            mixed = (CorrelationContainer*)file->Get(Form("%s/mixedEvent", flowDir.Data()));
        }
        if (!isMc) {
            trig = (THnSparseD*)file->Get(Form("%s/Trig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        } else {
            trig = (THnSparseD*)file->Get(Form("%s/MCTrue/MCTrig_hist_%s", flowDir.Data(), DihadronCorrTypeName[corrType].c_str()));
        }
        if (same && mixed && trig) {
            longRangeDir = flowDir;
        }
    }

    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    THnSparseF *sparSig = nullptr;
    THnSparseF *sparMix = nullptr;
    Int_t selectedStep = isMc ? CorrelationContainer::kCFStepAll : 6;
    sparSig = (THnSparseF*)same->getPairHist()->getTHn(selectedStep);
    sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(selectedStep);
    if (!sparSig || !sparMix || sparSig->GetNbins() <= 0 || sparMix->GetNbins() <= 0) {
        std::cerr << "Error: CF step " << selectedStep << " unavailable/empty (zero filled bins) for eta-diff " << fileNameSuffix << " in [" << minRange << ", " << maxRange << "] eta [" << etaMin << ", " << etaMax << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }
    std::cout << "Using CF step " << selectedStep << " for eta-diff " << fileNameSuffix << " " << DihadronCorrTypeName[corrType] << " [" << minRange << ", " << maxRange << "] eta [" << etaMin << ", " << etaMax << "]" << std::endl;

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Common axis settings for all samples
    // trigger particle: eta bins
    // associated particle: reference eta range
    sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaMin+0.0001, etaMax-0.0001);
    sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaMin+0.0001, etaMax-0.0001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    trig->GetAxis(trigAxis_pT)->SetRangeUser(etaMin+0.0001, etaMax-0.0001);

    // Create output file
    gSystem->mkdir("./ProcessOutput/EtaDiff", kTRUE);
    TFile* fout = TFile::Open(Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()), "RECREATE");
    if (!fout || !fout->IsOpen()) {
        std::cerr << "Error: Cannot create output file for eta differential analysis" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) {
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);

            // Normalize mixed event using peak-based method
            Double_t norm = 1.;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            
            // Find peak in eta distribution
            TH1D* hEtaNorm = hPhiEtaM->ProjectionY("hEtaNorm_temp");
            Int_t peakEtaBin = hEtaNorm->GetMaximumBin();
            Int_t binEta1 = std::max(1, peakEtaBin - 1);
            Int_t binEta2 = std::min(hPhiEtaM->GetNbinsY(), peakEtaBin + 1);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            delete hEtaNorm;
            
            // Fallback to full histogram if peak method fails
            if (!std::isfinite(norm) || norm <= 0.0) {
                const Int_t nAllBins = hPhiEtaM->GetNbinsX() * hPhiEtaM->GetNbinsY();
                const Double_t fullIntegral = hPhiEtaM->Integral(1, hPhiEtaM->GetNbinsX(), 1, hPhiEtaM->GetNbinsY());
                norm = (nAllBins > 0) ? (fullIntegral / nAllBins) : 1.0;
            }

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            Int_t dEtaNims = hPhiEtaM->GetYaxis()->GetNbins();
            TH2D* hPhiEtaMdiv = (TH2D*)hPhiEtaM->Clone(Form("dphideta_Mdiv_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            hPhiEtaMdiv->Rebin2D(1, dEtaNims);
            hPhiEtaSM->Rebin2D(1, dEtaNims);

            hPhiEtaSM->Divide(hPhiEtaMdiv);
            for (Int_t ix = 1; ix <= hPhiEtaSM->GetNbinsX(); ++ix) {
                for (Int_t iy = 1; iy <= hPhiEtaSM->GetNbinsY(); ++iy) {
                    const Double_t value = hPhiEtaSM->GetBinContent(ix, iy);
                    if (!std::isfinite(value)) {
                        hPhiEtaSM->SetBinContent(ix, iy, 0.0);
                        hPhiEtaSM->SetBinError(ix, iy, 0.0);
                    }
                }
            }

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }

            delete hPhiEtaSM;
            delete hPhiEtaMdiv;
            delete hPhiEtaS;
            delete hPhiEtaM;
            delete hTriggersS;
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");
        hPhiEtaSMsum->Rebin2D(1, 1);

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // set the ranges
        // hPhiEtaSMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaSsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Save histograms only (no dPhidEta canvas output in EtaDiff mode)
        fout->cd();

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        for (Int_t ibin = 1; ibin <= hPhiSameOverMixed_pos->GetNbinsX(); ++ibin) {
            const Double_t value = hPhiSameOverMixed_pos->GetBinContent(ibin);
            if (!std::isfinite(value)) {
                hPhiSameOverMixed_pos->SetBinContent(ibin, 0.0);
                hPhiSameOverMixed_pos->SetBinError(ibin, 0.0);
            }
        }
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete hEta;
        
    }

    // clean up
    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    std::cout << "Output file: " << Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()) << std::endl;
    std::cout << "Processing completed for all samples." << std::endl;
}

