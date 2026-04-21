/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 17:23:23 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 09:57:59
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TRandom3.h"
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "./include/BasicForDihadron.h"

struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange) {}
};

struct CorrelationInputSpec {
    std::string roleName;
    Int_t corrType;
    std::string fileNameSuffix;
    Bool_t isNch;
    Bool_t isEtadiff;
    std::vector<std::pair<Int_t, Int_t>> ranges;

    CorrelationInputSpec(std::string _roleName, Int_t _corrType, std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, std::vector<std::pair<Int_t, Int_t>> _ranges) :
        roleName(_roleName), corrType(_corrType), fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), ranges(_ranges) {}
};

std::vector<CorrelationInputSpec> BuildSystemInputSpecs(const std::string& systemName) {
    if (systemName == "O-O") {
        return {
            CorrelationInputSpec("FULL_LM_TPC_FT0A",  kTPCFT0A,  "LHC25ae_pass2_648788", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_MR_TPC_FT0C",  kTPCFT0C,  "LHC25ae_pass2_648788", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_LR_FT0A_FT0C", kFT0AFT0C, "LHC25ae_pass2_648798", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_638221", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_634099", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_634103", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_637591", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_648799", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_648800", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_653257", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_648788", kCent, kEtaDiffOff, {{0, 20}, {80, 100}})
        };
    }
    if (systemName == "O-O-DCAz05") {
        return {
            CorrelationInputSpec("FULL_LM_TPC_FT0A",  kTPCFT0A,  "LHC25ae_pass2_659116", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_MR_TPC_FT0C",  kTPCFT0C,  "LHC25ae_pass2_659116", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_LR_FT0A_FT0C", kFT0AFT0C, "LHC25ae_pass2_659116", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_659117", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_659117", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_659118", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_659118", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_659117", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_659120", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_659118", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_659119", kCent, kEtaDiffOff, {{0, 20}, {80, 100}})
        };
    }
    if (systemName == "Ne-Ne") {
        return {
            CorrelationInputSpec("FULL_LM_TPC_FT0A",  kTPCFT0A,  "LHC25af_pass2_645173",       kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_MR_TPC_FT0C",  kTPCFT0C,  "LHC25af_pass2_645173",       kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25af_pass2_637597",       kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25af_pass2_631290",       kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25af_pass2_637594",       kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_LR_FT0A_FT0C", kFT0AFT0C, "LHC25af_pass2_646139_id50585", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0A",     kFT0AFT0C, "LHC25af_pass2_646139_id50559", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("INNER_LR_FT0C",     kFT0AFT0C, "LHC25af_pass2_646139_id50560", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0A",     kFT0AFT0C, "LHC25af_pass2_646139_id50561", kCent, kEtaDiffOff, {{0, 20}, {80, 100}}),
            CorrelationInputSpec("OUTER_LR_FT0C",     kFT0AFT0C, "LHC25af_pass2_646139_id50562", kCent, kEtaDiffOff, {{0, 20}, {80, 100}})
        };
    }
    if (systemName == "p-O") {
        return {
            CorrelationInputSpec("FULL_LM_TPC_FT0A",  kTPCFT0A,  "LHC25ad_pass2_644389", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_MR_TPC_FT0C",  kTPCFT0C,  "LHC25ad_pass2_644389", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_LR_FT0A_FT0C", kFT0AFT0C, "LHC25ad_pass2_644389", kCent, kEtaDiffOff, {{0, 20}, {80, 100}})
        };
    }
    if (systemName == "p-O-Nch") {
        return {
            CorrelationInputSpec("NCH_LM_TPC_FT0A",   kTPCFT0A,  "LHC25ad_pass2_650295", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_MR_TPC_FT0C",   kTPCFT0C,  "LHC25ad_pass2_650297", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_LR_FT0A_FT0C",  kFT0AFT0C, "LHC25ad_pass2_650299", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ad_pass2_658082", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ad_pass2_658084", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0A",     kFT0AFT0C, "LHC25ad_pass2_658085", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0C",     kFT0AFT0C, "LHC25ad_pass2_658086", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0A",     kFT0AFT0C, "LHC25ad_pass2_658087", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0C",     kFT0AFT0C, "LHC25ad_pass2_658088", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ad_pass2_658089", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ad_pass2_658090", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}})
        };
    }
    if (systemName == "pp") {
        return {
            CorrelationInputSpec("FULL_LM_TPC_FT0A",  kTPCFT0A,  "LHC24af_pass1_644663", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_MR_TPC_FT0C",  kTPCFT0C,  "LHC24af_pass1_644663", kCent, kEtaDiffOn,  {{0, 20}, {80, 100}}),
            CorrelationInputSpec("FULL_LR_FT0A_FT0C", kFT0AFT0C, "LHC24af_pass1_644663", kCent, kEtaDiffOff, {{0, 20}, {80, 100}})
        };
    }
    if (systemName == "pp-Nch") {
        return {
            CorrelationInputSpec("NCH_LM_TPC_FT0A",   kTPCFT0A,  "LHC24af_pass1_650589", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_MR_TPC_FT0C",   kTPCFT0C,  "LHC24af_pass1_650589", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_LR_FT0A_FT0C",  kFT0AFT0C, "LHC24af_pass1_650588", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC24af_pass1_658242", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC24af_pass1_658244", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC24af_pass1_658243", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC24af_pass1_658245", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LR_FT0A",     kFT0AFT0C, "LHC24af_pass1_658238", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0C",     kFT0AFT0C, "LHC24af_pass1_658239", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0A",     kFT0AFT0C, "LHC24af_pass1_658240", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            // Temporary fallback: 658317 is the wrong pp outer-FT0C ring file; use the standard FT0AxFT0C case for now.
            CorrelationInputSpec("NCH_OUTER_LR_FT0C",     kFT0AFT0C, "LHC24af_pass1_650588", kNch, kEtaDiffOff, {{0, 10}, {10, 50}})
        };
    }
    if (systemName == "Ne-Ne-Nch") {
        return {
            CorrelationInputSpec("NCH_LM_TPC_FT0A",   kTPCFT0A,  "LHC25af_pass2_650316", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_MR_TPC_FT0C",   kTPCFT0C,  "LHC25af_pass2_650317", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_LR_FT0A_FT0C",  kFT0AFT0C, "LHC25af_pass2_650315", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC25af_pass2_660449", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25af_pass2_658367", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0A",     kFT0AFT0C, "LHC25af_pass2_660449", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0C",     kFT0AFT0C, "LHC25af_pass2_658367", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25af_pass2_658370", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25af_pass2_658371", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0A",     kFT0AFT0C, "LHC25af_pass2_658370", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0C",     kFT0AFT0C, "LHC25af_pass2_658371", kNch, kEtaDiffOff, {{0, 10}, {10, 50}})
        };
    }
    if (systemName == "O-O-Nch") {
        return {
            CorrelationInputSpec("NCH_LM_TPC_FT0A",   kTPCFT0A,  "LHC25ae_pass2_653254", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_MR_TPC_FT0C",   kTPCFT0C,  "LHC25ae_pass2_653254", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_LR_FT0A_FT0C",  kFT0AFT0C, "LHC25ae_pass2_653254", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_658351", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_658353", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LM_TPC_FT0A", kTPCFT0A,  "LHC25ae_pass2_658356", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_MR_TPC_FT0C", kTPCFT0C,  "LHC25ae_pass2_658357", kNch, kEtaDiffOn,  {{0, 10}, {10, 50}}),

            CorrelationInputSpec("NCH_INNER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_658351", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_INNER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_658353", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0A",     kFT0AFT0C, "LHC25ae_pass2_658356", kNch, kEtaDiffOff, {{0, 10}, {10, 50}}),
            CorrelationInputSpec("NCH_OUTER_LR_FT0C",     kFT0AFT0C, "LHC25ae_pass2_658357", kNch, kEtaDiffOff, {{0, 10}, {10, 50}})
        };
    }
    return {};
}

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange);
void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax);
void CreateAggregateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange);

void Process_CreateBootstrapSample(const std::string& systemOverride = "") {
    std::vector<InputUnit> inputList;
    // Options: "O-O", "O-O-DCAz05", "O-O-Nch", "Ne-Ne", "Ne-Ne-Nch", "p-O", "p-O-Nch", "pp", "pp-Nch"
    const std::string activeSystem = systemOverride.empty() ? "pp-Nch" : systemOverride;
    const auto specs = BuildSystemInputSpecs(activeSystem);

    std::cout << "[Config][Bootstrap] Active system: " << activeSystem << ", entries=" << specs.size() << std::endl;
    for (const auto& spec : specs) {
        std::cout << "  - " << spec.roleName << " | " << DihadronCorrTypeName[spec.corrType]
                  << " | dataset=" << spec.fileNameSuffix << std::endl;
        for (const auto& range : spec.ranges) {
            inputList.push_back(InputUnit(spec.fileNameSuffix, spec.corrType, spec.isNch, spec.isEtadiff, range.first, range.second));
        }
    }

    for (auto input : inputList) {
        if (input.isEtadiff) {
            std::cout << "Processing Bootstrap Sample eta diff: " << input.fileNameSuffix << std::endl;
            for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
                double etaMin = etaBins[iEta];
                double etaMax = etaBins[iEta + 1];
                CreateBootstrapSample_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax);
            }
            // After all EtaDiff samples are created, aggregate them into regular bootstrap files
            std::cout << "Creating aggregate Bootstrap Sample from EtaDiff: " << input.fileNameSuffix << std::endl;
            CreateAggregateBootstrapSample(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange);
        } else {
            std::cout << "Processing Bootstrap Sample: " << input.fileNameSuffix << std::endl;
            CreateBootstrapSample(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange);
        }
    }

}

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    TString inputPath = Form("./ProcessOutput/Mixed_%s_%s_%i_%i_%s.root",
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str());
    std::cout << "[CreateBootstrapSample] Attempting to open: " << inputPath.Data() << std::endl;
    
    // 打开输入文件
    TFile* file = TFile::Open(inputPath.Data(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图（clone后脱离输入文件，避免文件关闭后悬空指针）
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* hIn = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!hIn) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        TH1D* h = dynamic_cast<TH1D*>(hIn->Clone(Form("tmp_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample)));
        h->SetDirectory(nullptr);
        hists.push_back(h);
    }

    if (hists.empty()) {
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    // 创建输出文件
    TString outputPath = Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root",
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str());
    TFile* outFile = TFile::Open(outputPath.Data(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        for (auto* h : hists) delete h;
        delete file;
        return;
    }

    TH1D* hAllIn = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAllIn) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        outFile->Close();
        file->Close();
        for (auto* h : hists) delete h;
        delete outFile;
        delete file;
        return;
    }
    TH1D* hAll = dynamic_cast<TH1D*>(hAllIn->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    hAll->SetDirectory(nullptr);
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    outFile->cd();
    hAll->Write();
    delete hAll;

    // 初始化随机数生成器
    TRandom3 randGen;

    const Int_t availableSamples = static_cast<Int_t>(hists.size());
    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < availableSamples; ++i) {
            selectedIndices.push_back(randGen.Integer(availableSamples));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Int_t mergedCount = 0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                hmerge->SetDirectory(nullptr);
                mergedCount = 1;
            } else {       // 后续累加
                hmerge->Add(current);
                ++mergedCount;
            }
        }

        // 归一化并保存
        if (hmerge && mergedCount > 0) {
            hmerge->Scale(1.0 / mergedCount);
            outFile->cd();
            hmerge->Write();
            delete hmerge;
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    for (auto* h : hists) delete h;
    delete outFile;
    delete file;
    std::cout << "Created bootstrap file: " << outputPath.Data() << std::endl;
}

void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    gSystem->mkdir("./ProcessOutput/EtaDiff", kTRUE);
    
    TString inputPath = Form("./ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
    
    std::cout << "[CreateBootstrapSample_EtaDiff] Attempting to open: " << inputPath.Data() << std::endl;
    
    // 打开输入文件
    TFile* file = TFile::Open(inputPath.Data(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[CreateBootstrapSample_EtaDiff] Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图（clone后脱离输入文件，避免文件关闭后悬空指针）
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* hIn = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!hIn) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        TH1D* h = dynamic_cast<TH1D*>(hIn->Clone(Form("tmp_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample)));
        h->SetDirectory(nullptr);
        hists.push_back(h);
    }

    if (hists.empty()) {
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        return;
    }

    TH1D* hAllIn = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAllIn) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        outFile->Close();
        file->Close();
        for (auto* h : hists) delete h;
        delete outFile;
        delete file;
        return;
    }
    TH1D* hAll = dynamic_cast<TH1D*>(hAllIn->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    hAll->SetDirectory(nullptr);
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    outFile->cd();
    hAll->Write();
    delete hAll;

    // 初始化随机数生成器
    TRandom3 randGen;

    const Int_t availableSamples = static_cast<Int_t>(hists.size());
    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < availableSamples; ++i) {
            selectedIndices.push_back(randGen.Integer(availableSamples));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Int_t mergedCount = 0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                hmerge->SetDirectory(nullptr);
                mergedCount = 1;
            } else {       // 后续累加
                hmerge->Add(current);
                ++mergedCount;
            }
        }

        // 归一化并保存
        if (hmerge && mergedCount > 0) {
            hmerge->Scale(1.0 / mergedCount);
            outFile->cd();
            hmerge->Write();
            delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    for (auto* h : hists) delete h;
    delete outFile;
    delete file;
}

void CreateAggregateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // Create output directory if it doesn't exist
    gSystem->mkdir("./ProcessOutput", kTRUE);
    
    // Create output file for aggregate bootstrap samples
    TString outputPath = Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str());
    
    TFile* outFile = TFile::Open(outputPath.Data(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating aggregate output file!" << std::endl;
        return;
    }
    
    // Helper: open eta-diff bootstrap file, tolerant to formatting variants (0 vs 0.0)
    auto openEtaBootstrapFile = [&](double etaMin, double etaMax) -> TFile* {
        TString pathA = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root",
                            fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange,
                            etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
        TFile* file = TFile::Open(pathA.Data(), "READ");
        if (file && !file->IsZombie()) return file;
        if (file) {
            file->Close();
            delete file;
        }

        TString pathB = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%g_%g_%s.root",
                            fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange,
                            etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
        file = TFile::Open(pathB.Data(), "READ");
        if (file && !file->IsZombie()) return file;
        if (file) {
            file->Close();
            delete file;
        }
        return nullptr;
    };

    // First, sum the main histogram across all eta bins
    TH1D* hAggregate = nullptr;
    int nMainContributors = 0;
    for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
        double etaMin = etaBins[iEta];
        double etaMax = etaBins[iEta + 1];
        TFile* etaFile = openEtaBootstrapFile(etaMin, etaMax);
        if (!etaFile) {
            continue;
        }
        
        TH1D* hIn = dynamic_cast<TH1D*>(etaFile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
        if (!hIn) {
            etaFile->Close();
            delete etaFile;
            continue;
        }

        TH1D* h = dynamic_cast<TH1D*>(hIn->Clone(Form("tmp_main_%d_%d_%d", minRange, maxRange, iEta)));
        h->SetDirectory(nullptr);
        
        if (!hAggregate) {
            hAggregate = dynamic_cast<TH1D*>(h->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
            hAggregate->SetDirectory(nullptr);
        } else {
            hAggregate->Add(h);
        }
        nMainContributors++;
        delete h;
        
        etaFile->Close();
        delete etaFile;
    }
    
    if (hAggregate && nMainContributors > 0) {
        outFile->cd();
        hAggregate->Write();
        delete hAggregate;
    }
    
    // Now sum bootstrap samples across eta bins
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        TH1D* hBsSummed = nullptr;
        int nBsContributors = 0;
        
        for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
            double etaMin = etaBins[iEta];
            double etaMax = etaBins[iEta + 1];
            TFile* etaFile = openEtaBootstrapFile(etaMin, etaMax);
            if (!etaFile) {
                continue;
            }
            
            TH1D* hIn = dynamic_cast<TH1D*>(etaFile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)));
            if (!hIn) {
                etaFile->Close();
                delete etaFile;
                continue;
            }

            TH1D* h = dynamic_cast<TH1D*>(hIn->Clone(Form("tmp_bs_%d_%d_%d_%d", minRange, maxRange, bs, iEta)));
            h->SetDirectory(nullptr);
            
            if (!hBsSummed) {
                hBsSummed = dynamic_cast<TH1D*>(h->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)));
                hBsSummed->SetDirectory(nullptr);
            } else {
                hBsSummed->Add(h);
            }
            nBsContributors++;
            delete h;
            
            etaFile->Close();
            delete etaFile;
        }
        
        if (hBsSummed && nBsContributors > 0) {
            outFile->cd();
            hBsSummed->Write();
            delete hBsSummed;
        }
    }
    
    std::cout << "Created aggregate bootstrap file: " << outputPath.Data() << std::endl;
    
    outFile->Close();
    delete outFile;
}
