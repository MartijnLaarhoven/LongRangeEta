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

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange);
void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax);
void CreateAggregateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange);

void Process_CreateBootstrapSample() {
    std::vector<InputUnit> inputList;
    // p-O centrality-based datasets
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25ad_pass2_644389", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));

    // Ne-Ne centrality-based datasets
    // Dataset: LHC25af_pass2_632504 (TPC_FT0A), LHC25af_pass2_637596 (TPCFT0C), LHC25af_pass2_645746 (FT0AFT0C)
    inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_645746", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_645746", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));

    // O-O centrality-based datasets
    inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    inputList.push_back(InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kCent, kEtaDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kCent, kEtaDiffOff, 80, 100));

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
