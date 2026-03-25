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

    // Ne-Ne EtaDiff centrality-based
    // Dataset 1: LHC25af_pass2_632504 with TPC_FT0A
    inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_632504", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));
    // Dataset 2: LHC25af_pass2_637596 with TPC_FT0C
    inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_637596", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));
    // Dataset 3: LHC25af_pass2_640018 with FT0A_FT0C
    inputList.push_back(InputUnit("LHC25af_pass2_640018", kFT0AFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_640018", kFT0AFT0C, kCent, kEtaDiffOn, 80, 100));

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
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/Mixed_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图
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
        Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()), 
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
        for (auto* h : hists) delete h;
        outFile->Close();
        file->Close();
        delete outFile;
        delete file;
        return;
    }
    TH1D* hAll = dynamic_cast<TH1D*>(hAllIn->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    hAll->SetDirectory(nullptr);
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();
    delete hAll;

    // 初始化随机数生成器
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            // delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    for (auto* h : hists) delete h;
    delete outFile;
    delete file;
}

void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* h = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!h) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
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

    TH1D* hAll = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAll) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        return;
    }
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();

    // 初始化随机数生成器
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}

void CreateAggregateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // Create output directory if it doesn't exist
    gSystem->mkdir("./ProcessOutput", kTRUE);
    
    // Read all EtaDiff bootstrap files and sum histograms across eta bins
    std::vector<std::string> etaDiffFiles;
    for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
        double etaMin = etaBins[iEta];
        double etaMax = etaBins[iEta + 1];
        etaDiffFiles.push_back(
            Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", 
                 fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str())
        );
    }
    
    // Create output file for aggregate bootstrap samples
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating aggregate output file!" << std::endl;
        return;
    }
    
    // First, sum the main histogram across all eta bins
    TH1D* hAggregate = nullptr;
    for (const auto& filePath : etaDiffFiles) {
        TFile* etaFile = TFile::Open(filePath.c_str(), "READ");
        if (!etaFile || etaFile->IsZombie()) {
            std::cerr << "Cannot open eta-diff file: " << filePath << std::endl;
            continue;
        }
        
        TH1D* h = dynamic_cast<TH1D*>(etaFile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
        if (!h) {
            std::cerr << "Cannot find histogram in: " << filePath << std::endl;
            etaFile->Close();
            delete etaFile;
            continue;
        }
        
        if (!hAggregate) {
            hAggregate = dynamic_cast<TH1D*>(h->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange)));
            hAggregate->SetDirectory(nullptr);
        } else {
            hAggregate->Add(h);
        }
        
        etaFile->Close();
        delete etaFile;
    }
    
    if (hAggregate) {
        hAggregate->Write();
        delete hAggregate;
    }
    
    // Now sum bootstrap samples across eta bins
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        TH1D* hBsSummed = nullptr;
        
        for (const auto& filePath : etaDiffFiles) {
            TFile* etaFile = TFile::Open(filePath.c_str(), "READ");
            if (!etaFile || etaFile->IsZombie()) {
                continue;
            }
            
            TH1D* h = dynamic_cast<TH1D*>(etaFile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)));
            if (!h) {
                etaFile->Close();
                delete etaFile;
                continue;
            }
            
            if (!hBsSummed) {
                hBsSummed = dynamic_cast<TH1D*>(h->Clone(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)));
                hBsSummed->SetDirectory(nullptr);
            } else {
                hBsSummed->Add(h);
            }
            
            etaFile->Close();
            delete etaFile;
        }
        
        if (hBsSummed) {
            hBsSummed->Write();
            delete hBsSummed;
        }
    }
    
    std::cout << "Created aggregate bootstrap file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root", 
        fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()) << std::endl;
    
    outFile->Close();
    delete outFile;
}
