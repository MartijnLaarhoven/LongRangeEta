/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 13:51:01 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 17:03:59
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"
#include "./include/ErrorPropagation.h"
#include "./include/plotting.h"

// define struct
struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Int_t method;
    Int_t minRange;
    Int_t maxRange;
    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Int_t _method, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), method(_method), minRange(_minRange), maxRange(_maxRange) {}
};

struct ConfigUnit {
    Bool_t isNch;
    Bool_t isEtaDiff;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    Bool_t constructed;
    ConfigUnit(Bool_t _isNch, std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), isEtaDiff(false), dataList(_dataList), outputFileName(_outputFileName) {
            constructed = false;
            if (_dataList.size() != 3) {
                constructed = false;
                Printf("The input vector don't have a size of 3!");
                return;
            }
            if (!(_dataList[0].corrType == kTPCFT0A && _dataList[1].corrType == kTPCFT0C && _dataList[2].corrType == kFT0AFT0C)) {
                constructed = false;
                Printf("This calculation only valid for the sequence of kTPCFT0A, kTPCFT0C, kFT0AFT0C");
                return;
            }
            if (!(_dataList[0].minRange == _dataList[1].minRange && _dataList[0].minRange == _dataList[2].minRange && _dataList[0].maxRange == _dataList[1].maxRange && _dataList[0].maxRange == _dataList[2].maxRange)) {
                Printf("Warning: you are not using the same Cent/Nch range for one config:");
                Printf("kTPCFT0A %d %d, kTPCFT0C %d %d, kFT0AFT0C %d %d", _dataList[0].minRange, _dataList[0].maxRange, _dataList[1].minRange, _dataList[1].maxRange, _dataList[2].minRange, _dataList[2].maxRange);
            }
            constructed = true;
        }
};

struct VnUnit {
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err) {}
};

// 3times2PC processing
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* Vn3times2PC(Bool_t isNch, Bool_t isPtDiff, std::vector<InputUnit> dataList, std::vector<TFile*>& vnDeltaFiles, Double_t pTMin=0, Double_t pTMax=0);
VnUnit* GetResultsFromVnDeltaFiles(int isample, TFile* file, Bool_t isPtDiff, Int_t minRange, Int_t maxRange,  Double_t pTMin=0, Double_t pTMax=0);
double Get3times2PC(double vnDelta_LM, double vnDelta_MR, double vnDelta_LR);
double Get3times2PC_Error(double LM, double LM_err, double MR, double MR_err, double LR, double LR_err);
VnUnit* Get3times2PC(VnUnit* LMvalues, VnUnit* MRvalues, VnUnit* LRvalues);

std::string collisionSystemName = "Unknown";

void Process_3times2PC() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    
    // Create output directory
    gSystem->Exec("mkdir -p ./3times2PC");
    
    std::vector<ConfigUnit> configList;

    // LM, MR, LR
    // Ne-Ne datasets (af)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25af_pass2_632504", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_637596", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_642734", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25af_pass2_642734"));

    // O-O datasets (ae) - Disabled while focusing on af dataset debugging
    // configList.push_back(ConfigUnit(kCent,
    // {InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_645320", kFT0AFT0C, kTemplateFit, 0, 20)},
    // "LHC25ae_pass2_645320"));

    for (auto config : configList) {
        if (!config.constructed) continue;
        if (!config.dataList.empty()) {
            collisionSystemName = GetCollisionSystemNameFromDataset(config.dataList[0].fileNameSuffix);
        }
        ProcessConfig(config.isNch, config.dataList, config.outputFileName);
    }
}

void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    if (dataList.size() != 3) {
        Printf("The dataList don't have a size of 3!");
        return;
    }
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    auto dataLM = dataList[0];
    auto dataMR = dataList[1];
    auto dataLR = dataList[2];

    TFile* vnDeltaFile_LM = new TFile(Form("./%s/EtaDiff/VnDelta_%s_%s_%i_%i_%s.root",
        DihadronMethodName[dataLM.method].c_str(), dataLM.fileNameSuffix.c_str(), splitName.c_str(), dataLM.minRange, dataLM.maxRange, DihadronCorrTypeName[dataLM.corrType].c_str()), "READ");
    TFile* vnDeltaFile_MR = new TFile(Form("./%s/EtaDiff/VnDelta_%s_%s_%i_%i_%s.root",
        DihadronMethodName[dataMR.method].c_str(), dataMR.fileNameSuffix.c_str(), splitName.c_str(), dataMR.minRange, dataMR.maxRange, DihadronCorrTypeName[dataMR.corrType].c_str()), "READ");
    TFile* vnDeltaFile_LR = new TFile(Form("./%s/VnDelta_%s_%s_%s.root",
        DihadronMethodName[dataLR.method].c_str(), dataLR.fileNameSuffix.c_str(), splitName.c_str(), DihadronCorrTypeName[dataLR.corrType].c_str()), "READ");

    if (!vnDeltaFile_LM || !vnDeltaFile_LM->IsOpen()) {
        Printf("can not get LM EtaDiff file!");
        return;
    }
    if (!vnDeltaFile_MR || !vnDeltaFile_MR->IsOpen()) {
        Printf("can not get MR EtaDiff file!");
        vnDeltaFile_LM->Close();
        delete vnDeltaFile_LM;
        return;
    }
    if (!vnDeltaFile_LR || !vnDeltaFile_LR->IsOpen()) {
        Printf("can not get LR file!");
        vnDeltaFile_LM->Close();
        vnDeltaFile_MR->Close();
        delete vnDeltaFile_LM;
        delete vnDeltaFile_MR;
        return;
    }

    TH1D* hLM_v2 = (TH1D*)vnDeltaFile_LM->Get("hV2");
    TH1D* hLM_v3 = (TH1D*)vnDeltaFile_LM->Get("hV3");
    TH1D* hLM_v4 = (TH1D*)vnDeltaFile_LM->Get("hV4");
    TH1D* hMR_v2 = (TH1D*)vnDeltaFile_MR->Get("hV2");
    TH1D* hMR_v3 = (TH1D*)vnDeltaFile_MR->Get("hV3");
    TH1D* hMR_v4 = (TH1D*)vnDeltaFile_MR->Get("hV4");
    TH1D* hLR_v2 = (TH1D*)vnDeltaFile_LR->Get("hV2");
    TH1D* hLR_v3 = (TH1D*)vnDeltaFile_LR->Get("hV3");
    TH1D* hLR_v4 = (TH1D*)vnDeltaFile_LR->Get("hV4");

    if (!hLM_v2 || !hLM_v3 || !hLM_v4 || !hMR_v2 || !hMR_v3 || !hMR_v4 || !hLR_v2 || !hLR_v3 || !hLR_v4) {
        Printf("Can not get hV2/hV3/hV4 from one or more input files!");
        vnDeltaFile_LM->Close();
        vnDeltaFile_MR->Close();
        vnDeltaFile_LR->Close();
        delete vnDeltaFile_LM;
        delete vnDeltaFile_MR;
        delete vnDeltaFile_LR;
        return;
    }

    int nEtaBinsLM = hLM_v2->GetNbinsX();
    int nEtaBinsMR = hMR_v2->GetNbinsX();
    int nEtaBins = std::min(nEtaBinsLM, nEtaBinsMR);
    if (nEtaBins <= 0) {
        Printf("No eta bins found in EtaDiff files!");
        vnDeltaFile_LM->Close();
        vnDeltaFile_MR->Close();
        vnDeltaFile_LR->Close();
        delete vnDeltaFile_LM;
        delete vnDeltaFile_MR;
        delete vnDeltaFile_LR;
        return;
    }

    double lr_v2 = hLR_v2->GetBinContent(1);
    double lr_v2_err = hLR_v2->GetBinError(1);
    double lr_v3 = hLR_v3->GetBinContent(1);
    double lr_v3_err = hLR_v3->GetBinError(1);
    double lr_v4 = hLR_v4->GetBinContent(1);
    double lr_v4_err = hLR_v4->GetBinError(1);

    TFile outputFile(Form("./3times2PC/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), dataLR.minRange, dataLR.maxRange), "RECREATE");
    TH1D* hV2 = dynamic_cast<TH1D*>(hLM_v2->Clone("hV2"));
    TH1D* hV3 = dynamic_cast<TH1D*>(hLM_v3->Clone("hV3"));
    TH1D* hV4 = dynamic_cast<TH1D*>(hLM_v4->Clone("hV4"));
    hV2->Reset();
    hV3->Reset();
    hV4->Reset();
    hV2->SetTitle("v_{2}");
    hV3->SetTitle("v_{3}");
    hV4->SetTitle("v_{4}");

    std::cout << "[3times2PC] Processing eta-differential 3times2PC with " << nEtaBins << " bins" << std::endl;

    for (int ibin = 1; ibin <= nEtaBins; ++ibin) {
        double etaMin = hLM_v2->GetXaxis()->GetBinLowEdge(ibin);
        double etaMax = hLM_v2->GetXaxis()->GetBinUpEdge(ibin);

        double lm_v2 = hLM_v2->GetBinContent(ibin);
        double lm_v2_err = hLM_v2->GetBinError(ibin);
        double mr_v2 = hMR_v2->GetBinContent(ibin);
        double mr_v2_err = hMR_v2->GetBinError(ibin);

        double lm_v3 = hLM_v3->GetBinContent(ibin);
        double lm_v3_err = hLM_v3->GetBinError(ibin);
        double mr_v3 = hMR_v3->GetBinContent(ibin);
        double mr_v3_err = hMR_v3->GetBinError(ibin);

        double lm_v4 = hLM_v4->GetBinContent(ibin);
        double lm_v4_err = hLM_v4->GetBinError(ibin);
        double mr_v4 = hMR_v4->GetBinContent(ibin);
        double mr_v4_err = hMR_v4->GetBinError(ibin);

        double v2 = Get3times2PC(lm_v2, mr_v2, lr_v2);
        double v2_err = Get3times2PC_Error(lm_v2, lm_v2_err, mr_v2, mr_v2_err, lr_v2, lr_v2_err);
        double v3 = Get3times2PC(lm_v3, mr_v3, lr_v3);
        double v3_err = Get3times2PC_Error(lm_v3, lm_v3_err, mr_v3, mr_v3_err, lr_v3, lr_v3_err);
        double v4 = Get3times2PC(lm_v4, mr_v4, lr_v4);
        double v4_err = Get3times2PC_Error(lm_v4, lm_v4_err, mr_v4, mr_v4_err, lr_v4, lr_v4_err);

        hV2->SetBinContent(ibin, v2);
        hV2->SetBinError(ibin, v2_err);
        hV3->SetBinContent(ibin, v3);
        hV3->SetBinError(ibin, v3_err);
        hV4->SetBinContent(ibin, v4);
        hV4->SetBinError(ibin, v4_err);

        std::cout << Form("[3times2PC] Eta [%.1f, %.1f]: v2=%g +/- %g", etaMin, etaMax, v2, v2_err) << std::endl;
    }

    auto BuildSideSummary = [](TH1D* sourceHist, const char* name, const char* title) {
        TH1D* sideHist = new TH1D(name, title, 2, 0.5, 2.5);
        sideHist->GetXaxis()->SetBinLabel(1, "FT0C_left");
        sideHist->GetXaxis()->SetBinLabel(2, "FT0A_right");

        double leftSum = 0.0;
        double leftErr2 = 0.0;
        int leftCount = 0;
        double rightSum = 0.0;
        double rightErr2 = 0.0;
        int rightCount = 0;

        for (int ibin = 1; ibin <= sourceHist->GetNbinsX(); ++ibin) {
            double center = sourceHist->GetXaxis()->GetBinCenter(ibin);
            double value = sourceHist->GetBinContent(ibin);
            double error = sourceHist->GetBinError(ibin);

            if (!std::isfinite(value) || !std::isfinite(error) || value < 0.0 || error >= 9.9) {
                continue;
            }

            if (center < 0) {
                leftSum += value;
                leftErr2 += error * error;
                leftCount++;
            } else if (center > 0) {
                rightSum += value;
                rightErr2 += error * error;
                rightCount++;
            }
        }

        if (leftCount > 0) {
            sideHist->SetBinContent(1, leftSum / leftCount);
            sideHist->SetBinError(1, std::sqrt(leftErr2) / leftCount);
        } else {
            sideHist->SetBinContent(1, -1.0);
            sideHist->SetBinError(1, 10.0);
        }

        if (rightCount > 0) {
            sideHist->SetBinContent(2, rightSum / rightCount);
            sideHist->SetBinError(2, std::sqrt(rightErr2) / rightCount);
        } else {
            sideHist->SetBinContent(2, -1.0);
            sideHist->SetBinError(2, 10.0);
        }

        return sideHist;
    };

    TH1D* hV2_Sides = BuildSideSummary(hV2, "hV2_Sides", "v_{2};side;v_{2}");
    TH1D* hV3_Sides = BuildSideSummary(hV3, "hV3_Sides", "v_{3};side;v_{3}");
    TH1D* hV4_Sides = BuildSideSummary(hV4, "hV4_Sides", "v_{4};side;v_{4}");

    hV2->Write();
    hV3->Write();
    hV4->Write();
    hV2_Sides->Write();
    hV3_Sides->Write();
    hV4_Sides->Write();

    std::cout << "[3times2PC] Side summary v2: FT0C_left=" << hV2_Sides->GetBinContent(1)
              << " +/- " << hV2_Sides->GetBinError(1)
              << ", FT0A_right=" << hV2_Sides->GetBinContent(2)
              << " +/- " << hV2_Sides->GetBinError(2) << std::endl;

    outputFile.Close();

    vnDeltaFile_LM->Close();
    vnDeltaFile_MR->Close();
    vnDeltaFile_LR->Close();
    delete vnDeltaFile_LM;
    delete vnDeltaFile_MR;
    delete vnDeltaFile_LR;

    std::cout << "ALICE " << collisionSystemName << " 3times2PC (eta differential)" << std::endl;
    std::cout << "Output file: " << Form("./3times2PC/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), dataLR.minRange, dataLR.maxRange) << std::endl;
}

VnUnit* Vn3times2PC(Bool_t isNch, Bool_t isPtDiff, std::vector<InputUnit> dataList, std::vector<TFile*>& vnDeltaFiles, Double_t pTMin=0, Double_t pTMax=0) {
    if (dataList.size() != 3) {
        Printf("The dataList don't have a size of 3!");
        return nullptr;
    }
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // central values
    VnUnit* resultsLM = nullptr;
    VnUnit* resultsMR = nullptr;
    VnUnit* resultsLR = nullptr;
    resultsLM = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[0], isPtDiff, dataList[0].minRange, dataList[0].maxRange, pTMin, pTMax);
    resultsMR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[1], isPtDiff, dataList[1].minRange, dataList[1].maxRange, pTMin, pTMax);
    resultsLR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[2], kFALSE, dataList[2].minRange, dataList[2].maxRange, pTMin, pTMax);
    if (!resultsLM || !resultsMR || !resultsLR) {
        Printf("Can not get results from vnDelta files!");
        return nullptr;
    }
    std::cout << "[3times2PC Debug] LM (TPC_FT0A): v2=" << resultsLM->v2 << " +/- " << resultsLM->v2_err << std::endl;
    std::cout << "[3times2PC Debug] MR (TPC_FT0C): v2=" << resultsMR->v2 << " +/- " << resultsMR->v2_err << std::endl;
    std::cout << "[3times2PC Debug] LR (FT0_FT0C): v2=" << resultsLR->v2 << " +/- " << resultsLR->v2_err << std::endl;
    VnUnit* vnResult = Get3times2PC(resultsLM, resultsMR, resultsLR);

    // Store aggregated result errors for fallback
    double v2_err_agg = vnResult->v2_err;
    double v3_err_agg = vnResult->v3_err;
    double v4_err_agg = vnResult->v4_err;

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=3;//v22,v32,v42
    int NofSample = maxSample;
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    int validSampleCount = 0;
    for(int sample=0;sample<NofSample;sample++) {
        VnUnit* tempLM = nullptr;
        VnUnit* tempMR = nullptr;
        VnUnit* tempLR = nullptr;
        tempLM = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[0], isPtDiff, dataList[0].minRange, dataList[0].maxRange, pTMin, pTMax);
        tempMR = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[1], isPtDiff, dataList[1].minRange, dataList[1].maxRange, pTMin, pTMax);
        tempLR = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[2], kFALSE, dataList[2].minRange, dataList[2].maxRange, pTMin, pTMax);
        if (!tempLM || !tempMR || !tempLR) {
            Printf("Can not get subsample %d results from vnDelta files!", sample);
            return nullptr;
        }
        VnUnit* vnTemp = Get3times2PC(tempLM, tempMR, tempLR);
        // Count valid samples (not sentinel values v=-1, err=10)
        if (!(vnTemp->v2 == -1.0 && vnTemp->v2_err == 10.0)) {
            validSampleCount++;
        }
        ValueArray[0][sample][0] = vnTemp->v2;
        ValueErrorArray[0][sample][0] = vnTemp->v2_err;
        ValueArray[1][sample][0] = vnTemp->v3;
        ValueErrorArray[1][sample][0] = vnTemp->v3_err;
        ValueArray[2][sample][0] = vnTemp->v4;
        ValueErrorArray[2][sample][0] = vnTemp->v4_err;
        delete tempLM;
        delete tempMR;
        delete tempLR;
        delete vnTemp;
    }

    std::cout << "[3times2PC] Valid bootstrap samples: " << validSampleCount << "/" << NofSample << std::endl;

    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.);
    }

    // Check if bootstrap error calculation failed (< 3 valid samples) and fallback to aggregated result errors
    if (validSampleCount < 3) {
        std::cout << "[3times2PC] Too few valid bootstrap samples (" << validSampleCount << "), using aggregated result errors instead" << std::endl;
        vnResult->v2_err = v2_err_agg;
        vnResult->v3_err = v3_err_agg;
        vnResult->v4_err = v4_err_agg;
    } else {
        vnResult->v2_err = ErrorArray[0][0];
        vnResult->v3_err = ErrorArray[1][0];
        vnResult->v4_err = ErrorArray[2][0];
    }

    std::cout << "print result: " << dataList[2].fileNameSuffix << std::endl;
    std::cout << "v2: " << vnResult->v2 << " +/- " << vnResult->v2_err << std::endl;
    std::cout << "v3: " << vnResult->v3 << " +/- " << vnResult->v3_err << std::endl;
    std::cout << "v4: " << vnResult->v4 << " +/- " << vnResult->v4_err << std::endl;

    return vnResult;
}


VnUnit* GetResultsFromVnDeltaFiles(int isample, TFile* file, Bool_t isPtDiff, Int_t minRange, Int_t maxRange,  Double_t pTMin=0, Double_t pTMax=0) {
    TH1D* hv2 = nullptr;
    TH1D* hv3 = nullptr;
    TH1D* hv4 = nullptr;
    double v2 = -1;
    double v3 = -1;
    double v4 = -1;
    double v2_err = 10;
    double v3_err = 10;
    double v4_err = 10;
    if(isample > 0) {
        hv2 = (TH1D*)file->Get(Form("Subsamples/hV2_subsample_%d", isample));
        hv3 = (TH1D*)file->Get(Form("Subsamples/hV3_subsample_%d", isample));
        hv4 = (TH1D*)file->Get(Form("Subsamples/hV4_subsample_%d", isample));
    } else {
        hv2 = (TH1D*)file->Get("hV2");
        hv3 = (TH1D*)file->Get("hV3");
        hv4 = (TH1D*)file->Get("hV4");
    }
    if (!hv2 || !hv3 || !hv4) {
        Printf("Can not get histograms from file!");
        return nullptr;
    }
    if (isPtDiff) {
        // Get bincontent from the bin of [pTMin, pTMax] (exactly that one bin)
        int binIndex = hv2->GetXaxis()->FindBin((pTMin + pTMax) / 2.0);
        v2 = hv2->GetBinContent(binIndex);
        v2_err = hv2->GetBinError(binIndex);
        v3 = hv3->GetBinContent(binIndex);
        v3_err = hv3->GetBinError(binIndex);
        v4 = hv4->GetBinContent(binIndex);
        v4_err = hv4->GetBinError(binIndex);
    } else {
        // Check if histogram has multiple bins (eta-differential case)
        // If so, average all bins for the full eta range
        if (hv2->GetNbinsX() > 1) {
            // EtaDiff file: average across all eta bins for full eta-range result
            double sumV2 = 0, sumV2Err2 = 0;
            double sumV3 = 0, sumV3Err2 = 0;
            double sumV4 = 0, sumV4Err2 = 0;
            int nbins = hv2->GetNbinsX();
            for (int i = 1; i <= nbins; ++i) {
                double val = hv2->GetBinContent(i);
                double err = hv2->GetBinError(i);
                sumV2 += val;
                sumV2Err2 += err * err;
                
                val = hv3->GetBinContent(i);
                err = hv3->GetBinError(i);
                sumV3 += val;
                sumV3Err2 += err * err;
                
                val = hv4->GetBinContent(i);
                err = hv4->GetBinError(i);
                sumV4 += val;
                sumV4Err2 += err * err;
            }
            v2 = sumV2 / nbins;
            v2_err = std::sqrt(sumV2Err2) / nbins;  // Propagated error
            v3 = sumV3 / nbins;
            v3_err = std::sqrt(sumV3Err2) / nbins;
            v4 = sumV4 / nbins;
            v4_err = std::sqrt(sumV4Err2) / nbins;
        } else {
            // Single-bin histogram: use the single value directly
            v2 = hv2->GetBinContent(1);
            v2_err = hv2->GetBinError(1);
            v3 = hv3->GetBinContent(1);
            v3_err = hv3->GetBinError(1);
            v4 = hv4->GetBinContent(1);
            v4_err = hv4->GetBinError(1);
        }
    }
    VnUnit* vnResult = new VnUnit(v2, v2_err, v3, v3_err, v4, v4_err);
    return vnResult;
}

double Get3times2PC(double vnDelta_LM, double vnDelta_MR, double vnDelta_LR) {
    double vn = (vnDelta_LM * vnDelta_MR) / vnDelta_LR;
    if (vn >= 0) {
        return std::sqrt(vn);
    } else {
        return -1;
    }
}

double Get3times2PC_Error(double LM, double LM_err, double MR, double MR_err, double LR, double LR_err) {
    if (LM < 0 || MR < 0 || LR < 0) {
        return 10.;
    }
    double part1 = 0.5*pow(LM, -0.5)*pow(MR, 0.5)*pow(LR,-0.5)*LM_err;
    part1 = pow(part1, 2);
    double part2 = 0.5*pow(LM, 0.5)*pow(MR, -0.5)*pow(LR,-0.5)*MR_err;
    part2 = pow(part2, 2);
    double part3 = 0.5*pow(LM, 0.5)*pow(MR, 0.5)*pow(LR,-1.5)*LR_err;
    part3 = pow(part3, 2);
    return sqrt(part1 + part2 + part3);
}

VnUnit* Get3times2PC(VnUnit* LMvalues, VnUnit* MRvalues, VnUnit* LRvalues) {
    double v2 = Get3times2PC(LMvalues->v2, MRvalues->v2, LRvalues->v2);
    double v2_err = Get3times2PC_Error(LMvalues->v2, LMvalues->v2_err, MRvalues->v2, MRvalues->v2_err, LRvalues->v2, LRvalues->v2_err);
    double v3 = Get3times2PC(LMvalues->v3, MRvalues->v3, LRvalues->v3);
    double v3_err = Get3times2PC_Error(LMvalues->v3, LMvalues->v3_err, MRvalues->v3, MRvalues->v3_err, LRvalues->v3, LRvalues->v3_err);
    double v4 = Get3times2PC(LMvalues->v4, MRvalues->v4, LRvalues->v4);
    double v4_err = Get3times2PC_Error(LMvalues->v4, LMvalues->v4_err, MRvalues->v4, MRvalues->v4_err, LRvalues->v4, LRvalues->v4_err);
    VnUnit* vnResult = new VnUnit(v2, v2_err, v3, v3_err, v4, v4_err);
    return vnResult;
}
