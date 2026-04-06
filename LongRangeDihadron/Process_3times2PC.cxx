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
    // configList.push_back(ConfigUnit(kCent,
    // {InputUnit("LHC25af_pass2_637596", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_632504", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_642734", kFT0AFT0C, kTemplateFit, 0, 20)},
    // "LHC25af_pass2_642734"));

    // Old FT0A-FT0C ultra long range input (commented out)
    // configList.push_back(ConfigUnit(kCent,
    // {InputUnit("LHC25af_pass2_632504", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_637596", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_642734", kFT0AFT0C, kTemplateFit, 0, 20)},
    // "LHC25af_pass2_642734"));

    // Ne-Ne full-range combination (baseline)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25af_pass2_632504", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_637596", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_645746", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25af_pass2_645746"));

    // Ne-Ne ring combinations (TPC ring bins + side-specific FT0-FT0 ids)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25af_pass2_632504", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_631290", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_645746", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25af_pass2_innerRing"));

    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25af_pass2_637597", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_637594", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_645746", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25af_pass2_outerRing"));

    // O-O full-range combination (baseline)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25ae_pass2_644429", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_644429", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25ae_pass2_645657"));

    // O-O ring combinations
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25ae_pass2_638221", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_634099", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25ae_pass2_innerRing"));

    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25ae_pass2_634103", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_637591", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25ae_pass2_645657", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25ae_pass2_outerRing"));

    // p-O datasets (ad)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC25ad_pass2_644389", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25ad_pass2_644389", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25ad_pass2_644389", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC25ad_pass2_644389"));

    // p-p datasets (LHC24af)
    configList.push_back(ConfigUnit(kCent,
    {InputUnit("LHC24af_pass1_644663", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC24af_pass1_644663", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC24af_pass1_644663", kFT0AFT0C, kTemplateFit, 0, 20)},
    "LHC24af_pass1_644663"));

    // Ne-Ne Nch-dependent full-range combination (data: 10-50, template: 0-10)
    configList.push_back(ConfigUnit(kNch,
    {InputUnit("LHC25af_pass2_650316", kTPCFT0A, kTemplateFit, 10, 50), InputUnit("LHC25af_pass2_650317", kTPCFT0C, kTemplateFit, 10, 50), InputUnit("LHC25af_pass2_650315", kFT0AFT0C, kTemplateFit, 10, 50)},
    "LHC25af_pass2_650315_nch10_50"));

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

    // Default: use the same LR for both FT0 sides. Ring outputs can override per side.
    double lrC_v2 = lr_v2, lrC_v2_err = lr_v2_err;
    double lrA_v2 = lr_v2, lrA_v2_err = lr_v2_err;
    double lrC_v3 = lr_v3, lrC_v3_err = lr_v3_err;
    double lrA_v3 = lr_v3, lrA_v3_err = lr_v3_err;
    double lrC_v4 = lr_v4, lrC_v4_err = lr_v4_err;
    double lrA_v4 = lr_v4, lrA_v4_err = lr_v4_err;

    auto LoadLRFromDataset = [&](const std::string& datasetSuffix,
                                 double& outV2, double& outV2Err,
                                 double& outV3, double& outV3Err,
                                 double& outV4, double& outV4Err) {
        TFile* f = new TFile(Form("./TemplateFit/VnDelta_%s_%s_%s.root",
                                  datasetSuffix.c_str(), splitName.c_str(), DihadronCorrTypeName[kFT0AFT0C].c_str()), "READ");
        if (!f || !f->IsOpen()) {
            Printf("[3times2PC] Warning: cannot open side LR file for %s", datasetSuffix.c_str());
            if (f) {
                f->Close();
                delete f;
            }
            return false;
        }
        TH1D* hv2 = (TH1D*)f->Get("hV2");
        TH1D* hv3 = (TH1D*)f->Get("hV3");
        TH1D* hv4 = (TH1D*)f->Get("hV4");
        if (!hv2 || !hv3 || !hv4) {
            Printf("[3times2PC] Warning: missing hV2/hV3/hV4 in side LR file for %s", datasetSuffix.c_str());
            f->Close();
            delete f;
            return false;
        }

        const double v2 = hv2->GetBinContent(1);
        const double v2Err = hv2->GetBinError(1);
        const double v3 = hv3->GetBinContent(1);
        const double v3Err = hv3->GetBinError(1);
        const double v4 = hv4->GetBinContent(1);
        const double v4Err = hv4->GetBinError(1);

        auto isInvalidSideValue = [](double value, double error) {
            return (!std::isfinite(value) || !std::isfinite(error) || value < 0.0 || error >= 9.9);
        };
        if (isInvalidSideValue(v2, v2Err) || isInvalidSideValue(v3, v3Err) || isInvalidSideValue(v4, v4Err)) {
            Printf("[3times2PC] Warning: invalid side LR value(s) in %s (v2=%g+-%g, v3=%g+-%g, v4=%g+-%g)",
                   datasetSuffix.c_str(), v2, v2Err, v3, v3Err, v4, v4Err);
            f->Close();
            delete f;
            return false;
        }

        outV2 = v2;
        outV2Err = v2Err;
        outV3 = v3;
        outV3Err = v3Err;
        outV4 = v4;
        outV4Err = v4Err;
        f->Close();
        delete f;
        return true;
    };

    auto LoadLRFromDatasetWithOptionalId = [&](const std::string& datasetSuffix, const std::string& optionalId,
                                               double& outV2, double& outV2Err,
                                               double& outV3, double& outV3Err,
                                               double& outV4, double& outV4Err) {
        if (!optionalId.empty()) {
            if (LoadLRFromDataset(Form("%s_%s", datasetSuffix.c_str(), optionalId.c_str()),
                                  outV2, outV2Err, outV3, outV3Err, outV4, outV4Err)) {
                return true;
            }
        }
        return LoadLRFromDataset(datasetSuffix, outV2, outV2Err, outV3, outV3Err, outV4, outV4Err);
    };

    if (outputFileName.find("innerRing") != std::string::npos) {
        if (outputFileName.find("LHC25af_pass2") != std::string::npos) {
            // Ne-Ne inner ring: C from id50560 (reject C outside), A from id50559 (reject A outside)
            LoadLRFromDataset("LHC25af_pass2_646139_id50560", lrC_v2, lrC_v2_err, lrC_v3, lrC_v3_err, lrC_v4, lrC_v4_err);
            LoadLRFromDataset("LHC25af_pass2_646139_id50559", lrA_v2, lrA_v2_err, lrA_v3, lrA_v3_err, lrA_v4, lrA_v4_err);
        } else if (outputFileName.find("LHC25ae_pass2") != std::string::npos) {
            // O-O inner ring: C from 648800, A from 648799 (try id50564 variant first, then plain dataset)
            LoadLRFromDatasetWithOptionalId("LHC25ae_pass2_648800", "id50564", lrC_v2, lrC_v2_err, lrC_v3, lrC_v3_err, lrC_v4, lrC_v4_err);
            LoadLRFromDatasetWithOptionalId("LHC25ae_pass2_648799", "id50564", lrA_v2, lrA_v2_err, lrA_v3, lrA_v3_err, lrA_v4, lrA_v4_err);
        }
    } else if (outputFileName.find("outerRing") != std::string::npos) {
        if (outputFileName.find("LHC25af_pass2") != std::string::npos) {
            // Ne-Ne outer ring: C from id50562 (reject C inside), A from id50561 (reject A inside)
            LoadLRFromDataset("LHC25af_pass2_646139_id50562", lrC_v2, lrC_v2_err, lrC_v3, lrC_v3_err, lrC_v4, lrC_v4_err);
            LoadLRFromDataset("LHC25af_pass2_646139_id50561", lrA_v2, lrA_v2_err, lrA_v3, lrA_v3_err, lrA_v4, lrA_v4_err);
        } else if (outputFileName.find("LHC25ae_pass2") != std::string::npos) {
            // O-O outer ring: C from 648788, A from 644433 (try id50564 variant first, then plain dataset)
            LoadLRFromDatasetWithOptionalId("LHC25ae_pass2_648788", "id50564", lrC_v2, lrC_v2_err, lrC_v3, lrC_v3_err, lrC_v4, lrC_v4_err);
            LoadLRFromDatasetWithOptionalId("LHC25ae_pass2_644433", "id50564", lrA_v2, lrA_v2_err, lrA_v3, lrA_v3_err, lrA_v4, lrA_v4_err);
        }
    }

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

    auto ComputeHistogramAverage = [](TH1D* sourceHist) {
        double sum = 0.0;
        double sumErr2 = 0.0;
        int nValid = 0;
        if (!sourceHist) return std::make_tuple(-1.0, 10.0, 0);
        for (int ibin = 1; ibin <= sourceHist->GetNbinsX(); ++ibin) {
            double value = sourceHist->GetBinContent(ibin);
            double error = sourceHist->GetBinError(ibin);
            if (!std::isfinite(value) || !std::isfinite(error) || value < 0.0 || error >= 9.9) continue;
            sum += value;
            sumErr2 += error * error;
            ++nValid;
        }
        if (nValid <= 0) return std::make_tuple(-1.0, 10.0, 0);
        return std::make_tuple(sum / nValid, std::sqrt(sumErr2) / nValid, nValid);
    };

    auto BuildSideSummaryFromLMRLR = [&](TH1D* hLM, TH1D* hMR,
                                         double lrForC, double lrForCErr,
                                         double lrForA, double lrForAErr,
                                         const char* name, const char* title) {
        TH1D* sideHist = new TH1D(name, title, 2, 0.5, 2.5);
        sideHist->GetXaxis()->SetBinLabel(1, "FT0C [-3.3,-2.1]"); // Left
        sideHist->GetXaxis()->SetBinLabel(2, "FT0A [3.5,4.9]"); // Right

        double lmAvg = -1.0, lmAvgErr = 10.0;
        double mrAvg = -1.0, mrAvgErr = 10.0;
        int lmValidBins = 0;
        int mrValidBins = 0;
        std::tie(lmAvg, lmAvgErr, lmValidBins) = ComputeHistogramAverage(hLM);
        std::tie(mrAvg, mrAvgErr, mrValidBins) = ComputeHistogramAverage(hMR);

        // Calculate both values using the full TPC coverage on each side
        double vC = -1.0, vCerr = 10.0, vA = -1.0, vAerr = 10.0;
        if (lmValidBins > 0 && mrValidBins > 0) {
            // LM=t*a, MR=t*c, LR=a*c  => c=sqrt(MR*LR/LM), a=sqrt(LM*LR/MR)
            vC = Get3times2PC(mrAvg, lrForC, lmAvg);
            vCerr = Get3times2PC_Error(mrAvg, mrAvgErr, lrForC, lrForCErr, lmAvg, lmAvgErr);
            vA = Get3times2PC(lmAvg, lrForA, mrAvg);
            vAerr = Get3times2PC_Error(lmAvg, lmAvgErr, lrForA, lrForAErr, mrAvg, mrAvgErr);
        }
        // Keep physical mapping: bin 1 is FT0C, bin 2 is FT0A
        sideHist->SetBinContent(1, vC);
        sideHist->SetBinError(1, vCerr);
        sideHist->SetBinContent(2, vA);
        sideHist->SetBinError(2, vAerr);

        return sideHist;
    };

    TH1D* hV2_Sides = BuildSideSummaryFromLMRLR(hLM_v2, hMR_v2, lrC_v2, lrC_v2_err, lrA_v2, lrA_v2_err, "hV2_Sides", "v_{2};side;v_{2}");
    TH1D* hV3_Sides = BuildSideSummaryFromLMRLR(hLM_v3, hMR_v3, lrC_v3, lrC_v3_err, lrA_v3, lrA_v3_err, "hV3_Sides", "v_{3};side;v_{3}");
    TH1D* hV4_Sides = BuildSideSummaryFromLMRLR(hLM_v4, hMR_v4, lrC_v4, lrC_v4_err, lrA_v4, lrA_v4_err, "hV4_Sides", "v_{4};side;v_{4}");

    // Combined histogram-style output: FT0C + 16 eta bins + FT0A in one object
    int nCombinedBins = nEtaBins + 2;
    TH1D* hV2_Combined = new TH1D("hV2_Combined", "v_{2};bin;v_{2}", nCombinedBins, 0.5, nCombinedBins + 0.5);
    TH1D* hV3_Combined = new TH1D("hV3_Combined", "v_{3};bin;v_{3}", nCombinedBins, 0.5, nCombinedBins + 0.5);
    TH1D* hV4_Combined = new TH1D("hV4_Combined", "v_{4};bin;v_{4}", nCombinedBins, 0.5, nCombinedBins + 0.5);

    hV2_Combined->GetXaxis()->SetBinLabel(1, "FT0C[-3.3,-2.1]");
    hV3_Combined->GetXaxis()->SetBinLabel(1, "FT0C[-3.3,-2.1]");
    hV4_Combined->GetXaxis()->SetBinLabel(1, "FT0C[-3.3,-2.1]");
    hV2_Combined->SetBinContent(1, hV2_Sides->GetBinContent(1));
    hV2_Combined->SetBinError(1, hV2_Sides->GetBinError(1));
    hV3_Combined->SetBinContent(1, hV3_Sides->GetBinContent(1));
    hV3_Combined->SetBinError(1, hV3_Sides->GetBinError(1));
    hV4_Combined->SetBinContent(1, hV4_Sides->GetBinContent(1));
    hV4_Combined->SetBinError(1, hV4_Sides->GetBinError(1));

    for (int ibin = 1; ibin <= nEtaBins; ++ibin) {
        int combinedBin = ibin + 1;
        double etaLow = hV2->GetXaxis()->GetBinLowEdge(ibin);
        double etaUp = hV2->GetXaxis()->GetBinUpEdge(ibin);
        TString etaLabel = Form("[%.1f,%.1f]", etaLow, etaUp);
        hV2_Combined->GetXaxis()->SetBinLabel(combinedBin, etaLabel);
        hV3_Combined->GetXaxis()->SetBinLabel(combinedBin, etaLabel);
        hV4_Combined->GetXaxis()->SetBinLabel(combinedBin, etaLabel);

        hV2_Combined->SetBinContent(combinedBin, hV2->GetBinContent(ibin));
        hV2_Combined->SetBinError(combinedBin, hV2->GetBinError(ibin));
        hV3_Combined->SetBinContent(combinedBin, hV3->GetBinContent(ibin));
        hV3_Combined->SetBinError(combinedBin, hV3->GetBinError(ibin));
        hV4_Combined->SetBinContent(combinedBin, hV4->GetBinContent(ibin));
        hV4_Combined->SetBinError(combinedBin, hV4->GetBinError(ibin));
    }

    int lastCombinedBin = nCombinedBins;
    hV2_Combined->GetXaxis()->SetBinLabel(lastCombinedBin, "FT0A[3.5,4.9]");
    hV3_Combined->GetXaxis()->SetBinLabel(lastCombinedBin, "FT0A[3.5,4.9]");
    hV4_Combined->GetXaxis()->SetBinLabel(lastCombinedBin, "FT0A[3.5,4.9]");
    hV2_Combined->SetBinContent(lastCombinedBin, hV2_Sides->GetBinContent(2));
    hV2_Combined->SetBinError(lastCombinedBin, hV2_Sides->GetBinError(2));
    hV3_Combined->SetBinContent(lastCombinedBin, hV3_Sides->GetBinContent(2));
    hV3_Combined->SetBinError(lastCombinedBin, hV3_Sides->GetBinError(2));
    hV4_Combined->SetBinContent(lastCombinedBin, hV4_Sides->GetBinContent(2));
    hV4_Combined->SetBinError(lastCombinedBin, hV4_Sides->GetBinError(2));

    hV2->Write();
    hV3->Write();
    hV4->Write();
    hV2_Sides->Write();
    hV3_Sides->Write();
    hV4_Sides->Write();
    hV2_Combined->Write();
    hV3_Combined->Write();
    hV4_Combined->Write();

    std::cout << "[3times2PC] Side summary v2: FT0C_left=" << hV2_Sides->GetBinContent(1)
              << " +/- " << hV2_Sides->GetBinError(1)
              << ", FT0A_right=" << hV2_Sides->GetBinContent(2)
              << " +/- " << hV2_Sides->GetBinError(2) << std::endl;
    std::cout << "[3times2PC] Wrote combined histogram with FT0C + 16 TPC bins + FT0A" << std::endl;

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
    if(isample >= 0) {
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
            // EtaDiff file: use inverse-variance weighted average across eta bins.
            // This avoids giving noisy bins equal influence and helps recover system trends.
            double numV2 = 0.0, denV2 = 0.0;
            double numV3 = 0.0, denV3 = 0.0;
            double numV4 = 0.0, denV4 = 0.0;
            int validV2 = 0, validV3 = 0, validV4 = 0;
            int nbins = hv2->GetNbinsX();
            for (int i = 1; i <= nbins; ++i) {
                double val2 = hv2->GetBinContent(i);
                double err2 = hv2->GetBinError(i);
                if (std::isfinite(val2) && std::isfinite(err2) && err2 > 0.0 && err2 < 9.9 && val2 >= 0.0) {
                    double w2 = 1.0 / (err2 * err2);
                    numV2 += w2 * val2;
                    denV2 += w2;
                    ++validV2;
                }

                double val3 = hv3->GetBinContent(i);
                double err3 = hv3->GetBinError(i);
                if (std::isfinite(val3) && std::isfinite(err3) && err3 > 0.0 && err3 < 9.9 && val3 >= 0.0) {
                    double w3 = 1.0 / (err3 * err3);
                    numV3 += w3 * val3;
                    denV3 += w3;
                    ++validV3;
                }

                double val4 = hv4->GetBinContent(i);
                double err4 = hv4->GetBinError(i);
                if (std::isfinite(val4) && std::isfinite(err4) && err4 > 0.0 && err4 < 9.9 && val4 >= 0.0) {
                    double w4 = 1.0 / (err4 * err4);
                    numV4 += w4 * val4;
                    denV4 += w4;
                    ++validV4;
                }
            }

            if (denV2 > 0.0 && validV2 > 0) {
                v2 = numV2 / denV2;
                v2_err = std::sqrt(1.0 / denV2);
            }
            if (denV3 > 0.0 && validV3 > 0) {
                v3 = numV3 / denV3;
                v3_err = std::sqrt(1.0 / denV3);
            }
            if (denV4 > 0.0 && validV4 > 0) {
                v4 = numV4 / denV4;
                v4_err = std::sqrt(1.0 / denV4);
            }
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
