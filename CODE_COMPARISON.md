# Side-by-Side Code Comparison - Key Functions

## 1. VnUnit Structure Comparison

### REFERENCE
```cpp
struct VnUnit {
    Double_t v1;
    Double_t v1_err;
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    std::vector<Double_t> subsample_v1;
    std::vector<Double_t> subsample_v1_err;
    std::vector<Double_t> subsample_v2;
    std::vector<Double_t> subsample_v2_err;
    std::vector<Double_t> subsample_v3;
    std::vector<Double_t> subsample_v3_err;
    std::vector<Double_t> subsample_v4;
    std::vector<Double_t> subsample_v4_err;
    VnUnit(Double_t _v1, Double_t _v1_err, Double_t _v2, Double_t _v2_err, 
           Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) :
        v1(_v1), v1_err(_v1_err), v2(_v2), v2_err(_v2_err), 
        v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err) {}
    void Fillsample(int isample, Double_t _v1, Double_t _v1_err, 
                    Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, 
                    Double_t _v4, Double_t _v4_err) {
        subsample_v1[isample] = _v1;
        subsample_v1_err[isample] = _v1_err;
        subsample_v2[isample] = _v2;
        subsample_v2_err[isample] = _v2_err;
        subsample_v3[isample] = _v3;
        subsample_v3_err[isample] = _v3_err;
        subsample_v4[isample] = _v4;
        subsample_v4_err[isample] = _v4_err;
    }
};
```

### ETADIFF
```cpp
struct VnUnit {
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    std::vector<Double_t> subsample_v2;
    std::vector<Double_t> subsample_v2_err;
    std::vector<Double_t> subsample_v3;
    std::vector<Double_t> subsample_v3_err;
    std::vector<Double_t> subsample_v4;
    std::vector<Double_t> subsample_v4_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, 
           Double_t _v4, Double_t _v4_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), 
        v4(_v4), v4_err(_v4_err) {}
    void Fillsample(int isample, Double_t _v2, Double_t _v2_err, 
                    Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) {
        if (isample >= subsample_v2.size()) {
            Printf("You are filling out of range: isample %d, vector size: %zu", isample, subsample_v2.size());
            return;
        }
        subsample_v2[isample] = _v2;
        subsample_v2_err[isample] = _v2_err;
        subsample_v3[isample] = _v3;
        subsample_v3_err[isample] = _v3_err;
        subsample_v4[isample] = _v4;
        subsample_v4_err[isample] = _v4_err;
    }
};
```

**DIFFERENCE:** v1 members completely removed in EtaDiff; bounds checking added to Fillsample

---

## 2. FourierFit Function - Input File Handling

### REFERENCE
```cpp
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TFile* datafile = new TFile(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", 
             data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, 
             data.maxRange, DihadronCorrTypeName[data.corrType].c_str()), "READ");
    
    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "Cannot open input file: " << Form("./ProcessOutput/BootstrapSample_..." << std::endl;
        return nullptr;
    }
```

### ETADIFF
```cpp
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TFile* datafile = new TFile(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root",
             data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, 
             data.maxRange, DihadronCorrTypeName[data.corrType].c_str()), "READ");
    
    if (!datafile || !datafile->IsOpen()) {
        if (datafile) {
            datafile->Close();
            delete datafile;
            datafile = nullptr;
        }
        // FALLBACK mechanism
        datafile = new TFile(
            Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_0.0_0.1_%s.root",
                 data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, 
                 data.maxRange, DihadronCorrTypeName[data.corrType].c_str()), "READ");
        
        if (!datafile || !datafile->IsOpen()) {
            std::cerr << "Cannot open input file (default and EtaDiff fallback both failed) for " 
                      << data.fileNameSuffix << std::endl;
            exit(1);
        }
        std::cout << "Using EtaDiff fallback bootstrap file for FourierFit: " 
                  << data.fileNameSuffix << std::endl;
    }
```

**DIFFERENCE:** EtaDiff adds fallback to EtaDiff directory if primary path fails; changes error behavior from `return nullptr` to `exit(1)`

---

## 3. RooFourierFit - Initial Parameter Setup (CRITICAL!)

### REFERENCE - HistFFT-Based
```cpp
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr) {
    if (!hist) {
        Error("FourierFit", "Invalid histogram pointer!");
        return;
    }
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(5);
    fParamErr.resize(5);

    // Define Fourier function
    auto fourierFunc = [](double *x, double *p) {
        return p[0] + 2*p[1]*TMath::Cos(x[0]) 
                    + 2*p[2]*TMath::Cos(2*x[0]) 
                    + 2*p[3]*TMath::Cos(3*x[0])
                    + 2*p[4]*TMath::Cos(4*x[0]);
    };
    
    TF1 *fitFunc = new TF1("fourierFit", fourierFunc, 
                           hist->GetXaxis()->GetXmin(), 
                           hist->GetXaxis()->GetXmax(), 5);

    // *** CRITICAL DIFFERENCE: Use FFT ***
    std::vector<double> Coeff;
    HistFFT(hist, Coeff);
    std::cout << "Coeff: " << Coeff[0] << " " << Coeff[1] << " " 
              << Coeff[2] << " " << Coeff[3] << " " << Coeff[4] << " " 
              << Coeff[5] << std::endl;
    
    // Initialize from FFT coefficients
    fitFunc->SetParameter(0, Coeff[0]);
    for (int i = 1; i < 5; ++i) {
        fitFunc->SetParameter(i, Coeff[i]);
    }

    fitFunc->SetParName(0, "a0");
    fitFunc->SetParName(1, "a1");
    fitFunc->SetParName(2, "a2");
    fitFunc->SetParName(3, "a3");
    fitFunc->SetParName(4, "a4");

    TFitResultPtr fitResult = hist->Fit(fitFunc, "QN0");
    // ... rest same ...
```

### ETADIFF - Histogram Statistics-Based
```cpp
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr) {
    if (!hist) {
        Error("FourierFit", "Invalid histogram pointer!");
        return;
    }
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(5);
    fParamErr.resize(5);

    // Define same Fourier function
    auto fourierFunc = [](double *x, double *p) {
        return p[0] + 2*p[1]*TMath::Cos(x[0]) 
                    + 2*p[2]*TMath::Cos(2*x[0]) 
                    + 2*p[3]*TMath::Cos(3*x[0])
                    + 2*p[4]*TMath::Cos(4*x[0]);
    };

    TF1 *fitFunc = new TF1("fourierFit", fourierFunc, 
                           hist->GetXaxis()->GetXmin(), 
                           hist->GetXaxis()->GetXmax(), 5);

    // *** CRITICAL DIFFERENCE: Use histogram statistics ***
    // Use histogram statistics rather than potentially huge FFT coefficients
    double histMean = hist->GetMean();
    double histRMS = hist->GetRMS();
    
    fitFunc->SetParameter(0, histMean);          // a0 = Mean
    for (int i = 1; i < 5; ++i) {
        fitFunc->SetParameter(i, 0.1 * histRMS); // a{i} = 10% RMS
    }

    fitFunc->SetParName(0, "a0");
    fitFunc->SetParName(1, "a1");
    fitFunc->SetParName(2, "a2");
    fitFunc->SetParName(3, "a3");
    fitFunc->SetParName(4, "a4");

    TFitResultPtr fitResult = hist->Fit(fitFunc, "QN0");
    // ... rest same ...
```

**DIFFERENCE:** 
- Reference uses `HistFFT(hist, Coeff)` to compute frequency-domain coefficients
- EtaDiff uses `hist→GetMean()` and `hist→GetRMS()` for simpler, more conservative initial guesses
- EtaDiff comment explicitly states reason: "Use histogram statistics rather than potentially huge FFT coefficients"

---

## 4. fitSample Return Statement

### REFERENCE
```cpp
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0) {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!hm) {
        std::cerr << "[DEBUG] Cannot find histogram in file: " << data.fileNameSuffix << ...
        return 0;
    }
    RooFourierFit(hm, fParamVal, fParamErr);
    if (sample == -1) {
        PlotFitting(hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, data.corrType, pTMin, pTMax);
    }
    
    // fParamVal[0]=v2², [1]=v3², [2]=v4², [3]=v1², [4]=a0
    VnUnit* vnResult = new VnUnit(
        fParamVal[3], fParamErr[3],  // v1
        fParamVal[0], fParamErr[0],  // v2
        fParamVal[1], fParamErr[1],  // v3
        fParamVal[2], fParamErr[2]   // v4
    );
    return vnResult;
}
```

### ETADIFF
```cpp
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0) {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!hm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", ...) << std::endl;
        return 0;
    }
    RooFourierFit(hm, fParamVal, fParamErr);
    if (sample == -1) {
        PlotFitting(hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, data.corrType, pTMin, pTMax);
    }
    
    // fParamVal[0]=v2², [1]=v3², [2]=v4², [3]=v1², [4]=a0
    // BUT only use first 3 harmonics
    VnUnit* vnResult = new VnUnit(
        fParamVal[0], fParamErr[0],  // v2
        fParamVal[1], fParamErr[1],  // v3
        fParamVal[2], fParamErr[2]   // v4
    );
    return vnResult;
}
```

**DIFFERENCE:** Return order changes from (v1, v2, v3, v4) to (v2, v3, v4) only; histogram name error message differs

---

## 5. Bootstrap Error Calculation

### REFERENCE
```cpp
std::vector<std::vector<std::vector<double>>> ValueArray;
std::vector<std::vector<std::vector<double>>> ValueErrorArray;
std::vector<std::vector<double>> ErrorArray;
int Nobs = 3;  // v22, v32, v42
int NofSample = maxSample * maxSample;
int Nbin = 1;
ResizeValueArray(ValueArray, ValueErrorArray, ErrorArray, Nobs, NofSample, Nbin);
vnResult->ResizeSubsample(NofSample);
if (pTMin > 0 && pTMax > 0) vnResult_PtDiff->ResizeSubsample(NofSample);

for(int sample = 0; sample < NofSample; sample++) {
    VnUnit* vnTemp = fitSample(isNch, datafile, data, sample);
    if (!vnTemp) {
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << " sample: " << sample << std::endl;
        return nullptr;
    }
    ValueArray[0][sample][0] = vnTemp->v2;
    ValueErrorArray[0][sample][0] = vnTemp->v2_err;
    ValueArray[1][sample][0] = vnTemp->v3;
    ValueErrorArray[1][sample][0] = vnTemp->v3_err;
    ValueArray[2][sample][0] = vnTemp->v4;
    ValueErrorArray[2][sample][0] = vnTemp->v4_err;
    
    vnResult->Fillsample(sample, vnTemp->v1, vnTemp->v1_err, 
                         vnTemp->v2, vnTemp->v2_err, 
                         vnTemp->v3, vnTemp->v3_err, 
                         vnTemp->v4, vnTemp->v4_err);
    
    if (pTMin > 0 && pTMax > 0) {
        VnUnit* vnTemp_PtDiff = fitSample(isNch, datafile_PtDiff, data, sample, pTMin, pTMax);
        // ... similar processing ...
        vnResult_PtDiff->Fillsample(sample, vnTemp_PtDiff->v1, vnTemp_PtDiff->v1_err, ...);
    }
    delete vnTemp;
}

for(int iobs = 0; iobs < Nobs; iobs++) {
    CalculateBootstrapError(ValueArray[iobs], ValueErrorArray[iobs], ErrorArray[iobs], 1.);
}

vnResult->v2_err = ErrorArray[0][0];
vnResult->v3_err = ErrorArray[1][0];
vnResult->v4_err = ErrorArray[2][0];
if (pTMin > 0 && pTMax > 0) {
    vnResult_PtDiff->v2_err = ErrorArray[0][0];
    vnResult_PtDiff->v3_err = ErrorArray[1][0];
    vnResult_PtDiff->v4_err = ErrorArray[2][0];
}
```

### ETADIFF
```cpp
std::vector<std::vector<std::vector<double>>> ValueArray;
std::vector<std::vector<std::vector<double>>> ValueErrorArray;
std::vector<std::vector<double>> ErrorArray;
int Nobs = 3;  // v22, v32, v42
int NofSample = maxSample * maxSample;
int Nbin = 1;
ResizeValueArray(ValueArray, ValueErrorArray, ErrorArray, Nobs, NofSample, Nbin);
vnResult->ResizeSubsample(NofSample);
if (pTMin > 0 && pTMax > 0) vnResult_PtDiff->ResizeSubsample(NofSample);

for(int sample = 0; sample < NofSample; sample++) {
    VnUnit* vnTemp = fitSample(isNch, datafile, data, sample);
    if (!vnTemp) {
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << " sample: " << sample << std::endl;
        exit(1);  // Different: exit() instead of return nullptr
    }
    ValueArray[0][sample][0] = vnTemp->v2;
    ValueErrorArray[0][sample][0] = vnTemp->v2_err;
    ValueArray[1][sample][0] = vnTemp->v3;
    ValueErrorArray[1][sample][0] = vnTemp->v3_err;
    ValueArray[2][sample][0] = vnTemp->v4;
    ValueErrorArray[2][sample][0] = vnTemp->v4_err;
    
    vnResult->Fillsample(sample, vnTemp->v2, vnTemp->v2_err, 
                         vnTemp->v3, vnTemp->v3_err, 
                         vnTemp->v4, vnTemp->v4_err);  // NO v1!
    
    if (pTMin > 0 && pTMax > 0) {
        VnUnit* vnTemp_PtDiff = fitSample(isNch, datafile_PtDiff, data, sample, pTMin, pTMax);
        // ... similar processing ...
        vnResult_PtDiff->Fillsample(sample, vnTemp_PtDiff->v2, vnTemp_PtDiff->v2_err, ...);  // NO v1!
    }
    delete vnTemp;
}

for(int iobs = 0; iobs < Nobs; iobs++) {
    CalculateBootstrapError(ValueArray[iobs], ValueErrorArray[iobs], ErrorArray[iobs], 1.);
}

vnResult->v2_err = ErrorArray[0][0];
vnResult->v3_err = ErrorArray[1][0];
vnResult->v4_err = ErrorArray[2][0];
if (pTMin > 0 && pTMax > 0) {
    vnResult_PtDiff->v2_err = ErrorArray[0][0];
    vnResult_PtDiff->v3_err = ErrorArray[1][0];
    vnResult_PtDiff->v4_err = ErrorArray[2][0];
}
```

**DIFFERENCE:** Identical bootstrap calculation structure (Nobs=3), but Fillsample calls omit v1; error handling changes from return to exit

---

## 6. Output Histogram Creation

### REFERENCE
```cpp
if (isContinuous) {
    // Create variable bin width histograms
    TH1D* hV1 = new TH1D("hV1", "v_{1};Centrality;v_{1}", ...);  // v1 INCLUDED
    TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", ...);
    TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", ...);
    TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", ...);

    for (size_t i = 0; i < vnResults.size(); ++i) {
        hV1->SetBinContent(i+1, vnResults[i]->v1);
        hV1->SetBinError(i+1, vnResults[i]->v1_err);
        hV2->SetBinContent(i+1, vnResults[i]->v2);
        hV2->SetBinError(i+1, vnResults[i]->v2_err);
        hV3->SetBinContent(i+1, vnResults[i]->v3);
        hV3->SetBinError(i+1, vnResults[i]->v3_err);
        hV4->SetBinContent(i+1, vnResults[i]->v4);
        hV4->SetBinError(i+1, vnResults[i]->v4_err);
    }

    hV1->Write();
    hV2->Write();
    hV3->Write();
    hV4->Write();
```

### ETADIFF
```cpp
if (isContinuous) {
    // Create variable bin width histograms
    // NO hV1 created
    TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", ...);
    TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", ...);
    TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", ...);

    for (size_t i = 0; i < vnResults.size(); ++i) {
        // NO hV1 operations
        hV2->SetBinContent(i+1, vnResults[i]->v2);
        hV2->SetBinError(i+1, vnResults[i]->v2_err);
        hV3->SetBinContent(i+1, vnResults[i]->v3);
        hV3->SetBinError(i+1, vnResults[i]->v3_err);
        hV4->SetBinContent(i+1, vnResults[i]->v4);
        hV4->SetBinError(i+1, vnResults[i]->v4_err);
    }

    // Only write 3
    hV2->Write();
    hV3->Write();
    hV4->Write();
```

**DIFFERENCE:** Reference creates and writes 4 histograms (v1, v2, v3, v4); EtaDiff creates only 3 (v2, v3, v4)

---

## 7. PlotFitting Legend - Display Differences

### REFERENCE
```cpp
leg->AddEntry(hm, "Data", "lep");
leg->AddEntry(gCopy, "a_{0} + #Sigma_{n=1}^{4}2a_{n}cos(n#Delta#phi)", "l");
leg->AddEntry(gPeri, "Baseline a0", "l");
leg->AddEntry(fit_p1, Form("a_{0} + 2a_{1}cos(#Delta#phi), v1^{2} = %0.6f", v12), "l");
leg->AddEntry(fit_p2, Form("a_{0} + 2a_{2}cos(2#Delta#phi), v2^{2} = %0.6f", v21), "l");
leg->AddEntry(fit_p3, Form("a_{0} + 2a_{3}cos(3#Delta#phi), v3^{2} = %0.6f", v31), "l");
leg->AddEntry(fit_p4, Form("a_{0} + 2a_{4}cos(4#Delta#phi), v4^{2} = %0.6f", v41), "l");
```

### ETADIFF
```cpp
leg->AddEntry(hm, "Data", "lep");
leg->AddEntry(gCopy, "a_{0} + #Sigma_{n=1}^{4}2a_{n}cos(n#Delta#phi)", "l");
leg->AddEntry(gPeri, "Baseline a0", "l");
leg->AddEntry(fit_p1, Form("a_{0} + 2a_{1}cos(#Delta#phi), v1^{2}#times10^{3} = %0.2f", v12*1e3), "l");
leg->AddEntry(fit_p2, Form("a_{0} + 2a_{2}cos(2#Delta#phi), v2^{2}#times10^{3} = %0.2f", v21*1e3), "l");
leg->AddEntry(fit_p3, Form("a_{0} + 2a_{3}cos(3#Delta#phi), v3^{2}#times10^{3} = %0.2f", v31*1e3), "l");
leg->AddEntry(fit_p4, Form("a_{0} + 2a_{4}cos(4#Delta#phi), v4^{2}#times10^{3} = %0.2f", v41*1e3), "l");
```

**DIFFERENCE:** 
- EtaDiff scales all values by ×10³
- EtaDiff shows 2 decimal places (%0.2f) vs reference 6 (%0.6f)

---

## 8. Main Function - Dataset Configuration

### REFERENCE - 17 Datasets
```cpp
void Process_FourierFit() {
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;
    collisionSystemName = "Ne-Ne";
    kOutputVnDelta = true;
    
    // O-O MC centrality-based (3)
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, 
        {InputUnit("LHC25ae_pass2_638219", kFT0AFT0C, 0, 20)}, ...));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, 
        {InputUnit("LHC25ae_pass2_638219", kTPCFT0A, 0, 20)}, ...));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, 
        {InputUnit("LHC25ae_pass2_638219", kTPCFT0C, 0, 20)}, ...));
    
    // O-O PtDiff (2)
    configList.push_back(ConfigUnit(kCent, kPtDiffOn, ...));
    configList.push_back(ConfigUnit(kCent, kPtDiffOn, ...));
    
    // O-O Nch-based (5)
    configList.push_back(ConfigUnit(kNch, kPtDiffOff, ...));
    // ... (4 more Nch configs with PtDiff)
    
    // p-O Nch-based (5)
    configList.push_back(ConfigUnit(kNch, kPtDiffOff, ...));
    // ... (4 more)
    
    // p-p Nch-based (5)
    configList.push_back(ConfigUnit(kNch, kPtDiffOff, ...));
    // ... (4 more)
    
    // Ne-Ne Nch-based (5)
    configList.push_back(ConfigUnit(kNch, kPtDiffOff, ...));
    // ... (4 more)
    
    // MC Nch-based (5)
    configList.push_back(ConfigUnit(kNch, kPtDiffOff, ...));
    // ... (4 more)
    
    for (auto config : configList) {
        if (config.isPtDiff) {
            ProcessConfig_PtDiff(config.isNch, config.dataList, config.outputFileName);
        } else {
            ProcessConfig(config.isNch, config.dataList, config.outputFileName);
        }
    }
}
```

### ETADIFF - 3 Datasets Only
```cpp
void Process_FourierFit() {
    gROOT->SetBatch(kTRUE);
    
    std::vector<ConfigUnit> configList;
    collisionSystemName = "Unknown";
    kOutputVnDelta = true;
    
    // Ne-Ne inner ring datasets (Cent 0-20) - 3 ONLY
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
        {InputUnit("LHC25af_pass2_632504", kTPCFT0A, 0, 20)},
        "LHC25af_pass2_632504"));
    
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
        {InputUnit("LHC25af_pass2_637596", kTPCFT0C, 0, 20)},
        "LHC25af_pass2_637596"));
    
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
        {InputUnit("LHC25af_pass2_640018", kFT0AFT0C, 0, 20)},
        "LHC25af_pass2_640018"));
    
    for (auto config : configList) {
        if (!config.dataList.empty()) {
            collisionSystemName = GetCollisionSystemNameFromDataset(config.dataList[0].fileNameSuffix);
        }
        ProcessConfig(config.isNch, config.dataList, config.outputFileName);  // ONLY ProcessConfig
    }
}
```

**DIFFERENCE:** Reference operates on ~17 datasets with both ProcessConfig and ProcessConfig_PtDiff; EtaDiff operates on only 3 datasets with ProcessConfig only

---

## Summary Table

| Aspect | Reference | EtaDiff | Impact |
|--------|-----------|---------|--------|
| v1 members | Present | Removed | Smaller output files; physics focus shifts |
| fParamVal size | 5 | 5 | Same fit coefficients calculated |
| Initial parameters | HistFFT() | histogram Mean/RMS | Simpler but potentially less precise |
| File fallback | None | EtaDiff path | More robust to missing centrality data |
| Error behavior | return nullptr | exit(1) | Crashes on failure vs. graceful skip |
| Output histograms | 4 (v1,v2,v3,v4) | 3 (v2,v3,v4) | 25% smaller output |
| Datasets | 17 total | 3 total | Specialized vs. general purpose |
| PtDiff support | Yes (5+configs) | No | Different analysis scope |
| Scale in plots | Raw values | ×10³ | Different numerical display |

