# Detailed Comparison: Process_FourierFit.cxx

## Overview
This document compares the Reference version (`/home/martijn-laarhoven/Work/dihadronanalysis-master/LongRangeDihadron/Process_FourierFit.cxx`) with the EtaDiff version (`/home/martijn-laarhoven/Work/dihadronanalysis-master/LongRangeEta/dihadronanalysis/LongRangeDihadron/Process_FourierFit.cxx`).

---

## 1. VnUnit Structure - MAJOR DIFFERENCE

### Reference Version:
```cpp
struct VnUnit {
    Double_t v1;        // v1 harmonic
    Double_t v1_err;
    Double_t v2;        // v2 harmonic – OUTPUT PARAMETER
    Double_t v2_err;
    Double_t v3;        // v3 harmonic – OUTPUT PARAMETER
    Double_t v3_err;
    Double_t v4;        // v4 harmonic – OUTPUT PARAMETER
    Double_t v4_err;
    std::vector<Double_t> subsample_v1;     // Bootstrap subsamples include v1
    std::vector<Double_t> subsample_v1_err;
    std::vector<Double_t> subsample_v2;
    std::vector<Double_t> subsample_v2_err;
    std::vector<Double_t> subsample_v3;
    std::vector<Double_t> subsample_v3_err;
    std::vector<Double_t> subsample_v4;
    std::vector<Double_t> subsample_v4_err;
    
    // Constructor includes v1 parameters
    VnUnit(Double_t _v1, Double_t _v1_err, Double_t _v2, Double_t _v2_err, 
           Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err)
    
    // Fillsample includes v1
    void Fillsample(int isample, Double_t _v1, Double_t _v1_err, Double_t _v2, ...)
};
```

### EtaDiff Version:
```cpp
struct VnUnit {
    Double_t v2;        // v2 harmonic only
    Double_t v2_err;
    Double_t v3;        // v3 harmonic only
    Double_t v3_err;
    Double_t v4;        // v4 harmonic only
    Double_t v4_err;
    std::vector<Double_t> subsample_v2;     // NO v1 in subsamples
    std::vector<Double_t> subsample_v2_err;
    std::vector<Double_t> subsample_v3;
    std::vector<Double_t> subsample_v3_err;
    std::vector<Double_t> subsample_v4;
    std::vector<Double_t> subsample_v4_err;
    
    // Constructor omits v1 parameters
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, 
           Double_t _v4, Double_t _v4_err)
    
    // Fillsample omits v1
    void Fillsample(int isample, Double_t _v2, Double_t _v2_err, Double_t _v3, ...)
};
```

**IMPACT:** The EtaDiff version strips out the v1 (directed flow) parameter completely, focusing only on v2, v3, v4.

---

## 2. ProcessConfig Function - Input File Handling

### Reference Version:
```cpp
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Processes diverse datasets (O-O, p-O, p-p, Ne-Ne with multiple collision systems)
    // Handles failed fits by skipping datasets:
    
    for (const auto& data : dataList) {
        std::cout << "[DEBUG] Processing dataset: " << data.fileNameSuffix << ...
        VnUnit* result = FourierFit(isNch, data, (!kOutputVnDelta));
        if (!result) {
            std::cerr << "[WARN] Skipping dataset due to missing input or failed fit: ...
            continue;  // Skip failed entries
        }
        vnResults.push_back(result);
        validDataList.push_back(data);
    }
}
```

**Features:**
- Includes error handling for missing files or fit failures
- Maintains separate `validDataList` to track which datasets were successfully processed
- Creates output using only valid results

### EtaDiff Version:
```cpp
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Simpler processing for Ne-Ne (Cent 0-20) datasets only
    // No explicit error handling:
    
    for (const auto& data : dataList) {
        vnResults.push_back(FourierFit(isNch, data, (!kOutputVnDelta)));
    }
    
    // Uses original dataList directly for output
}
```

**Features:**
- No validation of fit results (assumes all succeed)
- Uses original `dataList` directly rather than creating a valid subset
- More vulnerable to crashes from missing input files

---

## 3. FourierFit Function - Input File Path Strategy

### Reference Version (Pure Centrality):
```cpp
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, ...) {
    TFile* datafile = new TFile(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root",
             data.fileNameSuffix.c_str(), splitName.c_str(), 
             data.minRange, data.maxRange, 
             DihadronCorrTypeName[data.corrType].c_str()), "READ");
    
    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "Cannot open input file: ...
        return nullptr;  // No fallback
    }
}
```

### EtaDiff Version (Fallback to EtaDiff):
```cpp
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, ...) {
    // FIRST try: default path
    TFile* datafile = new TFile(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", ...), "READ");
    
    if (!datafile || !datafile->IsOpen()) {
        if (datafile) {
            datafile->Close();
            delete datafile;
            datafile = nullptr;
        }
        
        // FALLBACK: try EtaDiff version
        datafile = new TFile(
            Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_0.0_0.1_%s.root",
                 data.fileNameSuffix.c_str(), splitName.c_str(), 
                 data.minRange, data.maxRange, 
                 DihadronCorrTypeName[data.corrType].c_str()), "READ");
        
        if (!datafile || !datafile->IsOpen()) {
            std::cerr << "Cannot open input file (default and EtaDiff fallback both failed) ...
            exit(1);
        }
        std::cout << "Using EtaDiff fallback bootstrap file for FourierFit: ...
    }
}
```

**CRITICAL DIFFERENCE:**
- Reference: Strict file path matching - fails if file not found
- EtaDiff: Implements fallback mechanism to use `./ProcessOutput/EtaDiff/BootstrapSample_*_Eta_0.0_0.1_*.root` if primary path fails
- EtaDiff also changes error handling from `return nullptr` to `exit(1)` for the fallback case

---

## 4. RooFourierFit Function - Initial Parameter Setup (MAJOR DIFFERENCE)

### Reference Version - Uses FFT:
```cpp
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, ...) {
    fParamVal.resize(5);
    fParamErr.resize(5);
    
    // Define Fourier function: a0 + 2*a1*cos(x) + 2*a2*cos(2x) + ... + 2*a4*cos(4x)
    
    // CRITICAL: Compute FFT coefficients
    std::vector<double> Coeff;
    HistFFT(hist, Coeff);
    std::cout << "Coeff: " << Coeff[0] << " " << Coeff[1] << " " << Coeff[2] 
              << " " << Coeff[3] << " " << Coeff[4] << " " << Coeff[5] << std::endl;
    
    // Initialize parameters using FFT coefficients
    fitFunc->SetParameter(0, Coeff[0]);      // a0 from FFT
    for (int i = 1; i < 5; ++i) {
        fitFunc->SetParameter(i, Coeff[i]);  // a{i} from FFT
    }
}
```

**Notes on HistFFT:**
- Function defined elsewhere (likely in BasicForDihadron.h or included header)
- Computes FFT coefficients directly from histogram
- Provides sophisticated initial parameter guesses based on frequency content
- Prints all 6 FFT coefficients for debugging

### EtaDiff Version - Uses Histogram Statistics:
```cpp
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, ...) {
    // NO HistFFT call - completely removed
    
    fParamVal.resize(5);
    fParamErr.resize(5);
    
    // Define same Fourier function
    
    // Simpler initialization using histogram statistics
    double histMean = hist->GetMean();
    double histRMS = hist->GetRMS();
    
    fitFunc->SetParameter(0, histMean);           // a0 = Mean of histogram
    for (int i = 1; i < 5; ++i) {
        fitFunc->SetParameter(i, 0.1 * histRMS); // a{i} = 10% of RMS
    }
}
```

**RATIONALE FOR CHANGE:**
- HistFFT potentially produces huge coefficients that may cause numerical issues
- Simpler histogram statistics provide conservative starting values
- Comment explicitly states the reason for the change: "Use histogram statistics rather than potentially huge FFT coefficients"

---

## 5. Nobs Value - CRITICAL DIFFERENCE

### Reference Version:
```cpp
int Nobs = 3;  // v22, v32, v42
int NofSample = maxSample * maxSample;
```

**Bootstrap Error Calculation loop:**
```cpp
ValueArray[0][sample][0] = vnTemp->v2;        // v2² (index 0)
ValueErrorArray[0][sample][0] = vnTemp->v2_err;
ValueArray[1][sample][0] = vnTemp->v3;        // v3² (index 1)
ValueErrorArray[1][sample][0] = vnTemp->v3_err;
ValueArray[2][sample][0] = vnTemp->v4;        // v4² (index 2)
ValueErrorArray[2][sample][0] = vnTemp->v4_err;

vnResult->Fillsample(sample, vnTemp->v1, vnTemp->v1_err, vnTemp->v2, ...);
```

### EtaDiff Version:
```cpp
int Nobs = 3;  // v22, v32, v42
int NofSample = maxSample * maxSample;
```

Same structure, but **Fillsample** call is different:
```cpp
ValueArray[0][sample][0] = vnTemp->v2;        // v2² only
ValueErrorArray[0][sample][0] = vnTemp->v2_err;
ValueArray[1][sample][0] = vnTemp->v3;        // v3² only
ValueErrorArray[1][sample][0] = vnTemp->v3_err;
ValueArray[2][sample][0] = vnTemp->v4;        // v4² only
ValueErrorArray[2][sample][0] = vnTemp->v4_err;

// NO v1 in the call
vnResult->Fillsample(sample, vnTemp->v2, vnTemp->v2_err, vnTemp->v3, vnTemp->v3_err, vnTemp->v4, vnTemp->v4_err);
```

**Note:** Both versions have the same `Nobs = 3` count, suggesting bootstrap errors are calculated identically (for v2², v3², v4² only).

---

## 6. Output File Naming and Structure - Output Handling

### Reference Version - Creates More Histograms:
```cpp
TH1D* hV1 = new TH1D("hV1", "v_{1};Centrality;v_{1}", ...);  // v1 INCLUDED
TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", ...);
TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", ...);
TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", ...);

// Fill all 4 harmonics
for (size_t i = 0; i < vnResults.size(); ++i) {
    hV1->SetBinContent(i+1, vnResults[i]->v1);
    hV1->SetBinError(i+1, vnResults[i]->v1_err);
    hV2->SetBinContent(i+1, vnResults[i]->v2);
    ...
}

// Write all 4
hV1->Write();
hV2->Write();
hV3->Write();
hV4->Write();

// Subsamples include v1
for (Int_t sample = 0; sample < vnResults[0]->subsample_v2.size(); sample++) {
    TH1D* hV1_sub = new TH1D(Form("hV1_subsample_%d", sample), ...);
    TH1D* hV2_sub = new TH1D(Form("hV2_subsample_%d", sample), ...);
    // Fill and write all 4 subsamples
}
```

### EtaDiff Version - Only v2, v3, v4:
```cpp
// NO hV1 created
TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", ...);
TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", ...);
TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", ...);

// Fill only 3 harmonics
for (size_t i = 0; i < vnResults.size(); ++i) {
    hV2->SetBinContent(i+1, vnResults[i]->v2);
    hV2->SetBinError(i+1, vnResults[i]->v2_err);
    hV3->SetBinContent(i+1, vnResults[i]->v3);
    ...
}

// Write only 3
hV2->Write();
hV3->Write();
hV4->Write();

// Subsamples only have v2, v3, v4
for (Int_t sample = 0; sample < vnResults[0]->subsample_v2.size(); sample++) {
    // NO hV1_sub created
    TH1D* hV2_sub = new TH1D(Form("hV2_subsample_%d", sample), ...);
    TH1D* hV3_sub = new TH1D(Form("hV3_subsample_%d", sample), ...);
    TH1D* hV4_sub = new TH1D(Form("hV4_subsample_%d", sample), ...);
    // Fill and write only 3 subsamples
}
```

---

## 7. fitSample Function - Histogram Name and Parameter Return Order

### Reference Version:
```cpp
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, ...) {
    TH1D* hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", ...));
    
    RooFourierFit(hm, fParamVal, fParamErr);
    
    // Return order: v1 (index 3), v1_err, v2 (index 0), v2_err, v3, v3_err, v4, v4_err
    VnUnit* vnResult = new VnUnit(
        fParamVal[3], fParamErr[3],  // v1 = a1/a0
        fParamVal[0], fParamErr[0],  // v2 = a2/a0
        fParamVal[1], fParamErr[1],  // v3 = a3/a0
        fParamVal[2], fParamErr[2]   // v4 = a4/a0
    );
}
```

### EtaDiff Version:
```cpp
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, ...) {
    TH1D* hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", ...));
    
    RooFourierFit(hm, fParamVal, fParamErr);
    
    // Return order: v2 (index 0), v2_err, v3 (index 1), v3_err, v4, v4_err
    // NO v1!
    VnUnit* vnResult = new VnUnit(
        fParamVal[0], fParamErr[0],  // v2 = a2/a0
        fParamVal[1], fParamErr[1],  // v3 = a3/a0
        fParamVal[2], fParamErr[2]   // v4 = a4/a0
    );
}
```

**Key Point:** The histogram names `bsSample_hPhiSameOverMixed_*` are identical; only the parameter order changes due to removed v1.

---

## 8. RooFourierFit Output Parameters - Array Structure

### Both Versions - Return Parameter Order:
```
fParamVal[0] = v22square  (v2²)
fParamVal[1] = v32square  (v3²)
fParamVal[2] = v42square  (v4²)
fParamVal[3] = v12square  (v1²)  -- Reference only
fParamVal[4] = a0         (baseline)
```

**Reference interpretation in fitSample:**
- Constructs VnUnit(fParamVal[3], ..., fParamVal[0], ...) = (v1, v2, v3, v4)

**EtaDiff interpretation in fitSample:**
- Constructs VnUnit(fParamVal[0], ..., fParamVal[1], ...) = (v2, v3, v4)
- Array elements remain but only subset is used

---

## 9. PlotFitting Function - Output Differences

### Reference Version:
```cpp
void PlotFitting(...) {
    // Display collision system dynamically
    std::string displayCollisionName = collisionSystemName;
    if (fileSuffix.find("LHC25af") == 0) displayCollisionName = "Ne-Ne";
    else if (fileSuffix.find("LHC25ae") == 0) displayCollisionName = "O-O";
    else if (fileSuffix.find("LHC25ac") == 0) displayCollisionName = "p-p";
    else if (fileSuffix.find("LHC25ad") == 0) displayCollisionName = "p-O";
    
    // Plot legend includes v1 and shows v1², v2², v3², v4² values (not 1e3 scaled)
    leg->AddEntry(fit_p1, Form("a_{0} + 2a_{1}cos(#Delta#phi), v1^{2} = %0.6f", v12), "l");
    // ... all 4 components plotted
}
```

### EtaDiff Version:
```cpp
void PlotFitting(...) {
    // Simpler collision name (uses global collisionSystemName)
    // NO dynamic collision system detection
    
    // Plot legend does NOT include v1, scales v1² and v2² by 1e3
    // Note: Code has "v1^{2}#times10^{3}" but no actual v1 is plotted
    leg->AddEntry(fit_p1, Form("... v1^{2}#times10^{3} = %0.2f", v12*1e3), "l");
    leg->AddEntry(fit_p2, Form("... v2^{2}#times10^{3} = %0.2f", v21*1e3), "l");
}
```

---

## 10. Dataset Configuration - Process_FourierFit() Main Function

### Reference Version:
Processes multiple collision systems and configurations:
```cpp
// 11 config entries total
// O-O datasets (3 configs)
configList.push_back(ConfigUnit(kCent, kPtDiffOff, 
    {InputUnit("LHC25ae_pass2_638219", kFT0AFT0C, 0, 20)}, ...));
configList.push_back(ConfigUnit(kCent, kPtDiffOff, 
    {InputUnit("LHC25ae_pass2_638219", kTPCFT0A, 0, 20)}, ...));
    
// O-O with PtDiff (2 configs)
configList.push_back(ConfigUnit(kCent, kPtDiffOn, ...));

// Nch-based for multiple systems (p-O, p-p, Ne-Ne, MC)
// O-O Nch (5 configs with PtDiff)
// p-O Nch (5 configs)
// p-p Nch (5 configs)
// Ne-Ne Nch (5 configs)
// MC Nch (5 configs)

// Call conditional processing
for (auto config : configList) {
    if (config.isPtDiff) {
        ProcessConfig_PtDiff(config.isNch, config.dataList, config.outputFileName);
    } else {
        ProcessConfig(config.isNch, config.dataList, config.outputFileName);
    }
}
```

### EtaDiff Version:
Simplified - only Ne-Ne inner ring datasets (Cent 0-20):
```cpp
// 3 config entries total
configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_632504", kTPCFT0A, 0, 20)},
    "LHC25af_pass2_632504"));
configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_637596", kTPCFT0C, 0, 20)},
    "LHC25af_pass2_637596"));
configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_640018", kFT0AFT0C, 0, 20)},
    "LHC25af_pass2_640018"));

// ONLY calls ProcessConfig, NOT ProcessConfig_PtDiff
for (auto config : configList) {
    if (!config.dataList.empty()) {
        collisionSystemName = GetCollisionSystemNameFromDataset(config.dataList[0].fileNameSuffix);
    }
    ProcessConfig(config.isNch, config.dataList, config.outputFileName);
}
```

**MAJOR CHANGE:**
- Reference: General-purpose script for multiple collision systems
- EtaDiff: Specialized only for Ne-Ne 0-20% centrality (3 detector configurations)
- EtaDiff: No PtDiff processing, simplified dataset handling

---

## Summary of Key Differences

| Aspect | Reference | EtaDiff |
|--------|-----------|---------|
| **v1 Harmonic** | Included throughout | REMOVED completely |
| **VnUnit Members** | v1, v2, v3, v4 (8 members) | v2, v3, v4 only (6 members) |
| **RooFourierFit Strategy** | Uses HistFFT for initial parameters | Uses histogram Mean/RMS statistics |
| **File Fallback** | No fallback, strict paths | EtaDiff fallback to `Eta_0.0_0.1` |
| **Error Handling** | Skips failed datasets | Exits on failure (fallback case) |
| **Output Histograms** | 4 harmonics (v1, v2, v3, v4) | 3 harmonics (v2, v3, v4) |
| **Dataset Scope** | Multiple collision systems | Ne-Ne only (3 detector configs) |
| **PtDiff Processing** | Implemented | Disabled/removed |
| **Plot Display** | Dynamic collision system detection | Static collision name |
| **Legend Scaling** | Raw v2 values (×10⁶ implicit) | Values scaled by ×10³ |

---

## Implications

1. **Physics Analysis:** EtaDiff focuses on eta-differential analysis, hence removal of v1 (which probes forward-backward asymmetry less relevant in eta-diff context)
2. **Initial Fit Parameters:** Switching from FFT to histogram statistics reduces potential numerical issues
3. **Robustness:** Fallback mechanism in EtaDiff provides more tolerant handling of missing intermediate output files
4. **Specialization:** EtaDiff version is optimized for specific Ne-Ne configurations rather than general use
5. **Bootstrap Error:** Both versions use identical bootstrap methodology (Nobs=3 for v2², v3², v4²)

