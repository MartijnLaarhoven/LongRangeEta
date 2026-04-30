# O-O vs Ne-Ne System Processing Analysis

## Executive Summary

The analysis shows that **O-O and Ne-Ne systems use identical mathematical formulas and processing code paths**, but are differentiated solely by **different datasets** and **physical properties** of the collision systems themselves. No explicit system-specific scaling, normalization factors, or weighting has been implemented in the code.

---

## 1. CORE CALCULATION: 3times2PC Formula

### Location: `Process_3times2PC.cxx` lines 813-825

```cpp
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
    double v2_err = Get3times2PC_Error(LMvalues->v2, LMvalues->v2_err, 
                                       MRvalues->v2, MRvalues->v2_err, 
                                       LRvalues->v2, LRvalues->v2_err);
    // ... same for v3 and v4 ...
    VnUnit* vnResult = new VnUnit(v2, v2_err, v3, v3_err, v4, v4_err);
    return vnResult;
}
```

**Formula:** 
$$v_n = \sqrt{\frac{v_{n,LM} \times v_{n,MR}}{v_{n,LR}}}$$

**Error Propagation:**
$$\sigma_{v_n} = \sqrt{\left(\frac{\partial v_n}{\partial LM} \sigma_{LM}\right)^2 + \left(\frac{\partial v_n}{\partial MR} \sigma_{MR}\right)^2 + \left(\frac{\partial v_n}{\partial LR} \sigma_{LR}\right)^2}$$

Where:
- **LM** = TPC × FT0A correlation (Left-Middle)
- **MR** = TPC × FT0C correlation (Middle-Right)  
- **LR** = FT0A × FT0C correlation (Long-Range)

---

## 2. v2 COEFFICIENT EXTRACTION

### TemplateFit Initialization: `Process_TemplateFit.cxx` lines 1224-1232

```cpp
if (!kRefit){
    ft->AddParameter("Fa","Fa",0.2,0,1000);        // Template amplitude factor
    ft->AddParameter("Ga","Ga",0.5,0,30000);       // Background constant
    ft->AddParameter("v2","v2",4e-3,-1.0,1.0);     // v2 initial: 0.004
    ft->AddParameter("v3","v3",6e-4,-1.0,1.0);     // v3 initial: 0.0006
    ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);   // v4 initial: 0.00018
}
```

**IDENTICAL for both O-O and Ne-Ne** - same initial parameter values and bounds.

### TemplateFit Core Function: `Process_TemplateFit.cxx` lines 1181-1234

```cpp
// Estimate Fa and Ga via a weighted linear regression: hm ~= Fa * lm + Ga
// ... (fitting implementation identical for all systems)
```

---

## 3. DATASET CONFIGURATIONS

### O-O System Configuration

**Location:** `Process_TemplateFit.cxx` lines 116-126

```cpp
if (systemName == "O-O") {
    return {
        NamedConfigUnit("FULL_LM_TPC_FT0A", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25ae_pass2_648788", kTPCFT0A, 80, 100), 
            {InputUnit("LHC25ae_pass2_648788", kTPCFT0A, 0, 20)}, 
            "LHC25ae_pass2_648788")),
        NamedConfigUnit("FULL_MR_TPC_FT0C", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25ae_pass2_648788", kTPCFT0C, 80, 100), 
            {InputUnit("LHC25ae_pass2_648788", kTPCFT0C, 0, 20)}, 
            "LHC25ae_pass2_648788")),
        NamedConfigUnit("FULL_LR_FT0A_FT0C", ConfigUnit(kCent, false, false, 
            InputUnit("LHC25ae_pass2_648798", kFT0AFT0C, 80, 100), 
            {InputUnit("LHC25ae_pass2_648798", kFT0AFT0C, 0, 20)}, 
            "LHC25ae_pass2_648798")),
        // Ring-specific configurations...
        NamedConfigUnit("INNER_LM_TPC_FT0A", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25ae_pass2_638221", kTPCFT0A, 80, 100), 
            {InputUnit("LHC25ae_pass2_638221", kTPCFT0A, 0, 20)}, 
            "LHC25ae_pass2_638221")),
        NamedConfigUnit("OUTER_LM_TPC_FT0A", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25ae_pass2_634103", kTPCFT0A, 80, 100), 
            {InputUnit("LHC25ae_pass2_634103", kTPCFT0A, 0, 20)}, 
            "LHC25ae_pass2_634103")),
        // ... more configurations
    };
}
```

### Ne-Ne System Configuration

**Location:** `Process_TemplateFit.cxx` lines 128-140

```cpp
if (systemName == "Ne-Ne") {
    return {
        NamedConfigUnit("FULL_LM_TPC_FT0A", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25af_pass2_645173", kTPCFT0A, 80, 100), 
            {InputUnit("LHC25af_pass2_645173", kTPCFT0A, 0, 20)}, 
            "LHC25af_pass2_645173")),
        NamedConfigUnit("FULL_MR_TPC_FT0C", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25af_pass2_645173", kTPCFT0C, 80, 100), 
            {InputUnit("LHC25af_pass2_645173", kTPCFT0C, 0, 20)}, 
            "LHC25af_pass2_645173")),  // ⚠️ SAME FILE AS FULL_LM!
        NamedConfigUnit("FULL_LR_FT0A_FT0C", ConfigUnit(kCent, false, false, 
            InputUnit("LHC25af_pass2_645746", kFT0AFT0C, 80, 100), 
            {InputUnit("LHC25af_pass2_645746", kFT0AFT0C, 0, 20)}, 
            "LHC25af_pass2_645746")),
        NamedConfigUnit("OUTER_LM_TPC_FT0A", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25af_pass2_637597", kTPCFT0A, 80, 100), 
            {InputUnit("LHC25af_pass2_637597", kTPCFT0A, 0, 20)}, 
            "LHC25af_pass2_637597")),
        NamedConfigUnit("INNER_MR_TPC_FT0C", ConfigUnit(kCent, false, true,  
            InputUnit("LHC25af_pass2_631290", kTPCFT0C, 80, 100), 
            {InputUnit("LHC25af_pass2_631290", kTPCFT0C, 0, 20)}, 
            "LHC25af_pass2_631290")),
        // ... more configurations
    };
}
```

**CRITICAL DIFFERENCE:** Ne-Ne FULL configuration uses **same file** (LHC25af_pass2_645173) for both LM and MR!

---

## 4. 3times2PC CORRELATION COMBINATION

### Location: `Process_3times2PC.cxx` lines 606-700

```cpp
VnUnit* Vn3times2PC(Bool_t isNch, Bool_t isPtDiff, std::vector<InputUnit> dataList, 
                    std::vector<TFile*>& vnDeltaFiles, Double_t pTMin=0, Double_t pTMax=0) {
    if (dataList.size() != 3) {
        Printf("The dataList don't have a size of 3!");
        return nullptr;
    }
    
    // Retrieve central values
    VnUnit* resultsLM = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[0], isPtDiff, 
                                                    dataList[0].minRange, dataList[0].maxRange, 
                                                    pTMin, pTMax);
    VnUnit* resultsMR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[1], isPtDiff, 
                                                    dataList[1].minRange, dataList[1].maxRange, 
                                                    pTMin, pTMax);
    VnUnit* resultsLR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[2], kFALSE, 
                                                    dataList[2].minRange, dataList[2].maxRange, 
                                                    pTMin, pTMax);
    
    std::cout << "[3times2PC Debug] LM (TPC_FT0A): v2=" << resultsLM->v2 << " +/- " 
              << resultsLM->v2_err << std::endl;
    std::cout << "[3times2PC Debug] MR (TPC_FT0C): v2=" << resultsMR->v2 << " +/- " 
              << resultsMR->v2_err << std::endl;
    std::cout << "[3times2PC Debug] LR (FT0_FT0C): v2=" << resultsLR->v2 << " +/- " 
              << resultsLR->v2_err << std::endl;
    
    // Combine using 3times2PC formula
    VnUnit* vnResult = Get3times2PC(resultsLM, resultsMR, resultsLR);
    
    // Bootstrap error calculation (IDENTICAL for all systems)
    // ... [same code path for O-O and Ne-Ne]
    
    return vnResult;
}
```

**IDENTICAL CODE PATH** - no system-specific branching or scaling applied.

### Histogram Value Retrieval with Inverse-Variance Weighting

**Location:** `Process_3times2PC.cxx` lines 735-780

```cpp
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
        if (std::isfinite(val2) && std::isfinite(err2) && err2 > 0.0 && 
            err2 < 9.9 && val2 >= 0.0) {
            double w2 = 1.0 / (err2 * err2);  // Inverse variance weighting
            numV2 += w2 * val2;
            denV2 += w2;
            ++validV2;
        }
    }
    
    if (denV2 > 0.0 && validV2 > 0) {
        v2 = numV2 / denV2;                  // Weighted average
        v2_err = std::sqrt(1.0 / denV2);     // Error from weights
    }
    // ... same for v3, v4
} else {
    // Single-bin histogram: use the single value directly
    v2 = hv2->GetBinContent(1);
    v2_err = hv2->GetBinError(1);
}
```

**IDENTICAL WEIGHTING SCHEME** for both O-O and Ne-Ne.

---

## 5. MIX-EVENT NORMALIZATION

### Location: `include/BasicForDihadron.h` lines 78-82

```cpp
std::map<int, float> MixEventNormalizationEta = {
    {kTPCFT0A, -4.25},    // TPC × FT0A: normalized at eta = -4.25
    {kTPCFT0C, 2.5},      // TPC × FT0C: normalized at eta = 2.5
    {kFT0AFT0C, 6.9}      // FT0A × FT0C: normalized at eta = 6.9
};
```

**SYSTEM-AGNOSTIC** - same normalization positions used for all systems (O-O, Ne-Ne, p-O, pp).

### Usage in Correlation Extraction

**Location:** `Process_dPhidEta.cxx` lines 765-766

```cpp
Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
```

No branching based on collision system.

---

## 6. PROCESSING CODE PATHS: IDENTICAL FOR BOTH SYSTEMS

### Main Processing Loop: `Process_3times2PC.cxx` lines 292-380

```cpp
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Open VnDelta files (same logic for all systems)
    TFile* vnDeltaFile_LM = new TFile(Form("./%s/EtaDiff/VnDelta_%s_%s_%i_%i_%s.root",
        "TemplateFit", outputFileName.c_str(), splitName.c_str(), 
        dataList[0].minRange, dataList[0].maxRange, 
        DihadronCorrTypeName[dataList[0].corrType].c_str()), "READ");
    
    TFile* vnDeltaFile_MR = new TFile(Form("./%s/EtaDiff/VnDelta_%s_%s_%i_%i_%s.root",
        "TemplateFit", outputFileName.c_str(), splitName.c_str(), 
        dataList[1].minRange, dataList[1].maxRange, 
        DihadronCorrTypeName[dataList[1].corrType].c_str()), "READ");
    
    TFile* vnDeltaFile_LR = new TFile(Form("./%s/VnDelta_%s_%s_%s.root",
        "TemplateFit", outputFileName.c_str(), outputSplitName.c_str(), 
        DihadronCorrTypeName[dataList[2].corrType].c_str()), "READ");
    
    // Process configurations with identical code path
    std::vector<TFile*> vnDeltaFiles = {vnDeltaFile_LM, vnDeltaFile_MR, vnDeltaFile_LR};
    
    VnUnit* vn = Vn3times2PC(isNch, kFALSE, dataList, vnDeltaFiles);
    // ... write results ...
}
```

**NO SYSTEM-SPECIFIC CONDITIONALS** - same processing regardless of O-O vs Ne-Ne.

---

## 7. SUMMARY TABLE: DIFFERENCES vs SIMILARITIES

| Aspect | O-O | Ne-Ne | Status |
|--------|-----|-------|--------|
| **Fitting Function** | `TemplateFitter` | `TemplateFitter` | ✅ IDENTICAL |
| **Initial v2 Parameter** | 0.004 | 0.004 | ✅ IDENTICAL |
| **3times2PC Formula** | √((LM×MR)/LR) | √((LM×MR)/LR) | ✅ IDENTICAL |
| **Error Propagation** | Quadrature sum | Quadrature sum | ✅ IDENTICAL |
| **Weighting Scheme** | Inverse-variance | Inverse-variance | ✅ IDENTICAL |
| **Normalization Eta** | -4.25, 2.5, 6.9 | -4.25, 2.5, 6.9 | ✅ IDENTICAL |
| **Centrality Ranges** | 0-20%, 80-100% | 0-20%, 80-100% | ✅ IDENTICAL |
| **Datasets** | LHC25ae_* | LHC25af_* | ❌ DIFFERENT |
| **Explicit Scaling** | None found | None found | ✅ NO SCALING |
| **System Weights** | None found | None found | ✅ NO WEIGHTS |

---

## 8. POSSIBLE SOURCES OF O-O > Ne-Ne DIFFERENCE

Since **no explicit mathematical difference** is coded, the higher O-O values likely come from:

1. **Different physical data** - O-O and Ne-Ne are distinct collision systems
2. **Dataset statistics variation** - Different pass2 productions may have different event yields
3. **Detector response asymmetry** - FT0 detector may respond differently to the neutron-rich Ne nucleus
4. **Event multiplicity distribution** - Different charged particle multiplicities between systems
5. **Tracking efficiency** - Different particle detection efficiency in each system

**Recommendation:** Check the input bootstrap sample statistics (event counts, multiplicity distributions) for O-O vs Ne-Ne to understand the source of the difference.

---

## 9. FILES ANALYZED

- `Process_TemplateFit.cxx` - v2 extraction via template fitting
- `Process_3times2PC.cxx` - v2 combination formula
- `Process_dPhidEta.cxx` - Correlation function input preparation
- `include/BasicForDihadron.h` - System constants and normalization

**No system-specific scaling or weighting code identified in any file.**
