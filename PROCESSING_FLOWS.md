# Key Processing Flows Comparison

## Reference Version - Full Data Processing Flow

```
Process_FourierFit() 
    ↓ (11 configs: Cent/Nch, multiple collision systems)
    ├─→ ProcessConfig() OR ProcessConfig_PtDiff()
    │
    └─→ For each dataset:
        ├─ FourierFit(isNch, data, cn2Tovn2)
        │  ├─ Open: ./ProcessOutput/BootstrapSample_*.root
        │  │         (STRICT path, fail if not found)
        │  │
        │  ├─ fitSample(sample=-1):
        │  │  ├─ Load: bsSample_hPhiSameOverMixed_*
        │  │  ├─ RooFourierFit(hist)
        │  │  │  ├─ HistFFT(hist, Coeff)  ← FFT-based initial guess
        │  │  │  │  ├─ a0 = Coeff[0]
        │  │  │  │  ├─ a1 = Coeff[1]
        │  │  │  │  ├─ a2 = Coeff[2]
        │  │  │  │  ├─ a3 = Coeff[3]
        │  │  │  │  └─ a4 = Coeff[4]
        │  │  │  ├─ Fit: a0 + 2a1*cos(x) + ... + 2a4*cos(4x)
        │  │  │  └─ Returns: v1², v2², v3², v4² (4 harmonics!)
        │  │  ├─ VnUnit(v1, v1_err, v2, v2_err, v3, v3_err, v4, v4_err)
        │  │  └─ PlotFitting() → PDFs
        │  │
        │  └─ Loop Nsamples times: fitSample(sample=i)
        │     └─ Accumulate bootstrap errors
        │
        └─ Output: ./FourierFit/VnDelta_*.root
           ├─ hV1, hV2, hV3, hV4  (4 histograms)
           └─ Subsamples/hV1_subsample_i, hV2_subsample_i, ...
```

## EtaDiff Version - Eta-Differential Processing Flow

```
Process_FourierFit()
    ↓ (3 configs: Ne-Ne Cent 0-20, 3 detector combinations)
    └─→ ProcessConfig() ONLY (PtDiff disabled)
        ↓
        └─→ For each dataset:
            ├─ FourierFit(isNch, data, cn2Tovn2)
            │  ├─ TRY: Open ./ProcessOutput/BootstrapSample_*.root
            │  │
            │  ├─ IF FAILS, FALLBACK:  ← NEW IN EtaDiff
            │  │  └─ Open ./ProcessOutput/EtaDiff/BootstrapSample_*_Eta_0.0_0.1_*.root
            │  │     (Allows use of eta-binned bootstrap samples)
            │  │
            │  ├─ fitSample(sample=-1):
            │  │  ├─ Load: bsSample_hPhiSameOverMixed_*
            │  │  ├─ RooFourierFit(hist)
            │  │  │  ├─ histMean = hist→GetMean()     ← Simpler init
            │  │  │  ├─ histRMS = hist→GetRMS()
            │  │  │  ├─ a0 = histMean
            │  │  │  ├─ a1 = 0.1 * histRMS
            │  │  │  ├─ a2 = 0.1 * histRMS
            │  │  │  ├─ a3 = 0.1 * histRMS
            │  │  │  └─ a4 = 0.1 * histRMS
            │  │  ├─ Fit: a0 + 2a1*cos(x) + ... + 2a4*cos(4x)
            │  │  └─ Returns: v2², v3², v4² only (v1 REMOVED!)
            │  │
            │  ├─ VnUnit(v2, v2_err, v3, v3_err, v4, v4_err)
            │  │  (Constructor simplified - no v1)
            │  │
            │  └─ Loop Nsamples times: fitSample(sample=i)
            │
            └─ Output: ./FourierFit/VnDelta_*.root
               ├─ hV2, hV3, hV4  (3 histograms only)
               └─ Subsamples/hV2_subsample_i, hV3_subsample_i, ...
```

## Critical Function: RooFourierFit Design Change

### Reference: FFT-Based Initial Parameters
```cpp
// Pro: Uses frequency content of data
// Con: FFT coefficients can be very large → numerical issues
// Suitable for: General centrality-based analysis

std::vector<double> Coeff;
HistFFT(hist, Coeff);           // Compute FFT coefficients
fitFunc→SetParameter(0, Coeff[0]);
for (int i = 1; i < 5; ++i) {
    fitFunc→SetParameter(i, Coeff[i]);  // Use FFT directly
}
```

### EtaDiff: Histogram Statistics-Based Parameters
```cpp
// Pro: Conservative, stable starting values
// Con: May miss initial frequency content
// Suitable for: Eta-differential analysis where fewer harmonics present

double histMean = hist→GetMean();
double histRMS = hist→GetRMS();
fitFunc→SetParameter(0, histMean);          // 0.1×RMS is very small
for (int i = 1; i < 5; ++i) {
    fitFunc→SetParameter(i, 0.1 * histRMS); // Conservative start
}
```

## Bootstrap Error Calculation - Both Versions

```
Same Nobs = 3 (v2², v3², v4²)

Reference:
  ValueArray[0] = v2² ↘
  ValueArray[1] = v3² ├→ CalculateBootstrapError() → ErrorArray
  ValueArray[2] = v4² ↗
  
EtaDiff (identical logic):
  ValueArray[0] = v2² ↘
  ValueArray[1] = v3² ├→ CalculateBootstrapError() → ErrorArray
  ValueArray[2] = v4² ↗

Result: Both versions calculate bootstrap errors identically
        (v1² excluded in EtaDiff from the start)
```

## Input File Organization

### Reference Pattern:
```
./ProcessOutput/
  ├─ BootstrapSample_LHC25ae_pass2_638219_Cent_0_20_*.root
  ├─ BootstrapSample_LHC25ac_pass2_639097_Cent_10_50_*.root
  ├─ BootstrapSample_LHC25ad_pass2_635471_Mult_10_50_*.root
  └─ ... (many collision systems)

./ProcessOutput/PtDiff/
  ├─ BootstrapSample_..._Pt_0.5_1.0_*.root
  ├─ BootstrapSample_..._Pt_1.0_1.5_*.root
  └─ ... (pT bins)
```

### EtaDiff Pattern:
```
./ProcessOutput/
  ├─ BootstrapSample_LHC25af_pass2_632504_Cent_0_20_*.root
  ├─ BootstrapSample_LHC25af_pass2_637596_Cent_0_20_*.root
  └─ BootstrapSample_LHC25af_pass2_640018_Cent_0_20_*.root

./ProcessOutput/EtaDiff/
  ├─ BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_0.0_0.1_*.root
  ├─ BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_-0.1_0.0_*.root
  ├─ BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_0.1_0.2_*.root
  ├─ BootstrapSample_LHC25af_pass2_637596_Cent_0_20_Eta_0.0_0.1_*.root
  └─ ... (eta bins for each dataset)

Note: No PtDiff/ directory used in EtaDiff version
```

## Output File Naming

### Reference:
```
VnDelta_LHC25ae_pass2_638219_Cent_FT0AFT0C.root
VnDelta_LHC25ae_pass2_638219_Mult_TPCFT0A.root  (for Nch-based)
PtDiff/VnDelta_LHC25ae_pass2_638219_Cent_0_20_TPCFT0A.root
```

### EtaDiff:
```
VnDelta_LHC25af_pass2_632504_Cent_TPCFT0A.root
VnDelta_LHC25af_pass2_637596_Cent_TPCFT0C.root
VnDelta_LHC25af_pass2_640018_Cent_FT0AFT0C.root
(Only Cent-based, no PtDiff variants)
```

## Configuration Management

### Reference main() - Comprehensive:
```cpp
Process_FourierFit() loads ~17 datasets:
  - O-O (kCent, 3 configs)
  - O-O with PtDiff (kCent, 2 configs)
  - O-O Nch (kNch, 5 configs)
  - p-O Nch (kNch, 5 configs)
  - p-p Nch (kNch, 5 configs)
  - Ne-Ne Nch (kNch, 5 configs)
  - MC Nch (kNch, 5 configs)

Most are commented out or active based on current needs
```

### EtaDiff main() - Focused:
```cpp
Process_FourierFit() loads 3 datasets only:
  - LHC25af_pass2_632504 (TPCFT0A)
  - LHC25af_pass2_637596 (TPCFT0C)
  - LHC25af_pass2_640018 (FT0AFT0C)

All are Ne-Ne, all are Cent 0-20, all are kPtDiffOff
Specialized for EtaDiff analysis only
```
