# Process_FourierFit.cxx Comparison - Complete Analysis Index

This directory contains a comprehensive comparison between two versions of `Process_FourierFit.cxx`:

- **Reference:** `/home/martijn-laarhoven/Work/dihadronanalysis-master/LongRangeDihadron/Process_FourierFit.cxx` (1224 lines)
- **EtaDiff:** `/home/martijn-laarhoven/Work/dihadronanalysis-master/LongRangeEta/dihadronanalysis/LongRangeDihadron/Process_FourierFit.cxx` (968 lines)

---

## Document Structure

### 1. [DETAILED_COMPARISON.md](DETAILED_COMPARISON.md) - START HERE
**Best for:** Understanding all differences comprehensively

**Contains:**
- Section 1: VnUnit structure differences (v1 removal)
- Section 2: ProcessConfig function - input file handling & error handling
- Section 3: FourierFit function - input file paths & fallback strategy
- Section 4: RooFourierFit - FFT vs. histogram statistics (CRITICAL CHANGE)
- Section 5: Nobs and bootstrap error calculation
- Section 6: Output file naming and structure
- Section 7: fitSample function - histogram names and parameter return order
- Section 8: RooFourierFit output parameters - array structure
- Section 9: PlotFitting function - output differences
- Section 10: Dataset configuration in main function
- Summary table comparing all aspects

### 2. [PROCESSING_FLOWS.md](PROCESSING_FLOWS.md) - FOR VISUALIZATION
**Best for:** Understanding workflow and data flow patterns

**Contains:**
- Reference version data processing flow diagram
- EtaDiff version data processing flow diagram
- Critical RooFourierFit design comparison (FFT-based vs. statistics-based)
- Bootstrap error calculation visualization
- Input file organization patterns
- Output file naming conventions
- Configuration management comparison

### 3. [CODE_COMPARISON.md](CODE_COMPARISON.md) - FOR IMPLEMENTATION DETAILS
**Best for:** Understanding exact code changes and line-by-line differences

**Contains:**
- Side-by-side code snippets for 8 key functions:
  1. VnUnit structure
  2. FourierFit input file handling
  3. RooFourierFit initial parameters (MOST IMPORTANT CHANGE)
  4. fitSample return statement
  5. Bootstrap error calculation
  6. Output histogram creation
  7. PlotFitting legend display
  8. Main function dataset configuration
- Summary table with key metrics

---

## Key Findings at a Glance

### CRITICAL CHANGES

**1. V1 Harmonic Removed Completely** (EtaDiff-specific)
```
Reference: Outputs v1, v2, v3, v4 (4 harmonics)
EtaDiff:   Outputs v2, v3, v4 only (3 harmonics)
```
- VnUnit structure changed (8 → 6 members)
- All subsamples omit v1
- Output files 25% smaller

**2. Initial Parameter Setup in RooFourierFit** (MOST IMPACTFUL)
```
Reference: Uses HistFFT(hist, Coeff) 
           - Computes FFT coefficients
           - Sophisticated frequency-domain initialization
           - Risk: Large coefficients cause numerical issues
           
EtaDiff:   Uses histogram->GetMean() and histogram->GetRMS()
           - Simple statistical initialization: a0=Mean, a{i}=0.1*RMS
           - Conservative, stable approach
           - Explicit comment: "Use histogram statistics rather than 
                               potentially huge FFT coefficients"
```

**3. Input File Fallback Mechanism** (Robustness feature)
```
Reference: ./ProcessOutput/BootstrapSample_*.root ONLY
           If missing → return nullptr
           
EtaDiff:   PRIMARY:  ./ProcessOutput/BootstrapSample_*.root
           FALLBACK: ./ProcessOutput/EtaDiff/BootstrapSample_*_Eta_0.0_0.1_*.root
           If both fail → exit(1)
```
- Allows graceful use of eta-binned bootstrap samples
- More robust to missing centrality-based files

### MAJOR STRUCTURAL CHANGES

**4. Scope and Specialization**
- Reference: General-purpose FourierFit processor
  - Supports 17 different dataset configurations
  - Multiple collision systems (O-O, p-O, p-p, Ne-Ne, MC)
  - Handles both Centrality and Nch-based analyses
  - PtDiff support enabled
  
- EtaDiff: Specialized eta-differential processor
  - Only 3 Ne-Ne configurations
  - Centrality-based only (0-20%)
  - PtDiff support DISABLED
  - Single collision system focus

**5. Error Handling Philosophy**
- Reference: Graceful degradation (skips failed datasets)
- EtaDiff: Fail-fast with explicit exits

### COMPUTATIONAL DIFFERENCES

**6. Bootstrap Error Calculation** (Identical)
- Both use Nobs=3 (v2², v3², v4² only)
- Both calculate errors identically
- Both return 3 ErrorArray values
- Difference: EtaDiff stores fewer members in structures

**7. Fitting Formula** (Identical)
```
Both: a0 + 2a1*cos(x) + 2a2*cos(2x) + 2a3*cos(3x) + 2a4*cos(4x)

Where:
- a0 = baseline constant
- a1 = v1²*baseline coefficient
- a2 = v2²*baseline coefficient  
- a3 = v3²*baseline coefficient
- a4 = v4²*baseline coefficient
```

---

## Physics Context

### Why v1 is Removed in EtaDiff
- **v1 (Directed Flow):** Sensitive to forward-backward asymmetry, bulk elliptic response
- **EtaDiff Analysis:** Focuses on η-dependent correlations at fixed η⊙Δη ranges
- **Rationale:** v1 signal is washed out in eta-differential analysis due to symmetry assumptions
- **Result:** EtaDiff optimizes for v2, v3, v4 only (stronger eta-differential signals)

### Why HistFFT is Replaced
- **HistFFT Approach:** Uses Fast Fourier Transform coefficients as initial parameters
  - Advantage: Leverages frequency content of data
  - Disadvantage: FFT amplifies noise → very large coefficients → numerical issues
  
- **Histogram Statistics Approach:** Uses simple mean/RMS statistics
  - Advantage: Conservative, stable, prevents numerical runaway
  - Disadvantage: Loses frequency-domain information
  - Trade-off: More robust but potentially slower convergence

---

## File I/O Organization

### Reference Input Pattern
```
./ProcessOutput/
├── BootstrapSample_LHC25ae_pass2_638219_Cent_0_20_FT0AFT0C.root
├── BootstrapSample_LHC25ae_pass2_638219_Cent_0_20_TPCFT0A.root
├── BootstrapSample_LHC25ae_pass2_640606_Mult_10_50_FT0AFT0C.root
└── [many more for multiple systems]

./ProcessOutput/PtDiff/
├── BootstrapSample_..._Pt_0.5_1.0_*.root
├── BootstrapSample_..._Pt_1.0_1.5_*.root
└── [many pT bins]
```

### EtaDiff Input Pattern with Fallback
```
PRIMARY PATH (Centrality):
./ProcessOutput/BootstrapSample_LHC25af_pass2_632504_Cent_0_20_*.root

FALLBACK PATH (Eta-differential bins):
./ProcessOutput/EtaDiff/BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_-0.8_-0.7_*.root
./ProcessOutput/EtaDiff/BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_-0.7_-0.6_*.root
./ProcessOutput/EtaDiff/BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_0.0_0.1_*.root
./ProcessOutput/EtaDiff/BootstrapSample_LHC25af_pass2_632504_Cent_0_20_Eta_0.1_0.2_*.root
[... many eta bins ...]
```

---

## Output Comparison

### Reference Output Files
```
./FourierFit/VnDelta_LHC25ae_pass2_638219_Cent_FT0AFT0C.root
├── hV1 (histogram with v1 values)
├── hV2 (histogram with v2 values)
├── hV3 (histogram with v3 values)
├── hV4 (histogram with v4 values)
└── Subsamples/
    ├── hV1_subsample_0, hV1_subsample_1, ...
    ├── hV2_subsample_0, hV2_subsample_1, ...
    ├── hV3_subsample_0, hV3_subsample_1, ...
    └── hV4_subsample_0, hV4_subsample_1, ...

./FourierFit/PtDiff/VnDelta_LHC25ae_pass2_638219_Cent_0_20_TPCFT0A.root
└── [similar structure with pT dependence]

./FourierFit/PDFs/
├── FourierFit_LHC25ae_pass2_638219_Cent_0_20_FT0AFT0C.pdf
├── FourierFit_LHC25ae_pass2_638219_Cent_0_20_TPCFT0A.pdf
└── [many other fit plots]
```

### EtaDiff Output Files
```
./FourierFit/VnDelta_LHC25af_pass2_632504_Cent_TPCFT0A.root
├── hV2 (no hV1!)
├── hV3
├── hV4
└── Subsamples/
    ├── hV2_subsample_0, hV2_subsample_1, ...
    ├── hV3_subsample_0, hV3_subsample_1, ...
    └── hV4_subsample_0, hV4_subsample_1, ...

./FourierFit/PDFs/
├── FourierFit_LHC25af_pass2_632504_Cent_0_20_TPCFT0A.pdf
├── FourierFit_LHC25af_pass2_637596_Cent_0_20_TPCFT0C.pdf
└── ... (only 3 Ne-Ne configs)

Note: No PtDiff/ directory created
```

---

## Line Count Comparison

| Aspect | Reference | EtaDiff | Difference |
|--------|-----------|---------|-----------|
| Total lines | 1224 | 968 | -256 lines (-21%) |
| VnUnit members | 16 (8×2) | 12 (6×2) | -4 members |
| Output histograms | 4 | 3 | -1 histogram |
| Dataset configs | 17 | 3 | -14 configs |
| Main function size | Large | Compact | Much simpler |

---

## Migration Path (if needed)

If you need to apply specific changes from EtaDiff back to Reference or vice versa:

### To Enable EtaDiff Features in Reference:
1. Remove v1 from VnUnit struct (requires careful index updates)
2. Replace HistFFT with histogram statistics in RooFourierFit
3. Add EtaDiff fallback path in FourierFit function
4. Update fitSample return to omit v1
5. Reduce output histograms from 4 to 3

### To Add Reference Features to EtaDiff:
1. Restore v1 member in VnUnit (full reimplementation)
2. Implement HistFFT-based parameter initialization
3. Enable ProcessConfig_PtDiff functionality
4. Expand dataset configuration in main()
5. Support multiple collision systems

---

## Recommendations

### Use Reference Version If:
- You need comprehensive multi-system analysis (O-O, p-O, p-p, Ne-Ne)
- pT-differential analysis is required
- You need v1 directed flow measurements
- Numerical stability with FFT is acceptable
- You want general-purpose framework

### Use EtaDiff Version If:
- Focus on eta-differential dihadron correlations
- Only v2, v3, v4 harmonics needed
- Robustness over sophistication (statistics vs. FFT)
- Simpler, more focused codebase preferred
- Memory/output size constraints exist
- Fallback file handling needed for missing centrality data

---

## Cross-Reference Guide

| Question | Document | Section |
|----------|----------|---------|
| What's the main difference? | DETAILED_COMPARISON | Summary table |
| How does file I/O work? | PROCESSING_FLOWS | Input file organization |
| What changed in RooFourierFit? | CODE_COMPARISON | Section 3 (CRITICAL) |
| How does bootstrap work? | CODE_COMPARISON | Section 5 |
| What datasets are processed? | PROCESSING_FLOWS | Configuration management |
| How are v1 values handled? | DETAILED_COMPARISON | Sections 1, 7 |
| Where do outputs go? | PROCESSING_FLOWS | Output file naming |
| Why remove HistFFT? | DETAILED_COMPARISON | Section 4 |
| What's the error handling? | DETAILED_COMPARISON | Section 2, 3 |
| How are plots generated? | CODE_COMPARISON | Section 7 |

---

**Last Updated:** 2025-03-25
**Comparison Scope:** Complete feature-by-feature analysis
**Total Documentation:** 3 detailed markdown files
