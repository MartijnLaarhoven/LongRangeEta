# Quick Reference: O-O vs Ne-Ne Processing

## 🎯 Key Finding
**No system-specific scaling or normalization factors** are applied in the code. Both systems use identical mathematical formulas, just different input datasets.

---

## ⚡ Core Formula

```cpp
// IDENTICAL for both O-O and Ne-Ne
v_n = √((v_n_LM × v_n_MR) / v_n_LR)

where:
  LM = TPC × FT0A correlation
  MR = TPC × FT0C correlation  
  LR = FT0A × FT0C correlation
```

---

## 📊 Dataset Differences

### O-O (Centrality 0-20%)
```
FULL_LM: LHC25ae_pass2_648788  (TPC×FT0A)
FULL_MR: LHC25ae_pass2_648788  (TPC×FT0C)
FULL_LR: LHC25ae_pass2_648798  (FT0A×FT0C)
```

### Ne-Ne (Centrality 0-20%)
```
FULL_LM: LHC25af_pass2_645173  (TPC×FT0A)
FULL_MR: LHC25af_pass2_645173  (TPC×FT0C)  ⚠️ SAME FILE!
FULL_LR: LHC25af_pass2_645746  (FT0A×FT0C)
```

---

## 🔬 Processing Steps (IDENTICAL CODE)

### 1️⃣ Template Fitting
```cpp
// Process_TemplateFit.cxx - SAME FOR BOTH SYSTEMS
ft->AddParameter("v2", "v2", 4e-3, -1.0, 1.0);    // Initial guess: 0.004
ft->AddParameter("v3", "v3", 6e-4, -1.0, 1.0);    // Initial guess: 0.0006
ft->AddParameter("v4", "v4", 1.8e-4, -1.0, 1.0);  // Initial guess: 0.00018
```

### 2️⃣ Inverse-Variance Weighting (for eta-differential data)
```cpp
// Process_3times2PC.cxx - SAME FOR BOTH SYSTEMS
for (int i = 1; i <= nbins; ++i) {
    double w = 1.0 / (err * err);    // Weight by inverse variance
    weighted_val += w * val;
    weight_sum += w;
}
result = weighted_val / weight_sum;
```

### 3️⃣ 3times2PC Combination
```cpp
// Process_3times2PC.cxx - SAME FOR BOTH SYSTEMS
double Get3times2PC(double LM, double MR, double LR) {
    double vn = (LM * MR) / LR;
    return (vn >= 0) ? sqrt(vn) : -1;
}
```

### 4️⃣ Error Propagation
```cpp
// Process_3times2PC.cxx - SAME FOR BOTH SYSTEMS
part1 = 0.5 * pow(LM, -0.5) * pow(MR, 0.5) * pow(LR,-0.5) * LM_err;
part2 = 0.5 * pow(LM, 0.5) * pow(MR, -0.5) * pow(LR,-0.5) * MR_err;
part3 = 0.5 * pow(LM, 0.5) * pow(MR, 0.5) * pow(LR,-1.5) * LR_err;
total_err = sqrt(part1² + part2² + part3²);
```

---

## ✅ Checked for System-Specific Variations

| Check | Result |
|-------|--------|
| Explicit system-specific scaling | ❌ NOT FOUND |
| System-dependent weights | ❌ NOT FOUND |
| Conditional branching for O-O vs Ne-Ne | ❌ NOT FOUND (except dataset selection) |
| Different fitting parameters | ❌ NOT FOUND |
| Different error propagation | ❌ NOT FOUND |
| Different normalization constants | ❌ NOT FOUND |

---

## 🔍 Why O-O > Ne-Ne?

**The difference is NOT in the code logic.** Likely sources:

1. **Different collision physics**
   - O-O: symmetric system (8 protons, 8 neutrons each)
   - Ne-Ne: symmetric system (10 protons, 10 neutrons each)
   - Different nuclear structure → different azimuthal anisotropy

2. **Different detector responses**
   - FT0 detector may respond differently to different nuclei
   - Different particle production patterns

3. **Event statistics**
   - Different event yields in different datasets
   - Different multiplicity distributions

4. **Tracking efficiency**
   - Different tracking performance for different systems

---

## 📁 Files with Key Code

| File | Function | Lines |
|------|----------|-------|
| `Process_TemplateFit.cxx` | v2 extraction | 1224-1232 |
| `Process_3times2PC.cxx` | 3times2PC formula | 813-825 |
| `Process_3times2PC.cxx` | Combination logic | 606-700 |
| `Process_3times2PC.cxx` | Weighting | 735-780 |
| `Process_dPhidEta.cxx` | Dataset config | 68-200 |
| `BasicForDihadron.h` | Normalization eta | 78-82 |

---

## 🔐 Conclusion

The O-O vs Ne-Ne difference in v₂ values is **NOT due to code-implemented scaling or system-specific weighting**. Both systems:

- Use identical fitting methodology
- Use identical 3times2PC formula
- Use identical error propagation
- Use identical weighting schemes

The difference emerges from **different physical properties of the collision systems themselves**, reflected in the different input datasets.
