void check_hists() {
    TFile f("ProcessOutput/EtaDiff/BootstrapSample_LHC25af_pass2_632504_Cent_80_100_Eta_-0.1_0.0_TPC_FT0A.root");
    auto h = f.Get("bsSample_hPhiSameOverMixed_80_100_0");
    std::cout << "Sample 0 hist: " << h << std::endl;
    auto h70 = f.Get("bsSample_hPhiSameOverMixed_80_100_70");
    std::cout << "Sample 70 hist: " << h70 << std::endl;
}
