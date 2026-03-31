void examine_input() {
    TFile* f = TFile::Open("../../../AnalysisResultsROOTFiles/LongRangeEta/AnalysisResults_LHC25af_pass2_640018.root");
    if (!f || f->IsZombie()) {
        cout << "ERROR: Cannot open input file!" << endl;
        return;
    }
    
    cout << "\n=== AVAILABLE FLOWS ===" << endl;
    f->ls();
    
    cout << "\n=== flow-decorrelation_Cent_0_20 ===" << endl;
    TDirectory* dir = (TDirectory*)f->Get("flow-decorrelation_Cent_0_20");
    if (dir) {
        dir->ls();
    } else {
        cout << "Not found!" << endl;
    }
    
    f->Close();
}
