void check_ft0a_FT0C() {
    TFile* f = TFile::Open("ProcessOutput/Mixed_LHC25af_pass2_640018_Cent_0_20_FT0A_FT0C.root");
    if (!f || f->IsZombie()) {
        cout << "ERROR: Cannot open file!" << endl;
        return;
    }
    
    cout << "\n=== FILE CONTENTS ===" << endl;
    f->ls();
    
    TH1D* h = (TH1D*)f->Get("hPhiSameOverMixed_0_20");
    if (h) {
        cout << "\n=== hPhiSameOverMixed_0_20 ===" << endl;
        cout << "Integral: " << h->Integral() << endl;
        cout << "Max: " << h->GetMaximum() << endl;
        cout << "Entries: " << h->GetEntries() << endl;
        cout << "Bins: " << h->GetNbinsX() << endl;
    } else {
        cout << "\nHistogram hPhiSameOverMixed_0_20 not found!" << endl;
    }
    
    TH2D* h2d = (TH2D*)f->Get("hPhiEtaSameOverMixed_0_20");
    if (h2d) {
        cout << "\n=== hPhiEtaSameOverMixed_0_20 2D ===" << endl;
        cout << "Integral: " << h2d->Integral() << endl;
        cout << "Max: " << h2d->GetMaximum() << endl;
        cout << "Entries: " << h2d->GetEntries() << endl;
    } else {
        cout << "\nhPhiEtaSameOverMixed_0_20 2D not found!" << endl;
    }
    
    f->Close();
}
