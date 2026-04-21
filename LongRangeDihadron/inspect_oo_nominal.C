#include <TFile.h>
#include <TH1D.h>
#include <iostream>

void printHist(const char* path, const char* histName) {
    TFile f(path, "READ");
    if (!f.IsOpen() || f.IsZombie()) {
        std::cout << path << " [open failed]" << std::endl;
        return;
    }
    TH1D* h = dynamic_cast<TH1D*>(f.Get(histName));
    if (!h) {
        std::cout << path << " missing " << histName << std::endl;
        return;
    }
    std::cout << path << " :: " << histName << std::endl;
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        std::cout << "  bin " << bin << " = " << h->GetBinContent(bin) << " +/- " << h->GetBinError(bin) << std::endl;
    }
}

void inspect_oo_nominal() {
    printHist("./TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_648788_Cent_0_20_TPC_FT0A.root", "hV2");
    printHist("./TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_648788_Cent_0_20_TPC_FT0C.root", "hV2");
    printHist("./TemplateFit/VnDelta_LHC25ae_pass2_648788_Cent_FT0A_FT0C.root", "hV2");
    printHist("./3times2PC/Vn_LHC25ae_pass2_648798_Cent_0_20.root", "hV2");
}
