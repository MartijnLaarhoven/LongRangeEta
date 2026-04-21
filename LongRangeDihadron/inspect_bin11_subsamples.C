#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TDirectory.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

void inspect_bin11_subsamples() {
    TFile f("./TemplateFit/EtaDiff/VnDelta_LHC25ae_pass2_648788_Cent_0_20_TPC_FT0A.root", "READ");
    if (!f.IsOpen() || f.IsZombie()) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }

    TDirectory* sub = dynamic_cast<TDirectory*>(f.Get("Subsamples"));
    if (!sub) {
        std::cout << "No Subsamples directory" << std::endl;
        return;
    }

    std::vector<double> vals;
    std::vector<double> errs;
    int sentinel = 0;
    int largeErr = 0;

    TIter next(sub->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next())) {
        std::string name = key->GetName();
        if (name.rfind("hV2_subsample_", 0) != 0) continue;
        TH1D* h = dynamic_cast<TH1D*>(sub->Get(name.c_str()));
        if (!h) continue;

        double v = h->GetBinContent(11);
        double e = h->GetBinError(11);
        vals.push_back(v);
        errs.push_back(e);

        if (std::fabs(v + 1.0) < 1e-9 && std::fabs(e - 10.0) < 1e-9) sentinel++;
        if (e > 0.01) largeErr++;
    }

    std::cout << "N subsamples: " << vals.size() << std::endl;
    std::cout << "Sentinel (-1,10) count: " << sentinel << std::endl;
    std::cout << "Large error (>0.01) count: " << largeErr << std::endl;

    if (vals.empty()) return;

    std::vector<double> sortedVals = vals;
    std::sort(sortedVals.begin(), sortedVals.end());
    double median = sortedVals[sortedVals.size()/2];

    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (double x : vals) var += (x - mean) * (x - mean);
    var /= vals.size();
    double rms = std::sqrt(var);

    auto minmaxV = std::minmax_element(vals.begin(), vals.end());
    auto minmaxE = std::minmax_element(errs.begin(), errs.end());

    std::cout << "v2 bin11 mean=" << mean << " median=" << median << " rms=" << rms
              << " min=" << *minmaxV.first << " max=" << *minmaxV.second << std::endl;
    std::cout << "err bin11 min=" << *minmaxE.first << " max=" << *minmaxE.second << std::endl;
}
