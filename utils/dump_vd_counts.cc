// dump_vd_counts.cc  (build: g++ -O2 dump_vd_counts.cc $(root-config --cflags --libs) -o dump_vd_counts)

#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <TStopwatch.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <dirent.h>
#include <cstring>
#include <cstdlib>

static const char* kDefaultBaseDir = "/eos/cms/store/user/asimsek/SPSH4Results/v9";

int run_dump_counts(int initial_energy_GeV,
                 std::string input_dir = "",
                 std::string out_csv  = "")
{
    // Resolve paths from energy if not explicitly given
    if (input_dir.empty()) {
        std::ostringstream p; p << kDefaultBaseDir << "/" << initial_energy_GeV << "GeV";
        input_dir = p.str();
    }
    if (out_csv.empty()) {
        std::ostringstream o; o << "output_" << initial_energy_GeV << "GeV.csv";
        out_csv = o.str();
    }

    std::cout << "[info] energy=" << initial_energy_GeV
              << " GeV, in=" << input_dir
              << ", out=" << out_csv << std::endl;

    std::ofstream outfile(out_csv);
    if (!outfile.is_open()) {
        std::cerr << "[error] cannot open output file: " << out_csv << std::endl;
        return 1;
    }
    outfile << "Tree,All,Positron,Electron,Photon,Muon,Proton\n";

    const char* trees_list[] = {
        "BeamStartPoint","BeforeColi1","BeforeQuad4",
        "BeforeColi2","AfterColi3","AfterB5h","BeforeTrim02",
        "BeforeB6a","AfterB7b","AfterB7c","BeforeTrim03","AfterTrim03",
        "BeforeB8a","AafterB8b","AfterB9b","AfterCol10",
        "AfterQ24","NA64","BeforeXTDV22601","Det022640","BeforeMCPs",
        "AfterMCPs","ECALFront"
    };
    const int total_trees = int(sizeof(trees_list)/sizeof(trees_list[0]));

    // Collect ROOT files in input_dir
    std::vector<std::string> root_files;
    DIR* dir = opendir(input_dir.c_str());
    if (!dir) {
        std::cerr << "[error] cannot open directory: " << input_dir << std::endl;
        return 2;
    }
    if (dirent* e; (e = readdir(dir)) != nullptr) {
        // do nothing; we just needed to ensure readdir is usable in this build
    }
    rewinddir(dir);
    while (auto* entry = readdir(dir)) {
        const char* name = entry->d_name;
        if (!name || name[0]=='.') continue;
        if (strstr(name, ".root")) {
            root_files.emplace_back(input_dir + "/" + name);
        }
    }
    closedir(dir);
    std::sort(root_files.begin(), root_files.end());
    if (root_files.empty()) {
        std::cerr << "[error] no .root files found in: " << input_dir << std::endl;
        return 3;
    }

    // Speed-up if possible
    ROOT::EnableImplicitMT();

    TStopwatch sw; sw.Start();
    for (int i = 0; i < total_trees; ++i) {
        const std::string tree_name = trees_list[i];
        TChain chain(("VirtualDetector/" + tree_name).c_str());
        for (const auto& f : root_files) chain.Add(f.c_str());

        ROOT::RDataFrame df(chain);
        auto nAll       = df.Count();
        auto nPositron  = df.Filter("PDGid == -11").Count();
        auto nElectron  = df.Filter("PDGid == 11").Count();
        auto nPhoton    = df.Filter("PDGid == 22").Count();
        auto nMuon      = df.Filter("abs(PDGid) == 13").Count();
        auto nProton    = df.Filter("PDGid == 2212").Count();

        outfile << tree_name << ","
                << *nAll << ","
                << *nPositron << ","
                << *nElectron << ","
                << *nPhoton << ","
                << *nMuon << ","
                << *nProton << "\n";

        double progress = 100.0 * (i + 1) / double(total_trees);
        std::cout << "\rProcessing... " << std::fixed << std::setprecision(1)
                  << progress << "% [" << (i+1) << "/" << total_trees << "]" << std::flush;
    }
    sw.Stop();
    std::cout << "\n";
    sw.Print();
    std::cout << "[done] wrote " << out_csv << std::endl;
    return 0;
}

// CLI parser: --energy/-E, --in, --out or env INITIAL_ENERGY
int main(int argc, char** argv) {
    // Defaults
    int energy = 250;
    std::string inDir, outCsv;

    if (const char* envE = std::getenv("INITIAL_ENERGY")) {
        energy = std::atoi(envE);
    }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](int &i)->std::string {
            if (i+1 >= argc) { std::cerr << "[error] missing value for " << a << "\n"; std::exit(64); }
            return std::string(argv[++i]);
        };
        if (a == "--energy" || a == "-E") { energy = std::stoi(next(i)); }
        else if (a.rfind("--energy=",0)==0) { energy = std::stoi(a.substr(9)); }
        else if (a == "--in") { inDir = next(i); }
        else if (a.rfind("--in=",0)==0) { inDir = a.substr(5); }
        else if (a == "--out") { outCsv = next(i); }
        else if (a.rfind("--out=",0)==0) { outCsv = a.substr(6); }
        else if (a == "-h" || a == "--help") {
            std::cout <<
              "Usage: ./dump_vd_counts [--energy GeV] [--in DIR] [--out FILE]\n"
              "  You can also set INITIAL_ENERGY env var.\n";
            return 0;
        }
    }
    return run_dump_counts(energy, inDir, outCsv);
}

