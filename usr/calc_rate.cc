#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TStyle.h"

struct ScrData {
    uint32_t header_magic = 0;
    uint32_t n_words = 0;
    uint32_t heartbeat_count = 0;
    uint32_t heartbeat_frame_number = 0;
    uint32_t real_time = 0;
    uint32_t daq_run_time = 0;
    std::vector<uint32_t> channels;
};

static std::vector<uint32_t> read_u32_file(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("failed to open: " + path);
    }

    ifs.seekg(0, std::ios::end);
    std::streamsize size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if (size < 0 || size % 4 != 0) {
        throw std::runtime_error("file size is not a multiple of 4 bytes: " + path);
    }

    std::vector<uint32_t> words(static_cast<size_t>(size) / 4);
    if (!ifs.read(reinterpret_cast<char*>(words.data()), size)) {
        throw std::runtime_error("failed to read: " + path);
    }
    return words;
}

static uint64_t diff_u32(uint32_t after, uint32_t before) {
    if (after >= before) {
        return static_cast<uint64_t>(after - before);
    }
    return static_cast<uint64_t>(after) + (1ULL << 32) - static_cast<uint64_t>(before);
}

static ScrData load_scr(const std::string& path) {
    auto w = read_u32_file(path);

    if (w.size() < 20) {
        throw std::runtime_error("too short: " + path);
    }

    ScrData d;
    d.header_magic = w[0];
    d.n_words = w[1];

    if (w.size() < static_cast<size_t>(2 + d.n_words)) {
        throw std::runtime_error("truncated payload: " + path);
    }
    if (d.n_words < 18) {
        throw std::runtime_error("invalid n_words (<18): " + path);
    }

    d.heartbeat_count        = w[2];
    d.heartbeat_frame_number = w[3];
    d.real_time              = w[4];
    d.daq_run_time           = w[5];

    const size_t n_channels = d.n_words - 18;
    d.channels.resize(n_channels);
    for (size_t ch = 0; ch < n_channels; ++ch) {
        d.channels[ch] = w[20 + ch];
    }

    return d;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr
            << "Usage: " << argv[0]
            << " outbase in0.dat in1.dat [in2.dat ...]\n";
        return 1;
    }

    const std::string outbase = argv[1];
    const std::string out_root = outbase + ".root";
    const std::string out_pdf  = outbase + ".pdf";

    constexpr double HB_TICK_SEC = 524288e-9; // 2^19 ns
    constexpr int INVALID_INT = -1000;
    constexpr double INVALID_DOUBLE = -1000.0;

    try {
        std::vector<ScrData> all_data;
        all_data.reserve(argc - 2);

        for (int i = 2; i < argc; ++i) {
            all_data.push_back(load_scr(argv[i]));
        }

        if (all_data.size() < 2) {
            throw std::runtime_error("need at least 2 input files");
        }

        const int n_files = static_cast<int>(all_data.size());
        const int n_ch = static_cast<int>(all_data[0].channels.size());

        for (int i = 1; i < n_files; ++i) {
            if (static_cast<int>(all_data[i].channels.size()) != n_ch) {
                throw std::runtime_error("channel count mismatch among input files");
            }
        }

        TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            throw std::runtime_error("failed to create ROOT file: " + out_root);
        }

        gStyle->SetOptStat(0);

        TTree* tree = new TTree("tree", "scaler time series");

        int entry_index = 0;
        int frame_number = 0;
        int frame_diff = INVALID_INT;
        double time_sec = 0.0;

        std::vector<int> count(n_ch, 0);
        std::vector<int> count_diff(n_ch, INVALID_INT);
        std::vector<double> rate_hz(n_ch, INVALID_DOUBLE);

        tree->Branch("entry_index", &entry_index, "entry_index/I");
        tree->Branch("frame_number", &frame_number, "frame_number/I");
        tree->Branch("frame_diff", &frame_diff, "frame_diff/I");
        tree->Branch("time_sec", &time_sec, "time_sec/D");
        tree->Branch("count", &count);
        tree->Branch("count_diff", &count_diff);
        tree->Branch("rate_hz", &rate_hz);

        std::vector<TH1D*> h_diff_list;
        std::vector<TH1D*> h_rate_list;
        h_diff_list.reserve(n_files - 1);
        h_rate_list.reserve(n_files - 1);

        for (int i = 1; i < n_files; ++i) {
            auto* hdiff = new TH1D(
                Form("h_count_diff_%d", i),
                Form("Count diff (%d -> %d);Channel;Count diff", i - 1, i),
                n_ch, -0.5, n_ch - 0.5
            );
            auto* hrate = new TH1D(
                Form("h_rate_%d", i),
                Form("Rate (%d -> %d);Channel;Rate [Hz]", i - 1, i),
                n_ch, -0.5, n_ch - 0.5
            );
            h_diff_list.push_back(hdiff);
            h_rate_list.push_back(hrate);
        }

        std::vector<std::vector<double>> graph_x(n_ch);
        std::vector<std::vector<double>> graph_y(n_ch);

        const uint32_t frame0 = all_data[0].heartbeat_frame_number;

        for (int i = 0; i < n_files; ++i) {
            const auto& cur = all_data[i];

            entry_index = i;
            frame_number = static_cast<int>(cur.heartbeat_frame_number);

            uint64_t t_from_start_ticks = diff_u32(cur.heartbeat_frame_number, frame0);
            time_sec = static_cast<double>(t_from_start_ticks) * HB_TICK_SEC;

            for (int ch = 0; ch < n_ch; ++ch) {
                count[ch] = static_cast<int>(cur.channels[ch]);
            }

            if (i == 0) {
                frame_diff = INVALID_INT;
                for (int ch = 0; ch < n_ch; ++ch) {
                    count_diff[ch] = INVALID_INT;
                    rate_hz[ch] = INVALID_DOUBLE;
                }
            } else {
                const auto& prev = all_data[i - 1];
                uint64_t df = diff_u32(cur.heartbeat_frame_number, prev.heartbeat_frame_number);
                if (df > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
                    throw std::runtime_error("frame diff too large for int");
                }
                frame_diff = static_cast<int>(df);

                const double dt_sec = static_cast<double>(df) * HB_TICK_SEC;
                if (dt_sec <= 0.0) {
                    throw std::runtime_error("non-positive delta t");
                }

                TH1D* hdiff = h_diff_list[i - 1];
                TH1D* hrate = h_rate_list[i - 1];

                for (int ch = 0; ch < n_ch; ++ch) {
                    uint64_t dc = diff_u32(cur.channels[ch], prev.channels[ch]);
                    if (dc > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
                        throw std::runtime_error("count diff too large for int");
                    }

                    count_diff[ch] = static_cast<int>(dc);
                    rate_hz[ch] = static_cast<double>(dc) / dt_sec;

                    hdiff->SetBinContent(ch + 1, count_diff[ch]);
                    hrate->SetBinContent(ch + 1, rate_hz[ch]);

                    graph_x[ch].push_back(time_sec);
                    graph_y[ch].push_back(rate_hz[ch]);
                }
            }

            tree->Fill();
        }

        std::vector<TGraph*> graphs;
        graphs.reserve(n_ch);
        for (int ch = 0; ch < n_ch; ++ch) {
            TGraph* g = new TGraph(
                static_cast<int>(graph_x[ch].size()),
                graph_x[ch].data(),
                graph_y[ch].data()
            );
            g->SetName(Form("g_rate_ch%d", ch));
            g->SetTitle(Form("Rate vs time ch%d;Time [s];Rate [Hz]", ch));
            graphs.push_back(g);
        }

        TCanvas* c = new TCanvas("c", "plots", 1200, 800);
        c->Print((out_pdf + "[").c_str());

        for (auto* h : h_diff_list) {
            h->Draw("HIST");
            c->Print(out_pdf.c_str());
        }

        for (auto* h : h_rate_list) {
            h->Draw("HIST");
            c->Print(out_pdf.c_str());
        }

        for (auto* g : graphs) {
            g->Draw("ALP");
            c->Print(out_pdf.c_str());
        }

        c->Print((out_pdf + "]").c_str());

        tree->Write();
        for (auto* h : h_diff_list) h->Write();
        for (auto* h : h_rate_list) h->Write();
        for (auto* g : graphs) g->Write();
        c->Write();

        fout->Close();

        std::cout << "Saved ROOT file : " << out_root << "\n";
        std::cout << "Saved PDF       : " << out_pdf << "\n";
        std::cout << "Entries         : " << n_files << "\n";
        std::cout << "Channels        : " << n_ch << "\n";

    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
