#include <iostream>
#include <exception>
#include "../include/skrf_cpp/touchstone.h"
#include "../include/skrf_cpp/vector_fitting_advanced.h"
#include <fstream>

int main(int argc, char **argv) {
    std::string path;
    if(argc > 1) path = argv[1];
    else path = "../skrf/data/ntwk1.s2p"; // default sample from repository

    try {
        auto net = skrf_cpp::Touchstone::load(path);
        net.printSummary();

        if(argc > 2 && std::string(argv[2]) == "fit") {
            if(net.n_ports < 2) {
                std::cerr << "Network has fewer than 2 ports; cannot fit S21\n";
                return 3;
            }
            skrf_cpp::VectorFittingAdvanced vf;
            skrf_cpp::VectorFittingAdvanced::FitOptions opts;
            opts.npoles = 8;
            // extract frequencies
            std::vector<double> freqs;
            for(const auto &f : net.freqs) freqs.push_back(f.hz);
            // fit S21 (row=1,col=0 assuming 0-based indexing)
            vf.fit_network(net, 1, 0, opts);

            // write fitted vs measured to CSV
            std::ofstream out("vf_fitted.csv");
            out << "freq_hz,measured,fit_real,fit_imag\n";
            for(size_t i=0;i<freqs.size();++i) {
                double fhz = freqs[i];
                auto flat = net.s_at_index(i);
                std::complex<double> Hmeas = flat[1 * net.n_ports + 0];
                auto Hfit = vf.eval(fhz);
                out << fhz << "," << Hmeas.real() << "+" << Hmeas.imag() << "j," << Hfit.real() << "," << Hfit.imag() << "\n";
            }
            out.close();
            std::cout << "Vector fitting complete, results written to vf_fitted.csv\n";
        }

    } catch(const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }

    return 0;
}
