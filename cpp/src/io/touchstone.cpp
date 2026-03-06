#include "../include/skrf_cpp/touchstone.h"
#include "../include/skrf_cpp/util.h"
#include "../include/skrf_cpp/mathFunctions.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <iomanip>

namespace skrf_cpp {

Network Touchstone::load(const std::string &path) {
    std::ifstream ifs(path);
    if(!ifs) throw std::runtime_error("cannot open file: " + path);

    Network net;
    std::string line;
    bool header_parsed = false;

    TouchstoneOptions opts;
    // scan for option line starting with '#'
    while(std::getline(ifs, line)) {
        std::string t = trim(line);
        if(t.empty()) continue;
        if(t[0] == '!') continue;
        if(t[0] == '#') {
            std::string rest = t.substr(1);
            opts = parse_touchstone_options(rest);
            net.z0 = opts.z0;
            break;
        }
        if(t[0] != '!') break;
    }
    TouchstoneOptions::DataFormat format = opts.format;
    double freq_multiplier = opts.freq_multiplier;

    ifs.clear();
    ifs.seekg(0);

    // deduce port count from first numeric line
    while(std::getline(ifs, line)) {
        line = trim(line);
        if(line.empty()) continue;
        if(line[0] == '!' || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<double> vals;
        double v;
        while(iss >> v) vals.push_back(v);
        if(vals.empty()) continue;

        int numD = static_cast<int>(vals.size()) - 1;
        if(numD <= 0) continue;
        int nParamPairs = numD/2;
        int n_ports = static_cast<int>(std::round(std::sqrt((double)nParamPairs)));
        if(n_ports*n_ports != nParamPairs) {
            throw std::runtime_error("unsupported touchstone data layout (non-square S matrix)");
        }

        net.n_ports = n_ports;

        // parse all numeric lines
        ifs.clear();
        ifs.seekg(0);
        while(std::getline(ifs, line)) {
            line = trim(line);
            if(line.empty()) continue;
            if(line[0] == '!' || line[0] == '#') continue;

            std::istringstream iss2(line);
            std::vector<double> row;
            double w;
            while(iss2 >> w) row.push_back(w);
            if(row.empty()) continue;

            double freq = row[0] * freq_multiplier;
            net.freqs.emplace_back(freq);

            int expected = 1 + 2 * n_ports * n_ports;
            if(static_cast<int>(row.size()) < expected) {
                throw std::runtime_error("insufficient columns in touchstone data line");
            }

            std::vector<std::complex<double>> srow(n_ports*n_ports);
            size_t idx = 1;
            for(int i=0;i<n_ports*n_ports;i++) {
                double a = row[idx++];
                double b = row[idx++];
                std::complex<double> val;
                if(format == DataFormat::MA) {
                    double mag = a;
                    double phase_deg = b;
                    double phase = phase_deg * M_PI / 180.0;
                    val = std::complex<double>(mag * std::cos(phase), mag * std::sin(phase));
                } else if(format == DataFormat::DB) {
                    double mag = std::pow(10.0, a/20.0);
                    double phase_deg = b;
                    double phase = phase_deg * M_PI / 180.0;
                    val = std::complex<double>(mag * std::cos(phase), mag * std::sin(phase));
                } else { // RI
                    val = std::complex<double>(a, b);
                }
                srow[i] = val;
            }
            net.sparams.push_back(std::move(srow));
        }

        header_parsed = true;
        break;
    }

    if(!header_parsed) throw std::runtime_error("no data parsed from touchstone file");
    return net;
}

// write network to touchstone file in MA/DB/RI formats (default MA)
void Touchstone::save(const std::string &path, const Network &net, const std::string &format) {
    if(net.freqs.empty() || net.sparams.empty()) throw std::runtime_error("empty network cannot be written");
    int n_ports = net.n_ports;
    if(n_ports <= 0) throw std::runtime_error("invalid port count");

    // decide frequency unit for human-readable file: use GHz if max >=1e9, else MHz if >=1e6, else Hz
    double maxf = 0.0;
    for(const auto &f : net.freqs) if(f.hz > maxf) maxf = f.hz;
    double divisor = 1.0; std::string unit = "HZ";
    if(maxf >= 1e9) { divisor = 1e9; unit = "GHZ"; }
    else if(maxf >= 1e6) { divisor = 1e6; unit = "MHZ"; }
    else if(maxf >= 1e3) { divisor = 1e3; unit = "KHZ"; }

    std::ofstream ofs(path);
    if(!ofs) throw std::runtime_error("cannot open output file: " + path);

    // header line (include reference Z0 if available)
    ofs << "# " << unit << " S " << format << " R " << net.z0 << "\n";
    ofs << std::fixed << std::setprecision(6);

    // write each frequency line
    for(size_t i=0;i<net.freqs.size();++i) {
        double fval = net.freqs[i].hz / divisor;
        ofs << fval;
        const auto &flat = net.sparams[i];
        if(static_cast<int>(flat.size()) != n_ports * n_ports) throw std::runtime_error("sparam length mismatch on write");
        for(int k=0;k<n_ports*n_ports;++k) {
            std::complex<double> z = flat[k];
            if(format == "RI") {
                ofs << " " << std::setprecision(6) << std::real(z) << " " << std::imag(z);
            } else if(format == "DB") {
                double mag_db = 20.0 * std::log10(std::abs(z));
                double phase_deg = std::arg(z) * 180.0 / M_PI;
                ofs << " " << std::setprecision(6) << mag_db << " " << phase_deg;
            } else { // MA
                double mag = std::abs(z);
                double phase_deg = std::arg(z) * 180.0 / M_PI;
                ofs << " " << std::setprecision(6) << mag << " " << phase_deg;
            }
        }
        ofs << "\n";
    }
    ofs.close();
}

}
