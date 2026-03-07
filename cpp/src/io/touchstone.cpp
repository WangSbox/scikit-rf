#include "../include/skrf_cpp/io/touchstone.h"
#include "../include/skrf_cpp/util.h"
#include "../include/skrf_cpp/mathFunctions.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <iomanip>

namespace skrf_cpp {

Network Touchstone::load(const std::string &path, bool lenient, std::vector<std::string> *warnings) {
    std::ifstream ifs(path);
    if(!ifs) throw std::runtime_error("cannot open file: " + path);

    Network net;
    std::string line;
    bool header_parsed = false;

    TouchstoneOptions opts;
    // read logical lines (support continuation with '\\' and inline '!' comments)
    int phys_line_no = 0;
    auto read_logical = [&](std::ifstream &ifs_in, std::string &out_line) -> bool {
        out_line.clear();
        std::string phy;
        bool any = false;
        while(std::getline(ifs_in, phy)) {
            ++phys_line_no;
            any = true;
            // strip inline comment marker '!'
            auto pc = phy.find('!');
            if(pc != std::string::npos) phy = phy.substr(0, pc);
            // trim right
            while(!phy.empty() && (phy.back()=='\r' || phy.back()=='\n' || phy.back()==' ' || phy.back()=='\t')) phy.pop_back();
            bool cont = false;
            if(!phy.empty() && phy.back() == '\\') { cont = true; phy.pop_back(); }
            if(out_line.empty()) out_line = phy; else out_line += " "; out_line += phy;
            if(!cont) break;
        }
        return any && !out_line.empty();
    };

    // first pass: parse header options
    ifs.clear(); ifs.seekg(0);
    while(read_logical(ifs, line)) {
        std::string t = trim(line);
        if(t.empty()) continue;
        if(t[0] == '#') {
            std::string rest = t.substr(1);
            opts = parse_touchstone_options(rest);
            net.z0 = opts.z0;
            break;
        }
    }
    TouchstoneOptions::DataFormat format = opts.format;
    double freq_multiplier = opts.freq_multiplier;

    // second pass: deduce port count and collect numeric rows using logical lines
    ifs.clear(); ifs.seekg(0);
    // reset phys_line_no for second pass
    phys_line_no = 0;
    // collect HFSS special comment blocks if present
    std::vector<std::vector<double>> hfss_gamma_vals;
    std::vector<std::vector<double>> hfss_impedance_vals;
    while(read_logical(ifs, line)) {
        std::string t = trim(line);
        if(t.empty()) continue;
        if(t[0] == '#' ) continue;
        // skip comment-only lines
        if(t[0] == '!') continue;

        // handle bracketed v2-style sections minimally: record and skip
        if(!t.empty() && t.front() == '[') {
            // basic support: allow files with [Version] / [Network Data] sections
            // we don't fully implement v2 parsing here but accept and continue
            continue;
        }

        // remove inline comments again and parse numbers
        auto posc = line.find('!');
        std::string commentpart = (posc==std::string::npos) ? std::string() : line.substr(posc);
        std::string dataline = (posc==std::string::npos) ? line : line.substr(0,posc);
        std::istringstream iss(dataline);
        std::vector<double> vals; double v;
        while(iss >> v) vals.push_back(v);
        if(vals.empty()) continue;

        int numD = static_cast<int>(vals.size()) - 1;
        if(numD <= 0) continue;
        int nParamPairs = numD/2;
        int n_ports = static_cast<int>(std::round(std::sqrt((double)nParamPairs)));
        if(n_ports * n_ports != nParamPairs) continue; // look for a proper numeric line

        net.n_ports = n_ports;

        // now collect numeric rows from start using logical reader
        ifs.clear(); ifs.seekg(0);
        while(read_logical(ifs, line)) {
            std::string tl = trim(line);
            if(tl.empty()) continue;
            if(tl[0] == '#') continue;
            // capture HFSS per-frequency comments like '! gamma' or '! port impedance'
            auto pc = tl.find('!');
            std::string comment = (pc != std::string::npos) ? tl.substr(pc) : std::string();
            if(!comment.empty()) {
                std::string low = comment;
                for(auto &c : low) c = static_cast<char>(std::tolower((unsigned char)c));
                if(low.rfind("! gamma", 0) == 0) {
                    // parse following numbers in comment part
                    std::istringstream issc(comment.substr(6)); double vv; std::vector<double> rowv;
                    while(issc >> vv) rowv.push_back(vv);
                    if(!rowv.empty()) hfss_gamma_vals.push_back(std::move(rowv));
                } else if(low.find("! port impedance") != std::string::npos) {
                    // parse following numbers
                    std::size_t pos = low.find("! port impedance");
                    std::istringstream issc(comment.substr(pos + std::string("! port impedance").size())); double vv; std::vector<double> rowv;
                    while(issc >> vv) rowv.push_back(vv);
                    if(!rowv.empty()) hfss_impedance_vals.push_back(std::move(rowv));
                }
            }
            if(pc != std::string::npos) tl = tl.substr(0, pc);
            std::istringstream iss2(tl);
            std::vector<double> row; double w;
            while(iss2 >> w) row.push_back(w);
            if(row.empty()) continue;

            double freq = row[0] * freq_multiplier;
            if(!net.freqs.empty() && std::abs(net.freqs.back().hz - freq) < 1e-12) continue;
            net.freqs.emplace_back(freq);

            int expected = 1 + 2 * n_ports * n_ports;
            if(static_cast<int>(row.size()) < expected) {
                if(lenient && warnings) {
                    std::ostringstream ss; ss << "skipping malformed row (expected " << expected << " cols) at freq " << row[0] << " (approx file line " << phys_line_no << ")";
                    warnings->push_back(ss.str());
                    continue;
                } else if(lenient) {
                    continue;
                } else {
                    std::ostringstream ss; ss << "insufficient columns in touchstone data line at approx file line " << phys_line_no;
                    throw std::runtime_error(ss.str());
                }
            }

            std::vector<std::complex<double>> srow(n_ports*n_ports);
            size_t idx = 1;
            for(int i=0;i<n_ports*n_ports;i++) {
                double a = row[idx++];
                double b = row[idx++];
                std::complex<double> val;
                if(format == DataFormat::MA) {
                    double mag = a; double phase_deg = b; double phase = phase_deg * M_PI / 180.0;
                    val = std::complex<double>(mag * std::cos(phase), mag * std::sin(phase));
                } else if(format == DataFormat::DB) {
                    double mag = std::pow(10.0, a/20.0); double phase_deg = b; double phase = phase_deg * M_PI / 180.0;
                    val = std::complex<double>(mag * std::cos(phase), mag * std::sin(phase));
                } else { val = std::complex<double>(a, b); }
                srow[i] = val;
            }
            net.sparams.push_back(std::move(srow));
        }

        header_parsed = true;
        break;
    }

    if(!header_parsed) {
        if(lenient && warnings) warnings->push_back(std::string("no data parsed from touchstone file: ") + path);
        if(lenient) return net;
        throw std::runtime_error(std::string("no data parsed from touchstone file: ") + path);
    }
    // If HFSS impedance comments were present, set net.z0 to first parsed value (best-effort)
    if(!hfss_impedance_vals.empty()) {
        const auto &v = hfss_impedance_vals.front();
        if(!v.empty()) {
            // take first numeric as reference resistance (real part)
            net.z0 = v.front();
            if(warnings) {
                std::ostringstream ss; ss << "parsed HFSS port impedance block: setting network z0 to " << net.z0;
                warnings->push_back(ss.str());
            }
        }
    }
    // If parsed per-frequency impedance/gamma blocks align with frequency count, attach them to network
    if(!hfss_impedance_vals.empty() && hfss_impedance_vals.size() == net.freqs.size()) {
        net.per_freq_z0.clear(); net.per_freq_z0.reserve(net.freqs.size());
        for(size_t i=0;i<net.freqs.size();++i) {
            const auto &row = hfss_impedance_vals[i];
            if(!row.empty()) net.per_freq_z0.push_back(row.front()); else net.per_freq_z0.push_back(net.z0);
        }
        if(warnings) warnings->push_back(std::string("filled per-frequency impedance vector from HFSS comments"));
    }
    if(!hfss_gamma_vals.empty() && hfss_gamma_vals.size() == net.freqs.size()) {
        net.per_freq_gamma.clear(); net.per_freq_gamma.reserve(net.freqs.size());
        for(size_t i=0;i<net.freqs.size();++i) {
            const auto &row = hfss_gamma_vals[i];
            if(row.size() >= 2) net.per_freq_gamma.emplace_back(row[0], row[1]);
            else if(row.size() == 1) net.per_freq_gamma.emplace_back(row[0], 0.0);
            else net.per_freq_gamma.emplace_back(0.0, 0.0);
        }
        if(warnings) warnings->push_back(std::string("filled per-frequency gamma vector from HFSS comments"));
    }
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
