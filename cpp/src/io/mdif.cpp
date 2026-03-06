#include "../../include/skrf_cpp/io/mdif.h"
#include "../../include/skrf_cpp/frequency.h"
#include "../../include/skrf_cpp/mathFunctions.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace skrf_cpp {
namespace io {

static std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if(a==std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

// Very small helper to split whitespace
static std::vector<std::string> split_ws(const std::string &s) {
    std::istringstream iss(s);
    std::vector<std::string> out; std::string tok;
    while(iss >> tok) out.push_back(tok);
    return out;
}

std::vector<Network> read_mdif(const std::string &path) {
    std::ifstream ifs(path);
    if(!ifs) throw std::runtime_error("cannot open mdif: " + path);
    std::string line;
    std::vector<std::string> block;
    std::vector<Network> nets;
    // read all lines and group between BEGIN...END
    bool in_block = false;
    std::vector<std::string> header_comments;
    ifs.clear(); ifs.seekg(0);
    while(std::getline(ifs, line)) {
        std::string t = trim(line);
        if(t.empty()) continue;
        if(!in_block) {
            if(t.size()>0 && t[0]=='!') header_comments.push_back(t.substr(1));
        }
        std::string lo = t;
        std::transform(lo.begin(), lo.end(), lo.begin(), [](unsigned char c){ return std::tolower(c); });
        if(lo.rfind("begin",0)==0) { in_block = true; block.clear(); continue; }
        if(lo.rfind("end",0)==0) {
            if(in_block) {
                // parse block into Network (best-effort)
                // find option line (# ...) and kinds (%) and data lines
                std::string opt_line = "";
                std::vector<std::string> kinds_lines;
                std::vector<std::vector<double>> data_lines;
                for(const auto &ln : block) {
                    if(ln.size()>0 && ln[0]=='#') opt_line = ln;
                    else if(ln.size()>0 && ln[0]=='%') kinds_lines.push_back(ln.substr(1));
                    else if(ln.size()>0 && ln[0] != '!') {
                        // numeric
                        std::istringstream iss(ln);
                        std::vector<double> row; double v;
                        while(iss >> v) row.push_back(v);
                        if(!row.empty()) data_lines.push_back(std::move(row));
                    }
                }
                if(data_lines.empty()) { in_block = false; continue; }
                // flatten kinds
                std::vector<std::string> kinds;
                for(auto &kln : kinds_lines) {
                    auto parts = split_ws(kln);
                    kinds.insert(kinds.end(), parts.begin(), parts.end());
                }
                // detect format: default RI
                std::string fmt = "ri";
                std::string unit = "hz";
                double z0 = 50.0;
                if(!opt_line.empty()) {
                    auto toks = split_ws(opt_line.substr(1));
                    if(toks.size()>=3) fmt = toks[2];
                    if(toks.size()>=1) unit = toks[0];
                    if(toks.size()>=5) {
                        try { z0 = std::stod(toks[4]); } catch(...) { z0 = 50.0; }
                    }
                }
                // group data lines per frequency: assume each line already corresponds to one freq
                size_t npts = data_lines.size();
                // build arrays of doubles
                // convert to complex values per pair
                // determine number of complex columns per row (excluding freq)
                size_t ncols = data_lines.front().size();
                if(ncols < 3) { in_block = false; continue; }
                // number of complex values = (ncols-1)/2
                size_t nvals = (ncols - 1) / 2;
                // build matrix of complex values [npts x nvals]
                std::vector<std::vector<std::complex<double>>> values(npts, std::vector<std::complex<double>>(nvals));
                std::vector<double> freqs; freqs.reserve(npts);
                for(size_t ri=0; ri<npts; ++ri) {
                    const auto &r = data_lines[ri];
                    freqs.push_back(r[0]);
                    for(size_t ci=0; ci<nvals; ++ci) {
                        double a = r[1 + 2*ci];
                        double b = r[1 + 2*ci + 1];
                        if(fmt == "ri") values[ri][ci] = std::complex<double>(a,b);
                        else if(fmt == "ma") {
                            double mag = a; double ph = b * M_PI / 180.0;
                            values[ri][ci] = std::complex<double>(mag * std::cos(ph), mag * std::sin(ph));
                        } else if(fmt == "db") {
                            double mag = std::pow(10.0, a/20.0); double ph = b * M_PI / 180.0;
                            values[ri][ci] = std::complex<double>(mag * std::cos(ph), mag * std::sin(ph));
                        } else values[ri][ci] = std::complex<double>(a,b);
                    }
                }
                // try to detect S-parameter layout: if kinds contain s[1,1]
                bool has_s_bracket = false;
                for(auto &k : kinds) {
                    std::string lk = k; std::transform(lk.begin(), lk.end(), lk.begin(), [](unsigned char c){ return std::tolower(c); });
                    if(lk.find("s[")!=std::string::npos) { has_s_bracket = true; break; }
                }
                Network net;
                if(has_s_bracket) {
                    // deduce rank
                    int rank = static_cast<int>(std::round(std::sqrt((double)std::count_if(kinds.begin(), kinds.end(), [](const std::string &s){ return s.find('s')!=std::string::npos; }))));
                    net.n_ports = rank;
                    net.z0 = z0;
                    net.freqs.clear(); net.sparams.clear();
                    for(size_t ri=0; ri<npts; ++ri) {
                        net.freqs.emplace_back(freqs[ri]);
                        std::vector<std::complex<double>> flat(rank*rank);
                        // naive mapping: take first rank*rank complex values
                        for(int idx=0; idx<rank*rank && idx < static_cast<int>(nvals); ++idx) flat[idx] = values[ri][idx];
                        net.sparams.push_back(flat);
                    }
                } else {
                    // fallback: attempt Z or Y similarly
                    net.n_ports = 1;
                    net.z0 = z0;
                    net.freqs.clear(); net.sparams.clear();
                    for(size_t ri=0; ri<npts; ++ri) {
                        net.freqs.emplace_back(freqs[ri]);
                        std::vector<std::complex<double>> flat(1);
                        flat[0] = values[ri][0];
                        net.sparams.push_back(flat);
                    }
                }
                nets.push_back(std::move(net));
            }
            in_block = false; continue;
        }
        if(in_block) block.push_back(line);
    }

    return nets;
}

} // namespace io
} // namespace skrf_cpp
