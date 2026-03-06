#include "../../include/skrf_cpp/io/citi.h"
#include "../../include/skrf_cpp/frequency.h"
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

std::vector<Network> read_citi(const std::string &path) {
    std::ifstream ifs(path);
    if(!ifs) throw std::runtime_error("cannot open cti: " + path);
    std::vector<std::string> lines;
    std::string line;
    while(std::getline(ifs, line)) lines.push_back(line);

    // Very small parser: collect VAR blocks and DATA blocks
    std::map<std::string, std::vector<double>> var_values;
    std::map<std::string, std::pair<std::string, std::vector<std::complex<double>>>> data_values;
    std::vector<Network> nets;

    size_t idx = 0;
    while(idx < lines.size()) {
        std::string ln = trim(lines[idx]);
        if(ln.empty()) { ++idx; continue; }
        std::string up = ln; std::transform(up.begin(), up.end(), up.begin(), [](unsigned char c){ return std::toupper(c); });
        if(up.rfind("VAR_LIST_BEGIN",0)==0) {
            // read next param's values - assume previous VAR line set name
            // Find last VAR before this
            size_t p = idx;
            while(p>0 && trim(lines[p]).rfind("VAR",0)!=0) --p;
            if(trim(lines[p]).rfind("VAR",0)==0) {
                std::istringstream iss(trim(lines[p])); std::string tag, name, formt; int occ;
                iss >> tag >> name >> formt >> occ;
                std::vector<double> vals;
                for(int k=0;k<occ && idx+1+k < lines.size(); ++k) {
                    std::string vln = trim(lines[idx+1+k]);
                    try { vals.push_back(std::stod(vln)); } catch(...) { vals.push_back(0.0); }
                }
                var_values[name] = std::move(vals);
            }
        }
        if(up.rfind("BEGIN",0)==0) {
            // read the data block until END
            std::vector<std::string> block;
            ++idx;
            while(idx < lines.size()) {
                std::string t = trim(lines[idx]);
                if(std::string(t).size()>0) {
                    std::string tu = t; std::transform(tu.begin(), tu.end(), tu.begin(), [](unsigned char c){ return std::toupper(c); });
                    if(tu.rfind("END",0)==0) break;
                }
                block.push_back(lines[idx]); ++idx;
            }
            // simplistic conversion: assume block rows are numeric comma separated
            std::vector<std::vector<double>> nums;
            for(auto &b : block) {
                std::istringstream iss(b);
                std::vector<double> row; double v; char c;
                while(iss >> v) {
                    row.push_back(v);
                    if(!(iss.peek() == ',' || std::isspace(iss.peek()))) break;
                    if(iss.peek()==',') iss.get(c);
                }
                if(!row.empty()) nums.push_back(row);
            }
            if(nums.empty()) { ++idx; continue; }
            // assume first column freq
            std::vector<double> freqs; freqs.reserve(nums.size());
            for(auto &r : nums) freqs.push_back(r[0]);
            // build network assuming 2-port simplified layout
            Network net; net.n_ports = 2; net.z0 = 50.0; net.freqs.clear(); net.sparams.clear();
            for(double f : freqs) net.freqs.emplace_back(f);
            for(auto &r : nums) {
                std::vector<std::complex<double>> flat(4);
                if(r.size()>=5) {
                    flat[0] = std::complex<double>(r[1], r.size()>2 ? r[2]:0.0);
                    flat[1] = std::complex<double>(r.size()>3 ? r[3]:0.0, r.size()>4 ? r[4]:0.0);
                    // remaining as zero or reuse
                    flat[2] = std::complex<double>(0,0);
                    flat[3] = std::complex<double>(0,0);
                } else {
                    flat[0] = std::complex<double>(0,0); flat[1]=flat[2]=flat[3]=std::complex<double>(0,0);
                }
                net.sparams.push_back(flat);
            }
            nets.push_back(std::move(net));
        }
        ++idx;
    }

    return nets;
}

} // namespace io
} // namespace skrf_cpp
