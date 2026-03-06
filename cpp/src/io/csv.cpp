#include "../../include/skrf_cpp/io/csv.h"
#include "../../include/skrf_cpp/frequency.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace skrf_cpp {
namespace io {

static std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if(a==std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

CSVReadResult read_pna_csv(const std::string &path) {
    std::ifstream ifs(path);
    if(!ifs) throw std::runtime_error("cannot open csv: " + path);
    CSVReadResult out;
    std::string line;
    int begin_line = -2;
    int end_line = -1;
    int n_END = 0;
    std::string comments;
    int k = 0;
    std::vector<std::string> lines;
    while(std::getline(ifs, line)) {
        lines.push_back(line);
    }
    for(size_t i=0;i<lines.size();++i) {
        const auto &ln = lines[i];
        if(ln.size()>0 && ln[0]=='!') comments += ln.substr(1) + "\n";
        else if(ln.rfind("BEGIN",0)==0 && n_END==0) begin_line = static_cast<int>(i);
        else if(ln.rfind("END",0)==0) {
            if(n_END==0) end_line = static_cast<int>(i);
            ++n_END;
        }
        if(static_cast<int>(i) == begin_line+1) out.header = ln;
    }
    out.comments = comments;
    int footer = static_cast<int>(lines.size()) - end_line;
    // parse numeric block after begin_line+2 until before end_line
    int start = std::max(0, begin_line+2);
    int stop = (end_line>=0) ? end_line : static_cast<int>(lines.size());
    for(int i=start;i<stop;++i) {
        std::istringstream iss(lines[i]);
        std::vector<double> row;
        double v;
        while(iss >> v) row.push_back(v);
        if(!row.empty()) out.data.push_back(std::move(row));
    }
    return out;
}

Network pna_csv_to_network(const std::string &path) {
    CSVReadResult r = read_pna_csv(path);
    if(r.data.empty()) throw std::runtime_error("no numeric data in csv");
    // try to detect DB/deg style from header string
    std::string header = r.header;
    std::string lower = header;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c){ return std::tolower(c); });
    bool is_db_deg = (lower.find("db")!=std::string::npos && lower.find("deg")!=std::string::npos);
    // assume 2-port layout if possible: find columns corresponding to S11,S21,S12,S22 in header or fallback
    // simplistic strategy: if 7 or more columns, interpret as freq + pairs for S11,S21,S12,S22
    int ncols = static_cast<int>(r.data.front().size());
    if(ncols < 3) throw std::runtime_error("csv: insufficient columns");
    // frequency is col 0
    std::vector<double> freqs; freqs.reserve(r.data.size());
    for(const auto &row : r.data) freqs.push_back(row[0]);

    // build S array shape (n_points,2,2)
    size_t npts = r.data.size();
    Network net;
    net.n_ports = 2;
    net.freqs.clear();
    for(double f : freqs) net.freqs.emplace_back(f);
    net.sparams.reserve(npts);

    for(const auto &row : r.data) {
        if(static_cast<int>(row.size()) < 5) throw std::runtime_error("csv row too short for S-params");
        std::vector<std::complex<double>> flat(4);
        if(is_db_deg) {
            // heuristic: assume next columns are S11_mag_dB, S11_phase_deg, S21_mag_dB, S21_phase_deg, S12..., S22...
            // find groups by scanning header tokens; fallback to fixed ordering
            double s11_db = row.size()>1 ? row[1]:0.0;
            double s11_ph = row.size()>2 ? row[2]:0.0;
            double s21_db = row.size()>3 ? row[3]:0.0;
            double s21_ph = row.size()>4 ? row[4]:0.0;
            // best-effort mapping
            auto dbdeg_to_c = [](double db, double deg){ double mag = std::pow(10.0, db/20.0); double ph = deg * M_PI / 180.0; return std::complex<double>(mag*cos(ph), mag*sin(ph)); };
            flat[0] = dbdeg_to_c(s11_db, s11_ph);
            flat[2] = dbdeg_to_c(s21_db, s21_ph);
            // attempt remaining
            if(row.size() >= 7) flat[1] = dbdeg_to_c(row[5], row[6]); else flat[1] = std::complex<double>(0,0);
            if(row.size() >= 9) flat[3] = dbdeg_to_c(row[7], row[8]); else flat[3] = std::complex<double>(0,0);
        } else {
            // assume columns are freq, re1, im1, re2, im2, ... assign first two pairs
            flat[0] = std::complex<double>(row[1], row.size()>2 ? row[2] : 0.0);
            flat[1] = std::complex<double>(row.size()>3 ? row[3] : 0.0, row.size()>4 ? row[4] : 0.0);
            flat[2] = std::complex<double>(row.size()>5 ? row[5] : 0.0, row.size()>6 ? row[6] : 0.0);
            flat[3] = std::complex<double>(row.size()>7 ? row[7] : 0.0, row.size()>8 ? row[8] : 0.0);
        }
        net.sparams.push_back(std::move(flat));
    }
    return net;
}

} // namespace io
} // namespace skrf_cpp
