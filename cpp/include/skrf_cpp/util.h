#pragma once

#include <string>
#include <algorithm>

namespace skrf_cpp {

inline std::string trim(const std::string &s) {
    auto a = s.find_first_not_of(" \t\r\n");
    if(a==std::string::npos) return "";
    auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

// choose a human-friendly frequency unit and divisor based on maximum Hz
inline std::pair<double,std::string> freq_scale(double hz) {
    if(hz >= 1e9) return {1e9, "GHZ"};
    if(hz >= 1e6) return {1e6, "MHZ"};
    if(hz >= 1e3) return {1e3, "KHZ"};
    return {1.0, "HZ"};
}

// Touchstone option parsing results
struct TouchstoneOptions {
    enum class DataFormat {RI, MA, DB};
    double freq_multiplier = 1.0;
    DataFormat format = DataFormat::MA;
    double z0 = 50.0;
};

// Parse the option line (contents after '#') of a Touchstone file.
// Recognizes frequency units (HZ/KHZ/MHZ/GHZ), data formats (RI/MA/DB), and reference R (e.g. 'R 50' or 'R=50').
inline TouchstoneOptions parse_touchstone_options(const std::string &optline) {
    TouchstoneOptions out;
    std::istringstream iss(optline);
    std::string token;
    while(iss >> token) {
        std::string up = token;
        for(auto &c: up) c = static_cast<char>(std::toupper((unsigned char)c));
        if(up == "GHZ") out.freq_multiplier = 1e9;
        else if(up == "MHZ") out.freq_multiplier = 1e6;
        else if(up == "KHZ") out.freq_multiplier = 1e3;
        else if(up == "HZ") out.freq_multiplier = 1.0;
        else if(up == "RI") out.format = TouchstoneOptions::DataFormat::RI;
        else if(up == "MA") out.format = TouchstoneOptions::DataFormat::MA;
        else if(up == "DB") out.format = TouchstoneOptions::DataFormat::DB;
        else if(up == "R" || up == "R=" || up == "Z0" || up == "Z0=") {
            // try to read next numeric token or parse after '='
            if(up.find('=') != std::string::npos) {
                auto pos = token.find('=');
                if(pos != std::string::npos && pos + 1 < token.size()) {
                    std::string nv = token.substr(pos+1);
                    try { out.z0 = std::stod(nv); } catch(...) {}
                }
            } else {
                std::string v;
                if(iss >> v) {
                    try { out.z0 = std::stod(v); } catch(...) {}
                }
            }
        } else {
            // also accept forms like 'R=50' without splitting
            auto eq = up.find('R=');
            if(eq != std::string::npos) {
                try { out.z0 = std::stod(up.substr(eq+2)); } catch(...) {}
            }
        }
    }
    return out;
}

}
