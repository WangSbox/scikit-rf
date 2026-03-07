#pragma once

#include <vector>

namespace skrf_cpp {

class FrequencySeries {
public:
    std::vector<double> hz;

    FrequencySeries() = default;
    explicit FrequencySeries(std::vector<double> v): hz(std::move(v)) {}
    // construct from values with unit string (hz, khz, mhz, ghz)
    FrequencySeries(std::vector<double> v, const std::string &unit) {
        double mul = unit_multiplier(unit);
        hz.reserve(v.size());
        for(double x : v) hz.push_back(x * mul);
    }

    size_t size() const { return hz.size(); }
    double operator[](size_t i) const { return hz[i]; }
    // find index of exact match within tol, else return npos
    size_t index_of(double f, double tol = 1e-9) const {
        for(size_t i=0;i<hz.size();++i) if(std::abs(hz[i] - f) <= tol) return i;
        return static_cast<size_t>(-1);
    }
    // return index of closest frequency
    size_t closest_index(double f) const {
        if(hz.empty()) return static_cast<size_t>(-1);
        size_t best = 0; double bestd = std::abs(hz[0]-f);
        for(size_t i=1;i<hz.size();++i) { double d = std::abs(hz[i]-f); if(d < bestd) { bestd = d; best = i; } }
        return best;
    }

    // check monotonic increasing (strict if strict=true)
    bool is_monotonic(bool strict = true) const {
        for(size_t i=1;i<hz.size();++i) {
            if(strict) { if(!(hz[i] > hz[i-1])) return false; }
            else { if(!(hz[i] >= hz[i-1])) return false; }
        }
        return true;
    }

    // helper: map unit string to multiplier
    static double unit_multiplier(const std::string &unit) {
        std::string u = unit;
        for(auto &c : u) c = static_cast<char>(std::tolower((unsigned char)c));
        if(u == "hz") return 1.0;
        if(u == "khz") return 1e3;
        if(u == "mhz") return 1e6;
        if(u == "ghz") return 1e9;
        throw std::runtime_error(std::string("unknown frequency unit: ") + unit);
    }
};

}
#pragma once

namespace skrf_cpp {
struct Frequency {
    double hz{0.0};
    Frequency() = default;
    explicit Frequency(double f): hz(f) {}
};
}
