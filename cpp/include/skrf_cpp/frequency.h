#pragma once

#include <vector>

namespace skrf_cpp {

class FrequencySeries {
public:
    std::vector<double> hz;

    FrequencySeries() = default;
    explicit FrequencySeries(std::vector<double> v): hz(std::move(v)) {}

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
