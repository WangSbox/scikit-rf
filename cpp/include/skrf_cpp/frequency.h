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
