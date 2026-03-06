#pragma once

#include <vector>
#include <string>
#include "frequency.h"
#include "constants.h"

namespace skrf_cpp {

inline std::vector<double> linspace(double start, double stop, size_t n) {
    std::vector<double> out;
    if (n == 0) return out;
    if (n == 1) { out.push_back(start); return out; }
    out.resize(n);
    double step = (stop - start) / static_cast<double>(n - 1);
    for (size_t i = 0; i < n; ++i) out[i] = start + step * static_cast<double>(i);
    return out;
}

class StaticInstances {
public:
    static FrequencySeries f_wr10() { return FrequencySeries(linspace(75.0 * FREQ_UNIT_GHZ, 110.0 * FREQ_UNIT_GHZ, 1001)); }
    static FrequencySeries f_wr3()  { return FrequencySeries(linspace(220.0 * FREQ_UNIT_GHZ, 325.0 * FREQ_UNIT_GHZ, 1001)); }
    static FrequencySeries f_wr2p2(){ return FrequencySeries(linspace(330.0 * FREQ_UNIT_GHZ, 500.0 * FREQ_UNIT_GHZ, 1001)); }
    static FrequencySeries f_wr1p5(){ return FrequencySeries(linspace(500.0 * FREQ_UNIT_GHZ, 750.0 * FREQ_UNIT_GHZ, 1001)); }
    static FrequencySeries f_wr1()  { return FrequencySeries(linspace(750.0 * FREQ_UNIT_GHZ, 1100.0 * FREQ_UNIT_GHZ, 1001)); }

    // additional aliases
    static FrequencySeries f_wm1295() { return f_wr5p1(); }

    // Note: RectangularWaveguide and Freespace classes are not yet ported.
};

} // namespace skrf_cpp
