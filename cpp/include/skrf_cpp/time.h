#pragma once

#include <vector>

namespace skrf_cpp {

class TimeSeries {
public:
    std::vector<double> t;
    TimeSeries() = default;
};

} // namespace skrf_cpp
