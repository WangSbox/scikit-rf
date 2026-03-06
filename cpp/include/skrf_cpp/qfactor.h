#pragma once

#include "network.h"

namespace skrf_cpp {

class Qfactor {
public:
    Network ntwk;
    double Q_L{0.0};
    double f_L{0.0};
    bool fitted{false};

    Qfactor(const Network &n): ntwk(n) {}
    // placeholder fit method
    void fit(const std::string &method = "NLQFIT6") { (void)method; fitted = true; }
};

} // namespace skrf_cpp
