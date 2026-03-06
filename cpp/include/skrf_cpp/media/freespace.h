#pragma once

#include <complex>
#include <cmath>
#include "media_base.h"

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;
static const double C0 = 1.0/std::sqrt(MU0*EPS0);
static const double ETA0 = std::sqrt(MU0/EPS0);

class FreeSpace : public MediaModel {
public:
    // characteristic impedance in free space (approx 376.730313668)
    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        return std::complex<double>(ETA0, 0.0);
    }

    // propagation constant gamma = j*omega*sqrt(mu*eps)
    std::complex<double> propagation_constant(double freq_hz, double eps_r = 1.0, double mu_r = 1.0) const override {
        double omega = 2.0*M_PI*freq_hz;
        std::complex<double> val = std::complex<double>(0.0, omega * std::sqrt(MU0*mu_r * EPS0*eps_r));
        return val;
    }
};

} // namespace media
} // namespace skrf_cpp
