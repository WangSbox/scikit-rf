#pragma once

#include <complex>
#include <cmath>
#include <stdexcept>
#include "media_base.h"

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;

class CoaxialFull : public MediaModel {
public:
    double a; // inner radius (m)
    double b; // outer radius (m)
    double eps_r;
    double mu_r;

    CoaxialFull(double a_, double b_, double eps_r_=1.0, double mu_r_=1.0)
        : a(a_), b(b_), eps_r(eps_r_), mu_r(mu_r_) {
        if(a <= 0 || b <= 0 || b <= a) throw std::runtime_error("CoaxialFull: invalid radii");
    }

    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        double mu = MU0 * mu_r;
        double eps = EPS0 * eps_r;
        double Z0 = (1.0/(2.0*M_PI)) * std::sqrt(mu/eps) * std::log(b / a);
        return std::complex<double>(Z0, 0.0);
    }

    std::complex<double> propagation_constant(double freq_hz) const override {
        double omega = 2.0*M_PI*freq_hz;
        std::complex<double> val = std::complex<double>(0.0, omega * std::sqrt(MU0*mu_r * EPS0*eps_r));
        return val;
    }
};

} // namespace media
} // namespace skrf_cpp
