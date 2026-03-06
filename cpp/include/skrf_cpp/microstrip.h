#pragma once

#include <complex>
#include <cmath>
#include "media_base.h"

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;

class Microstrip : public MediaModel {
public:
    // w = conductor width, h = substrate thickness, eps_r = substrate relative permittivity
    double w, h, eps_r;

    Microstrip(double w_, double h_, double eps_r_=4.4) : w(w_), h(h_), eps_r(eps_r_) {}

    double effective_eps() const {
        double wh = w / h;
        double eps_eff = (eps_r + 1.0)/2.0 + (eps_r - 1.0)/2.0 * 1.0/std::sqrt(1.0 + 12.0/wh);
        return eps_eff;
    }

    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        double wh = w / h;
        double eps_eff = effective_eps();
        double Z0;
        if(wh <= 1.0) {
            Z0 = (60.0 / std::sqrt(eps_eff)) * std::log(8.0/wh + 0.25*wh);
        } else {
            Z0 = (120.0*M_PI) / (std::sqrt(eps_eff) * (wh + 1.393 + 0.667*std::log(wh + 1.444)));
        }
        return std::complex<double>(Z0, 0.0);
    }

    std::complex<double> propagation_constant(double freq_hz) const override {
        double eps_eff = effective_eps();
        double omega = 2.0*M_PI*freq_hz;
        double beta = omega * std::sqrt(MU0 * EPS0 * eps_eff);
        return std::complex<double>(0.0, beta);
    }
};

} // namespace media
} // namespace skrf_cpp
