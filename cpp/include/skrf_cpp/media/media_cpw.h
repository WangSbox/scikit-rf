#pragma once

#include <complex>
#include <cmath>
#include <stdexcept>
#include "media_base.h"

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;

// Complete elliptic integral of the first kind via AGM
static double ellipticK(double k) {
    if(k < 0.0 || k >= 1.0) {
        if(k == 0.0) return M_PI_2;
        throw std::runtime_error("ellipticK: modulus out of range");
    }
    double a = 1.0;
    double b = std::sqrt(1.0 - k*k);
    const double eps = 1e-12;
    while(std::abs(a - b) > eps) {
        double an = 0.5*(a + b);
        b = std::sqrt(a*b);
        a = an;
    }
    return M_PI/(2.0*a);
}

// Coplanar waveguide (CPW) quasi-static model (approximate)
class CPW : public MediaModel {
public:
    // w = center conductor width, s = gap, t = thickness, h = substrate height
    double w, s, t, h;
    double eps_r;

    CPW(double w_, double s_, double t_, double h_, double eps_r_=1.0)
        : w(w_), s(s_), t(t_), h(h_), eps_r(eps_r_) {}

    std::complex<double> characteristic_impedance(double freq_hz) const override {
        // quasi-static effective permittivity approximation
        double eps_eff = (eps_r + 1.0)/2.0;
        // modulus for conformal mapping
        double k = w/(w + 2.0*s);
        if(k <= 0.0 || k >= 1.0) throw std::runtime_error("CPW: invalid geometry k");
        double Kk = ellipticK(k);
        double Kkprime = ellipticK(std::sqrt(1.0 - k*k));
        double Z0 = (30.0*M_PI)/std::sqrt(eps_eff) * (Kkprime / Kk);
        return std::complex<double>(Z0, 0.0);
    }

    std::complex<double> propagation_constant(double freq_hz) const override {
        double eps_eff = (eps_r + 1.0)/2.0;
        double omega = 2.0*M_PI*freq_hz;
        double beta = omega * std::sqrt(MU0 * EPS0 * eps_eff);
        return std::complex<double>(0.0, beta);
    }
};

} // namespace media
} // namespace skrf_cpp
