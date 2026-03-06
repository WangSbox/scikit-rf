#pragma once

#include <complex>
#include <cmath>
#include <stdexcept>
#include "media_base.h"

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;
static const double C0 = 1.0/std::sqrt(MU0*EPS0);
static const double ETA0 = std::sqrt(MU0/EPS0);

// Rectangular waveguide (dominant TE10 mode) simple model
class RectangularWaveguide : public MediaModel {
public:
    // a = broad dimension (m), b = narrow dimension (m)
    double a, b;
    double eps_r;
    double mu_r;

    RectangularWaveguide(double a_, double b_, double eps_r_=1.0, double mu_r_=1.0)
        : a(a_), b(b_), eps_r(eps_r_), mu_r(mu_r_) {}

    // cutoff frequency for TE_mn: fc_mn = c/(2*sqrt(eps_r))*sqrt((m/a)^2 + (n/b)^2)
    double cutoff_freq(int m, int n) const {
        double term = (m==0?0.0: (m/(2.0*a)))*(m==0?0.0: (m/(2.0*a)));
        double term2 = (n==0?0.0: (n/(2.0*b)))*(n==0?0.0: (n/(2.0*b)));
        double fc = C0/std::sqrt(eps_r) * std::sqrt(term + term2);
        return fc;
    }

    std::complex<double> propagation_constant(double freq_hz, int m=1, int n=0) const override {
        double fc = cutoff_freq(m,n);
        double omega = 2.0*M_PI*freq_hz;
        double k0 = omega * std::sqrt(MU0 * EPS0 * mu_r * eps_r);
        if(freq_hz <= fc) {
            double alpha = std::sqrt((M_PI*0.0) + 1.0) * (fc - freq_hz); // simple positive real attenuation
            return std::complex<double>(alpha, 0.0);
        }
        double beta = k0 * std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(0.0, beta);
    }

    std::complex<double> characteristic_impedance(double freq_hz, int m=1, int n=0) const override {
        double fc = cutoff_freq(m,n);
        double eta = ETA0 * std::sqrt(mu_r/eps_r);
        if(freq_hz <= fc) {
            return std::complex<double>(INFINITY, 0.0);
        }
        double factor = 1.0 / std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        // For TE modes, Z_TE = eta * factor
        return std::complex<double>(eta * factor, 0.0);
    }
};

} // namespace media
} // namespace skrf_cpp
