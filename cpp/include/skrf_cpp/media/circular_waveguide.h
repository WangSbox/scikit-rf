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

// Circular waveguide simplified model (dominant TE11)
class CircularWaveguide : public MediaModel {
public:
    // a = radius (m)
    double a;
    double eps_r;
    double mu_r;

    CircularWaveguide(double a_, double eps_r_ = 1.0, double mu_r_ = 1.0)
        : a(a_), eps_r(eps_r_), mu_r(mu_r_) {
        if(a <= 0) throw std::runtime_error("CircularWaveguide: invalid radius");
    }

    // approximate cutoff frequency for TE11 mode using first root ~1.84118
    double cutoff_TE11() const {
        double uc = 1.841183781; // first root for TE11
        return (uc * C0) / (2.0 * M_PI * a * std::sqrt(eps_r));
    }

    std::complex<double> propagation_constant(double freq_hz) const override {
        double fc = cutoff_TE11();
        double omega = 2.0*M_PI*freq_hz;
        double k0 = omega * std::sqrt(MU0 * EPS0 * mu_r * eps_r);
        if(freq_hz <= fc) {
            // below cutoff: evanescent, real positive attenuation
            double alpha = std::sqrt((fc - freq_hz) * (fc - freq_hz));
            return std::complex<double>(alpha, 0.0);
        }
        double beta = k0 * std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(0.0, beta);
    }

    std::complex<double> characteristic_impedance(double freq_hz) const override {
        double fc = cutoff_TE11();
        double eta = std::sqrt(MU0/EPS0) * std::sqrt(mu_r/eps_r);
        if(freq_hz <= fc) return std::complex<double>(INFINITY, 0.0);
        double factor = 1.0 / std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(eta * factor, 0.0);
    }
};

} // namespace media
} // namespace skrf_cpp
