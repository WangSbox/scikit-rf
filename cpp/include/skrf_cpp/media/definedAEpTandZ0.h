#pragma once

#include <complex>
#include <tuple>
#include <cmath>

namespace skrf_cpp {
namespace media {

// Helper that returns (alpha, epsilon_eff, Z0) for a homogeneous medium
static std::tuple<double, double, std::complex<double>> definedAEpTandZ0(double eps_r, double mu_r=1.0, double freq_hz=1.0) {
    // alpha (attenuation) default 0 for lossless, eps_eff ~= eps_r
    double alpha = 0.0;
    double eps_eff = eps_r;
    double eta = std::sqrt((mu_r*4.0*M_PI*1e-7) / (eps_r*8.854187817e-12));
    std::complex<double> Z0(eta, 0.0);
    return std::make_tuple(alpha, eps_eff, Z0);
}

} // namespace media
} // namespace skrf_cpp
