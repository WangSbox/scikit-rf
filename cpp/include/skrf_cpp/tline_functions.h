#pragma once

#include <Eigen/Dense>
#include <complex>
#include <cmath>

namespace skrf_cpp {

using Matrix2cd = Eigen::Matrix2cd;

// Transmission-line ABCD calculators
// Lossy/general: gamma = alpha + j*beta (per meter), Z0 (ohm)
inline Matrix2cd abcd_of_transmission_line_complex(std::complex<double> gamma, double length_m, std::complex<double> Z0) {
    std::complex<double> gl = gamma * length_m;
    std::complex<double> c = std::cosh(gl);
    std::complex<double> s = std::sinh(gl);
    Matrix2cd M;
    M(0,0) = c; M(0,1) = Z0 * s;
    M(1,0) = s / Z0; M(1,1) = c;
    return M;
}

// Frequency-dependent lossless transmission line ABCD: uses effective permittivity eps_eff
// freq_hz: frequency, length_m: length, Z0: characteristic impedance, eps_eff: effective dielectric constant
inline Matrix2cd abcd_of_transmission_line(double freq_hz, double length_m, double Z0, double eps_eff = 1.0) {
    const double c0 = 299792458.0;
    double beta = 2.0 * M_PI * freq_hz / c0 * std::sqrt(eps_eff);
    std::complex<double> gamma(0.0, beta);
    return abcd_of_transmission_line_complex(gamma, length_m, std::complex<double>(Z0, 0.0));
}

// Frequency-dependent transmission line allowing user-supplied loss function alpha(f)
inline Matrix2cd abcd_of_transmission_line_with_alpha(double freq_hz, double length_m, double Z0, std::function<double(double)> alpha_func, double eps_eff = 1.0) {
    const double c0 = 299792458.0;
    double beta = 2.0 * M_PI * freq_hz / c0 * std::sqrt(eps_eff);
    double alpha = alpha_func ? alpha_func(freq_hz) : 0.0;
    std::complex<double> gamma(alpha, beta);
    return abcd_of_transmission_line_complex(gamma, length_m, std::complex<double>(Z0, 0.0));
}

inline Matrix2cd cascade_abcd(const Matrix2cd &A, const Matrix2cd &B) {
    return A * B;
}

// Convert series connection of multiple two-port ABCD matrices
inline Matrix2cd cascade_list(const std::vector<Matrix2cd> &list) {
    Matrix2cd acc = Matrix2cd::Identity();
    for(const auto &m : list) acc = acc * m;
    return acc;
}

} // namespace skrf_cpp
