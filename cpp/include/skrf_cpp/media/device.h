#pragma once

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "mline.h"

namespace skrf_cpp {
namespace media {

using namespace Eigen;

// Generic Device: represents a 2-port linear device by its frequency-dependent
// S-parameters stored per-frequency. Provides simple helpers to compute
// ABCD and convert representations.
class Device {
public:
    std::vector<double> freqs; // Hz
    std::vector<Matrix<std::complex<double>,2,2>> s_params; // per-frequency 2x2 S
    double z0 = 50.0;

    Device() {}

    // construct from frequency list and corresponding S matrices
    Device(const std::vector<double> &f, const std::vector<Matrix<std::complex<double>,2,2>> &S, double z0_=50.0)
        : freqs(f), s_params(S), z0(z0_) {
        if(freqs.size() != s_params.size()) throw std::runtime_error("Device: freqs and s_params size mismatch");
    }

    // return ABCD for given index
    Matrix<std::complex<double>,2,2> abcd_at(size_t idx) const {
        if(idx >= s_params.size()) throw std::out_of_range("Device: index out of range");
        const auto &S = s_params[idx];
        // convert S to ABCD for 2-port
        // Using standard formula assuming reference impedance z0
        std::complex<double> s11 = S(0,0), s12 = S(0,1), s21 = S(1,0), s22 = S(1,1);
        std::complex<double> delta = (1.0 - s11)*(1.0 - s22) - s12*s21;
        Matrix<std::complex<double>,2,2> A;
        A(0,0) = ((1.0 + s11)*(1.0 - s22) + s12*s21) / (2.0*s21);
        A(0,1) = z0 * ((1.0 + s11)*(1.0 + s22) - s12*s21) / (2.0*s21);
        A(1,0) = (1.0 / z0) * ((1.0 - s11)*(1.0 - s22) - s12*s21) / (2.0*s21);
        A(1,1) = ((1.0 - s11)*(1.0 + s22) + s12*s21) / (2.0*s21);
        return A;
    }
};

} // namespace media
} // namespace skrf_cpp
