#pragma once

#include <complex>
#include <vector>
#include <Eigen/Dense>

namespace skrf_cpp {
namespace tline {

using namespace Eigen;

// Simple transmission line segment model that can return ABCD matrix
class MLine {
public:
    double z0; // characteristic impedance (real or complex)
    std::complex<double> gamma; // propagation constant
    double length; // meters

    MLine(double z0_, std::complex<double> gamma_, double length_)
        : z0(z0_), gamma(gamma_), length(length_) {}

    // Return 2x2 ABCD matrix as Eigen::Matrix2cd
    Matrix<std::complex<double>,2,2> abcd() const {
        using cd = std::complex<double>;
        cd coshgl = std::cosh(gamma * length);
        cd sinhgl = std::sinh(gamma * length);
        Matrix<std::complex<double>,2,2> M;
        M(0,0) = coshgl;
        M(0,1) = z0 * sinhgl;
        M(1,0) = (1.0 / z0) * sinhgl;
        M(1,1) = coshgl;
        return M;
    }
};

} // namespace tline
} // namespace skrf_cpp
