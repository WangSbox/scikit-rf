#pragma once

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "mline.h"

namespace skrf_cpp {
namespace media {

using namespace Eigen;

// DistributedCircuit: sequence of transmission-line segments
// Each segment defined by z0 (double), gamma (complex), length (double)
class DistributedCircuit {
public:
    struct Segment { double z0; std::complex<double> gamma; double length; };
    std::vector<Segment> segments;

    DistributedCircuit() {}

    void add_segment(double z0, std::complex<double> gamma, double length) {
        segments.push_back({z0,gamma,length});
    }

    // compute overall ABCD matrix (2x2) by cascading each segment's ABCD
    Matrix<std::complex<double>,2,2> overall_abcd() const {
        Matrix<std::complex<double>,2,2> M = Matrix<std::complex<double>,2,2>::Identity();
        for(const auto &seg : segments) {
            tline::MLine ml(seg.z0, seg.gamma, seg.length);
            M = M * ml.abcd();
        }
        return M;
    }

    // convert overall ABCD to S for given reference Z0
    Matrix<std::complex<double>,2,2> overall_s(double ref_z0 = 50.0) const {
        auto A = overall_abcd();
        using cd = std::complex<double>;
        cd a=A(0,0), b=A(0,1), c=A(1,0), d=A(1,1);
        cd denom = (a + b/ref_z0 + c*ref_z0 + d);
        Matrix<std::complex<double>,2,2> S;
        S(0,0) = (a + b/ref_z0 - c*ref_z0 - d)/denom;
        S(0,1) = (2.0*(a*d - b*c))/denom;
        S(1,0) = 2.0/denom;
        S(1,1) = (-a + b/ref_z0 - c*ref_z0 + d)/denom;
        return S;
    }
};

} // namespace media
} // namespace skrf_cpp
