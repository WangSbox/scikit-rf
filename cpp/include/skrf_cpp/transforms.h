#pragma once

#include <Eigen/Dense>

namespace skrf_cpp {

using MatrixXcd = Eigen::MatrixXcd;

// Convert S to Z: Z = Z0 * (I + S) * inv(I - S)
inline MatrixXcd s_to_z(const MatrixXcd &S, const MatrixXcd &Z0) {
    MatrixXcd I = MatrixXcd::Identity(S.rows(), S.cols());
    return Z0 * (I + S) * (I - S).inverse();
}

// convenience overload: scalar Z0
inline MatrixXcd s_to_z(const MatrixXcd &S, double Z0_scalar) {
    MatrixXcd Z0 = MatrixXcd::Identity(S.rows(), S.cols()) * Z0_scalar;
    return s_to_z(S, Z0);
}

// 2x2 overload
inline Eigen::Matrix2cd s_to_z(const Eigen::Matrix2cd &S, double Z0_scalar) {
    Eigen::Matrix2cd I = Eigen::Matrix2cd::Identity();
    return Z0_scalar * (I + S) * (I - S).inverse();
}

// Convert Z to S: S = (Z - Z0) * inv(Z + Z0)
inline MatrixXcd z_to_s(const MatrixXcd &Z, const MatrixXcd &Z0) {
    return (Z - Z0) * (Z + Z0).inverse();
}

// convenience overload: scalar Z0
inline MatrixXcd z_to_s(const MatrixXcd &Z, double Z0_scalar) {
    MatrixXcd Z0 = MatrixXcd::Identity(Z.rows(), Z.cols()) * Z0_scalar;
    return z_to_s(Z, Z0);
}

// 2x2 overload
inline Eigen::Matrix2cd z_to_s(const Eigen::Matrix2cd &Z, double Z0_scalar) {
    Eigen::Matrix2cd I = Eigen::Matrix2cd::Identity();
    return (Z - Z0_scalar * I) * (Z + Z0_scalar * I).inverse();
}

// 2-port S <-> ABCD conversions (Z0 scalar)
inline Eigen::Matrix2cd s_to_abcd_2port(const Eigen::Matrix2cd &S, double Z0) {
    using cplx = std::complex<double>;
    cplx S11 = S(0,0), S12 = S(0,1), S21 = S(1,0), S22 = S(1,1);
    cplx denom = 2.0 * S21;
    Eigen::Matrix2cd M;
    M(0,0) = ((1.0 + S11) * (1.0 - S22) + S12 * S21) / denom; // A
    M(0,1) = Z0 * ((1.0 + S11) * (1.0 + S22) - S12 * S21) / denom; // B
    M(1,0) = (1.0 / Z0) * ((1.0 - S11) * (1.0 - S22) - S12 * S21) / denom; // C
    M(1,1) = ((1.0 - S11) * (1.0 + S22) + S12 * S21) / denom; // D
    return M;
}

inline Eigen::Matrix2cd abcd_to_s_2port(const Eigen::Matrix2cd &ABCD, double Z0) {
    using cplx = std::complex<double>;
    cplx A = ABCD(0,0), B = ABCD(0,1), C = ABCD(1,0), D = ABCD(1,1);
    cplx denom = A + B / Z0 + C * Z0 + D;
    Eigen::Matrix2cd S;
    S(0,0) = (A + B / Z0 - C * Z0 - D) / denom;
    S(0,1) = 2.0 * (A * D - B * C) / denom;
    S(1,0) = 2.0 / denom;
    S(1,1) = (-A + B / Z0 - C * Z0 + D) / denom;
    return S;
}

}
