#pragma once

#include <vector>
#include <complex>
#include <Eigen/Dense>

namespace skrf_cpp {

// Lightweight vector fitting: simple pole-residue fit using fixed initial poles.
// This is not a full implementation of the Vector Fitting algorithm, but
// provides a practical rational fit interface for small-order models.
class VectorFitting {
public:
    VectorFitting() = default;

    // Fit rational model H(s) ≈ sum_k r_k / (s - p_k) + d + s*e
    // freqs: Hz, vals: complex-valued responses at freqs
    // n_poles: number of poles to use (default 6)
    void fit(const std::vector<double>& freqs,
             const std::vector<std::complex<double>>& vals,
             int n_poles = 6);

    // evaluate fitted model at complex frequency s (rad/s)
    std::complex<double> eval(std::complex<double> s) const;

    const std::vector<std::complex<double>>& poles() const { return poles_; }
    const std::vector<std::complex<double>>& residues() const { return residues_; }
    std::complex<double> constant_term() const { return d_; }
    std::complex<double> proportional_term() const { return e_; }

private:
    std::vector<std::complex<double>> poles_;
    std::vector<std::complex<double>> residues_;
    std::complex<double> d_{0.0};
    std::complex<double> e_{0.0};
};

} // namespace skrf_cpp
