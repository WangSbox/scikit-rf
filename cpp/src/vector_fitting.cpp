#include "../include/skrf_cpp/vector_fitting.h"
#include <complex>
#include <cmath>

namespace skrf_cpp {

void VectorFitting::fit(const std::vector<double>& freqs,
                        const std::vector<std::complex<double>>& vals,
                        int n_poles) {
    size_t N = freqs.size();
    if(N == 0 || vals.size() != N) throw std::runtime_error("invalid input to VectorFitting::fit");
    if(n_poles <= 0) n_poles = 6;

    // initial poles: spread on imaginary axis (rad/s)
    poles_.clear(); residues_.clear();
    double fmax = freqs.back();
    double wmax = 2.0 * M_PI * fmax;
    for(int k=0;k<n_poles;++k) {
        double frac = double(k+1) / double(n_poles+1);
        double omega = frac * wmax;
        // place poles slightly in LHP for stability
        std::complex<double> p = std::complex<double>(-omega*0.01, omega);
        poles_.push_back(p);
    }

    // Build least squares system: A x = y, x = [r_1..r_n, d, e]
    int ncols = n_poles + 2;
    Eigen::MatrixXcd A(N, ncols);
    Eigen::VectorXcd y(N);
    for(size_t i=0;i<N;++i) {
        std::complex<double> s = std::complex<double>(0.0, 2.0 * M_PI * freqs[i]);
        for(int k=0;k<n_poles;++k) {
            A(i,k) = 1.0 / (s - poles_[k]);
        }
        A(i, n_poles) = 1.0;      // d
        A(i, n_poles+1) = s;      // e * s
        y(i) = vals[i];
    }

    // Solve least squares
    Eigen::VectorXcd x = A.colPivHouseholderQr().solve(y);

    residues_.resize(n_poles);
    for(int k=0;k<n_poles;++k) residues_[k] = x(k);
    d_ = x(n_poles);
    e_ = x(n_poles+1);
}

std::complex<double> VectorFitting::eval(std::complex<double> s) const {
    std::complex<double> y = d_ + e_ * s;
    for(size_t k=0;k<poles_.size();++k) {
        y += residues_[k] / (s - poles_[k]);
    }
    return y;
}

} // namespace skrf_cpp
