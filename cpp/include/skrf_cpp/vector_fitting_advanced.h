#pragma once

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "network.h"
#include "network_eigen.h"

namespace skrf_cpp {

class VectorFittingAdvanced {
public:
    std::vector<std::complex<double>> poles; // current poles
    std::vector<std::complex<double>> residues;
    std::complex<double> d = 0.0;
    std::complex<double> e = 0.0;

    struct FitOptions {
        int npoles = 8;
        int max_iter = 12;
        double tol = 1e-6;
        double regularization = 0.0; // Tikhonov regularization lambda
        // optional weights per sample (same length as freqs/H). If empty, uniform weights used.
        std::vector<double> weights;
        // If true, adapt weights each iteration inversely proportional to residual magnitude
        bool adaptive_weights = false;
        // small floor added when adapting weights
        double adapt_eps = 1e-8;
    };

    // Fit H(freq) (complex) sampled at freqs (Hz) with options controlling poles, weights and regularization
    void fit(const std::vector<double> &freqs_hz,
             const std::vector<std::complex<double>> &H,
             const FitOptions &opts = FitOptions());

    // Convenience: fit from a single-column S-parameter vector (alias)
    void fit_sparameter(const std::vector<double> &freqs_hz,
                        const std::vector<std::complex<double>> &S,
                        const FitOptions &opts = FitOptions()) {
        fit(freqs_hz, S, opts);
    }

    // Fit from a `Network` (flat) selecting S(row,col)
    void fit_network(const skrf_cpp::Network &net, int row, int col, const FitOptions &opts = FitOptions());

    // Fit from a `NetworkEigen` selecting S(row,col)
    void fit_network_eigen(const skrf_cpp::NetworkEigen &net, int row, int col, const FitOptions &opts = FitOptions());

    // Evaluate fitted rational model at s = j*2pi*f
    std::complex<double> eval(double freq_hz) const;

    // Evaluate at multiple frequencies
    std::vector<std::complex<double>> eval_vector(const std::vector<double> &freqs_hz) const {
        std::vector<std::complex<double>> out; out.reserve(freqs_hz.size());
        for(double f : freqs_hz) out.push_back(eval(f));
        return out;
    }

    // Accessors
    int get_model_order() const { return static_cast<int>(poles.size()); }
    const std::vector<std::complex<double>>& get_poles() const { return poles; }
    const std::vector<std::complex<double>>& get_residues() const { return residues; }

    // Simple serialization: save poles/residues/d/e to a text file (one value per line)
    bool save_to_file(const std::string &path) const;
    bool load_from_file(const std::string &path);
};

} // namespace skrf_cpp
