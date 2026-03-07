#pragma once

#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include "frequency.h"
#include "network_eigen.h"
#include <stdexcept>
#include <algorithm>

namespace skrf_cpp {

class Network {
public:
    int n_ports{0};
    double z0{50.0};
    std::vector<Frequency> freqs;
    // optional per-frequency reference impedances (if present in source, e.g., HFSS per-frequency port impedance)
    std::vector<double> per_freq_z0;
    // optional per-frequency gamma/comments parsed from Touchstone (best-effort)
    std::vector<std::complex<double>> per_freq_gamma;
    // For each frequency, a flat vector of size n_ports*n_ports (row-major)
    std::vector<std::vector<std::complex<double>>> sparams;

    Network() = default;

    void printSummary() const {
        std::cout << "Network: ports=" << n_ports << " points=" << freqs.size() << "\n";
        if(!freqs.empty()) {
            std::cout << "  first freq: " << freqs.front().hz << " Hz\n";
            if(!sparams.empty() && !sparams.front().empty()) {
                std::cout << "  first S[0,0] = " << sparams.front()[0] << "\n";
            }
        }
    }

    // Return flat S vector at index (copy)
    std::vector<std::complex<double>> s_at_index(size_t idx) const {
        if(idx >= sparams.size()) throw std::out_of_range("freq index");
        return sparams[idx];
    }

    // Simple linear interpolation of S-parameters between nearest freq points
    std::vector<std::complex<double>> s_interp(double freq_hz, const std::string &method = "linear") const {
        if(freqs.empty()) throw std::runtime_error("no frequency points");
        if(freq_hz <= freqs.front().hz) return sparams.front();
        if(freq_hz >= freqs.back().hz) return sparams.back();
        if(method == "linear") {
            // find upper index
            size_t i = 1;
            while(i < freqs.size() && freqs[i].hz < freq_hz) ++i;
            double f0 = freqs[i-1].hz;
            double f1 = freqs[i].hz;
            double t = (freq_hz - f0) / (f1 - f0);

            const auto &s0 = sparams[i-1];
            const auto &s1 = sparams[i];
            if(s0.size() != s1.size()) throw std::runtime_error("sparam size mismatch");
            std::vector<std::complex<double>> out(s0.size());
            for(size_t k=0;k<s0.size();++k) out[k] = s0[k] * (1.0 - t) + s1[k] * t;
            return out;
        } else if(method == "cubic") {
            // natural cubic spline for each s-parameter component (real and imag separately)
            size_t nfreq = freqs.size();
            size_t ncomp = sparams.front().size();
            std::vector<double> xs(nfreq);
            for(size_t i=0;i<nfreq;++i) xs[i] = freqs[i].hz;
            std::vector<std::complex<double>> out(ncomp, std::complex<double>(0,0));
            // For each component build y vector and compute second derivs
            for(size_t k=0;k<ncomp;++k) {
                std::vector<double> yreal(nfreq), yimag(nfreq);
                for(size_t i=0;i<nfreq;++i) { yreal[i] = std::real(sparams[i][k]); yimag[i] = std::imag(sparams[i][k]); }
                auto y2r = compute_natural_cubic_second_derivs(xs, yreal);
                auto y2i = compute_natural_cubic_second_derivs(xs, yimag);
                double yr = spline_eval(xs, yreal, y2r, freq_hz);
                double yi = spline_eval(xs, yimag, y2i, freq_hz);
                out[k] = std::complex<double>(yr, yi);
            }
            return out;
        } else {
            throw std::runtime_error("s_interp: unknown method");
        }
    }

private:
    // compute second derivatives for natural cubic spline (vector y over x)
    static std::vector<double> compute_natural_cubic_second_derivs(const std::vector<double> &x, const std::vector<double> &y) {
        size_t n = x.size();
        if(n < 2) return std::vector<double>(n, 0.0);
        std::vector<double> u(n-1), y2(n);
        y2[0] = 0.0; u[0] = 0.0;
        for(size_t i=1;i<n-1;++i) {
            double sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
            double p = sig * y2[i-1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            double d1 = (y[i+1] - y[i]) / (x[i+1] - x[i]);
            double d0 = (y[i] - y[i-1]) / (x[i] - x[i-1]);
            u[i] = (6.0 * (d1 - d0) / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
        }
        y2[n-1] = 0.0;
        for(size_t k = n-1; k-- > 0;) {
            y2[k] = y2[k] * y2[k+1] + u[k];
        }
        return y2;
    }

    // evaluate natural cubic spline at xq given x,y and second derivs y2
    static double spline_eval(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &y2, double xq) {
        size_t n = x.size();
        if(n == 0) return 0.0;
        if(n == 1) return y[0];
        // find right place
        size_t klo = 0, khi = n-1;
        while(khi - klo > 1) {
            size_t k = (khi + klo) >> 1;
            if(x[k] > xq) khi = k; else klo = k;
        }
        double h = x[khi] - x[klo];
        if(h == 0.0) return y[klo];
        double a = (x[khi] - xq) / h;
        double b = (xq - x[klo]) / h;
        double val = a * y[klo] + b * y[khi] + ((a*a*a - a) * y2[klo] + (b*b*b - b) * y2[khi]) * (h*h) / 6.0;
        return val;
    }

    // Return bracketing indices and interpolation factor t in [0,1]
    // If freq_hz is below/above range, both indices will be equal to 0 or last respectively and t=0.
    std::tuple<size_t,size_t,double> freq_bracket(double freq_hz, double tol = 1e-12) const {
        if(freqs.empty()) throw std::runtime_error("no frequency points");
        if(freq_hz <= freqs.front().hz + tol) return {0, 0, 0.0};
        if(freq_hz >= freqs.back().hz - tol) {
            size_t last = freqs.size() - 1; return {last, last, 0.0};
        }
        size_t i = 1;
        while(i < freqs.size() && freqs[i].hz < freq_hz) ++i;
        double f0 = freqs[i-1].hz; double f1 = freqs[i].hz;
        double t = (f1==f0) ? 0.0 : ((freq_hz - f0) / (f1 - f0));
        return {i-1, i, t};
    }

    // Convenience alias for to_network_eigen()
    NetworkEigen as_eigen() const { return to_network_eigen(); }

    // Convert this flat-format Network to a Matrix-based NetworkEigen
    NetworkEigen to_network_eigen() const {
        NetworkEigen out;
        out.n_ports = n_ports;
        out.z0 = z0;
        out.per_freq_z0 = per_freq_z0;
        out.per_freq_gamma = per_freq_gamma;
        out.freqs = freqs;
        out.s_params.clear();
        out.s_params.reserve(sparams.size());
        for(const auto &flat : sparams) {
            if(static_cast<int>(flat.size()) != n_ports * n_ports) throw std::runtime_error("sparam length mismatch");
            MatrixXcd M = MatrixXcd::Zero(n_ports, n_ports);
            for(int r=0;r<n_ports;++r) {
                for(int c=0;c<n_ports;++c) {
                    M(r,c) = flat[r * n_ports + c];
                }
            }
            out.s_params.push_back(std::move(M));
        }
        return out;
    }

    // Resample Network to new frequency points (Hz)
    Network resample(const std::vector<double> &new_freqs) const {
        NetworkEigen ne = to_network_eigen();
        NetworkEigen re = ne.resample(new_freqs);
        return re.to_flat_network();
    }

    // Extract subset of ports (by index) into a new flat Network
    Network extract_ports(const std::vector<int> &ports_idx) const {
        NetworkEigen ne = to_network_eigen();
        NetworkEigen sub = ne.extract_ports(ports_idx);
        return sub.to_flat_network();
    }

    // Merge/connect ports according to pairs and return a new Network with fewer ports
    Network merge_ports(const std::vector<std::pair<int,int>> &pairs) const {
        NetworkEigen ne = to_network_eigen();
        NetworkEigen merged = ne.merge_ports(pairs);
        return merged.to_flat_network();
    }

    // Convenience: convert to Z-parameters (returns Matrix-based NetworkEigen)
    NetworkEigen to_eigen_z(double z0_scalar = 50.0) const {
        NetworkEigen ne = to_network_eigen();
        return ne.to_z(z0_scalar);
    }
};

}
