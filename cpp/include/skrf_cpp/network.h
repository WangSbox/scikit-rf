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
    std::vector<std::complex<double>> s_interp(double freq_hz) const {
        if(freqs.empty()) throw std::runtime_error("no frequency points");
        if(freq_hz <= freqs.front().hz) return sparams.front();
        if(freq_hz >= freqs.back().hz) return sparams.back();

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
    }

    // Convert this flat-format Network to a Matrix-based NetworkEigen
    NetworkEigen to_network_eigen() const {
        NetworkEigen out;
        out.n_ports = n_ports;
        out.z0 = z0;
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
