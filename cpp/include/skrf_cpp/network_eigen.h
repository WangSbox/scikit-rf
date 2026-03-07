#pragma once

#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "network.h"
#include "frequency.h"
#include "mathFunctions.h"
#include <algorithm>
#include "transforms.h"

namespace skrf_cpp {

using MatrixXcd = Eigen::MatrixXcd;

class NetworkEigen {
public:
    int n_ports{0};
    double z0{50.0};
    std::vector<Frequency> freqs;
    std::vector<MatrixXcd> s_params; // per-frequency S matrix

    NetworkEigen() = default;

    void resize(int ports, size_t n_freqs) {
        n_ports = ports;
        freqs.resize(n_freqs);
        s_params.resize(n_freqs);
        for(auto &m : s_params) m = MatrixXcd::Zero(ports, ports);
    }

    void printSummary() const {
        std::cout << "NetworkEigen: ports=" << n_ports << " points=" << freqs.size() << "\n";
        if(!freqs.empty() && !s_params.empty()) {
            std::cout << "  first freq: " << freqs.front().hz << " Hz\n";
            if(n_ports>0) std::cout << "  first S(0,0) = " << s_params.front()(0,0) << "\n";
        }
    }

    // return S matrix at index (copy)
    MatrixXcd s_at_index(size_t idx) const {
        if(idx >= s_params.size()) throw std::out_of_range("freq index");
        return s_params[idx];
    }

    // simple linear interpolation of S matrices between nearest freq points
    MatrixXcd s_interp(double freq_hz) const {
        if(freqs.empty()) throw std::runtime_error("no frequency points");
        if(freq_hz <= freqs.front().hz) return s_params.front();
        if(freq_hz >= freqs.back().hz) return s_params.back();

        // find upper index
        size_t i = 1;
        while(i < freqs.size() && freqs[i].hz < freq_hz) ++i;
        double f0 = freqs[i-1].hz;
        double f1 = freqs[i].hz;
        double t = (freq_hz - f0) / (f1 - f0);
        return s_params[i-1] * (1.0 - t) + s_params[i] * t;
    }

    // convert S-parameters to Z-parameters assuming scalar Z0 (same on all ports)
    // returns a new NetworkEigen with Z matrices
    NetworkEigen to_z(double z0_scalar = 50.0) const {
        NetworkEigen out;
        out.n_ports = n_ports;
        out.freqs = freqs;
        out.s_params.clear();
        MatrixXcd Z0 = MatrixXcd::Identity(n_ports, n_ports) * z0_scalar;
        for(const auto &S : s_params) {
            MatrixXcd Z = s_to_z(S, Z0);
            out.s_params.push_back(Z);
        }
        return out;
    }

    // convert S to Y (admittance): first to Z then invert
    NetworkEigen to_y(double z0_scalar = 50.0) const {
        NetworkEigen out = to_z(z0_scalar);
        for(auto &Z : out.s_params) {
            Z = Z.inverse();
        }
        return out;
    }

    // convert S to ABCD for 2-port networks (Z0 scalar). For n_ports != 2 throws.
    NetworkEigen to_abcd(double z0_scalar = 50.0) const {
        if(n_ports != 2) throw std::runtime_error("to_abcd implemented only for 2-port networks; for n-port consider extracting a 2-port subnetwork or implementing a custom conversion");
        NetworkEigen out;
        out.n_ports = n_ports;
        out.freqs = freqs;
        out.s_params.clear();
        for(const auto &S : s_params) {
            Eigen::Matrix2cd S2 = S.block<2,2>(0,0);
            Eigen::Matrix2cd ABCD = s_to_abcd_2port(S2, z0_scalar);
            MatrixXcd M = MatrixXcd::Zero(2,2);
            M(0,0) = ABCD(0,0);
            M(0,1) = ABCD(0,1);
            M(1,0) = ABCD(1,0);
            M(1,1) = ABCD(1,1);
            out.s_params.push_back(M);
        }
        return out;
    }

    // Convert matrix-based NetworkEigen back to flat-format Network
    Network to_flat_network() const {
        Network out;
        out.n_ports = n_ports;
        out.z0 = z0;
        out.freqs = freqs;
        out.sparams.clear();
        out.sparams.reserve(s_params.size());
        for(const auto &M : s_params) {
            if(M.rows() != static_cast<int>(n_ports) || M.cols() != static_cast<int>(n_ports))
                throw std::runtime_error("matrix size mismatch when flattening network");
            std::vector<std::complex<double>> flat(n_ports * n_ports);
            for(int r=0;r<static_cast<int>(n_ports);++r) {
                for(int c=0;c<static_cast<int>(n_ports);++c) {
                    flat[r * n_ports + c] = M(r,c);
                }
            }
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

    // Cascade this 2-port network with another (this followed by other).
    // Both networks must have same frequency points and be 2-port.
    NetworkEigen cascade_with(const NetworkEigen &other, double z0_scalar = 50.0) const {
        if(n_ports != 2 || other.n_ports != 2) throw std::runtime_error("cascade_with supports only 2-port networks; for n-port workflows use merge_block_diag/connect_and_merge_indices or extract_ports to create 2-port chains");
        if(freqs.size() != other.freqs.size()) throw std::runtime_error("frequency point mismatch for cascade");
        NetworkEigen out;
        out.n_ports = 2;
        out.z0 = z0_scalar;
        out.freqs = freqs;
        out.s_params.clear();
        for(size_t i=0;i<s_params.size();++i) {
            Eigen::Matrix2cd S1 = s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd S2 = other.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd A1 = s_to_abcd_2port(S1, z0_scalar);
            Eigen::Matrix2cd A2 = s_to_abcd_2port(S2, z0_scalar);
            Eigen::Matrix2cd Atot = A1 * A2;
            Eigen::Matrix2cd Stot = abcd_to_s_2port(Atot, z0_scalar);
            Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2,2);
            M(0,0) = Stot(0,0); M(0,1) = Stot(0,1); M(1,0) = Stot(1,0); M(1,1) = Stot(1,1);
            out.s_params.push_back(M);
        }
        return out;
    }

    // series-cascade a list of 2-port NetworkEigen objects (A1 * A2 * ...)
    static NetworkEigen cascade_list(const std::vector<NetworkEigen> &list, double z0_scalar = 50.0) {
        if(list.empty()) return NetworkEigen();
        // all must be 2-port and same freq size
        size_t nfreq = list.front().freqs.size();
        for(const auto &n : list) {
            if(n.n_ports != 2) throw std::runtime_error("cascade_list supports only 2-port networks; provide only 2-port elements or pre-extract 2-port chains");
            if(n.freqs.size() != nfreq) throw std::runtime_error("frequency point mismatch in cascade_list");
        }
        NetworkEigen out;
        out.n_ports = 2;
        out.z0 = z0_scalar;
        out.freqs = list.front().freqs;
        out.s_params.clear();
        for(size_t fi=0; fi<nfreq; ++fi) {
            Eigen::Matrix2cd Aacc = Eigen::Matrix2cd::Identity();
            for(const auto &n : list) {
                Eigen::Matrix2cd S = n.s_params[fi].block<2,2>(0,0);
                Eigen::Matrix2cd A = s_to_abcd_2port(S, z0_scalar);
                Aacc = Aacc * A;
            }
            Eigen::Matrix2cd Stot = abcd_to_s_2port(Aacc, z0_scalar);
            MatrixXcd M = MatrixXcd::Zero(2,2);
            M(0,0)=Stot(0,0); M(0,1)=Stot(0,1); M(1,0)=Stot(1,0); M(1,1)=Stot(1,1);
            out.s_params.push_back(M);
        }
        return out;
    }

    // Parallel (shunt) connection: Y_total = Y1 + Y2 ; returns combined S
    NetworkEigen parallel_with(const NetworkEigen &other, double z0_scalar = 50.0) const {
        if(n_ports != other.n_ports) throw std::runtime_error("parallel networks must have same port count");
        if(freqs.size() != other.freqs.size()) throw std::runtime_error("frequency point mismatch for parallel");
        NetworkEigen out;
        out.n_ports = n_ports;
        out.freqs = freqs;
        out.s_params.clear();
        for(size_t i=0;i<s_params.size();++i) {
            const MatrixXcd &S1 = s_params[i];
            const MatrixXcd &S2 = other.s_params[i];
            MatrixXcd Z1 = s_to_z(S1, z0_scalar);
            MatrixXcd Z2 = s_to_z(S2, z0_scalar);
            MatrixXcd Y1 = Z1.inverse();
            MatrixXcd Y2 = Z2.inverse();
            MatrixXcd Y = Y1 + Y2;
            MatrixXcd Z = Y.inverse();
            MatrixXcd S = z_to_s(Z, z0_scalar);
            out.s_params.push_back(S);
        }
        return out;
    }

    // Reorder ports according to 'order' vector (order.size()==n_ports).
    // Semantics: new(i,j) = old(order[i], order[j])
    NetworkEigen reorder_ports(const std::vector<int> &order) const {
        if(order.size() != static_cast<size_t>(n_ports)) throw std::runtime_error("order size mismatch");
        NetworkEigen out;
        out.n_ports = n_ports;
        out.freqs = freqs;
        out.s_params.clear();
        for(const auto &M : s_params) {
            MatrixXcd N = MatrixXcd::Zero(n_ports, n_ports);
            for(int i=0;i<static_cast<int>(n_ports);++i) {
                for(int j=0;j<static_cast<int>(n_ports);++j) {
                    int oi = order[i];
                    int oj = order[j];
                    if(oi < 0 || oi >= static_cast<int>(n_ports) || oj < 0 || oj >= static_cast<int>(n_ports))
                        throw std::runtime_error("invalid port index in reorder_ports");
                    N(i,j) = M(oi, oj);
                }
            }
            out.s_params.push_back(std::move(N));
        }
        return out;
    }

    // Extract a subset of ports (by index), preserving order
    NetworkEigen extract_ports(const std::vector<int> &ports_idx) const {
        NetworkEigen out;
        if(ports_idx.empty()) return out;
        int newp = static_cast<int>(ports_idx.size());
        out.n_ports = newp;
        out.z0 = z0;
        out.freqs = freqs;
        out.s_params.clear();
        for(const auto &M : s_params) {
            MatrixXcd N = MatrixXcd::Zero(newp, newp);
            for(int i=0;i<newp;++i) {
                for(int j=0;j<newp;++j) {
                    int oi = ports_idx[i];
                    int oj = ports_idx[j];
                    if(oi < 0 || oi >= n_ports || oj < 0 || oj >= n_ports) throw std::runtime_error("invalid port index in extract_ports");
                    N(i,j) = M(oi, oj);
                }
            }
            out.s_params.push_back(std::move(N));
        }
        return out;
    }

    // Resample/interpolate this NetworkEigen to a new frequency vector (Hz)
    NetworkEigen resample(const std::vector<double> &new_freqs) const {
        NetworkEigen out;
        out.n_ports = n_ports;
        out.z0 = z0;
        out.freqs.clear();
        out.s_params.clear();
        for(double f : new_freqs) out.freqs.emplace_back(f);
        for(double f : new_freqs) {
            MatrixXcd S = s_interp(f);
            out.s_params.push_back(S);
        }
        return out;
    }

    // Merge pairs of ports: each pair (i,j) will be connected (zero-impedance) and collapsed into one port.
    // Pairs are applied sequentially; indices refer to the current numbering and must be valid.
    NetworkEigen merge_ports(const std::vector<std::pair<int,int>> &pairs) const {
        NetworkEigen out = *this;
        for(const auto &pr : pairs) {
            int i = pr.first;
            int j = pr.second;
            if(i == j) continue;
            if(i < 0 || j < 0 || i >= out.n_ports || j >= out.n_ports) throw std::runtime_error("merge_ports: invalid port index");
            // ensure i < j for removal simplicity
            if(i > j) std::swap(i,j);
            int n = out.n_ports;
            int n2 = n - 1;
            NetworkEigen next;
            next.n_ports = n2;
            next.z0 = out.z0;
            next.freqs = out.freqs;
            next.s_params.clear();
            for(const auto &S : out.s_params) {
                // convert S to Z then to Y
                MatrixXcd Z = s_to_z(S, out.z0);
                MatrixXcd Y = Z.inverse();
                MatrixXcd Y2 = MatrixXcd::Zero(n2, n2);
                std::vector<int> old_to_new(n, -1);
                int idx = 0;
                for(int old=0; old<n; ++old) {
                    if(old == j) continue;
                    old_to_new[old] = idx++;
                }

                for(int a=0; a<n; ++a) {
                    if(a == j) continue;
                    for(int b=0; b<n; ++b) {
                        if(b == j) continue;
                        if(a==i && b==i) {
                            // merged diagonal
                            Y2(old_to_new[a], old_to_new[b]) = Y(i,i) + Y(j,j) + Y(i,j) + Y(j,i);
                        } else if(a==i && b!=i) {
                            Y2(old_to_new[a], old_to_new[b]) = Y(i,b) + Y(j,b);
                        } else if(a!=i && b==i) {
                            Y2(old_to_new[a], old_to_new[b]) = Y(a,i) + Y(a,j);
                        } else {
                            Y2(old_to_new[a], old_to_new[b]) = Y(a,b);
                        }
                    }
                }
                // convert back to Z then to S
                MatrixXcd Z2 = Y2.inverse();
                MatrixXcd S2 = z_to_s(Z2, out.z0);
                next.s_params.push_back(S2);
            }
            out = std::move(next);
        }
        return out;
    }
};

}
