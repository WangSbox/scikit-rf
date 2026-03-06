#pragma once

#include <vector>
#include <string>
#include <functional>
#include <filesystem>
#include <algorithm>
#include <map>
#include "network.h"
#include "io/general.h"
#include <numeric>
#include <cmath>

namespace skrf_cpp {

class NetworkSet {
public:
    std::vector<Network> ntwk_set;
    std::string name;

    NetworkSet() = default;
    explicit NetworkSet(const std::vector<Network>& v): ntwk_set(v) {}

    size_t size() const { return ntwk_set.size(); }

    void append(const Network &n) { ntwk_set.push_back(n); }

    void remove_at(size_t idx) {
        if(idx >= ntwk_set.size()) throw std::out_of_range("remove_at index");
        ntwk_set.erase(ntwk_set.begin() + idx);
    }

    // Minimal API stubs mirroring Python NetworkSet. Implementations
    // can be expanded later to compute mean/std/plotting etc.
    Network mean_s() const {
        Network out;
        if(ntwk_set.empty()) return out;
        // assume all networks have same shape
        const Network &ref = ntwk_set.front();
        out.n_ports = ref.n_ports;
        out.freqs = ref.freqs;
        size_t nfreq = ref.freqs.size();
        out.sparams.resize(nfreq);
        size_t nnet = ntwk_set.size();

        for(size_t fi=0; fi<nfreq; ++fi) {
            const auto &ref_s = ref.sparams[fi];
            std::vector<std::complex<double>> accum(ref_s.size(), std::complex<double>(0,0));
            for(const auto &net : ntwk_set) {
                const auto &s = net.sparams[fi];
                for(size_t k=0;k<s.size();++k) accum[k] += s[k];
            }
            for(auto &c : accum) c /= static_cast<double>(nnet);
            out.sparams[fi] = std::move(accum);
        }
        return out;
    }

    Network std_s() const {
        Network out;
        if(ntwk_set.empty()) return out;
        const Network &ref = ntwk_set.front();
        out.n_ports = ref.n_ports;
        out.freqs = ref.freqs;
        size_t nfreq = ref.freqs.size();
        out.sparams.resize(nfreq);
        size_t nnet = ntwk_set.size();

        // first compute mean
        Network mean = mean_s();

        for(size_t fi=0; fi<nfreq; ++fi) {
            const auto &mean_s = mean.sparams[fi];
            std::vector<std::complex<double>> var(mean_s.size(), std::complex<double>(0,0));
            for(const auto &net : ntwk_set) {
                const auto &s = net.sparams[fi];
                for(size_t k=0;k<s.size();++k) {
                    std::complex<double> d = s[k] - mean_s[k];
                    var[k] += d * std::conj(d);
                }
            }
            // std = sqrt( var / (n-1) ) ; take magnitude
            std::vector<std::complex<double>> out_s(mean_s.size());
            for(size_t k=0;k<mean_s.size();++k) {
                double v = std::real(var[k]) / static_cast<double>(std::max<size_t>(1, nnet-1));
                out_s[k] = std::complex<double>(std::sqrt(v), 0.0);
            }
            out.sparams[fi] = std::move(out_s);
        }
        return out;
    }

    // Create a NetworkSet from all networks in a directory using C++ IO helpers
    static NetworkSet from_dir(const std::string &dir) {
        NetworkSet out;
        auto pairs = skrf_cpp::io::read_all_networks(dir);
        for(const auto &p : pairs) out.ntwk_set.push_back(p.second);
        return out;
    }

    // Create from a map of name -> flat S arrays. Frequency must be provided (Hz)
    static NetworkSet from_s_dict(const std::map<std::string, std::vector<std::complex<double>>> &d,
                                  const std::vector<Frequency> &frequency) {
        NetworkSet out;
        for(const auto &kv : d) {
            Network n;
            n.freqs = frequency;
            size_t nfreq = frequency.size();
            // infer ports from vector length
            size_t flat_len = kv.second.size();
            int nports = static_cast<int>(std::round(std::sqrt(flat_len)));
            if(static_cast<size_t>(nports * nports) != flat_len) throw std::runtime_error("from_s_dict: invalid flat size");
            n.n_ports = nports;
            n.sparams.resize(nfreq);
            // assume kv.second is arranged as (f0 flat, f1 flat, ...)
            if(flat_len % nfreq != 0) throw std::runtime_error("from_s_dict: cannot chunk data by frequency");
            size_t per_f = flat_len / nfreq;
            if(per_f != static_cast<size_t>(nports * nports)) throw std::runtime_error("from_s_dict: unexpected per-frequency length");
            for(size_t fi=0; fi<nfreq; ++fi) {
                std::vector<std::complex<double>> flat(per_f);
                for(size_t k=0;k<per_f;++k) flat[k] = kv.second[fi*per_f + k];
                n.sparams[fi] = std::move(flat);
            }
            out.ntwk_set.push_back(std::move(n));
        }
        return out;
    }

    // Create from MDIF file (uses C++ mdif reader)
    static NetworkSet from_mdif(const std::string &path) {
        NetworkSet out;
        auto nets = skrf_cpp::io::read_mdif(path);
        out.ntwk_set = std::move(nets);
        return out;
    }

    // Create from CITI file
    static NetworkSet from_citi(const std::string &path) {
        NetworkSet out;
        auto nets = skrf_cpp::io::read_citi(path);
        out.ntwk_set = std::move(nets);
        return out;
    }

    // Apply a function to each Network and collect results
    NetworkSet apply_each(const std::function<Network(const Network&)> &fn) const {
        NetworkSet out;
        out.ntwk_set.reserve(ntwk_set.size());
        for(const auto &n : ntwk_set) out.ntwk_set.push_back(fn(n));
        return out;
    }

    // Return scalar matrix (freq, obs, port_flat) similar to Python scalar_mat
    std::vector<std::vector<std::vector<std::complex<double>>>> scalar_mat(const std::string &param = "s") const {
        // returns 3D vector: [freq][obs][flat_ports]
        if(ntwk_set.empty()) return {};
        size_t nfreq = ntwk_set.front().freqs.size();
        size_t napp = ntwk_set.size();
        size_t flat_len = ntwk_set.front().sparams.front().size();
        std::vector<std::vector<std::vector<std::complex<double>>>> out(nfreq, std::vector<std::vector<std::complex<double>>>(napp, std::vector<std::complex<double>>(flat_len)));
        for(size_t fi=0; fi<nfreq; ++fi) {
            for(size_t ni=0; ni<napp; ++ni) out[fi][ni] = ntwk_set[ni].sparams[fi];
        }
        return out;
    }

    // covariance over observations per frequency
    std::vector<std::vector<std::vector<double>>> cov() const {
        // returns per-frequency covariance (flat_len x flat_len) as doubles
        auto sm = scalar_mat();
        if(sm.empty()) return {};
        size_t nfreq = sm.size();
        size_t napp = sm.front().size();
        size_t dim = sm.front().front().size();
        std::vector<std::vector<std::vector<double>>> out(nfreq, std::vector<std::vector<double>>(dim, std::vector<double>(dim, 0.0)));
        for(size_t fi=0; fi<nfreq; ++fi) {
            // build matrix observations x dim (real only using magnitude)
            std::vector<std::vector<double>> mat(napp, std::vector<double>(dim));
            for(size_t ni=0; ni<napp; ++ni) for(size_t k=0;k<dim;++k) mat[ni][k] = std::abs(sm[fi][ni][k]);
            // compute covariance (naive)
            for(size_t i=0;i<dim;++i) for(size_t j=0;j<dim;++j) {
                double mean_i=0, mean_j=0;
                for(size_t r=0;r<napp;++r) { mean_i += mat[r][i]; mean_j += mat[r][j]; }
                mean_i /= static_cast<double>(napp); mean_j /= static_cast<double>(napp);
                double covv = 0.0;
                for(size_t r=0;r<napp;++r) covv += (mat[r][i]-mean_i)*(mat[r][j]-mean_j);
                covv /= static_cast<double>(std::max<size_t>(1, napp-1));
                out[fi][i][j] = covv;
            }
        }
        return out;
    }

    Network mean_s_db() const {
        Network mean = mean_s();
        for(auto &flat : mean.sparams) {
            for(auto &c : flat) {
                double mag = std::abs(c);
                double db = 20.0 * std::log10(std::max(1e-300, mag));
                c = std::complex<double>(db, 0.0);
            }
        }
        return mean;
    }

    Network std_s_db() const {
        Network s = std_s();
        for(auto &flat : s.sparams) {
            for(auto &c : flat) {
                double mag = std::abs(c);
                double db = 20.0 * std::log10(std::max(1e-300, mag));
                c = std::complex<double>(db, 0.0);
            }
        }
        return s;
    }

    NetworkSet inv() const {
        NetworkSet out;
        out.ntwk_set.reserve(ntwk_set.size());
        for(const auto &n : ntwk_set) {
            NetworkEigen ne = n.to_network_eigen();
            // invert S by converting to Z then inverting (approx)
            NetworkEigen ze = ne.to_z(ne.z0);
            for(auto &Z : ze.s_params) Z = Z.inverse();
            out.ntwk_set.push_back(ze.to_flat_network());
        }
        return out;
    }

    // Interpolate each network to new frequency grid (freqs_hz in Hz)
    NetworkSet interpolate_frequency(const std::vector<double> &freqs_hz, const std::string &basis = "s") const {
        NetworkSet out;
        out.ntwk_set.reserve(ntwk_set.size());
        for(const auto &n : ntwk_set) out.ntwk_set.push_back(n.resample(freqs_hz));
        return out;
    }

    // Interpolate across networks: ntw_param length must equal number of networks.
    Network interpolate_from_network(const std::vector<double> &ntw_param, double x, const std::string &interp_kind = "linear") const {
        if(ntwk_set.empty()) return Network();
        if(ntw_param.size() != ntwk_set.size()) throw std::runtime_error("interpolate_from_network: parameter length mismatch");
        // simple linear interpolation per-frequency per-element
        const Network &ref = ntwk_set.front();
        size_t nfreq = ref.freqs.size();
        int nports = ref.n_ports;
        Network out;
        out.n_ports = nports;
        out.freqs = ref.freqs;
        out.sparams.resize(nfreq);
        for(size_t fi=0; fi<nfreq; ++fi) {
            std::vector<std::complex<double>> vals(ntwk_set.size()*nports*nports);
            for(size_t ni=0; ni<ntwk_set.size(); ++ni) {
                const auto &flat = ntwk_set[ni].sparams[fi];
                for(size_t k=0;k<flat.size();++k) vals[ni*flat.size()+k] = flat[k];
            }
            // for each component, interpolate simple linear between nearest
            std::vector<std::complex<double>> outflat(nports*nports);
            for(size_t k=0;k<nports*nports;++k) {
                // build vector of samples
                std::vector<std::pair<double,std::complex<double>>> samples;
                for(size_t ni=0; ni<ntwk_set.size(); ++ni) samples.emplace_back(ntw_param[ni], vals[ni*(nports*nports)+k]);
                // find segment
                size_t idx = 1;
                while(idx < samples.size() && samples[idx].first < x) ++idx;
                if(idx==0) outflat[k] = samples.front().second;
                else if(idx>=samples.size()) outflat[k] = samples.back().second;
                else {
                    double x0 = samples[idx-1].first; double x1 = samples[idx].first;
                    double t = (x - x0) / (x1 - x0);
                    outflat[k] = samples[idx-1].second * (1.0 - t) + samples[idx].second * t;
                }
            }
            out.sparams[fi] = std::move(outflat);
        }
        return out;
    }

    // Interpolate based on named params not implemented in this C++ port
    Network interpolate_from_params(const std::string &param, double x, const std::map<std::string, double> &sub_params = {}) const {
        throw std::runtime_error("interpolate_from_params not supported in this C++ port (no params stored)");
    }

    // Merge networks into a single block-diagonal Network (ports concatenated)
    Network merge_block_diag(const std::vector<size_t> &indices) const {
        if(indices.empty()) return Network();
        // validate indices
        for(size_t id : indices) if(id >= ntwk_set.size()) throw std::out_of_range("merge_block_diag index out of range");
        // use first network's frequency vector as reference
        const Network &ref = ntwk_set[indices.front()];
        size_t nfreq = ref.freqs.size();
        for(size_t id : indices) {
            if(ntwk_set[id].freqs.size() != nfreq) throw std::runtime_error("frequency point mismatch in merge_block_diag");
        }
        // compute total ports
        int total_ports = 0;
        std::vector<int> ports_per;
        for(size_t id : indices) { ports_per.push_back(ntwk_set[id].n_ports); total_ports += ntwk_set[id].n_ports; }

        Network out;
        out.n_ports = total_ports;
        out.freqs = ref.freqs;
        out.sparams.resize(nfreq);

        for(size_t fi=0; fi<nfreq; ++fi) {
            std::vector<std::complex<double>> flat(total_ports * total_ports, std::complex<double>(0.0,0.0));
            int offset = 0;
            for(size_t k=0;k<indices.size();++k) {
                const auto &src = ntwk_set[indices[k]].sparams[fi];
                int p = ports_per[k];
                for(int r=0;r<p;++r) for(int c=0;c<p;++c) {
                    flat[(offset + r) * total_ports + (offset + c)] = src[r * p + c];
                }
                offset += p;
            }
            out.sparams[fi] = std::move(flat);
        }
        return out;
    }

    // Cascade a list of 2-port networks (by indices) into one 2-port Network
    Network cascade_indices(const std::vector<size_t> &indices) const {
        if(indices.empty()) return Network();
        std::vector<NetworkEigen> list;
        for(size_t id : indices) {
            if(id >= ntwk_set.size()) throw std::out_of_range("cascade_indices index out of range");
            list.push_back(ntwk_set[id].to_network_eigen());
        }
        auto result = NetworkEigen::cascade_list(list, list.front().z0);
        return result.to_flat_network();
    }

    // Parallel combine networks that have the same port count and frequency grid
    Network parallel_indices(const std::vector<size_t> &indices) const {
        if(indices.empty()) return Network();
        std::vector<NetworkEigen> list;
        for(size_t id : indices) {
            if(id >= ntwk_set.size()) throw std::out_of_range("parallel_indices index out of range");
            list.push_back(ntwk_set[id].to_network_eigen());
        }
        // reduce by pairwise parallel
        NetworkEigen acc = list.front();
        for(size_t i=1;i<list.size();++i) acc = acc.parallel_with(list[i], acc.z0);
        return acc.to_flat_network();
    }

    // Connect/merge port pairs inside a given network (by index) and replace that network
    // pairs are given as vector of (i,j) on the current port numbering
    void connect_pairs_in_network(size_t net_index, const std::vector<std::pair<int,int>> &pairs) {
        if(net_index >= ntwk_set.size()) throw std::out_of_range("connect_pairs_in_network index");
        Network &n = ntwk_set[net_index];
        Network merged = n.merge_ports(pairs);
        ntwk_set[net_index] = std::move(merged);
    }

    // More general: select several networks (by indices), merge them block-diagonally,
    // then perform cross-network connections described as groups. Each group is a vector
    // of (local_net_idx, port_idx) pairs where local_net_idx is the index into the
    // provided 'indices' vector (0-based) and port_idx is the port number within that network.
    // All ports in a group will be collapsed into a single port (via sequential merging).
    // The merged result replaces the first network in 'indices' and the other networks
    // are removed from the set.
    void connect_and_merge_indices(const std::vector<size_t> &indices,
                                   const std::vector<std::vector<std::pair<int,int>>> &groups) {
        if(indices.empty()) return;
        // validate indices
        for(size_t id : indices) if(id >= ntwk_set.size()) throw std::out_of_range("connect_and_merge_indices index out of range");

        // create block-diagonal merged network
        Network merged = merge_block_diag(indices);

        // compute offsets for each selected network
        std::vector<int> offsets;
        offsets.reserve(indices.size());
        int acc = 0;
        for(size_t id : indices) {
            offsets.push_back(acc);
            acc += ntwk_set[id].n_ports;
        }

        // build global pairs for all groups
        std::vector<std::pair<int,int>> global_pairs;
        for(const auto &grp : groups) {
            if(grp.empty()) continue;
            // map to global port indices
            std::vector<int> globals;
            globals.reserve(grp.size());
            for(const auto &lp : grp) {
                int local_net_idx = lp.first;
                int port_idx = lp.second;
                if(local_net_idx < 0 || static_cast<size_t>(local_net_idx) >= indices.size()) throw std::out_of_range("group local_net_idx out of range");
                size_t sel_net = indices[local_net_idx];
                if(port_idx < 0 || port_idx >= ntwk_set[sel_net].n_ports) throw std::out_of_range("group port_idx out of range");
                int g = offsets[local_net_idx] + port_idx;
                globals.push_back(g);
            }
            if(globals.size() <= 1) continue;
            int base = globals.front();
            for(size_t k=1;k<globals.size();++k) {
                if(globals[k] == base) continue;
                global_pairs.emplace_back(base, globals[k]);
            }
        }

        // apply merges on the merged network
        if(!global_pairs.empty()) {
            Network after = merged.merge_ports(global_pairs);
            // replace first index with result and erase the others
            size_t first_idx = indices.front();
            ntwk_set[first_idx] = std::move(after);
            // erase remaining indices (descending order to keep positions valid)
            std::vector<size_t> to_erase(indices.begin()+1, indices.end());
            std::sort(to_erase.begin(), to_erase.end(), std::greater<size_t>());
            for(size_t e : to_erase) ntwk_set.erase(ntwk_set.begin() + e);
        } else {
            // no merges to do: still replace first with block-diag and remove others
            size_t first_idx = indices.front();
            ntwk_set[first_idx] = std::move(merged);
            std::vector<size_t> to_erase(indices.begin()+1, indices.end());
            std::sort(to_erase.begin(), to_erase.end(), std::greater<size_t>());
            for(size_t e : to_erase) ntwk_set.erase(ntwk_set.begin() + e);
        }
    }
};

} // namespace skrf_cpp
