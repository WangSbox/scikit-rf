#pragma once

#include <vector>
#include <string>
#include <functional>
#include <filesystem>
#include <sstream>
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
    // Optional named parameters associated with each network. If non-empty,
    // params[i] holds a map of parameter-name -> value for ntwk_set[i].
    std::vector<std::map<std::string, double>> params;

    NetworkSet() = default;
    explicit NetworkSet(const std::vector<Network>& v): ntwk_set(v) {}

    // convenient constructor that also accepts parameter maps per-network
    NetworkSet(const std::vector<Network>& v, const std::vector<std::map<std::string,double>> &p): ntwk_set(v), params(p) {}

    size_t size() const { return ntwk_set.size(); }

    void append(const Network &n) { ntwk_set.push_back(n); }

    // append with optional params map
    void append(const Network &n, const std::map<std::string,double> &p) { ntwk_set.push_back(n); params.push_back(p); }

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
        // Implementation: look up networks that have the named parameter and
        // perform a linear interpolation across that parameter using
        // `interpolate_from_network` helper. If no networks contain the
        // parameter, raise with guidance.
        if(ntwk_set.empty()) return Network();
        // collect samples where param exists
        std::vector<double> samples_val;
        std::vector<size_t> samples_idx;
        for(size_t i=0;i<params.size() && i<ntwk_set.size(); ++i) {
            auto it = params[i].find(param);
            if(it != params[i].end()) { samples_val.push_back(it->second); samples_idx.push_back(i); }
        }
        if(samples_idx.empty()) {
            throw std::runtime_error(std::string("interpolate_from_params: parameter '") + param + "' not found in any stored network params; provide params when constructing NetworkSet or use interpolate_from_network.");
        }

        // build ntw_param vector and subset networks
        std::vector<double> ntw_param;
        std::vector<Network> subset;
        ntw_param.reserve(samples_idx.size()); subset.reserve(samples_idx.size());
        for(size_t k=0;k<samples_idx.size();++k) {
            ntw_param.push_back(samples_val[k]);
            subset.push_back(ntwk_set[samples_idx[k]]);
        }

        // If the provided x is outside sample range, clamp to min/max (linear extrapolation not supported)
        double xmin = *std::min_element(ntw_param.begin(), ntw_param.end());
        double xmax = *std::max_element(ntw_param.begin(), ntw_param.end());
        double xclamped = x;
        if(x < xmin) xclamped = xmin;
        if(x > xmax) xclamped = xmax;

        // use existing interpolate_from_network by delegating to a temporary NetworkSet
        NetworkSet tmp(subset);
        return tmp.interpolate_from_network(ntw_param, xclamped, "linear");
    }

    // Multi-parameter interpolation (separable, iterative along params order).
    // `param_names` and `xvals` must have same length. `sub_params` are fixed
    // filters applied to select candidate networks (must match exactly).
    Network interpolate_from_params(const std::vector<std::string> &param_names,
                                    const std::vector<double> &xvals,
                                    const std::map<std::string, double> &sub_params = {}) const {
        if(param_names.size() != xvals.size()) throw std::runtime_error("interpolate_from_params: param_names and xvals size mismatch");
        if(ntwk_set.empty()) return Network();

        // ensure networks share a common frequency grid for meaningful interpolation
        // create a mutable copy and align frequencies (intersection) before interpolating
        NetworkSet tmp = *this;
        tmp.align_frequency_grid("intersection");

        // initial candidate indices: those matching sub_params (if any)
        std::vector<size_t> candidates;
        for(size_t i=0;i<tmp.ntwk_set.size();++i) {
            bool ok = true;
            for(const auto &kv : sub_params) {
                if(i >= params.size()) { ok = false; break; }
                auto it = params[i].find(kv.first);
                if(it == params[i].end() || it->second != kv.second) { ok = false; break; }
            }
            if(ok) candidates.push_back(i);
        }
        if(candidates.empty()) throw std::runtime_error("interpolate_from_params: no networks match sub_params filters");

        // Start with subset networks corresponding to candidates (from tmp aligned set)
        std::vector<Network> current_set;
        std::vector<std::map<std::string,double>> current_params;
        for(size_t idx : candidates) {
            current_set.push_back(tmp.ntwk_set[idx]);
            if(idx < tmp.params.size()) current_params.push_back(tmp.params[idx]); else current_params.emplace_back();
        }

        // Iteratively interpolate along each parameter
        for(size_t pi=0; pi<param_names.size(); ++pi) {
            const std::string &pname = param_names[pi];
            double x = xvals[pi];

            // build samples for this parameter grouped by other params values
            // grouping key: string of other params
            std::map<std::string, std::vector<std::pair<double, Network>>> groups;
            for(size_t si=0; si<current_set.size(); ++si) {
                const auto &pm = current_params[si];
                auto it = pm.find(pname);
                if(it == pm.end()) continue; // skip networks without this param
                // build key from other parameters (excluding pname)
                std::ostringstream keyss;
                for(const auto &kv : pm) {
                    if(kv.first == pname) continue;
                    keyss << kv.first << "=" << kv.second << ";";
                }
                groups[keyss.str()].emplace_back(it->second, current_set[si]);
            }

            if(groups.empty()) throw std::runtime_error(std::string("interpolate_from_params: parameter '") + pname + " not present in any candidate");

            // For each group, perform 1D interpolation along this parameter
            std::vector<Network> next_set;
            std::vector<std::map<std::string,double>> next_params;
            for(const auto &grp : groups) {
                const auto &samples = grp.second;
                if(samples.empty()) continue;
                // build temporary NetworkSet from samples and param vector
                std::vector<double> sample_vals;
                std::vector<Network> sample_nets;
                for(const auto &s : samples) { sample_vals.push_back(s.first); sample_nets.push_back(s.second); }
                NetworkSet tmp(sample_nets);
                // clamp x to sample range
                double xmin = *std::min_element(sample_vals.begin(), sample_vals.end());
                double xmax = *std::max_element(sample_vals.begin(), sample_vals.end());
                double xcl = x;
                if(xcl < xmin) xcl = xmin; if(xcl > xmax) xcl = xmax;
                Network out = tmp.interpolate_from_network(sample_vals, xcl, "linear");
                next_set.push_back(out);
                // build representative params map: set pname to xcl and keep other keys from first sample
                std::map<std::string,double> rep = current_params[0];
                rep[pname] = xcl;
                next_params.push_back(rep);
            }

            // replace current set
            current_set = std::move(next_set);
            current_params = std::move(next_params);
        }

        // After all parameters interpolated, if multiple resulting networks exist,
        // return the first (users can reduce groups via sub_params to get unique result).
        if(current_set.empty()) return Network();
        return current_set.front();
    }

    // Align frequency grids of all networks in-place.
    // mode: "union" (all unique points) or "intersection" (points present in all networks within tol)
    void align_frequency_grid(const std::string &mode = "union", double tol = 1e-6) {
        if(ntwk_set.empty()) return;
        // collect all frequency Hz arrays
        std::vector<std::vector<double>> all_freqs;
        all_freqs.reserve(ntwk_set.size());
        for(const auto &n : ntwk_set) {
            std::vector<double> v; v.reserve(n.freqs.size());
            for(const auto &f : n.freqs) v.push_back(f.hz);
            all_freqs.push_back(std::move(v));
        }

        std::vector<double> target;
        if(mode == "union") {
            std::vector<double> tmp;
            for(const auto &v : all_freqs) for(double x : v) tmp.push_back(x);
            std::sort(tmp.begin(), tmp.end());
            // unique with tolerance
            for(double x : tmp) {
                if(target.empty() || std::abs(target.back() - x) > tol) target.push_back(x);
            }
        } else { // intersection
            // start from first list and retain points that exist in all others within tol
            target = all_freqs.front();
            std::sort(target.begin(), target.end());
            std::vector<double> keep;
            for(double x : target) {
                bool present = true;
                for(size_t i=1;i<all_freqs.size();++i) {
                    const auto &v = all_freqs[i];
                    bool found = false;
                    for(double y : v) { if(std::abs(y - x) <= tol) { found = true; break; } }
                    if(!found) { present = false; break; }
                }
                if(present) keep.push_back(x);
            }
            target = std::move(keep);
        }

        if(target.empty()) throw std::runtime_error("align_frequency_grid resulted in empty target grid");

        // resample each network to target
        for(auto &n : ntwk_set) {
            std::vector<double> hz; hz.reserve(target.size());
            for(double x : target) hz.push_back(x);
            n = n.resample(hz);
        }
    }

    // Numeric check helper: verify that cascading a list of 2-port networks via
    // `cascade_list` gives same result as doing block-diagonal merge then
    // connecting ports via `connect_and_merge_indices` (used as equivalence test).
    // Returns relative max error across S elements.
    double numeric_verify_cascade_equivalence(const std::vector<size_t> &indices) const {
        if(indices.empty()) return 0.0;
        // check all indices valid and networks are 2-port
        for(size_t id : indices) if(id >= ntwk_set.size()) throw std::out_of_range("numeric_verify_cascade_equivalence index out of range");
        std::vector<NetworkEigen> list;
        for(size_t id : indices) list.push_back(ntwk_set[id].to_network_eigen());
        // method A: cascade_list (existing optimized 2-port method)
        NetworkEigen A = NetworkEigen::cascade_list(list, list.front().z0);

        // method B: build a temporary NetworkSet with the listed networks in order,
        // then perform connect_and_merge_indices to collapse the cascade into one network.
        NetworkSet tmp;
        // copy the networks indicated by indices into tmp in the same order
        for(size_t id : indices) tmp.ntwk_set.push_back(ntwk_set[id]);
        // build groups: each connection between net k port 1 and net k+1 port 0 becomes a group
        std::vector<std::vector<std::pair<int,int>>> groups;
        for(size_t k=0;k+1<indices.size();++k) {
            std::vector<std::pair<int,int>> grp;
            grp.emplace_back(static_cast<int>(k), 1);   // (local_net_idx, port_idx)
            grp.emplace_back(static_cast<int>(k+1), 0);
            groups.push_back(std::move(grp));
        }
        // indices into tmp are 0..n-1
        std::vector<size_t> tmp_indices(indices.size());
        for(size_t i=0;i<tmp_indices.size();++i) tmp_indices[i] = i;
        // perform connect_and_merge_indices on tmp
        tmp.connect_and_merge_indices(tmp_indices, groups);
        // after operation the merged network is at tmp.ntwk_set[tmp_indices.front()]
        NetworkEigen B = tmp.ntwk_set.empty() ? NetworkEigen() : tmp.ntwk_set[tmp_indices.front()].to_network_eigen();

        // Now compare A and B on frequency grid and S elements (only when sizes match)
        if(A.n_ports != B.n_ports || A.freqs.size() != B.freqs.size()) {
            return std::numeric_limits<double>::infinity();
        }
        double max_rel_err = 0.0;
        for(size_t fi=0; fi<A.freqs.size(); ++fi) {
            const MatrixXcd &Sa = A.s_params[fi];
            const MatrixXcd &Sb = B.s_params[fi];
            for(int r=0;r<A.n_ports;++r) for(int c=0;c<A.n_ports;++c) {
                double na = std::abs(Sa(r,c));
                double nb = std::abs(Sb(r,c));
                double denom = std::max(1e-12, na);
                double rel = std::abs(na - nb) / denom;
                if(rel > max_rel_err) max_rel_err = rel;
            }
        }
        return max_rel_err;
    }

    // Multilinear interpolation for arbitrary D (requires samples on a complete Cartesian grid)
    // The implementation finds for each dimension the surrounding grid points
    // and performs multilinear interpolation using the 2^D corner values.
    Network interpolate_from_params_multilinear(const std::vector<std::string> &param_names,
                                               const std::vector<double> &xvals,
                                               const std::map<std::string,double> &sub_params = {}) const {
        if(param_names.size() != xvals.size()) throw std::runtime_error("interpolate_from_params_multilinear: sizes mismatch");
        size_t D = param_names.size();
        if(D == 0) throw std::runtime_error("interpolate_from_params_multilinear: empty param_names");
        if(ntwk_set.empty()) return Network();

        // collect candidate samples matching sub_params
        struct Sample { std::map<std::string,double> pm; Network net; };
        std::vector<Sample> samples;
        for(size_t i=0;i<ntwk_set.size();++i) {
            bool ok = true;
            for(const auto &kv: sub_params) {
                if(i >= params.size()) { ok = false; break; }
                auto it = params[i].find(kv.first);
                if(it == params[i].end() || it->second != kv.second) { ok = false; break; }
            }
            if(!ok) continue;
            Sample s; if(i < params.size()) s.pm = params[i]; s.net = ntwk_set[i]; samples.push_back(std::move(s));
        }
        if(samples.empty()) throw std::runtime_error("interpolate_from_params_multilinear: no candidate networks after applying sub_params");

        // build sorted unique grid vectors for each dimension
        std::vector<std::vector<double>> grid_vals(D);
        for(const auto &s : samples) {
            for(size_t d=0; d<D; ++d) {
                const std::string &p = param_names[d];
                auto it = s.pm.find(p);
                if(it == s.pm.end()) continue; // sample may not have this param
                grid_vals[d].push_back(it->second);
            }
        }
        for(size_t d=0; d<D; ++d) {
            auto &gv = grid_vals[d];
            if(gv.empty()) throw std::runtime_error(std::string("interpolate_from_params_multilinear: parameter '") + param_names[d] + " not present in any candidate");
            std::sort(gv.begin(), gv.end());
            gv.erase(std::unique(gv.begin(), gv.end()), gv.end());
        }

        // verify complete Cartesian grid: product of sizes should equal number of samples that have all params
        size_t expected = 1;
        for(size_t d=0; d<D; ++d) expected *= grid_vals[d].size();
        // count samples that have all params
        size_t full_samples = 0;
        for(const auto &s : samples) {
            bool ok = true;
            for(size_t d=0; d<D; ++d) if(s.pm.find(param_names[d])==s.pm.end()) { ok = false; break; }
            if(ok) ++full_samples;
        }
        if(full_samples != expected) throw std::runtime_error("interpolate_from_params_multilinear: samples do not form a complete Cartesian grid for requested params");

        // map from coordinate key to Network
        auto make_key = [&](const std::vector<double> &coords)->std::string{
            std::ostringstream os; os.setf(std::ios::fmtflags(0), std::ios::floatfield);
            for(size_t i=0;i<coords.size();++i) {
                if(i) os << '|'; os << std::setprecision(17) << coords[i];
            }
            return os.str();
        };
        std::map<std::string, Network> grid_map;
        for(const auto &s : samples) {
            std::vector<double> coords; coords.reserve(D);
            bool ok = true;
            for(size_t d=0; d<D; ++d) {
                auto it = s.pm.find(param_names[d]); if(it==s.pm.end()) { ok=false; break; } coords.push_back(it->second);
            }
            if(!ok) continue;
            grid_map[make_key(coords)] = s.net;
        }

        // for each dim, find lower/upper indices and t
        std::vector<size_t> li(D), ui(D);
        std::vector<double> t(D);
        for(size_t d=0; d<D; ++d) {
            const auto &vec = grid_vals[d]; double x = xvals[d];
            if(x <= vec.front()) { li[d]=ui[d]=0; t[d]=0.0; }
            else if(x >= vec.back()) { li[d]=ui[d]=vec.size()-1; t[d]=0.0; }
            else {
                size_t idx = 1; while(idx < vec.size() && vec[idx] < x) ++idx; li[d]=idx-1; ui[d]=idx;
                double v0 = vec[li[d]]; double v1 = vec[ui[d]]; t[d] = (v1==v0)?0.0:((x - v0)/(v1 - v0));
            }
        }

        // choose reference network to validate shape
        Network ref = grid_map.begin()->second;
        size_t nfreq = ref.freqs.size(); int nports = ref.n_ports;

        // prepare output network
        Network out; out.n_ports = nports; out.freqs = ref.freqs; out.sparams.assign(nfreq, std::vector<std::complex<double>>(nports*nports, {0,0}));

        // iterate over all 2^D corners
        size_t corners = 1ull << D;
        std::vector<double> corner_coords(D);
        for(size_t c=0; c<corners; ++c) {
            // build corner coordinate values and compute weight
            double weight = 1.0;
            for(size_t d=0; d<D; ++d) {
                bool bit = ((c >> d) & 1u);
                size_t idx = bit ? ui[d] : li[d];
                double val = grid_vals[d][idx]; corner_coords[d] = val;
                double td = t[d];
                weight *= (bit ? td : (1.0 - td));
            }
            std::string key = make_key(corner_coords);
            auto git = grid_map.find(key);
            if(git == grid_map.end()) throw std::runtime_error("interpolate_from_params_multilinear: missing corner sample in grid_map");
            const Network &corner_net = git->second;
            if(corner_net.n_ports != static_cast<int>(nports) || corner_net.freqs.size() != nfreq) throw std::runtime_error("interpolate_from_params_multilinear: inconsistent network shapes among grid samples");

            // accumulate weighted corner values
            for(size_t fi=0; fi<nfreq; ++fi) {
                const auto &corner_flat = corner_net.sparams[fi];
                auto &outflat = out.sparams[fi];
                for(size_t k=0;k<outflat.size();++k) outflat[k] += corner_flat[k] * weight;
            }
        }

        return out;
    }

    // Automated n-port cascade helper: given ordered `indices` of networks to
    // cascade and `adjacent_pairs` describing how ports between successive
    // networks are connected. `adjacent_pairs` must have length indices.size()-1,
    // and each element is a vector of pairs (port_in_k, port_in_k+1) referring to
    // port indices within the respective networks (0-based). The helper builds
    // a block-diagonal merge then applies all specified connections (merging
    // corresponding ports) and returns the resulting Network.
    Network cascade_with_connections(const std::vector<size_t> &indices,
                                     const std::vector<std::vector<std::pair<int,int>>> &adjacent_pairs) const {
        if(indices.empty()) return Network();
        if(adjacent_pairs.size() != (indices.size()>=1?indices.size()-1:0)) throw std::runtime_error("cascade_with_connections: adjacent_pairs size mismatch");
        for(size_t id : indices) if(id >= ntwk_set.size()) throw std::out_of_range("cascade_with_connections index out of range");
        // merge block-diagonal
        Network merged = merge_block_diag(indices);
        // compute offsets
        std::vector<int> offsets; offsets.reserve(indices.size()); int acc=0; for(size_t id: indices) { offsets.push_back(acc); acc += ntwk_set[id].n_ports; }
        // build global pairs to merge
        std::vector<std::pair<int,int>> global_pairs;
        for(size_t k=0;k+1<indices.size();++k) {
            const auto &pairs = adjacent_pairs[k];
            for(const auto &pp : pairs) {
                int p_k = pp.first; int p_k1 = pp.second;
                if(p_k < 0 || p_k >= ntwk_set[indices[k]].n_ports) throw std::out_of_range("cascade_with_connections: port index out of range");
                if(p_k1 < 0 || p_k1 >= ntwk_set[indices[k+1]].n_ports) throw std::out_of_range("cascade_with_connections: port index out of range");
                int g1 = offsets[k] + p_k; int g2 = offsets[k+1] + p_k1;
                global_pairs.emplace_back(g1, g2);
            }
        }
        if(global_pairs.empty()) return merged;
        Network out = merged.merge_ports(global_pairs);
        return out;
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
