#pragma once

#include "network_eigen.h"
#include "touchstone.h"
#include "transforms.h"

namespace skrf_cpp {

// Minimal calibration helper that holds an error network (2-port)
// and can de-embed measurements by removing the error network via ABCD math.
class Calibration {
public:
    Calibration() = default;

    // load an error network from a touchstone file
    static Calibration load_from_touchstone(const std::string &path) {
        Calibration cal;
        Network net = Touchstone::load(path);
        cal.error_ = net.to_network_eigen();
        return cal;
    }

    explicit Calibration(const NetworkEigen &err): error_(err) {}

    // De-embed a 2-port Network (measurement) using the stored error network.
    // Returns a new Network with the same freq axis containing the de-embedded S-parameters.
    Network deembed(const Network &meas) const;

    // access to error network
    const NetworkEigen &error_network() const { return error_; }

private:
    NetworkEigen error_; // expected to be 2-port
};

// Deembedding base class and common deembedding methods ported from
// skrf.calibration.deembedding (subset). Implementations operate on
// Matrix-based NetworkEigen objects and return a flattened Network result.
class Deembedding {
public:
    Deembedding() = default;
    explicit Deembedding(const std::vector<NetworkEigen> &dummies)
        : dummies_(dummies) {}

    virtual ~Deembedding() = default;

    // Apply de-embedding to a measurement Network and return corrected Network
    virtual Network deembed(const Network &meas) const = 0;

protected:
    // ensure dummy networks are resampled to measurement frequencies when needed
    static std::vector<NetworkEigen> ensure_resampled(const std::vector<NetworkEigen>& dummies,
                                                      const NetworkEigen &meas_ne) {
        std::vector<NetworkEigen> out = dummies;
        // target freqs
        std::vector<double> target_freqs;
        for(const auto &f: meas_ne.freqs) target_freqs.push_back(f.hz);
        for(auto &d : out) {
            if(d.freqs.size() != meas_ne.freqs.size()) {
                d = d.resample(target_freqs);
            } else {
                // quick value check, if mismatch resample
                bool same = true;
                for(size_t i=0;i<d.freqs.size();++i) {
                    if(std::abs(d.freqs[i].hz - meas_ne.freqs[i].hz) > 1e-9) { same = false; break; }
                }
                if(!same) d = d.resample(target_freqs);
            }
        }
        return out;
    }

    std::vector<NetworkEigen> dummies_;
};

// OpenShort: remove open parasitics followed by short parasitics.
class OpenShort : public Deembedding {
public:
    OpenShort(const Network &dummy_open, const Network &dummy_short) {
        dummies_.push_back(dummy_open.to_network_eigen());
        dummies_.push_back(dummy_short.to_network_eigen());
    }

    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &op = d[0];
        const NetworkEigen &sh = d[1];

        NetworkEigen Y_meas = meas_ne.to_y(meas_ne.z0);
        NetworkEigen Y_op = op.to_y(op.z0);
        NetworkEigen Y_sh = sh.to_y(sh.z0);

        Network out;
        out.n_ports = meas.n_ports;
        out.z0 = meas.z0;
        out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());

        for(size_t i=0;i<Y_meas.s_params.size();++i) {
            MatrixXcd Ydsh = Y_sh.s_params[i] - Y_op.s_params[i];
            MatrixXcd Ycaled = Y_meas.s_params[i] - Y_op.s_params[i];
            MatrixXcd Zdsh = Ydsh.inverse();
            MatrixXcd Zcaled = Ycaled.inverse();
            MatrixXcd Sfinal = z_to_s(Zcaled - Zdsh, meas_ne.z0);
            // flatten
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sfinal(0,0); flat[1]=Sfinal(0,1); flat[2]=Sfinal(1,0); flat[3]=Sfinal(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// Open: remove open (parallel) parasitics only
class Open : public Deembedding {
public:
    explicit Open(const Network &dummy_open) { dummies_.push_back(dummy_open.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &op = d[0];
        NetworkEigen Y_meas = meas_ne.to_y(meas_ne.z0);
        NetworkEigen Y_op = op.to_y(op.z0);

        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());
        for(size_t i=0;i<Y_meas.s_params.size();++i) {
            MatrixXcd Ycaled = Y_meas.s_params[i] - Y_op.s_params[i];
            MatrixXcd Zcaled = Ycaled.inverse();
            MatrixXcd Sfinal = z_to_s(Zcaled, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sfinal(0,0); flat[1]=Sfinal(0,1); flat[2]=Sfinal(1,0); flat[3]=Sfinal(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// ShortOpen: short then open
class ShortOpen : public Deembedding {
public:
    ShortOpen(const Network &dummy_short, const Network &dummy_open) {
        dummies_.push_back(dummy_open.to_network_eigen());
        dummies_.push_back(dummy_short.to_network_eigen());
    }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &op = d[0];
        const NetworkEigen &sh = d[1];

        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());

        NetworkEigen Z_meas = meas_ne.to_z(meas_ne.z0);
        NetworkEigen Z_op = op.to_z(op.z0);
        NetworkEigen Z_sh = sh.to_z(sh.z0);

        for(size_t i=0;i<Z_meas.s_params.size();++i) {
            MatrixXcd Zdeop = Z_op.s_params[i] - Z_sh.s_params[i];
            MatrixXcd Ydeop = Zdeop.inverse();
            MatrixXcd Zcaled = Z_meas.s_params[i] - Z_sh.s_params[i];
            MatrixXcd Ycaled = Zcaled.inverse();
            MatrixXcd Sfinal = z_to_s(Ycaled.inverse(), meas_ne.z0); // Ycaled -> Z -> S
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sfinal(0,0); flat[1]=Sfinal(0,1); flat[2]=Sfinal(1,0); flat[3]=Sfinal(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// Short: remove series parasitics only
class Short : public Deembedding {
public:
    explicit Short(const Network &dummy_short) { dummies_.push_back(dummy_short.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &sh = d[0];
        NetworkEigen Z_meas = meas_ne.to_z(meas_ne.z0);
        NetworkEigen Z_sh = sh.to_z(sh.z0);
        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());
        for(size_t i=0;i<Z_meas.s_params.size();++i) {
            MatrixXcd Zcaled = Z_meas.s_params[i] - Z_sh.s_params[i];
            MatrixXcd Sfinal = z_to_s(Zcaled, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sfinal(0,0); flat[1]=Sfinal(0,1); flat[2]=Sfinal(1,0); flat[3]=Sfinal(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// SplitPi: see Python implementation; builds left/right from thru Y and applies
class SplitPi : public Deembedding {
public:
    explicit SplitPi(const Network &dummy_thru) { dummies_.push_back(dummy_thru.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &thru = d[0];

        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());

        NetworkEigen Y_thru = thru.to_y(thru.z0);
        for(size_t i=0;i<Y_thru.s_params.size();++i) {
            MatrixXcd Yt = Y_thru.s_params[i];
            MatrixXcd left_y = MatrixXcd::Zero(2,2);
            left_y(0,0) = (Yt(0,0) - Yt(1,0) + Yt(1,1) - Yt(0,1)) / 2.0;
            left_y(0,1) = Yt(1,0) + Yt(0,1);
            left_y(1,0) = left_y(0,1);
            left_y(1,1) = - Yt(1,0) - Yt(0,1);
            MatrixXcd Zleft = left_y.inverse();
            MatrixXcd Sleft = z_to_s(Zleft, meas_ne.z0);
            // right is flipped
            MatrixXcd Sright = MatrixXcd::Zero(2,2);
            Sright(0,0)=Sleft(1,1); Sright(0,1)=Sleft(1,0);
            Sright(1,0)=Sleft(0,1); Sright(1,1)=Sleft(0,0);

            Eigen::Matrix2cd Aleft = s_to_abcd_2port(Sleft, meas_ne.z0);
            Eigen::Matrix2cd Aright = s_to_abcd_2port(Sright, meas_ne.z0);
            Eigen::Matrix2cd Ameas = s_to_abcd_2port(meas_ne.s_params[i].block<2,2>(0,0), meas_ne.z0);

            Eigen::Matrix2cd A_corr = Aleft.inverse() * Ameas * Aright.inverse();
            Eigen::Matrix2cd S_corr = abcd_to_s_2port(A_corr, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=S_corr(0,0); flat[1]=S_corr(0,1); flat[2]=S_corr(1,0); flat[3]=S_corr(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// SplitTee: similar to SplitPi but builds left_z from thru.z
class SplitTee : public Deembedding {
public:
    explicit SplitTee(const Network &dummy_thru) { dummies_.push_back(dummy_thru.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &thru = d[0];

        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());

        NetworkEigen Z_thru = thru.to_z(thru.z0);
        for(size_t i=0;i<Z_thru.s_params.size();++i) {
            MatrixXcd Zt = Z_thru.s_params[i];
            MatrixXcd left_z = MatrixXcd::Zero(2,2);
            left_z(0,0) = (Zt(0,0) + Zt(1,0) + Zt(1,1) + Zt(0,1)) / 2.0;
            left_z(0,1) = Zt(1,0) + Zt(0,1);
            left_z(1,0) = left_z(0,1);
            left_z(1,1) = Zt(1,0) + Zt(0,1);
            MatrixXcd Sleft = z_to_s(left_z, meas_ne.z0);
            MatrixXcd Sright = MatrixXcd::Zero(2,2);
            Sright(0,0)=Sleft(1,1); Sright(0,1)=Sleft(1,0);
            Sright(1,0)=Sleft(0,1); Sright(1,1)=Sleft(0,0);

            Eigen::Matrix2cd Aleft = s_to_abcd_2port(Sleft, meas_ne.z0);
            Eigen::Matrix2cd Aright = s_to_abcd_2port(Sright, meas_ne.z0);
            Eigen::Matrix2cd Ameas = s_to_abcd_2port(meas_ne.s_params[i].block<2,2>(0,0), meas_ne.z0);

            Eigen::Matrix2cd A_corr = Aleft.inverse() * Ameas * Aright.inverse();
            Eigen::Matrix2cd S_corr = abcd_to_s_2port(A_corr, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=S_corr(0,0); flat[1]=S_corr(0,1); flat[2]=S_corr(1,0); flat[3]=S_corr(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// AdmittanceCancel and ImpedanceCancel: left-right mirroring
class AdmittanceCancel : public Deembedding {
public:
    explicit AdmittanceCancel(const Network &dummy_thru) { dummies_.push_back(dummy_thru.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &thru = d[0];
        // compute h = ntwk ** thru.inv  -> in ABCD domain
        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());
        for(size_t i=0;i<meas_ne.s_params.size();++i) {
            Eigen::Matrix2cd Smeas = meas_ne.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Athru = s_to_abcd_2port(thru.s_params[i].block<2,2>(0,0), thru.z0);
            Eigen::Matrix2cd Ameas = s_to_abcd_2port(Smeas, meas_ne.z0);
            Eigen::Matrix2cd h = Ameas * Athru.inverse();
            // h_ is flipped (swap ports)
            Eigen::Matrix2cd S_h = abcd_to_s_2port(h, meas_ne.z0);
            Eigen::Matrix2cd S_h_flipped;
            S_h_flipped(0,0)=S_h(1,1); S_h_flipped(0,1)=S_h(1,0); S_h_flipped(1,0)=S_h(0,1); S_h_flipped(1,1)=S_h(0,0);
            // compute average in Y domain
            MatrixXcd Yh = s_to_z(S_h, meas_ne.z0).inverse();
            MatrixXcd Yhf = s_to_z(S_h_flipped, meas_ne.z0).inverse();
            MatrixXcd Yavg = (Yh + Yhf) / 2.0;
            MatrixXcd Zout = Yavg.inverse();
            MatrixXcd Sout = z_to_s(Zout, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sout(0,0); flat[1]=Sout(0,1); flat[2]=Sout(1,0); flat[3]=Sout(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

class ImpedanceCancel : public Deembedding {
public:
    explicit ImpedanceCancel(const Network &dummy_thru) { dummies_.push_back(dummy_thru.to_network_eigen()); }
    Network deembed(const Network &meas) const override {
        NetworkEigen meas_ne = meas.to_network_eigen();
        auto d = ensure_resampled(dummies_, meas_ne);
        const NetworkEigen &thru = d[0];
        Network out; out.n_ports = meas.n_ports; out.z0 = meas.z0; out.freqs = meas.freqs;
        out.sparams.reserve(meas.sparams.size());
        for(size_t i=0;i<meas_ne.s_params.size();++i) {
            Eigen::Matrix2cd Smeas = meas_ne.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Athru = s_to_abcd_2port(thru.s_params[i].block<2,2>(0,0), thru.z0);
            Eigen::Matrix2cd Ameas = s_to_abcd_2port(Smeas, meas_ne.z0);
            Eigen::Matrix2cd h = Ameas * Athru.inverse();
            Eigen::Matrix2cd S_h = abcd_to_s_2port(h, meas_ne.z0);
            Eigen::Matrix2cd S_h_flipped;
            S_h_flipped(0,0)=S_h(1,1); S_h_flipped(0,1)=S_h(1,0); S_h_flipped(1,0)=S_h(0,1); S_h_flipped(1,1)=S_h(0,0);
            MatrixXcd Zh = s_to_z(S_h, meas_ne.z0);
            MatrixXcd Zhf = s_to_z(S_h_flipped, meas_ne.z0);
            MatrixXcd Zavg = (Zh + Zhf) / 2.0;
            MatrixXcd Sout = z_to_s(Zavg, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=Sout(0,0); flat[1]=Sout(0,1); flat[2]=Sout(1,0); flat[3]=Sout(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }
};

// Simplified IEEEP370 helper utilities (partial port).
// Provides basic DC extrapolation and thru generation used by some
// IEEEP370-based de-embedding algorithms. This is a pragmatic, numerically
// conservative implementation that avoids external FFT dependencies.
class IEEEP370 {
public:
    // Simple cubic-like interpolation to estimate DC from first few points.
    // This implements a conservative approximation of Python's dc_interp.
    static double dc_interp(const std::vector<std::complex<double>> &svec,
                             const std::vector<double> &freqs) {
        size_t n = std::min<size_t>(9, svec.size());
        if(n < 2) return std::real(svec.front());
        // use least-squares linear fit on real part vs freq for small window
        double sumx=0, sumy=0, sumxx=0, sumxy=0;
        for(size_t i=0;i<n;++i) {
            double x = freqs[i];
            double y = std::real(svec[i]);
            sumx += x; sumy += y; sumxx += x*x; sumxy += x*y;
        }
        double denom = n*sumxx - sumx*sumx;
        double intercept = (denom==0) ? sumy/n : (sumy*sumxx - sumx*sumxy)/denom;
        return intercept;
    }

    // Add DC point by simple interpolation of each S entry (1 or 2 port assumed)
    static Network add_dc(const Network &ntwk) {
        if(ntwk.freqs.empty()) return ntwk;
        size_t nfreq = ntwk.freqs.size();
        int ports = ntwk.n_ports;
        Network out;
        out.n_ports = ports;
        out.z0 = ntwk.z0;
        // build new freqs (prepend 0)
        out.freqs.push_back(Frequency(0.0));
        for(const auto &f: ntwk.freqs) out.freqs.push_back(f);
        out.sparams.reserve(nfreq+1);

        // extract per-frequency matrices and compute DC per-element
        // flatten format: row-major r * n + c
        std::vector<std::complex<double>> first = ntwk.sparams[0];
        // compute DC for each element using linear extrapolation using first two points
        for(size_t fi=0; fi<1; ++fi) {
            (void)fi;
        }
        // compute DC values per matrix element
        std::vector<std::complex<double>> dc_flat(ports*ports);
        if(nfreq >= 2) {
            // linear extrapolation from first two points
            for(int r=0;r<ports;++r) for(int c=0;c<ports;++c) {
                std::complex<double> s0 = ntwk.sparams[0][r*ports + c];
                std::complex<double> s1 = ntwk.sparams[1][r*ports + c];
                // slope = s1-s0 over f1-f0
                double f0 = ntwk.freqs[0].hz;
                double f1 = ntwk.freqs[1].hz;
                if(std::abs(f1 - f0) < 1e-12) dc_flat[r*ports + c] = s0;
                else {
                    std::complex<double> slope = (s1 - s0) / (f1 - f0);
                    dc_flat[r*ports + c] = s0 - slope * f0;
                }
            }
        } else {
            // fallback: use first point
            for(int i=0;i<ports*ports;++i) dc_flat[i] = ntwk.sparams[0][i];
        }

        // push DC then all original sparams
        out.sparams.push_back(dc_flat);
        for(const auto &flat : ntwk.sparams) out.sparams.push_back(flat);
        return out;
    }

    // Create a lossless perfect thru 2-port network matching ntwk frequency/z0
    static Network thru(const Network &ntwk) {
        Network out = ntwk;
        if(out.n_ports < 2) return out;
        for(auto &flat : out.sparams) {
            if(static_cast<int>(flat.size()) >= 4) {
                flat[0] = std::complex<double>(0.0,0.0);
                flat[1] = std::complex<double>(1.0,0.0);
                flat[2] = std::complex<double>(1.0,0.0);
                flat[3] = std::complex<double>(0.0,0.0);
            }
        }
        return out;
    }
};


} // namespace skrf_cpp
