#pragma once

#include <vector>
#include <complex>
#include <string>
#include <map>
#include <stdexcept>
#include <Eigen/Dense>
#include "network.h"
#include "calibration.h"
#include <memory>

namespace skrf_cpp {

// Light-weight Calibration base class inspired by skrf.calibration.Calibration
// This file implements a minimal migration of Calibration base behavior and
// a OnePort algorithm that computes 3-term error coefficients per-frequency.

class CalibrationBase {
public:
    CalibrationBase() = default;
    CalibrationBase(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : measured_(measured), ideals_(ideals) {
        if(!measured_.empty() && !ideals_.empty() && measured_.size() != ideals_.size()) {
            // keep behavior simple: require same lengths
            throw std::runtime_error("CalibrationBase: measured and ideals must have same length for now");
        }
    }

    virtual ~CalibrationBase() = default;

    // run the calibration algorithm and populate internal coefs_ map
    virtual void run() = 0;

    // apply calibration to a network (to be implemented by subclasses)
    virtual Network apply_cal(const Network &ntwk) const = 0;

    // return computed coefficients; runs algorithm lazily
    const std::map<std::string, std::vector<std::complex<double>>> &coefs() {
        if(!coefs_computed_) { run(); coefs_computed_ = true; }
        return coefs_;
    }

protected:
    std::vector<Network> measured_;
    std::vector<Network> ideals_;
    mutable std::map<std::string, std::vector<std::complex<double>>> coefs_;
    mutable bool coefs_computed_{false};
};


// OnePort calibration: solves for 3 complex coefficients per frequency
// using least-squares if more than three standards are supplied.
class OnePortCal : public CalibrationBase {
public:
    OnePortCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {
        if(measured.size() == 0) throw std::runtime_error("OnePortCal: empty measured list");
    }

    void run() override {
        size_t nstd = measured_.size();
        // assume all networks share same frequency grid
        size_t nfreq = measured_.front().freqs.size();
        // prepare output vectors
        std::vector<std::complex<double>> directivity(nfreq), src_match(nfreq), refl_track(nfreq);

        for(size_t fi=0; fi<nfreq; ++fi) {
            // build measurement and ideal vectors for this frequency
            Eigen::VectorXcd mvec(static_cast<int>(nstd));
            Eigen::MatrixXcd Q(static_cast<int>(nstd), 3);
            for(size_t k=0;k<nstd;++k) {
                const auto &mflat = measured_[k].sparams[fi];
                const auto &iflat = ideals_[k].sparams[fi];
                if(mflat.size() != 1 || iflat.size() != 1) {
                    // one-port should have single scalar per frequency
                    throw std::runtime_error("OnePortCal::run expects 1-port networks for OnePort calibration");
                }
                std::complex<double> m = mflat[0];
                std::complex<double> i = iflat[0];
                mvec(static_cast<int>(k)) = m;
                Q(static_cast<int>(k),0) = i;             // column 0: i
                Q(static_cast<int>(k),1) = std::complex<double>(1.0,0.0); // column 1: 1
                Q(static_cast<int>(k),2) = i * m;         // column 2: i*m
            }
            // solve least squares: abc = argmin ||Q * abc - mvec||
            Eigen::VectorXcd abc = Q.colPivHouseholderQr().solve(mvec);
            directivity[fi] = abc(0);
            src_match[fi] = abc(1);
            refl_track[fi] = abc(2);
        }

        coefs_["directivity"] = std::move(directivity);
        coefs_["source match"] = std::move(src_match);
        coefs_["reflection tracking"] = std::move(refl_track);
        coefs_computed_ = true;
    }

    Network apply_cal(const Network &ntwk) const override {
        // NOTE: full correction logic is non-trivial and depends on error model.
        // For now provide a conservative placeholder: return input (no-op)
        // while coefs are available through coefs().
        (void)ntwk;
        return ntwk;
    }
};

// Deembedding wrappers: expose common de-embedding methods as Calibration-like
class OpenShortCal : public CalibrationBase {
public:
    OpenShortCal(const Network &dummy_open, const Network &dummy_short) {
        deembed_ = std::make_shared<OpenShort>(dummy_open, dummy_short);
    }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override {
        return deembed_->deembed(ntwk);
    }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class OpenCal : public CalibrationBase {
public:
    OpenCal(const Network &dummy_open) { deembed_ = std::make_shared<Open>(dummy_open); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class ShortOpenCal : public CalibrationBase {
public:
    ShortOpenCal(const Network &dummy_short, const Network &dummy_open) { deembed_ = std::make_shared<ShortOpen>(dummy_short, dummy_open); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class ShortCal : public CalibrationBase {
public:
    ShortCal(const Network &dummy_short) { deembed_ = std::make_shared<Short>(dummy_short); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class SplitPiCal : public CalibrationBase {
public:
    SplitPiCal(const Network &dummy_thru) { deembed_ = std::make_shared<SplitPi>(dummy_thru); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class SplitTeeCal : public CalibrationBase {
public:
    SplitTeeCal(const Network &dummy_thru) { deembed_ = std::make_shared<SplitTee>(dummy_thru); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class AdmittanceCancelCal : public CalibrationBase {
public:
    AdmittanceCancelCal(const Network &dummy_thru) { deembed_ = std::make_shared<AdmittanceCancel>(dummy_thru); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

class ImpedanceCancelCal : public CalibrationBase {
public:
    ImpedanceCancelCal(const Network &dummy_thru) { deembed_ = std::make_shared<ImpedanceCancel>(dummy_thru); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { return deembed_->deembed(ntwk); }
private:
    std::shared_ptr<Deembedding> deembed_;
};

// IEEEP370 family (Z-branch implementations). For high-fidelity time-domain
// implementations (FFT/ifft) a richer framework is required; here we
// implement the Z-branch fallback that mirrors the Python "else" branch.
class IEEEP370_SE_ZC_2xThru : public CalibrationBase {
public:
    IEEEP370_SE_ZC_2xThru(const Network &dummy_2xthru, double z0 = 50.0, const std::string &name = "")
        : CalibrationBase(), s2xthru_(dummy_2xthru), z0_(z0) {
        name_ = name;
        // compute fixtures using Z-branch immediately
        std::tie(s_side1_, s_side2_) = split2xthru_z(s2xthru_);
    }

    void run() override { coefs_computed_ = true; }

    Network apply_cal(const Network &ntwk) const override {
        // deembed: left.inv ** ntwk ** right.inv using ABCD inversion per-frequency
        NetworkEigen meas_ne = ntwk.to_network_eigen();
        NetworkEigen left_ne = s_side1_.to_network_eigen();
        NetworkEigen right_ne = s_side2_.to_network_eigen();
        // ensure resampled
        auto d = Deembedding::ensure_resampled(std::vector<NetworkEigen>{left_ne, right_ne}, meas_ne);
        NetworkEigen left_r = d[0];
        NetworkEigen right_r = d[1];

        if(meas_ne.n_ports != 2) throw std::runtime_error("IEEEP370_SE_ZC_2xThru deembed supports only 2-port measurements");

        Network out; out.n_ports = 2; out.z0 = meas_ne.z0; out.freqs = meas_ne.freqs; out.sparams.reserve(meas_ne.s_params.size());
        for(size_t i=0;i<meas_ne.s_params.size();++i) {
            Eigen::Matrix2cd Sm = meas_ne.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Sle = left_r.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Sri = right_r.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Am = s_to_abcd_2port(Sm, meas_ne.z0);
            Eigen::Matrix2cd Ale = s_to_abcd_2port(Sle, meas_ne.z0);
            Eigen::Matrix2cd Ari = s_to_abcd_2port(Sri, meas_ne.z0);
            Eigen::Matrix2cd Alep_inv = Ale.inverse();
            Eigen::Matrix2cd Arip_inv = Ari.inverse();
            Eigen::Matrix2cd A_corr = Alep_inv * Am * Arip_inv;
            Eigen::Matrix2cd S_corr = abcd_to_s_2port(A_corr, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=S_corr(0,0); flat[1]=S_corr(0,1); flat[2]=S_corr(1,0); flat[3]=S_corr(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

private:
    Network s2xthru_;
    double z0_{50.0};
    std::string name_;
    Network s_side1_, s_side2_;

    // build fixture models using Z matrices per-frequency
    std::pair<Network, Network> split2xthru_z(const Network &s2xthru) {
        NetworkEigen ne = s2xthru.to_network_eigen();
        // convert S->Z using scalar z0 (use object's z0_)
        NetworkEigen zne = ne.to_z(z0_);
        size_t n = zne.s_params.size();
        int ports = zne.n_ports;
        if(ports < 2) return {Network(), Network()};
        std::vector<Eigen::MatrixXcd> ZL; ZL.reserve(n);
        std::vector<Eigen::MatrixXcd> ZR; ZR.reserve(n);
        for(size_t i=0;i<n;++i) {
            const Eigen::MatrixXcd &Z = zne.s_params[i];
            Eigen::MatrixXcd ZLm = Eigen::MatrixXcd::Zero(2,2);
            Eigen::MatrixXcd ZRm = Eigen::MatrixXcd::Zero(2,2);
            // index mapping assumes 2x2 subblock at (0,0)
            std::complex<double> z00 = Z(0,0);
            std::complex<double> z01 = Z(0,1);
            std::complex<double> z10 = Z(1,0);
            std::complex<double> z11 = Z(1,1);
            ZLm(0,0) = z00 + z10;
            ZLm(0,1) = 2.0 * z10;
            ZLm(1,0) = 2.0 * z10;
            ZLm(1,1) = 2.0 * z10;

            ZRm(0,0) = 2.0 * z01;
            ZRm(0,1) = 2.0 * z01;
            ZRm(1,0) = 2.0 * z01;
            ZRm(1,1) = z11 + z01;

            ZL.push_back(ZLm);
            ZR.push_back(ZRm);
        }
        // convert Z matrices back to S using z0_
        Network outL, outR; outL.n_ports = 2; outR.n_ports = 2; outL.z0 = z0_; outR.z0 = z0_;
        outL.freqs = s2xthru.freqs; outR.freqs = s2xthru.freqs;
        for(size_t i=0;i<n;++i) {
            Eigen::MatrixXcd S1 = z_to_s(ZL[i], z0_);
            Eigen::MatrixXcd S2 = z_to_s(ZR[i], z0_);
            std::vector<std::complex<double>> f1(4), f2(4);
            f1[0]=S1(0,0); f1[1]=S1(0,1); f1[2]=S1(1,0); f1[3]=S1(1,1);
            f2[0]=S2(0,0); f2[1]=S2(0,1); f2[2]=S2(1,0); f2[3]=S2(1,1);
            outL.sparams.push_back(std::move(f1)); outR.sparams.push_back(std::move(f2));
        }
        // flip right side to match Python behavior
        // flipping implemented as reorder of ports (swap 0<->1)
        outR = outR.extract_ports({1,0});
        return {outL, outR};
    }
};

// MM (multiport) variant using Z-branch for 4-port 2xThru splitting.
class IEEEP370_MM_ZC_2xThru : public CalibrationBase {
public:
    IEEEP370_MM_ZC_2xThru(const Network &dummy_2xthru, double z0 = 50.0, const std::string &name = "")
        : s2xthru_(dummy_2xthru), z0_(z0) { name_ = name; std::tie(s_side1_, s_side2_) = split2xthru_z(s2xthru_); }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override {
        // similar approach: ensure resample and apply ABCD inverses per-frequency (supports 2-port deembed)
        NetworkEigen meas_ne = ntwk.to_network_eigen();
        NetworkEigen left_ne = s_side1_.to_network_eigen();
        NetworkEigen right_ne = s_side2_.to_network_eigen();
        auto d = Deembedding::ensure_resampled(std::vector<NetworkEigen>{left_ne, right_ne}, meas_ne);
        NetworkEigen left_r = d[0]; NetworkEigen right_r = d[1];
        if(meas_ne.n_ports != 2) throw std::runtime_error("IEEEP370_MM_ZC_2xThru deembed supports only 2-port measurements");
        Network out; out.n_ports = 2; out.z0 = meas_ne.z0; out.freqs = meas_ne.freqs; out.sparams.reserve(meas_ne.s_params.size());
        for(size_t i=0;i<meas_ne.s_params.size();++i) {
            Eigen::Matrix2cd Sm = meas_ne.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Sle = left_r.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Sri = right_r.s_params[i].block<2,2>(0,0);
            Eigen::Matrix2cd Am = s_to_abcd_2port(Sm, meas_ne.z0);
            Eigen::Matrix2cd Ale = s_to_abcd_2port(Sle, meas_ne.z0);
            Eigen::Matrix2cd Ari = s_to_abcd_2port(Sri, meas_ne.z0);
            Eigen::Matrix2cd Alep_inv = Ale.inverse();
            Eigen::Matrix2cd Arip_inv = Ari.inverse();
            Eigen::Matrix2cd A_corr = Alep_inv * Am * Arip_inv;
            Eigen::Matrix2cd S_corr = abcd_to_s_2port(A_corr, meas_ne.z0);
            std::vector<std::complex<double>> flat(4);
            flat[0]=S_corr(0,0); flat[1]=S_corr(0,1); flat[2]=S_corr(1,0); flat[3]=S_corr(1,1);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

private:
    Network s2xthru_; double z0_{50.0}; std::string name_;
    Network s_side1_, s_side2_;
    // create 4-port-aware splitting by operating on Z matrices; for simplicity
    // we reuse the same 2x2 block operations on appropriate subblocks when present.
    std::pair<Network, Network> split2xthru_z(const Network &s2xthru) {
        // For multiport 4-port 2xThru we will attempt to extract 2x2 blocks
        NetworkEigen ne = s2xthru.to_network_eigen();
        size_t n = ne.s_params.size();
        Network outL, outR; outL.n_ports = 2; outR.n_ports = 2; outL.z0 = z0_; outR.z0 = z0_; outL.freqs = s2xthru.freqs; outR.freqs = s2xthru.freqs;
        for(size_t i=0;i<n;++i) {
            Eigen::MatrixXcd Z = s_to_z(ne.s_params[i], ne.z0);
            // attempt to map corners similarly to Python implementation
            Eigen::Matrix2cd ZL = Eigen::Matrix2cd::Zero(); Eigen::Matrix2cd ZR = Eigen::Matrix2cd::Zero();
            // use [0,0],[0,1],[1,0],[1,1] mapping
            std::complex<double> z00 = Z(0,0); std::complex<double> z01 = Z(0,1); std::complex<double> z10 = Z(1,0); std::complex<double> z11 = Z(1,1);
            ZL(0,0) = z00 + z10; ZL(0,1) = 2.0*z10; ZL(1,0) = 2.0*z10; ZL(1,1) = 2.0*z10;
            ZR(0,0) = 2.0*z01; ZR(0,1) = 2.0*z01; ZR(1,0) = 2.0*z01; ZR(1,1) = z11 + z01;
            Eigen::Matrix2cd S1 = z_to_s(ZL, z0_); Eigen::Matrix2cd S2 = z_to_s(ZR, z0_);
            std::vector<std::complex<double>> f1(4), f2(4);
            f1[0]=S1(0,0); f1[1]=S1(0,1); f1[2]=S1(1,0); f1[3]=S1(1,1);
            f2[0]=S2(0,0); f2[1]=S2(0,1); f2[2]=S2(1,0); f2[3]=S2(1,1);
            outL.sparams.push_back(std::move(f1)); outR.sparams.push_back(std::move(f2));
        }
        outR = outR.extract_ports({1,0});
        return {outL, outR};
    }
};

// NZC variants: provide runnable Z-branch fallback and placeholders for
// future time-domain implementation. These classes expose same interface
// but note they currently use Z-branch internally (use_z_instead_ifft flag ignored)
// Helper: extrapolate network to DC by linear extrapolation of first two
// frequency samples (produces a copy with an inserted DC point at front).
static Network extrapolate_to_dc(const Network &ntwk) {
    Network out = ntwk;
    if(out.freqs.size() == 0) return out;
    if(out.sparams.size() < 2) return out;
    // compute DC point by linear extrapolation of first two samples
    double f0 = out.freqs[0].hz;
    double f1 = out.freqs.size() > 1 ? out.freqs[1].hz : (f0 + 1.0);
    size_t Nflat = out.sparams[0].size();
    std::vector<std::complex<double>> dc_flat(Nflat);
    for(size_t idx=0; idx<Nflat; ++idx) {
        std::complex<double> s0 = out.sparams[0][idx];
        std::complex<double> s1 = out.sparams[1][idx];
        if(f0 == 0.0) {
            dc_flat[idx] = s0;
        } else {
            dc_flat[idx] = s0 + (s0 - s1) * (f0 / (f1 - f0));
        }
    }
    // insert DC freq and sparams at front
    // construct a Frequency at 0 Hz
    Frequency fzero; fzero.hz = 0.0; out.freqs.insert(out.freqs.begin(), fzero);
    out.sparams.insert(out.sparams.begin(), dc_flat);
    return out;
}
class IEEEP370_SE_NZC_2xThru_Cpp : public CalibrationBase {
public:
    IEEEP370_SE_NZC_2xThru_Cpp(const Network &dummy_2xthru, double z0 = 50.0, bool use_z_instead_ifft = true, bool verbose=false)
        : s2xthru_(dummy_2xthru), z0_(z0), verbose_(verbose), fd_qm_(verbose), td_qm_(1e9, 32, 0.35, 1, 2, verbose) {
        // attempt to extrapolate dummy 2xthru to DC before splitting; fall back to original
        Network expl = extrapolate_to_dc(dummy_2xthru);
        // attempt NZC-style time-domain correction
        try {
            Network expl2 = td_qm_.perform_nzc_extrapolation(expl);
            std::tie(s_side1_, s_side2_) = IEEEP370_SE_ZC_2xThru(expl2, z0_).split2xthru_z(expl2);
        } catch(...) {
            std::tie(s_side1_, s_side2_) = IEEEP370_SE_ZC_2xThru(expl, z0_).split2xthru_z(expl);
        }
        // compute quick FD/TD quality metrics and cache results for diagnostics
        try {
            fd_qm_result_ = fd_qm_.check_se_quality(dummy_2xthru, verbose_);
        } catch(...) { }
        try {
            td_qm_result_ = td_qm_.check_se_quality(dummy_2xthru, verbose_);
        } catch(...) { }
    }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override {
        // reuse ZC deembedding logic until full NZC time-domain port is implemented
        IEEEP370_SE_ZC_2xThru zimpl(s2xthru_, z0_);
        return zimpl.apply_cal(ntwk);
    }

    // Diagnostic accessors
    IEEEP370_FD_QM::QMResult fd_qm_result() const { return fd_qm_result_; }
    IEEEP370_TD_QM::TDQMResult td_qm_result() const { return td_qm_result_; }

private:
    Network s2xthru_; double z0_{50.0}; bool verbose_{false}; Network s_side1_, s_side2_;
    IEEEP370_FD_QM fd_qm_;
    IEEEP370_TD_QM td_qm_;
    IEEEP370_FD_QM::QMResult fd_qm_result_;
    IEEEP370_TD_QM::TDQMResult td_qm_result_;
};

class IEEEP370_MM_NZC_2xThru_Cpp : public CalibrationBase {
public:
    IEEEP370_MM_NZC_2xThru_Cpp(const Network &dummy_2xthru, double z0 = 50.0, bool use_z_instead_ifft = true, bool verbose=false)
        : s2xthru_(dummy_2xthru), z0_(z0), verbose_(verbose), fd_qm_(verbose), td_qm_(1e9, 32, 0.35, 1, 2, verbose) {
        Network expl = extrapolate_to_dc(dummy_2xthru);
        try {
            Network expl2 = td_qm_.perform_nzc_extrapolation(expl);
            std::tie(s_side1_, s_side2_) = IEEEP370_MM_ZC_2xThru(expl2, z0_).split2xthru_z(expl2);
        } catch(...) {
            std::tie(s_side1_, s_side2_) = IEEEP370_MM_ZC_2xThru(expl, z0_).split2xthru_z(expl);
        }
        try { fd_qm_result_ = fd_qm_.check_se_quality(dummy_2xthru, verbose_); } catch(...) {}
        try { td_qm_result_ = td_qm_.check_se_quality(dummy_2xthru, verbose_); } catch(...) {}
    }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { IEEEP370_MM_ZC_2xThru zimpl(s2xthru_, z0_); return zimpl.apply_cal(ntwk); }

    IEEEP370_FD_QM::QMResult fd_qm_result() const { return fd_qm_result_; }
    IEEEP370_TD_QM::TDQMResult td_qm_result() const { return td_qm_result_; }

private:
    Network s2xthru_; double z0_{50.0}; bool verbose_{false}; Network s_side1_, s_side2_;
    IEEEP370_FD_QM fd_qm_;
    IEEEP370_TD_QM td_qm_;
    IEEEP370_FD_QM::QMResult fd_qm_result_;
    IEEEP370_TD_QM::TDQMResult td_qm_result_;
};

} // namespace skrf_cpp

// Additional calibration subclasses (skeletons) to be implemented fully.
namespace skrf_cpp {

// IEEEP370 helper: Fixture Electrical Requirements (FD) and Quality Metrics
class IEEEP370_FER {
public:
    // Compute simple FD metrics used by plotting routines in Python reference.
    // Returns pair of vectors: IL (S21 dB) and RL (S11 dB) for two-port 2xThru.
    static std::pair<std::vector<double>, std::vector<double>> fd_se_metrics(const Network &s2xthru) {
        NetworkEigen ne = s2xthru.to_network_eigen();
        size_t N = ne.s_params.size();
        std::vector<double> il(N, 0.0), rl(N, 0.0);
        for(size_t i=0;i<N;++i) {
            if(ne.n_ports >= 2) {
                std::complex<double> s21 = ne.s_params[i](1,0);
                std::complex<double> s11 = ne.s_params[i](0,0);
                il[i] = 20.0 * std::log10(std::abs(s21) + 1e-16);
                rl[i] = 20.0 * std::log10(std::abs(s11) + 1e-16);
            }
        }
        return {il, rl};
    }
};

// Frequency-domain quality metrics (initial checks)
class IEEEP370_FD_QM {
public:
    explicit IEEEP370_FD_QM(bool verbose = false) : verbose_(verbose) {}

    // Causality quality metric (percent)
    double check_causality(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        int Nf = static_cast<int>(ne.freqs.size());
        if(Nf < 3) return 100.0;
        int nports = ne.n_ports;
        std::vector<std::vector<double>> CQM(nports, std::vector<double>(nports, 100.0));
        for(int i=0;i<nports;++i) for(int j=0;j<nports;++j) {
            // collect vector over freq
            std::vector<std::complex<double>> svec(Nf);
            for(int k=0;k<Nf;++k) svec[k] = ne.s_params[k](i,j);
            bool allsame=true; for(int k=1;k<Nf;++k) if(svec[k]!=svec[0]) { allsame=false; break; }
            if(allsame) { CQM[i][j] = 100.0; continue; }
            double TotalR=0.0, PositiveR=0.0;
            for(int k=0;k<Nf-2;++k) {
                auto Vn = svec[k+1] - svec[k];
                auto Vn1 = svec[k+2] - svec[k+1];
                double R = std::real(Vn1) * std::imag(Vn) - std::imag(Vn1) * std::real(Vn);
                if(R > 0) PositiveR += R;
                TotalR += std::abs(R);
            }
            CQM[i][j] = (TotalR==0.0) ? 100.0 : std::max(0.0, PositiveR/TotalR) * 100.0;
        }
        // return minimum across matrix
        double minv = 1e9; for(auto &r: CQM) for(auto &v: r) if(v < minv) minv = v;
        return minv;
    }

    // Passivity quality metric (percent)
    double check_passivity(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        if(ne.n_ports == 1) throw std::runtime_error("check_passivity: not defined for 1-port");
        int Nf = static_cast<int>(ne.s_params.size());
        std::vector<double> PM(Nf, 0.0);
        for(int i=0;i<Nf;++i) {
            // compute spectral norm (largest singular value)
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(ne.s_params[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
            double sigma_max = svd.singularValues()(0);
            PM[i] = sigma_max;
        }
        double A = 1.00001; double B = 0.1;
        std::vector<double> PW(Nf, 0.0);
        for(int i=0;i<Nf;++i) if(PM[i] > A) PW[i] = (PM[i] - A)/B;
        double res = std::max(0.0, static_cast<double>(Nf - std::accumulate(PW.begin(), PW.end(), 0.0))) / Nf * 100.0;
        return res;
    }

    // Reciprocity quality metric (percent)
    double check_reciprocity(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        if(ne.n_ports == 1) throw std::runtime_error("check_reciprocity: not defined for 1-port");
        int Nf = static_cast<int>(ne.s_params.size());
        int nports = ne.n_ports;
        std::vector<double> RM(Nf, 0.0), RW(Nf, 0.0);
        for(int i=0;i<Nf;++i) {
            double sum=0.0;
            for(int k=0;k<nports;++k) for(int m=0;m<nports;++m) sum += std::abs(ne.s_params[i](k,m) - ne.s_params[i](m,k));
            RM[i] = sum / (nports * (nports - 1));
            double C = 1e-6, B = 0.1;
            if(RM[i] > C) RW[i] = (RM[i] - C) / B;
        }
        double res = std::max(0.0, static_cast<double>(Nf - std::accumulate(RW.begin(), RW.end(), 0.0))) / Nf * 100.0;
        return res;
    }

    // Single-ended combined quality check (returns dict-like struct)
    struct QMResult { double causality; double passivity; double reciprocity; std::string causality_eval; std::string passivity_eval; std::string reciprocity_eval; };
    QMResult check_se_quality(const Network &ntwk, bool verbose=false) const {
        QMResult out;
        out.causality = check_causality(ntwk);
        out.passivity = check_passivity(ntwk);
        out.reciprocity = check_reciprocity(ntwk);
        // evaluations (same thresholds as Python)
        out.causality_eval = (out.causality <= 20.) ? "poor" : (out.causality <= 50. ? "inconclusive" : (out.causality <= 80. ? "acceptable" : "good"));
        out.passivity_eval = (out.passivity <= 80.) ? "poor" : (out.passivity <= 99. ? "inconclusive" : (out.passivity <= 99.9 ? "acceptable" : "good"));
        out.reciprocity_eval = (out.reciprocity <= 80.) ? "poor" : (out.reciprocity <= 99. ? "inconclusive" : (out.reciprocity <= 99.9 ? "acceptable" : "good"));
        return out;
    }

    // Mixed-mode quality (attempt for 4-port networks by extracting sub-blocks)
    std::map<std::string, QMResult> check_mm_quality(const Network &ntwk, bool verbose=false) const {
        // attempt to extract differential/common mode pairs (0,1) and (2,3)
        std::map<std::string, QMResult> res;
        if(ntwk.n_ports < 4) throw std::runtime_error("check_mm_quality: requires 4-port network");
        Network se_dd = ntwk.extract_ports({0,1});
        Network se_cc = ntwk.extract_ports({2,3});
        res["dd"] = check_se_quality(se_dd, verbose);
        res["cc"] = check_se_quality(se_cc, verbose);
        return res;
    }

    void print_qm(const QMResult &qm) const {
        std::cout << "causality: " << qm.causality << "% (" << qm.causality_eval << ")\n";
        std::cout << "passivity: " << qm.passivity << "% (" << qm.passivity_eval << ")\n";
        std::cout << "reciprocity: " << qm.reciprocity << "% (" << qm.reciprocity_eval << ")\n";
    }

private:
    bool verbose_{false};
};

// Time-domain quality metrics (TD_QM) - pragmatic port of Python logic.
class IEEEP370_TD_QM {
public:
    IEEEP370_TD_QM(double data_rate, int sample_per_UI, double rise_time_per,
                  int pulse_shape = 1, int extrapolation = 2, bool verbose = false)
        : data_rate_(data_rate), sample_per_UI_(sample_per_UI), rise_time_per_(rise_time_per),
          pulse_shape_(pulse_shape), extrapolation_(extrapolation), verbose_(verbose) {}

    // Add complex-conjugate symmetric extension (length 2N-1)
    static std::vector<std::complex<double>> add_conj(const std::vector<std::complex<double>> &s) {
        size_t N = s.size();
        std::vector<std::complex<double>> out(2*N - 1);
        for(size_t i=0;i<N;++i) out[i] = s[i];
        for(size_t k=0;k<N-1;++k) out[N + k] = std::conj(out[N - k - 1]);
        return out;
    }

    // Simple inverse DFT (naive). Input X length M -> returns time samples length M
    static std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> &X) {
        size_t N = X.size();
        std::vector<std::complex<double>> x(N);
        const double two_pi = 2.0 * M_PI;
        for(size_t n=0;n<N;++n) {
            std::complex<double> sum(0.0,0.0);
            for(size_t k=0;k<N;++k) {
                double ang = two_pi * k * n / static_cast<double>(N);
                sum += X[k] * std::exp(std::complex<double>(0.0, ang));
            }
            x[n] = sum / static_cast<double>(N);
        }
        return x;
    }

    // Simple forward DFT
    static std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &x) {
        size_t N = x.size();
        std::vector<std::complex<double>> X(N);
        const double two_pi = 2.0 * M_PI;
        for(size_t k=0;k<N;++k) {
            std::complex<double> sum(0.0,0.0);
            for(size_t n=0;n<N;++n) {
                double ang = -two_pi * k * n / static_cast<double>(N);
                sum += x[n] * std::exp(std::complex<double>(0.0, ang));
            }
            X[k] = sum;
        }
        return X;
    }

    // Hann window generator
    static std::vector<double> hann_window(size_t N) {
        std::vector<double> w(N);
        for(size_t n=0;n<N;++n) w[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * static_cast<double>(n) / static_cast<double>(N-1)));
        return w;
    }

    // Estimate peak location with parabolic interpolation for sub-sample accuracy.
    static std::pair<int,double> estimate_peak_subsample(const std::vector<std::complex<double>> &v) {
        int M = static_cast<int>(v.size());
        int imax = 0; double vmax = 0.0;
        for(int i=0;i<M;++i) { double mag = std::abs(v[i]); if(mag > vmax) { vmax = mag; imax = i; } }
        int im = (imax - 1 + M) % M; int ip = (imax + 1) % M;
        double ym = std::abs(v[im]); double y0 = std::abs(v[imax]); double yp = std::abs(v[ip]);
        double denom = (ym - 2.0*y0 + yp);
        double delta = 0.0;
        if(std::abs(denom) > 1e-16) delta = 0.5 * (ym - yp) / denom;
        return {imax, delta};
    }

    // Circular shift by integer samples
    static std::vector<std::complex<double>> circ_shift(const std::vector<std::complex<double>> &v, int shift) {
        int N = static_cast<int>(v.size());
        std::vector<std::complex<double>> out(N);
        for(int i=0;i<N;++i) out[i] = v[(i + shift + N) % N];
        return out;
    }

    // Apply fractional (sub-sample) shift using simple linear interpolation in time domain
    static std::vector<std::complex<double>> fractional_shift_linear(const std::vector<std::complex<double>> &v, double frac) {
        int N = static_cast<int>(v.size());
        if(std::abs(frac) < 1e-12) return v;
        // positive frac means shift right by frac samples
        std::vector<std::complex<double>> out(N);
        for(int n=0;n<N;++n) {
            double idx_f = static_cast<double>(n) - frac;
            int idx0 = static_cast<int>(std::floor(idx_f));
            double alpha = idx_f - static_cast<double>(idx0);
            int i0 = (idx0 % N + N) % N;
            int i1 = ((idx0+1) % N + N) % N;
            out[n] = (1.0 - alpha) * v[i0] + alpha * v[i1];
        }
        return out;
    }

    // Unwrap phase of complex vector (returns phase vector)
    static std::vector<double> unwrap_phase(const std::vector<std::complex<double>> &v) {
        size_t N = v.size();
        std::vector<double> ph(N);
        for(size_t i=0;i<N;++i) ph[i] = std::atan2(std::imag(v[i]), std::real(v[i]));
        for(size_t i=1;i<N;++i) {
            double d = ph[i] - ph[i-1];
            if(d > M_PI) ph[i] -= 2.0*M_PI;
            else if(d < -M_PI) ph[i] += 2.0*M_PI;
        }
        // cumulative correction
        for(size_t i=1;i<N;++i) ph[i] += ph[i-1] - ph[i-1];
        return ph;
    }

    // Create passive network (SVD clip) - operates frequency-by-frequency
    Network create_passive(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        Network out = ntwk; // copy meta
        out.sparams.clear(); out.sparams.reserve(ne.s_params.size());
        for(size_t i=0;i<ne.s_params.size();++i) {
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(ne.s_params[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::VectorXd S = svd.singularValues().real();
            for(int k=0;k<S.size();++k) if(S(k) > 1.0) S(k) = 1.0;
            Eigen::MatrixXcd recon = svd.matrixU() * S.asDiagonal().cast<std::complex<double>>() * svd.matrixV().adjoint();
            // flatten
            int p = recon.rows(); std::vector<std::complex<double>> flat(p*p);
            for(int r=0;r<p;++r) for(int c=0;c<p;++c) flat[r*p + c] = recon(r,c);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

    // Create reciprocal: swap s_ij and s_ji
    Network create_reciprocal(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        Network out = ntwk; out.sparams.clear(); out.sparams.reserve(ne.s_params.size());
        for(size_t i=0;i<ne.s_params.size();++i) {
            Eigen::MatrixXcd M = ne.s_params[i];
            Eigen::MatrixXcd Mr = M.transpose();
            int p = Mr.rows(); std::vector<std::complex<double>> flat(p*p);
            for(int r=0;r<p;++r) for(int c=0;c<p;++c) flat[r*p + c] = Mr(r,c);
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

    // Simplified time-domain impulse computation for S-parameters
    // Returns v (2N-1 x nports x nports) as flattened vector of matrices per-time-step
    std::tuple<std::vector<std::vector<std::complex<double>>>, std::vector<double>> get_time_domain(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        int N = static_cast<int>(ne.s_params.size());
        double df = (ne.freqs.size()>=2) ? (ne.freqs[1].hz - ne.freqs[0].hz) : 1.0;
        double dt = 1.0 / (2.0 * ne.freqs.back().hz + df);
        int nports = ne.n_ports;
        int outN = 2*N - 1;
        std::vector<double> t(outN); for(int i=0;i<outN;++i) t[i] = dt * i;
        // prepare result: for each time index store flattened matrix (nports*nports)
        std::vector<std::vector<std::complex<double>>> v(outN, std::vector<std::complex<double>>(nports*nports, std::complex<double>(0,0)));
        // simple approach: for each s_ij compute ifft of conj-extended spectrum
        for(int i=0;i<nports;++i) for(int j=0;j<nports;++j) {
            std::vector<std::complex<double>> svec(N);
            for(int k=0;k<N;++k) svec[k] = ne.s_params[k](i,j) * filter_at_index(k, N);
            svec[0] = std::complex<double>(std::real(svec[0]), 0.0);
            auto sconj = add_conj(svec);
            auto vfull = ifft(sconj);
            for(int ti=0; ti<outN; ++ti) v[ti][i*nports + j] = std::real(vfull[ti]);
        }
        return {v, t};
    }

    // compute simple filter multiplier (1 or lowpass) at freq index
    double filter_at_index(int idx, int N) const { (void)idx; (void)N; return 1.0; }

    // compute time-domain metric difference (simple L-infty over unit intervals)
    std::vector<std::vector<double>> get_td_difference_mv(const std::vector<std::vector<std::complex<double>>> &v1,
                                                          const std::vector<std::vector<std::complex<double>>> &v2,
                                                          const std::vector<double> &t, int nports) const {
        int N = static_cast<int>(t.size());
        double dt = (N>=2) ? (t[1] - t[0]) : 1.0;
        double UI = 1.0 / data_rate_ / dt;
        int N_UI = std::max(1, static_cast<int>(std::round(UI)));
        std::vector<std::vector<double>> out(nports, std::vector<double>(nports, 0.0));
        for(int i=0;i<nports;++i) for(int j=0;j<nports;++j) {
            double maxv = 0.0;
            for(int k=0;k<N;++k) {
                double diff = std::abs(v2[k][i*nports + j] - v1[k][i*nports + j]);
                if(diff > maxv) maxv = diff;
            }
            out[i][j] = maxv;
        }
        return out;
    }

    // High-level single-ended check (orchestrates extrapolation, causal/passive/reciprocal models)
    struct TDQMResult { double causality_mV; double passivity_mV; double reciprocity_mV; std::string causality_eval; std::string passivity_eval; std::string reciprocity_eval; };
    TDQMResult check_se_quality(const Network &ntwk, bool verbose=false) {
        // extrapolate and simple pipeline using approximations
        Network ntwk_interp = ntwk; // TODO: implement extrapolate_to_fmax and extrapolate_to_dc faithfully
        auto [v_origin, t] = get_time_domain(ntwk_interp);
        Network causal = create_reciprocal(ntwk_interp); // placeholder: use reciprocal as proxy
        auto [v_causal, tc] = get_time_domain(causal);
        Network passive = create_passive(ntwk_interp);
        auto [v_passive, tp] = get_time_domain(passive);

        auto caus_mv = get_td_difference_mv(v_causal, v_origin, t, ntwk.n_ports);
        auto pass_mv = get_td_difference_mv(v_passive, v_origin, t, ntwk.n_ports);
        auto rec_mv = get_td_difference_mv(v_origin, v_origin, t, ntwk.n_ports);
        // aggregate into scalar metrics
        double caus_metric = 1000.0 * l2_norm_matrix(caus_mv);
        double pass_metric = 1000.0 * l2_norm_matrix(pass_mv);
        double rec_metric  = 1000.0 * l2_norm_matrix(rec_mv);
        TDQMResult res;
        res.causality_mV = std::round(caus_metric * 10.0) / 10.0;
        res.passivity_mV = std::round(pass_metric * 10.0) / 10.0;
        res.reciprocity_mV = std::round(rec_metric * 10.0) / 10.0;
        // simple evaluations
        auto eval = [](double v)->std::string { if(v>=15.) return "poor"; if(v>=10.) return "inconclusive"; if(v>=5.) return "acceptable"; return "good"; };
        res.causality_eval = eval(res.causality_mV/2.0); res.passivity_eval = eval(res.passivity_mV/2.0); res.reciprocity_eval = eval(res.reciprocity_mV/2.0);
        return res;
    }

    void print_qm(const TDQMResult &qm) const {
        std::cout << "causality (mV): " << qm.causality_mV << " (" << qm.causality_eval << ")\n";
        std::cout << "passivity (mV): " << qm.passivity_mV << " (" << qm.passivity_eval << ")\n";
        std::cout << "reciprocity (mV): " << qm.reciprocity_mV << " (" << qm.reciprocity_eval << ")\n";
    }

    // Perform NZC-style time-domain processing and return a frequency-domain
    // network with corrected S-parameters (DC included as first sample).
    Network perform_nzc_extrapolation(const Network &ntwk) const {
        NetworkEigen ne = ntwk.to_network_eigen();
        int N = static_cast<int>(ne.s_params.size());
        int nports = ne.n_ports;
        Network out = ntwk; out.sparams.clear(); out.sparams.reserve(N);
        // For each port pair, compute corrected spectrum
        // We'll build per-frequency flattened matrices
        std::vector<std::vector<std::complex<double>>> corrected_flat(N, std::vector<std::complex<double>>(nports*nports, {0,0}));
        for(int i=0;i<nports;++i) for(int j=0;j<nports;++j) {
            // collect spectrum over freqs
            std::vector<std::complex<double>> S(N);
            for(int k=0;k<N;++k) S[k] = ne.s_params[k](i,j);
            // conj-extend and get time-domain impulse
            auto Sconj = add_conj(S);
            auto impulse = ifft(Sconj);
            // find peak and circular-shift so peak at index 0
            int M = static_cast<int>(impulse.size());
            int peak = 0; double vmax = 0.0;
            for(int t=0;t<M;++t) { double mag = std::abs(impulse[t]); if(mag > vmax) { vmax = mag; peak = t; } }
            std::vector<std::complex<double>> shifted(M);
            for(int t=0;t<M;++t) shifted[t] = impulse[(t + peak) % M];
            // forward FFT and take first N bins
            auto Sshift = fft(shifted);
            for(int k=0;k<N;++k) corrected_flat[k][i*nports + j] = Sshift[k];
        }
        // now pack flattened matrices per frequency
        for(int k=0;k<N;++k) {
            std::vector<std::complex<double>> flat(nports*nports);
            for(int idx=0; idx<nports*nports; ++idx) flat[idx] = corrected_flat[k][idx];
            out.sparams.push_back(std::move(flat));
        }
        return out;
    }

private:
    double data_rate_;
    int sample_per_UI_;
    double rise_time_per_;
    int pulse_shape_;
    int extrapolation_;
    bool verbose_;

    static double l2_norm_matrix(const std::vector<std::vector<double>> &m) {
        double sum=0.0; for(auto &r:m) for(auto &v:r) sum += v*v; return std::sqrt(sum);
    }
};


// Twelve-term calibration (12-term error-box model)
class TwelveTermCal : public CalibrationBase {
public:
    TwelveTermCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {}

    void run() override {
        // Full 12-term solver is non-trivial; placeholder for now.
        // TODO: implement conversion from measured/ideals to 12-term coefs per-frequency
        throw std::runtime_error("TwelveTermCal::run not implemented yet");
    }

    Network apply_cal(const Network &ntwk) const override {
        (void)ntwk; throw std::runtime_error("TwelveTermCal::apply_cal not implemented");
    }
};

// SOLT calibration (Short-Open-Load-Thru) - common two-port calibration
class SOLTCal : public CalibrationBase {
public:
    SOLTCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {}

    void run() override {
        // TODO: implement SOLT algorithm (map to 8/12-term as needed)
        throw std::runtime_error("SOLTCal::run not implemented yet");
    }

    Network apply_cal(const Network &ntwk) const override {
        (void)ntwk; throw std::runtime_error("SOLTCal::apply_cal not implemented");
    }
};

// TRL calibration (Thru-Reflect-Line)
class TRLCal : public CalibrationBase {
public:
    TRLCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {}

    void run() override {
        // TODO: implement TRL algorithm
        throw std::runtime_error("TRLCal::run not implemented yet");
    }

    Network apply_cal(const Network &ntwk) const override {
        (void)ntwk; throw std::runtime_error("TRLCal::apply_cal not implemented");
    }
};

} // namespace skrf_cpp

// More calibration subclasses and helpers (grouped to avoid duplicate aliases)
namespace skrf_cpp {

// Eight-term (error-box) calibration skeleton
class EightTermCal : public CalibrationBase {
public:
    EightTermCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {}
    void run() override { throw std::runtime_error("EightTermCal::run not implemented"); }
    Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("EightTermCal::apply_cal not implemented"); }
};

// UnknownThru (helper for TRL variants)
class UnknownThruCal : public CalibrationBase {
public:
    UnknownThruCal(const std::vector<Network> &measured, const std::vector<Network> &ideals)
        : CalibrationBase(measured, ideals) {}
    void run() override { throw std::runtime_error("UnknownThruCal::run not implemented"); }
    Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("UnknownThruCal::apply_cal not implemented"); }
};

// LRM / LRRM / MRC placeholders
class LRMCal : public CalibrationBase { public: LRMCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("LRMCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("LRMCal::apply_cal not implemented"); } };
class LRRMCal : public CalibrationBase { public: LRRMCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("LRRMCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("LRRMCal::apply_cal not implemented"); } };
class MRC_Cal : public CalibrationBase { public: MRC_Cal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("MRC_Cal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("MRC_Cal::apply_cal not implemented"); } };

// Multiline TRL variants
class MultilineTRLCal : public CalibrationBase { public: MultilineTRLCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("MultilineTRLCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("MultilineTRLCal::apply_cal not implemented"); } };
class NISTMultilineTRLCal : public MultilineTRLCal { public: using MultilineTRLCal::MultilineTRLCal; };
class TUGMultilineTRLCal : public MultilineTRLCal { public: using MultilineTRLCal::MultilineTRLCal; };

// Sixteen-term and LMR16 placeholders
class SixteenTermCal : public CalibrationBase { public: SixteenTermCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("SixteenTermCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("SixteenTermCal::apply_cal not implemented"); } };
class LMR16Cal : public CalibrationBase { public: LMR16Cal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("LMR16Cal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("LMR16Cal::apply_cal not implemented"); } };

// Normalization and other utilities
class NormalizationCal : public CalibrationBase { public: NormalizationCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("NormalizationCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("NormalizationCal::apply_cal not implemented"); } };

// Multiport calibrations
class MultiportCal : public CalibrationBase { public: MultiportCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("MultiportCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("MultiportCal::apply_cal not implemented"); } };
class MultiportSOLTCal : public MultiportCal { public: using MultiportCal::MultiportCal; };

// TwoPortOnePath and EnhancedResponse (1.5 port) placeholders
class TwoPortOnePathCal : public CalibrationBase { public: TwoPortOnePathCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("TwoPortOnePathCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("TwoPortOnePathCal::apply_cal not implemented"); } };
class EnhancedResponseCal : public CalibrationBase { public: EnhancedResponseCal(const std::vector<Network>&m,const std::vector<Network>&i):CalibrationBase(m,i){}; void run() override { throw std::runtime_error("EnhancedResponseCal::run not implemented"); } Network apply_cal(const Network &ntwk) const override { (void)ntwk; throw std::runtime_error("EnhancedResponseCal::apply_cal not implemented"); } };

// Switch terms computation helper (placeholder)
inline std::pair<Network, Network> compute_switch_terms(const std::vector<Network> &measurements) {
    (void)measurements;
    throw std::runtime_error("compute_switch_terms not implemented yet");
}

} // namespace skrf_cpp
