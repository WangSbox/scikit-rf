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
class IEEEP370_SE_NZC_2xThru_Cpp : public CalibrationBase {
public:
    IEEEP370_SE_NZC_2xThru_Cpp(const Network &dummy_2xthru, double z0 = 50.0, bool use_z_instead_ifft = true, bool verbose=false)
        : s2xthru_(dummy_2xthru), z0_(z0), verbose_(verbose) {
        // fallback to Z-branch implementation for now
        std::tie(s_side1_, s_side2_) = IEEEP370_SE_ZC_2xThru(dummy_2xthru, z0_).split2xthru_z(dummy_2xthru);
    }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override {
        // reuse ZC deembedding logic
        IEEEP370_SE_ZC_2xThru zimpl(s2xthru_, z0_);
        return zimpl.apply_cal(ntwk);
    }
private:
    Network s2xthru_; double z0_{50.0}; bool verbose_{false}; Network s_side1_, s_side2_;
};

class IEEEP370_MM_NZC_2xThru_Cpp : public CalibrationBase {
public:
    IEEEP370_MM_NZC_2xThru_Cpp(const Network &dummy_2xthru, double z0 = 50.0, bool use_z_instead_ifft = true, bool verbose=false)
        : s2xthru_(dummy_2xthru), z0_(z0), verbose_(verbose) {
        std::tie(s_side1_, s_side2_) = IEEEP370_MM_ZC_2xThru(dummy_2xthru, z0_).split2xthru_z(dummy_2xthru);
    }
    void run() override { coefs_computed_ = true; }
    Network apply_cal(const Network &ntwk) const override { IEEEP370_MM_ZC_2xThru zimpl(s2xthru_, z0_); return zimpl.apply_cal(ntwk); }
private:
    Network s2xthru_; double z0_{50.0}; bool verbose_{false}; Network s_side1_, s_side2_;
};

} // namespace skrf_cpp

// Additional calibration subclasses (skeletons) to be implemented fully.
namespace skrf_cpp {

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
