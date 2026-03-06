#pragma once

#include <complex>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <Eigen/Dense>

namespace skrf_cpp {
namespace media {

static const double MU0 = 4.0*M_PI*1e-7;
static const double EPS0 = 8.854187817e-12;
static const double C0 = 1.0/std::sqrt(MU0*EPS0);
static const double ETA0 = std::sqrt(MU0/EPS0);

// Base class for media models
class MediaModel {
public:
    virtual ~MediaModel() {}
    virtual std::complex<double> characteristic_impedance(double freq_hz) const = 0;
    virtual std::complex<double> propagation_constant(double freq_hz) const = 0;
};

// Complete elliptic integral of the first kind via AGM
static double ellipticK(double k) {
    if(k < 0.0 || k >= 1.0) {
        if(k == 0.0) return M_PI_2;
        throw std::runtime_error("ellipticK: modulus out of range");
    }
    double a = 1.0;
    double b = std::sqrt(1.0 - k*k);
    const double eps = 1e-12;
    while(std::abs(a - b) > eps) {
        double an = 0.5*(a + b);
        b = std::sqrt(a*b);
        a = an;
    }
    return M_PI/(2.0*a);
}

// Free space
class FreeSpace : public MediaModel {
public:
    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        return std::complex<double>(ETA0, 0.0);
    }
    std::complex<double> propagation_constant(double freq_hz) const override {
        double omega = 2.0*M_PI*freq_hz;
        return std::complex<double>(0.0, omega * std::sqrt(MU0 * EPS0));
    }
};

// Coaxial (TEM) full model
class CoaxialFull : public MediaModel {
public:
    double a; // inner radius
    double b; // outer radius
    double eps_r;
    double mu_r;
    CoaxialFull(double a_, double b_, double eps_r_=1.0, double mu_r_=1.0)
        : a(a_), b(b_), eps_r(eps_r_), mu_r(mu_r_) {
        if(a <= 0 || b <= 0 || b <= a) throw std::runtime_error("CoaxialFull: invalid radii");
    }
    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        double mu = MU0 * mu_r;
        double eps = EPS0 * eps_r;
        double Z0 = (1.0/(2.0*M_PI)) * std::sqrt(mu/eps) * std::log(b / a);
        return std::complex<double>(Z0, 0.0);
    }
    std::complex<double> propagation_constant(double freq_hz) const override {
        double omega = 2.0*M_PI*freq_hz;
        return std::complex<double>(0.0, omega * std::sqrt(MU0*mu_r * EPS0*eps_r));
    }
};

// Microstrip simple model
class Microstrip : public MediaModel {
public:
    double w, h, eps_r;
    Microstrip(double w_, double h_, double eps_r_=4.4) : w(w_), h(h_), eps_r(eps_r_) {}
    double effective_eps() const {
        double wh = w / h;
        double eps_eff = (eps_r + 1.0)/2.0 + (eps_r - 1.0)/2.0 * 1.0/std::sqrt(1.0 + 12.0/wh);
        return eps_eff;
    }
    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        double wh = w / h;
        double eps_eff = effective_eps();
        double Z0;
        if(wh <= 1.0) {
            Z0 = (60.0 / std::sqrt(eps_eff)) * std::log(8.0/wh + 0.25*wh);
        } else {
            Z0 = (120.0*M_PI) / (std::sqrt(eps_eff) * (wh + 1.393 + 0.667*std::log(wh + 1.444)));
        }
        return std::complex<double>(Z0, 0.0);
    }
    std::complex<double> propagation_constant(double freq_hz) const override {
        double eps_eff = effective_eps();
        double omega = 2.0*M_PI*freq_hz;
        double beta = omega * std::sqrt(MU0 * EPS0 * eps_eff);
        return std::complex<double>(0.0, beta);
    }
};

// Coplanar waveguide (CPW)
class CPW : public MediaModel {
public:
    double w, s, t, h, eps_r;
    CPW(double w_, double s_, double t_, double h_, double eps_r_=1.0) : w(w_), s(s_), t(t_), h(h_), eps_r(eps_r_) {}
    std::complex<double> characteristic_impedance(double /*freq_hz*/) const override {
        double eps_eff = (eps_r + 1.0)/2.0;
        double k = w/(w + 2.0*s);
        if(k <= 0.0 || k >= 1.0) throw std::runtime_error("CPW: invalid geometry k");
        double Kk = ellipticK(k);
        double Kkprime = ellipticK(std::sqrt(1.0 - k*k));
        double Z0 = (30.0*M_PI)/std::sqrt(eps_eff) * (Kkprime / Kk);
        return std::complex<double>(Z0, 0.0);
    }
    std::complex<double> propagation_constant(double freq_hz) const override {
        double eps_eff = (eps_r + 1.0)/2.0;
        double omega = 2.0*M_PI*freq_hz;
        double beta = omega * std::sqrt(MU0 * EPS0 * eps_eff);
        return std::complex<double>(0.0, beta);
    }
};

// Rectangular waveguide (simplified TE10)
class RectangularWaveguide : public MediaModel {
public:
    double a, b, eps_r, mu_r;
    RectangularWaveguide(double a_, double b_, double eps_r_=1.0, double mu_r_=1.0)
        : a(a_), b(b_), eps_r(eps_r_), mu_r(mu_r_) {}
    double cutoff_freq(int m, int n) const {
        double term = (m==0?0.0: (m/(2.0*a)))*(m==0?0.0: (m/(2.0*a)));
        double term2 = (n==0?0.0: (n/(2.0*b)))*(n==0?0.0: (n/(2.0*b)));
        double fc = C0/std::sqrt(eps_r) * std::sqrt(term + term2);
        return fc;
    }
    std::complex<double> propagation_constant(double freq_hz, int m=1, int n=0) const override {
        double fc = cutoff_freq(m,n);
        double omega = 2.0*M_PI*freq_hz;
        double k0 = omega * std::sqrt(MU0 * EPS0 * mu_r * eps_r);
        if(freq_hz <= fc) {
            double alpha = std::sqrt((M_PI*0.0) + 1.0) * (fc - freq_hz);
            return std::complex<double>(alpha, 0.0);
        }
        double beta = k0 * std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(0.0, beta);
    }
    std::complex<double> characteristic_impedance(double freq_hz, int m=1, int n=0) const override {
        double fc = cutoff_freq(m,n);
        double eta = ETA0 * std::sqrt(mu_r/eps_r);
        if(freq_hz <= fc) return std::complex<double>(INFINITY, 0.0);
        double factor = 1.0 / std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(eta * factor, 0.0);
    }
};

// Circular waveguide (TE11 simplified)
class CircularWaveguide : public MediaModel {
public:
    double a, eps_r, mu_r;
    CircularWaveguide(double a_, double eps_r_ = 1.0, double mu_r_ = 1.0)
        : a(a_), eps_r(eps_r_), mu_r(mu_r_) {
        if(a <= 0) throw std::runtime_error("CircularWaveguide: invalid radius");
    }
    double cutoff_TE11() const {
        double uc = 1.841183781;
        return (uc * C0) / (2.0 * M_PI * a * std::sqrt(eps_r));
    }
    std::complex<double> propagation_constant(double freq_hz) const override {
        double fc = cutoff_TE11();
        double omega = 2.0*M_PI*freq_hz;
        double k0 = omega * std::sqrt(MU0 * EPS0 * mu_r * eps_r);
        if(freq_hz <= fc) {
            double alpha = std::sqrt((fc - freq_hz) * (fc - freq_hz));
            return std::complex<double>(alpha, 0.0);
        }
        double beta = k0 * std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(0.0, beta);
    }
    std::complex<double> characteristic_impedance(double freq_hz) const override {
        double fc = cutoff_TE11();
        double eta = std::sqrt(MU0/EPS0) * std::sqrt(mu_r/eps_r);
        if(freq_hz <= fc) return std::complex<double>(INFINITY, 0.0);
        double factor = 1.0 / std::sqrt(1.0 - (fc*fc)/(freq_hz*freq_hz));
        return std::complex<double>(eta * factor, 0.0);
    }
};

// helper: definedAEpTandZ0
static std::tuple<double, double, std::complex<double>> definedAEpTandZ0(double eps_r, double mu_r=1.0, double /*freq_hz*/=1.0) {
    double alpha = 0.0;
    double eps_eff = eps_r;
    double eta = std::sqrt((mu_r*MU0) / (eps_r*EPS0));
    std::complex<double> Z0(eta, 0.0);
    return std::make_tuple(alpha, eps_eff, Z0);
}

} // namespace media

// Transmission line / device helpers in top-level skrf_cpp::tline and media-device in skrf_cpp::device
namespace tline {
    using namespace Eigen;
    class MLine {
    public:
        double z0;
        std::complex<double> gamma;
        double length;
        MLine(double z0_, std::complex<double> gamma_, double length_)
            : z0(z0_), gamma(gamma_), length(length_) {}
        Matrix<std::complex<double>,2,2> abcd() const {
            using cd = std::complex<double>;
            cd coshgl = std::cosh(gamma * length);
            cd sinhgl = std::sinh(gamma * length);
            Matrix<std::complex<double>,2,2> M;
            M(0,0) = coshgl;
            M(0,1) = z0 * sinhgl;
            M(1,0) = (1.0 / z0) * sinhgl;
            M(1,1) = coshgl;
            return M;
        }
    };
}

namespace device {
    using namespace Eigen;
    class Device {
    public:
        std::vector<double> freqs; // Hz
        std::vector<Matrix<std::complex<double>,2,2>> s_params;
        double z0 = 50.0;
        Device() {}
        Device(const std::vector<double> &f, const std::vector<Matrix<std::complex<double>,2,2>> &S, double z0_=50.0)
            : freqs(f), s_params(S), z0(z0_) {
            if(freqs.size() != s_params.size()) throw std::runtime_error("Device: freqs and s_params size mismatch");
        }
        Matrix<std::complex<double>,2,2> abcd_at(size_t idx) const {
            if(idx >= s_params.size()) throw std::out_of_range("Device: index out of range");
            const auto &S = s_params[idx];
            std::complex<double> s11 = S(0,0), s12 = S(0,1), s21 = S(1,0), s22 = S(1,1);
            Matrix<std::complex<double>,2,2> A;
            A(0,0) = ((1.0 + s11)*(1.0 - s22) + s12*s21) / (2.0*s21);
            A(0,1) = z0 * ((1.0 + s11)*(1.0 + s22) - s12*s21) / (2.0*s21);
            A(1,0) = (1.0 / z0) * ((1.0 - s11)*(1.0 - s22) - s12*s21) / (2.0*s21);
            A(1,1) = ((1.0 - s11)*(1.0 + s22) + s12*s21) / (2.0*s21);
            return A;
        }
    };

    class DistributedCircuit {
    public:
        struct Segment { double z0; std::complex<double> gamma; double length; };
        std::vector<Segment> segments;
        void add_segment(double z0, std::complex<double> gamma, double length) {
            segments.push_back({z0,gamma,length});
        }
        Matrix<std::complex<double>,2,2> overall_abcd() const {
            Matrix<std::complex<double>,2,2> M = Matrix<std::complex<double>,2,2>::Identity();
            for(const auto &seg : segments) {
                tline::MLine ml(seg.z0, seg.gamma, seg.length);
                M = M * ml.abcd();
            }
            return M;
        }
        Matrix<std::complex<double>,2,2> overall_s(double ref_z0 = 50.0) const {
            auto A = overall_abcd();
            using cd = std::complex<double>;
            cd a=A(0,0), b=A(0,1), c=A(1,0), d=A(1,1);
            cd denom = (a + b/ref_z0 + c*ref_z0 + d);
            Matrix<std::complex<double>,2,2> S;
            S(0,0) = (a + b/ref_z0 - c*ref_z0 - d)/denom;
            S(0,1) = (2.0*(a*d - b*c))/denom;
            S(1,0) = 2.0/denom;
            S(1,1) = (-a + b/ref_z0 - c*ref_z0 + d)/denom;
            return S;
        }
    };
}

} // namespace skrf_cpp
#pragma once

#include "frequency.h"

namespace skrf_cpp {

class Freespace {
public:
    double z0_override{50.0};
    Freespace() = default;
    explicit Freespace(double z0_override_): z0_override(z0_override_) {}
};

class RectangularWaveguide {
public:
    FrequencySeries freq;
    double a{0.0}, b{0.0};
    double z0_override{50.0};

    RectangularWaveguide() = default;
    RectangularWaveguide(const FrequencySeries &f, double a_, double b_, double z0_override_=50.0)
        : freq(f), a(a_), b(b_), z0_override(z0_override_) {}

    // TE10 cutoff frequency (Hz) for rectangular waveguide with width a (meters)
    double fc_te10() const {
        if (a <= 0.0) return 0.0;
        return c0 / (2.0 * a);
    }
    double characteristic_impedance() const {
        // placeholder: return override if provided
        return z0_override;
    }
};

class Coaxial {
public:
    double a{0.0}; // inner conductor radius (m)
    double b{0.0}; // outer conductor inner radius (m)
    double eps_r{1.0};
    double mu_r{1.0};

    Coaxial() = default;
    Coaxial(double a_, double b_, double eps_r_ = 1.0, double mu_r_ = 1.0)
        : a(a_), b(b_), eps_r(eps_r_), mu_r(mu_r_) {}

    // Characteristic impedance (TEM approximation)
    double characteristic_impedance() const {
        if(a <= 0.0 || b <= 0.0 || b <= a) return 0.0;
        return (60.0 / std::sqrt(eps_r)) * std::log(b / a);
    }

    // Propagation constant ~ j*omega*sqrt(mu*eps) (lossless)
    std::complex<double> propagation_constant(double freq_hz) const {
        using std::sqrt;
        std::complex<double> j(0.0,1.0);
        double eps = eps0 * eps_r;
        double mu = mu0 * mu_r;
        return j * 2.0 * M_PI * freq_hz * sqrt(mu * eps);
    }
};

} // namespace skrf_cpp
