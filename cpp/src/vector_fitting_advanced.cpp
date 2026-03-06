#include "../include/skrf_cpp/vector_fitting_advanced.h"
#include <Eigen/Eigenvalues>
#include <stdexcept>
#include <fstream>

using namespace Eigen;
using namespace std;

namespace skrf_cpp {

// helper to extract H from Network flat format
static std::vector<std::complex<double>> extract_H_from_network(const skrf_cpp::Network &net, int row, int col) {
    std::vector<std::complex<double>> H;
    int n = net.n_ports;
    if(row < 0 || row >= n || col < 0 || col >= n) throw std::runtime_error("fit_network: invalid port indices");
    H.reserve(net.freqs.size());
    for(const auto &flat : net.sparams) {
        if(static_cast<int>(flat.size()) != n*n) throw std::runtime_error("fit_network: sparam size mismatch");
        H.push_back(flat[row * n + col]);
    }
    return H;
}

// helper to extract H from NetworkEigen
static std::vector<std::complex<double>> extract_H_from_network_eigen(const skrf_cpp::NetworkEigen &net, int row, int col) {
    std::vector<std::complex<double>> H;
    int n = net.n_ports;
    if(row < 0 || row >= n || col < 0 || col >= n) throw std::runtime_error("fit_network_eigen: invalid port indices");
    H.reserve(net.freqs.size());
    for(const auto &M : net.s_params) {
        if(static_cast<int>(M.rows()) != n || static_cast<int>(M.cols()) != n) throw std::runtime_error("fit_network_eigen: matrix size mismatch");
        H.push_back(M(row, col));
    }
    return H;
}


static std::vector<complex<double>> poly_multiply(const std::vector<complex<double>> &a, const std::vector<complex<double>> &b) {
    std::vector<complex<double>> out(a.size() + b.size() - 1, complex<double>(0.0,0.0));
    for(size_t i=0;i<a.size();++i) for(size_t j=0;j<b.size();++j) out[i+j] += a[i]*b[j];
    return out;
}

static std::vector<complex<double>> poly_divide_by_linear(const std::vector<complex<double>> &poly, complex<double> root) {
    // divide poly (degree n) by (s - root) returning quotient degree n-1
    size_t n = poly.size()-1;
    std::vector<complex<double>> q(n, complex<double>(0.0,0.0));
    complex<double> acc = poly[0];
    q[0] = acc;
    for(size_t i=1;i<n;i++) {
        acc = poly[i] + acc * root;
        q[i] = acc;
    }
    return q;
}


void VectorFittingAdvanced::fit(const std::vector<double> &freqs_hz,
                                const std::vector<complex<double>> &H,
                                const FitOptions &opts) {
    if(freqs_hz.size() != H.size()) throw std::runtime_error("freqs and H size mismatch");
    int m = static_cast<int>(H.size());
    int npoles = opts.npoles;

    // initial pole placement: place conjugate pairs on left-half plane scaled by max frequency
    double fmax = *std::max_element(freqs_hz.begin(), freqs_hz.end());
    poles.clear(); poles.resize(npoles);
    for(int k=0;k<npoles/2;++k) {
        double re = -2.0*M_PI*fmax * (0.05 + 0.05 * k);
        double im = 2.0*M_PI*fmax * (0.2 + 0.2 * k);
        poles[2*k] = complex<double>(re, im);
        if(2*k+1 < npoles) poles[2*k+1] = std::conj(poles[2*k]);
    }
    if(npoles % 2 == 1) poles[npoles-1] = complex<double>(-2.0*M_PI*fmax*0.1, 0.0);

    // prepare weights (sqrt) for weighting rows
    std::vector<double> w(m, 1.0);
    if(!opts.weights.empty()) {
        if(static_cast<int>(opts.weights.size()) != m) throw std::runtime_error("weights size mismatch");
        for(int i=0;i<m;++i) w[i] = opts.weights[i];
    }

    double lambda = opts.regularization;

    // iterative loop
    for(int iter=0; iter<opts.max_iter; ++iter) {
        int cols = npoles + 2 + npoles; // r_k, d, e, a_k
        MatrixXcd A = MatrixXcd::Zero(m, cols);
        VectorXcd b(m);
        for(int i=0;i<m;++i) {
            complex<double> s = complex<double>(0.0, 2.0*M_PI*freqs_hz[i]);
            double wi = std::sqrt(std::max(0.0, w[i]));
            for(int k=0;k<npoles;++k) {
                complex<double> denom = s - poles[k];
                A(i, k) = (1.0 / denom) * wi; // r_k columns
                A(i, npoles + 2 + k) = (-H[i] / denom) * wi; // a_k columns
            }
            A(i, npoles) = complex<double>(1.0,0.0) * wi; // d
            A(i, npoles+1) = s * wi; // e
            b(i) = H[i] * wi;
        }

        // Solve weighted least squares using SVD with Tikhonov regularization for stability
        // x = V * diag(s/(s^2+lambda)) * U^H * b
        Eigen::JacobiSVD<MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        VectorXcd s = svd.singularValues();
        MatrixXcd U = svd.matrixU();
        MatrixXcd V = svd.matrixV();
        VectorXcd Ut_b = U.adjoint() * b;
        VectorXcd y = VectorXcd::Zero(s.size());
        for(int ii=0; ii<s.size(); ++ii) {
            std::complex<double> denom = s(ii) * s(ii) + std::complex<double>(lambda, 0.0);
            if(std::abs(denom) == 0.0) y(ii) = std::complex<double>(0.0,0.0);
            else y(ii) = (s(ii) * Ut_b(ii)) / denom;
        }
        VectorXcd x = V * y;

        // extract r, d, e, a
        residues.assign(npoles, complex<double>(0.0,0.0));
        std::vector<complex<double>> acoef(npoles);
        for(int k=0;k<npoles;++k) residues[k] = x(k);
        d = x(npoles);
        e = x(npoles+1);
        for(int k=0;k<npoles;++k) acoef[k] = x(npoles + 2 + k);

        // build polynomial: poly = prod_j (s - p_j)
        std::vector<complex<double>> poly = {complex<double>(1.0,0.0)};
        for(int k=0;k<npoles;++k) {
            std::vector<complex<double>> factor = {complex<double>(1.0,0.0), -poles[k]};
            poly = poly_multiply(poly, factor);
        }

        // compute P(s) = poly + sum_k a_k * poly_excl_k
        std::vector<complex<double>> P = poly;
        for(int k=0;k<npoles;++k) {
            std::vector<complex<double>> poly_excl = poly_divide_by_linear(poly, poles[k]);
            for(size_t t=0;t<poly_excl.size();++t) {
                P[t] += acoef[k] * poly_excl[t];
            }
        }

        // companion matrix roots
        int deg = static_cast<int>(P.size()) - 1;
        if(deg != npoles) throw std::runtime_error("polynomial degree mismatch in VF");
        VectorXcd coeffs(deg+1);
        for(int i=0;i<=deg;++i) coeffs(i) = P[deg - i];
        complex<double> leading = coeffs(0);
        if(std::abs(leading) == 0.0) throw std::runtime_error("leading coefficient zero in VF");
        coeffs /= leading;
        MatrixXcd companion = MatrixXcd::Zero(deg, deg);
        for(int i=1;i<deg;++i) companion(i, i-1) = complex<double>(1.0,0.0);
        for(int j=0;j<deg;++j) companion(0, j) = -coeffs(j+1);
        Eigen::EigenSolver<MatrixXcd> es(companion);
        VectorXcd roots = es.eigenvalues();

        // update poles, enforce stability by mirroring RHP poles
        double maxchg = 0.0;
        for(int k=0;k<npoles;++k) {
            complex<double> newp = roots(k);
            if(newp.real() > 0.0) newp = complex<double>(-std::abs(newp.real()), newp.imag());
            double chg = std::abs(newp - poles[k]);
            if(chg > maxchg) maxchg = chg;
            poles[k] = newp;
        }

        // adapt weights if requested: compute residuals r = A * x - b (weighted system)
        if(opts.adaptive_weights) {
            VectorXcd resid = A * x - b;
            double meanw = 0.0;
            for(int i=0;i<m;++i) {
                double ri = std::abs(resid(i));
                // inverse magnitude weighting
                double neww = 1.0 / (ri + opts.adapt_eps);
                w[i] = neww;
                meanw += neww;
            }
            meanw /= static_cast<double>(m);
            // normalize weights to have mean 1 to avoid scaling issues
            for(int i=0;i<m;++i) w[i] /= meanw;
        }

        if(maxchg < opts.tol) break;
    }

    // final residues recompute with final poles (with weighting)
    int m2 = static_cast<int>(H.size());
    int cols2 = static_cast<int>(poles.size()) + 2;
    MatrixXcd A2 = MatrixXcd::Zero(m2, cols2);
    VectorXcd b2(m2);
    for(int i=0;i<m2;++i) {
        complex<double> s = complex<double>(0.0, 2.0*M_PI*freqs_hz[i]);
        for(int k=0;k<static_cast<int>(poles.size());++k) A2(i,k) = 1.0 / (s - poles[k]);
        A2(i, static_cast<int>(poles.size())) = 1.0;
        A2(i, static_cast<int>(poles.size())+1) = s;
        b2(i) = H[i];
    }
    // apply weights if provided
    if(!opts.weights.empty()) {
        for(int i=0;i<m2;++i) {
            double wi = std::sqrt(std::max(0.0, opts.weights[i]));
            A2.row(i) *= wi;
            b2(i) *= wi;
        }
    }
    // final solve using SVD regularized
    Eigen::JacobiSVD<MatrixXcd> svd2(A2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXcd s2 = svd2.singularValues();
    MatrixXcd U2 = svd2.matrixU();
    MatrixXcd V2 = svd2.matrixV();
    VectorXcd Ut_b2 = U2.adjoint() * b2;
    VectorXcd y2 = VectorXcd::Zero(s2.size());
    for(int ii=0; ii<s2.size(); ++ii) {
        std::complex<double> denom = s2(ii) * s2(ii) + std::complex<double>(lambda,0.0);
        if(std::abs(denom) == 0.0) y2(ii) = std::complex<double>(0.0,0.0);
        else y2(ii) = (s2(ii) * Ut_b2(ii)) / denom;
    }
    VectorXcd x2 = V2 * y2;
    residues.resize(poles.size());
    for(int k=0;k<static_cast<int>(poles.size());++k) residues[k] = x2(k);
    d = x2(static_cast<int>(poles.size()));
    e = x2(static_cast<int>(poles.size())+1);
}

std::complex<double> VectorFittingAdvanced::eval(double freq_hz) const {
    complex<double> s(0.0, 2.0*M_PI*freq_hz);
    complex<double> y = d + e * s;
    for(size_t k=0;k<poles.size();++k) y += residues[k] / (s - poles[k]);
    return y;
}

void VectorFittingAdvanced::fit_network(const skrf_cpp::Network &net, int row, int col, const FitOptions &opts) {
    if(net.freqs.empty()) throw std::runtime_error("fit_network: network has no frequencies");
    std::vector<double> freqs; freqs.reserve(net.freqs.size());
    for(const auto &f : net.freqs) freqs.push_back(f.hz);
    auto H = extract_H_from_network(net, row, col);
    fit(freqs, H, opts);
}

void VectorFittingAdvanced::fit_network_eigen(const skrf_cpp::NetworkEigen &net, int row, int col, const FitOptions &opts) {
    if(net.freqs.empty()) throw std::runtime_error("fit_network_eigen: network has no frequencies");
    std::vector<double> freqs; freqs.reserve(net.freqs.size());
    for(const auto &f : net.freqs) freqs.push_back(f.hz);
    auto H = extract_H_from_network_eigen(net, row, col);
    fit(freqs, H, opts);
}

bool VectorFittingAdvanced::save_to_file(const std::string &path) const {
    std::ofstream ofs(path);
    if(!ofs) return false;
    ofs << "POLES " << poles.size() << "\n";
    for(const auto &p : poles) ofs << p.real() << " " << p.imag() << "\n";
    ofs << "RESIDUES " << residues.size() << "\n";
    for(const auto &r : residues) ofs << r.real() << " " << r.imag() << "\n";
    ofs << "D " << d.real() << " " << d.imag() << "\n";
    ofs << "E " << e.real() << " " << e.imag() << "\n";
    return true;
}

bool VectorFittingAdvanced::load_from_file(const std::string &path) {
    std::ifstream ifs(path);
    if(!ifs) return false;
    std::string tag;
    size_t n;
    ifs >> tag >> n;
    if(tag != "POLES") return false;
    poles.clear(); poles.reserve(n);
    for(size_t i=0;i<n;++i) { double re,im; ifs >> re >> im; poles.emplace_back(re,im); }
    ifs >> tag >> n; if(tag != "RESIDUES") return false;
    residues.clear(); residues.reserve(n);
    for(size_t i=0;i<n;++i) { double re,im; ifs >> re >> im; residues.emplace_back(re,im); }
    ifs >> tag; if(tag != "D") return false; double dre,dim; ifs >> dre >> dim; d = std::complex<double>(dre,dim);
    ifs >> tag; if(tag != "E") return false; double ere,eim; ifs >> ere >> eim; e = std::complex<double>(ere,eim);
    return true;
}

// Note: eval(Vector) implemented inline in header as eval_vector wrapper

} // namespace skrf_cpp
