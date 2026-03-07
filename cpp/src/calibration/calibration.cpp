#include "../include/skrf_cpp/calibration/calibration.h"
#include "../include/skrf_cpp/mathFunctions.h"
#include <stdexcept>

namespace skrf_cpp {

Network Calibration::deembed(const Network &meas) const {
    if(error_.n_ports != 2) throw std::runtime_error("Calibration error network must be 2-port");
    if(meas.n_ports != 2) throw std::runtime_error("deembed currently supports only 2-port measurements");
    if(meas.freqs.size() != error_.freqs.size()) throw std::runtime_error("frequency axis mismatch between measurement and error network");

    Network out;
    out.n_ports = meas.n_ports;
    out.z0 = error_.z0;
    out.freqs = meas.freqs;
    out.sparams.reserve(meas.sparams.size());

    for(size_t i=0;i<meas.sparams.size();++i) {
        // measurement S and error S at this freq
        const auto &s_flat = meas.sparams[i];
        const auto &SerrM = error_.s_params[i];
        // build 2x2 S matrices
        Eigen::Matrix2cd Sm;
        Sm(0,0) = s_flat[0]; Sm(0,1) = s_flat[1];
        Sm(1,0) = s_flat[2]; Sm(1,1) = s_flat[3];

        Eigen::Matrix2cd Se = SerrM.block<2,2>(0,0);

        // convert to ABCD
        Eigen::Matrix2cd Am = s_to_abcd_2port(Sm, error_.z0);
        Eigen::Matrix2cd Ae = s_to_abcd_2port(Se, error_.z0);

        // invert error ABCD and apply: A_corrected = A_e^{-1} * A_m
        Eigen::Matrix2cd Ae_inv = Ae.inverse();
        Eigen::Matrix2cd A_corr = Ae_inv * Am;

        // back to S
        Eigen::Matrix2cd S_corr = abcd_to_s_2port(A_corr, error_.z0);

        std::vector<std::complex<double>> flat(4);
        flat[0] = S_corr(0,0); flat[1] = S_corr(0,1);
        flat[2] = S_corr(1,0); flat[3] = S_corr(1,1);
        out.sparams.push_back(std::move(flat));
    }

    return out;
}

} // namespace skrf_cpp
