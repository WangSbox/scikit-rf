#pragma once

#include <complex>
#include <cmath>

namespace skrf_cpp {

inline double db2linear(double db) { return std::pow(10.0, db/20.0); }
inline double linear2db(double lin) { return 20.0 * std::log10(std::max(lin, 1e-300)); }

inline std::complex<double> polar_to_complex(double mag, double phase_deg) {
    double phase = phase_deg * M_PI / 180.0;
    return std::complex<double>(mag * std::cos(phase), mag * std::sin(phase));
}

inline void complex_to_mag_phase(const std::complex<double>& z, double &mag, double &phase_deg) {
    mag = std::abs(z);
    phase_deg = std::arg(z) * 180.0 / M_PI;
}

}
