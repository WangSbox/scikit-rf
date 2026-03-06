#pragma once

namespace skrf_cpp {

constexpr double pi = 3.14159265358979323846;
constexpr double c0 = 299792458.0; // speed of light m/s
constexpr double mu0 = 4.0 * pi * 1e-7; // H/m
constexpr double eps0 = 1.0 / (mu0 * c0 * c0); // F/m

// Additional numerical and physical constants
constexpr double INF = 1e99;
constexpr double ALMOST_ZERO = 1e-12;
constexpr double ZERO = 1e-4;
constexpr double ONE = 1.0 + 1.0/1e14;
constexpr double LOG_OF_NEG = -100.0;
constexpr double K_BOLTZMANN = 1.38064852e-23;
constexpr double T0 = 290.0;
constexpr double EIG_COND = 1e-9;
constexpr double EIG_MIN = 1e-12;

// Minimal frequency unit multipliers
// Use string keys in code where needed (e.g., "Hz","kHz","MHz","GHz","THz").
inline const double FREQ_UNIT_HZ = 1.0;
inline const double FREQ_UNIT_KHZ = 1e3;
inline const double FREQ_UNIT_MHZ = 1e6;
inline const double FREQ_UNIT_GHZ = 1e9;
inline const double FREQ_UNIT_THz = 1e12;

// Convert value d in unit to meters. For time-like units (s, us, ns, ps)
// the caller may pass a group velocity v_g (default: speed of light).
inline double to_meters(double d, const char *unit = "m", double v_g = c0) {
	// quick checks for common units (case-insensitive by first char)
	if (!unit) return d;
	char c = unit[0];
	switch (c) {
		case 'm': // m, mm, mil, m... check second char
			if (unit[1] == 'm') return d * 1e-3; // mm
			return d * 1.0; // m
		case 'c': // cm
			return d * 1e-2;
		case 'u': // um or us
			if (unit[1] == 'm') return d * 1e-6; // um
			return d * 1e-6 * v_g; // us
		case 'i': // in
			return d * 0.0254;
		case 's': // s
			return d * v_g;
		case 'n': // ns
			return d * 1e-9 * v_g;
		case 'p': // ps
			return d * 1e-12 * v_g;
		default:
			return d;
	}
}

} // namespace skrf_cpp
