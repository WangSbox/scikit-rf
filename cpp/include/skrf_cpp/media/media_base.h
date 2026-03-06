#pragma once

#include <complex>

namespace skrf_cpp {
namespace media {

class MediaModel {
public:
    virtual ~MediaModel() {}
    virtual std::complex<double> characteristic_impedance(double freq_hz) const = 0;
    virtual std::complex<double> propagation_constant(double freq_hz) const = 0;
};

} // namespace media
} // namespace skrf_cpp
