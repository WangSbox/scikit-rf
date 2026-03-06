#pragma once

#include <vector>
#include <complex>

namespace skrf_cpp {

class CircuitComponent {
public:
    std::string name;
    virtual ~CircuitComponent() = default;
};

// A simple distributed circuit placeholder
class DistributedCircuit : public CircuitComponent {
public:
    DistributedCircuit() = default;
};

} // namespace skrf_cpp
