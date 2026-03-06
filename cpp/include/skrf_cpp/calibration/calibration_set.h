#pragma once

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include "calibration.h"
#include "network.h"

namespace skrf_cpp {

class CalibrationSet {
public:
    CalibrationSet() {}

    void add(const std::string &name, const Calibration &cal) {
        cals[name] = cal;
    }

    bool has(const std::string &name) const {
        return cals.find(name) != cals.end();
    }

    Calibration get(const std::string &name) const {
        auto it = cals.find(name);
        if(it == cals.end()) throw std::runtime_error("CalibrationSet: calibration not found");
        return it->second;
    }

    // Apply named calibration de-embedding to a single Network
    Network deembed(const Network &net, const std::string &name) const {
        auto cal = get(name);
        return cal.deembed(net);
    }

    // Batch deembed: apply named calibration to a list of networks
    std::vector<Network> deembed_all(const std::vector<Network> &nets, const std::string &name) const {
        std::vector<Network> out;
        out.reserve(nets.size());
        for(const auto &n : nets) out.push_back(deembed(n, name));
        return out;
    }

    void remove(const std::string &name) {
        cals.erase(name);
    }

private:
    std::map<std::string, Calibration> cals;
};

} // namespace skrf_cpp
