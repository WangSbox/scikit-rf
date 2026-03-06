#pragma once

#include <string>
#include <vector>
#include "../network.h"

namespace skrf_cpp {
namespace io {

// Minimal MDIF reader: parses MDIF file and returns a vector of Networks
// (each Network corresponds to one data block). This is a simplified
// but practical implementation covering common MDIF exports.
std::vector<Network> read_mdif(const std::string &path);

} // namespace io
} // namespace skrf_cpp
