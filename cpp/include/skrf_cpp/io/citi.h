#pragma once

#include <string>
#include <vector>
#include "../network.h"

namespace skrf_cpp {
namespace io {

// Minimal CTI reader: parses basic CITI (.cti) files and returns a list of Networks
std::vector<Network> read_citi(const std::string &path);

} // namespace io
} // namespace skrf_cpp
