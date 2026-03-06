#pragma once

#include <string>
#include <vector>
#include "../network.h"

namespace skrf_cpp {
namespace io {

// Write a Network to disk. If filename has no extension, append ".s2p".
void write_network(const std::string &path, const Network &net, bool overwrite = true);

// Read a Network from a file: supports Touchstone .sNp files.
Network read_network(const std::string &path);

// Read all networks in a directory (searches for *.s*p). Returns map as vector of pairs (name, Network).
std::vector<std::pair<std::string, Network>> read_all_networks(const std::string &dir);

} // namespace io
} // namespace skrf_cpp
