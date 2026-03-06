#pragma once

#include <string>
#include <vector>
#include <complex>
#include "../network.h"

namespace skrf_cpp {
namespace io {

struct CSVReadResult {
    std::string header;
    std::string comments;
    std::vector<std::vector<double>> data; // rows x cols
};

// Read an Agilent/PN A style CSV header/data block. Returns header, comments and numeric data.
CSVReadResult read_pna_csv(const std::string &path);

// Attempt to convert a PNA-style CSV (DB/deg) into a 2-port Network.
// If parsing fails throws std::runtime_error.
Network pna_csv_to_network(const std::string &path);

} // namespace io
} // namespace skrf_cpp
