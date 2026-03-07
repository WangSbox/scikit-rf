#pragma once

#include "network.h"
#include <string>

namespace skrf_cpp {

class Touchstone {
public:
    // Read a touchstone .sNp file and return a Network
    // If `lenient` is true, parser will skip malformed rows and collect
    // textual warnings in `warnings` (if provided) instead of throwing
    // immediately. Default `lenient=false`.
    static Network load(const std::string &path, bool lenient = false, std::vector<std::string> *warnings = nullptr);
    // Write a touchstone .sNp file (default MA format)
    static void save(const std::string &path, const Network &net, const std::string &format = "MA");
};

}
