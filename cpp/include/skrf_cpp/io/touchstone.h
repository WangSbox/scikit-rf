#pragma once

#include "network.h"
#include <string>

namespace skrf_cpp {

class Touchstone {
public:
    // Read a touchstone .sNp file and return a Network
    static Network load(const std::string &path);
    // Write a touchstone .sNp file (default MA format)
    static void save(const std::string &path, const Network &net, const std::string &format = "MA");
};

}
