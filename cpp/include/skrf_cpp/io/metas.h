#pragma once

#include <string>
#include "../network_set.h"

namespace skrf_cpp {
namespace io {

// Write a minimal SDATCV file from a NetworkSet (METAS format). This is
// a pragmatic implementation that writes frequency and mean S values.
void ns_to_sdatcv(const NetworkSet &ns, const std::string &fname, bool polar=false);

} // namespace io
} // namespace skrf_cpp
