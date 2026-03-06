#pragma once

// Single-include umbrella header for the C++ skrf port.
// Include this to pull in core APIs (Network, NetworkEigen, Touchstone,
// Media models, VectorFitting, NetworkSet, Calibration utilities, and helpers).

#include "network.h"
#include "network_eigen.h"
#include "network_set.h"
#include "touchstone.h"
#include "media.h"
#include "util.h"
#include "transforms.h"
#include "calibration.h"
#include "calibration_set.h"
#include "vector_fitting_advanced.h"
// I/O helpers
#include "io/general.h"
#include "io/csv.h"
#include "io/mdif.h"
#include "io/citi.h"
#include "io/metas.h"

// Convenience namespace aliasing not needed; users can include this single file.
