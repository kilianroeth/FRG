#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846;
#endif
#include <algorithm>
#include <omp.h>
#include <string>

#include "grid.hpp"
#include "integrator.hpp"
#include "utils.hpp"

/**
 * @brief quark_meson.hpp
 *
 * We want to solve the flow equation for the effective potential for the quark meson model in
 * vacuum on the LPA approximation
 *
 */

namespace QM {

// Paramters ---------------------------

struct Params {
    double m2 = -1. / (9 * M_PI * M_PI); // UV scale mass parameter
    double lambda = 1.0;                 // UV scale coupling constant
    double h = 0.5;                      // Yukawa coupling
    int N = 1;                           // Number of fields
    int Nc = 3;                          // Number of colors (SU(NC), Nc=3 for QCD)
    double d = 3;                        // spacetime dimension
    double t_start = 0.0;
    double t_end = -15.0;
    Grid grid;
    int warning_level = 0;
};

// Compute RHS ------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p);

} // namespace QM