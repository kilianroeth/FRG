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
#include <chrono>
#include <omp.h>
#include <string>

#include "grid.hpp"
#include "integrator.hpp"
#include "phi4_lpa.hpp"
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

// QM observables
struct Observables {
    double rho0 = 0.;
    double m2_sigma = 0.;
    double m2_pi = 0.;
};

// classical potential -----------------

std::vector<double> V_classical(const Params& p);
double V_min_classical(const Params& p);

// Compute RHS ------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p);

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Params& p);
void save_all(const std::vector<std::vector<double>>& snapshots,
              const std::vector<std::vector<double>>& rhs_snapshots,
              const std::vector<double>& k_values, const Params& p, const std::string& filename);
void save_dt_hist(const std::vector<double>& dt_values, const std::vector<double>& k_values,
                  const std::string& filename);

// Adaptive integrator (RK4 with step-doubling error estimate)
// dt init must be negative
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p,
                             const StepperConfig& cfg,
                             const std::string& filename = "results/QM/flow_adaptive.csv",
                             int n_snapshots = 100);

// Compute observables of the QM model
Observables compute_observables(std::vector<double>& V, const Params& p);

// Sweeps UV parameters m²_Λ, λ_Λ, h_Λ and compute ρ at k=0
void sweep_params(const std::vector<double>& m2, const std::vector<double>& lambda,
                  const std::vector<double>& h, const Params& p, const StepperConfig& cfg,
                  const std::string& filename = "results/QM/UV_params.csv");

} // namespace QM