#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
    constexpr double M_PI = 3.14159265358979323846;
#endif
#include <string>
#include <algorithm>
#include <omp.h>

#include "utils.hpp"
#include "grid.hpp"
#include "integrator.hpp"

/**
 * @brief phi4_lpa.hpp
 * 
 * We want to solve the Wetterich equation for φ4 theory in the LPA approximation, i.e. in a zeroth order derivative expansion
 * 
 * The wetterich equation, i.e. the flow of the 1PI-effective action for a constant field in d=3 dimensions where φ4 theory is well defined is
 * 
 *     ∂_t V_k(ρ) = 1/6π² k³/(1 + (V_k'+2ρV_k'')/k²)
 * 
 * where we use the O(N) invariant ρ = 1/2φ²
 * The classical poential at k=Λ is given by
 * 
 *      V_classical(ρ) = m_Λ² ρ + λ_φ,Λ/3!ρ²
 * 
 */

// Paramters ---------------------------

struct Params {
    double m2               = -1./(9*M_PI*M_PI);                // UV scale mass parameter
    double lambda           = 1.0;                              // UV scale coupling constant
    size_t N                = 1;                                // Number of fields
    size_t d                = 3;                                // spacetime dimension
    double t_start          = 0.0;
    double t_end            = -15.0;
    Grid grid;
};


// classical potential -----------------

std::vector<double> V_classical(const Params& p);
double V_min_classical(const Params& p);

// dimensionless classical potential
std::vector<double> u_classical(const Params& p);
double u_min_classical(const Params& p);

// Compute RHS -------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p);
std::vector<double> RHS_dimless(const std::vector<double>& u, double k, const Params& p);

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Params& p);
void save_all(const std::vector<std::vector<double>>& snapshots, const std::vector<std::vector<double>>& rhs_snapshots, const std::vector<double>& k_values, const Params& p, const std::string& filename);
void save_dt_hist(const std::vector<double>& dt_values, const std::vector<double>& k_values, const std::string& filename);


// ----------------------------------
// Refactor the following
// ----------------------------------
 

// Integrate RG time step --------------

// Simple Euler time step
std::vector<double> step(const std::vector<double>& V, double k, double dt, const Params& p);
// RK4 step used by adaptive integrator
std::vector<double> step_rk4(const std::vector<double>& V, double k, double dt, const Params& p);

// Integrate complete RG flow ----------

// Simple dt forward step integrator
void integrate_flow(const std::vector<double>& V_init, double dt, const Params& p, const std::string& filename = "results/flow.csv", int n_snapshots = 100);
// Adaptive integrator (RK4 with step-doubling error estimate)
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p, const std::string& filename = "results/flow_adaptive.csv", int n_snapshots = 100, double absolute_tolerance = 1e-8, double relative_tolerance = 1e-8, bool show_progress_bar = true);