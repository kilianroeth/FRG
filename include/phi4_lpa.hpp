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
    double m2               = -1./(9*M_PI*M_PI);
    double lambda           = 1.0;
    double t_start          = 0.0;
    double t_end            = -15.0;
    size_t n_rho            = 1000;
    double rho_max          = 0.5;
    double drho() const     { return rho_max / (n_rho - 1); } 
    double rho_at(size_t i) const { return i * drho(); }
};


// Finite difference derivatives -------

double dV(const std::vector<double>& V, size_t i, const Params& p);
double ddV(const std::vector<double>& V, size_t i, const Params& p);

// classical potential -----------------

std::vector<double> V_classical(const Params& p);
double V_min_classical(const Params& p);

// Compute RHS -------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p);

// Integrate RG time step --------------

std::vector<double> step(const std::vector<double>& V, double k, double dt, const Params& p);

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Params& p);
void save_all(const std::vector<std::vector<double>>& snapshots, const std::vector<std::vector<double>>& rhs_snapshots, const std::vector<double>& k_values, const Params& p, const std::string& filename);

// Integrate complete RG flow ----------

void integrate_flow(const std::vector<double>& V_init, double dt, const Params& p, const std::string& filename = "results/flow.csv", int n_snapshots = 100);
