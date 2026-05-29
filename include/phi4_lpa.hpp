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
 * @brief phi4_lpa.cpp
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
    double m2;
    double λ;

    double t_start;
    double t_end;

    double n_rho;
    double rho_max;
};

static double M2            = -1/(9*M_PI*M_PI);     // classical mass paramters
static double LAMBDA        = 1.0;      // classical coupling paramter  

// RG times t = ln(k/Λ)
static double T_START       = 0.0;
static double T_END         = -15.0;

// rho grid
static size_t N_RHO         = 1000;                     // number of rho grid points
static double RHO_MAX       = 0.5;                      // maximal value for rho
static double D_RHO         = RHO_MAX / (N_RHO - 1);    // rho spacing

class Solver {
    private:

    public:
        Params m_params;
};

// Finite difference derivatives -------

double dV(const std::vector<double>& V, size_t i);
double ddV(const std::vector<double>& V, size_t i);

// classical potential -----------------

std::vector<double> V_classical();

// Compute RHS -------------------------

std::vector<double> RHS(const std::vector<double>& V, double k);

// Integrate RG time step --------------

std::vector<double> step(const std::vector<double>& V, double k, double dt);

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string filename);

// Integrate complete RG flow ----------

void integrate_flow(const std::vector<double>& V_init, double dt);
