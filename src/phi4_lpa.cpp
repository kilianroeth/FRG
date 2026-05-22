#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

/**
 * @brief phi4_lpa.cpp
 * 
 * 
 * 
 */

// --- Paramters ---------------------------

static double M2        = -1.0;     // classical mass paramters
static double LAMBDA    = 0.5;      // classical coupling paramter  

// RG times t = ln(k/Λ)
static double T_START   = 0.0;
static double T_END     = -10.0;

// rho grid
static int N_RHO        = 500;      // number of rho grid points
static int RHO_MAX      = 10.0;     // maximal value for rho



// --- Finite difference derivatives -------

static double d_rho;

inline double dV(std::vector<double> V, int i) {
    if (i == 0) { return (V[1] - V[0])/d_rho; }
    if (i == N_RHO) { return (V[N_RHO] - V[N_RHO-1])/d_rho; }
    return (V[i+1] - V[i-1])/(2*d_rho);
}

inline double ddV(std::vector<double> V, int i) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(d_rho*d_rho); }
    if (i == N_RHO) { return (V[N_RHO] - 2.0*V[N_RHO-1] + V[N_RHO-2])/(d_rho*d_rho); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(d_rho*d_rho);
}


// --- Compute RHS -------------------------

inline std::vector<double> V_classical() {
    std::vector<double> V;
    V.reserve(N_RHO);

    for (size_t i = 0; i < N_RHO; ++i) {
        V.push_back(M2*i/RHO_MAX + LAMBDA/6.0*i/RHO_MAX*i/RHO_MAX);
    }

    return V;
}


// --- Integrate RG time step --------------