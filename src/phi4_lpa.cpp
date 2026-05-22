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

static double M2        = -0.5;     // classical mass paramters
static double LAMBDA    = 1.0;      // classical coupling paramter  

// RG times t = ln(k/Λ)
static double T_START   = 0.0;
static double T_END     = -10.0;

// rho grid
static int N_RHO        = 500;      // number of rho grid points
static int RHO_MAX      = 10.0;     // maximal value for rho

// --- Finite difference derivatives -------

double d_rho;

double dV(std::vector<double> V, int i) {
    if (i == 0) { return (V[1] - V[0])/d_rho; }
    if (i = N_RHO) { return (V[N_RHO] - V[N_RHO-1])/d_rho; }
    return (V[i+1] - V[i-1])/(2*d_rho);
}

double ddV(std::vector<double> V, int i) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(d_rho*d_rho); }
    if (i == N_RHO) { return (V[N_RHO] - 2.0*V[N_RHO-1] + V[N_RHO-2])/(d_rho*d_rho); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(d_rho*d_rho);
}


// --- Compute RHS -------------------------


// --- Integrate RG time step --------------