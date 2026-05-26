#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

/**
 * @brief phi4_lpa.cpp
 * 
 * We want to solve the Wetterich equation for φ4 theory in the LPA approximation, i.e. in a zeroth order derivative expansion
 * 
 * The wetterich equation, i.e. the flow of the 1PI-effective action for a constant field in d=3 dimensions where φ4 theory is well defined is
 * 
 *     ∂_t V_k(ρ) = 1/6π k³/(1 + (V_k'+2ρV_k'')/k²)
 * 
 * where we use the O(N) invariant ρ = 1/2φ²
 * The classical poential at k=Λ is given by
 * 
 *      V_classical(ρ) = m_Λ² ρ + λ_φ,Λ/3!ρ²
 * 
 */

// Paramters ---------------------------

static double M2            = -1.0;     // classical mass paramters
static double LAMBDA        = 2.0;      // classical coupling paramter  

// RG times t = ln(k/Λ)
static double T_START       = 0.0;
static double T_END         = -15.0;

// rho grid
static size_t N_RHO         = 500;      // number of rho grid points
static double RHO_MAX       = 10.0;  // maximal value for rho



// Finite difference derivatives -------

static double d_rho = RHO_MAX / (N_RHO - 1);

inline double dV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[1] - V[0])/d_rho; }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - V[N_RHO-2])/d_rho; }
    return (V[i+1] - V[i-1])/(2*d_rho);
}

inline double ddV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(d_rho*d_rho); }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - 2.0*V[N_RHO-2] + V[N_RHO-3])/(d_rho*d_rho); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(d_rho*d_rho);
}


// Compute RHS -------------------------

inline std::vector<double> V_classical() {
    std::vector<double> V;
    V.reserve(N_RHO);
    double rho;

    for (size_t i = 0; i < N_RHO; ++i) {
        rho = i * d_rho;
        V.push_back(M2*rho + LAMBDA/6.0*rho*rho);
    }

    return V;
}

inline std::vector<double> RHS(std::vector<double> V, double k) {
    std::vector<double> RHS(N_RHO);

    double prefactor = std::pow(k,3.0)/(6.0*M_PI*M_PI);
    double rho;

    for (size_t i = 0; i < N_RHO; ++i) {
        rho = i * d_rho;
        RHS[i] = prefactor * 1.0 / (1.0 + (dV(V, i) + 2.0*rho*ddV(V, i))/(k*k));
    }

    return RHS;
}

// Integrate RG time step --------------

inline std::vector<double> integrate(const std::vector<double>& V, double k, double dt) {
    std::vector<double> V_next(N_RHO);

    std::vector<double> RHS_vals = RHS(V, k);
    for (size_t i = 0; i < N_RHO; ++i) {
        V_next[i] = V[i] +  dt * RHS_vals[i];
    }
    return V_next;
}

// save current potential --------------

inline void save_V(const std::vector<double>& V, const std::string filename) {
    std::ofstream file(filename);

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << i*d_rho << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

// Integrate complete RG flow ----------

inline void integrate_flow(const std::vector<double>& V_init, double dt) {
    size_t N = static_cast<int>( std::ceil(T_START - T_END)/ std::abs(dt));
    std::vector<double> V = V_init;
    for (size_t i = 0; i < N; ++i) {
        double t = T_START + i*(T_END - T_START)/(N - 1);
        double k = exp(t);
        V = integrate(V, k, dt);
        if (i % 10 == 0) {
            save_V(V, "results/V_" + std::to_string(std::abs(t)) + ".txt");
        }
    }
}