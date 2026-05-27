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

static double M2            = -1.0;     // classical mass paramters
static double LAMBDA        = 2.0;      // classical coupling paramter  

// RG times t = ln(k/Λ)
static double T_START       = 0.0;
static double T_END         = -7.5;

// rho grid
static size_t N_RHO         = 1000;                     // number of rho grid points
static double RHO_MAX       = 7.5;                      // maximal value for rho
static double D_RHO         = RHO_MAX / (N_RHO - 1);    // rho spacing


// Finite difference derivatives -------

inline double dV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[1] - V[0])/D_RHO; }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - V[N_RHO-2])/D_RHO; }
    return (V[i+1] - V[i-1])/(2*D_RHO);
}

inline double ddV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(D_RHO*D_RHO); }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - 2.0*V[N_RHO-2] + V[N_RHO-3])/(D_RHO*D_RHO); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(D_RHO*D_RHO);
}

// classical potential -----------------
inline std::vector<double> V_classical() {
    std::vector<double> V;
    V.reserve(N_RHO);
    double rho;

    for (size_t i = 0; i < N_RHO; ++i) {
        rho = i * D_RHO;
        V.push_back(M2*rho + LAMBDA/6.0*rho*rho);
    }

    return V;
}

// Compute RHS -------------------------

inline std::vector<double> RHS(const std::vector<double>& V, double k) {
    std::vector<double> RHS(N_RHO);

    double prefactor = k*k*k/(6.0*M_PI*M_PI);
    // std::cout << "k = " << k << std::endl;

    for (size_t i = 0; i < N_RHO; ++i) {
        const double rho = i * D_RHO;
        double denom = (1.0 + (dV(V, i) + 2.0*rho*ddV(V, i))/(k*k));
        if (k == 0) {
            denom = 1.0;
        }
        if (std::abs(denom) < 1e-5) {
            std::cerr << "[WARNING] |denom| = " << denom << "rho = " << rho << std::endl;
            denom = 1e-5;
        }
        RHS[i] = prefactor / denom;
    }

    return RHS;
}

// Integrate RG time step --------------

inline std::vector<double> integrate(const std::vector<double>& V, double k, double dt) {
    std::vector<double> V_next(N_RHO);
    if (dt >= 0) {
        std::cerr << "[ERROR] dt must be negative" << std::endl;
        return V_next;
    }

    std::vector<double> RHS_vals = RHS(V, k);
    // std::cout << "RHS_0 = " << RHS_vals[0] << "\nRHS_15 = " << RHS_vals[15] <<  std::endl;
    for (size_t i = 0; i < N_RHO; ++i) {
        V_next[i] = V[i] + dt * RHS_vals[i];
    }
    return V_next;
}

// save current potential --------------

inline void save_V(const std::vector<double>& V, const std::string filename) {
    std::ofstream file(filename);

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << i*D_RHO << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

// Integrate complete RG flow ----------

inline void integrate_flow(const std::vector<double>& V_init, double dt) {
    if (dt >= 0) {
        std::cerr << "[ERROR] dt must be negative" << std::endl;
        return;
    }
    const double total_t = T_START - T_END;
    size_t N = static_cast<size_t>(std::ceil(total_t / std::abs(dt))) + 1;
    double dt_t = -total_t / (N - 1);
    std::vector<double> V = V_init;
    for (size_t i = 0; i < N; ++i) {
        double t = T_START + i*(T_END - T_START)/(N - 1);
        double k = exp(t);

        // if (i < 200) {
        // double denom0 = 1.0 + (dV(V,0) + 0.0*ddV(V,0)) / (k*k);
        // std::cerr << "i=" << i << " k=" << k 
        //           << " denom0=" << denom0 
        //           << " V[0]=" << V[0] 
        //           << " V[1]=" << V[1] << "\n";
        // }

        if (i % 100 == 0) {
            save_V(V, "results/V_" + std::to_string(std::abs(t)) + ".txt");
        }
        if (i + 1 < N) {
            V = integrate(V, k, dt_t);
        }
    }
}



namespace HEATEQ
{

/**
 * @brief 
 * 
 * We want to solve the 1+1 dimemsional Heat equation
 * 
 * (∂t - ∂x²) ϕ = 0
 * 
 * With initial condition
 * 
 * ϕ(t=0) = cos(x)
 * 
 */

// constants -------------------------

static double T_I       = 0;                // initial time
static double T_F       = 20;               // final time

static double X_MIN     = -2*M_PI;              
static double X_MAX     = 2*M_PI;
static double N_X       = 100;
static double DX        = (X_MAX - X_MIN)/(N_X - 1);

static std::vector<double> X_GRID = []() {
    std::vector<double> x_grid(N_X);
    for (size_t i = 0; i < N_X; ++i) {
        x_grid[i] = X_MIN + i * DX;
    }
    return x_grid;
}();



// initial ϕ(t=0)
inline std::vector<double> phi0() {
    std::vector<double> phi(N_X);
    for (size_t i = 0; i < N_X; ++i) {
        const double x = X_GRID[i];
        phi[i] = sin(x);
        // phi[i] = exp(-x*x);
    }
    return phi;
}

// Finite difference derivatives -------

inline double ddphi(const std::vector<double>& phi, size_t i) {
    if (i == 0) { return (phi[2] - 2.0*phi[1] + phi[0])/(DX*DX); }
    if (i == N_X - 1) { return (phi[N_X-1] - 2.0*phi[N_X-2] + phi[N_X-3])/(DX*DX); }
    return (phi[i+1] - 2.0*phi[i] + phi[i-1])/(DX*DX);
}

// Integrate time step -----------------

inline std::vector<double> integrate_phi(const std::vector<double>& phi, double dt) {
    std::vector<double> phi_next(N_X);
    if (dt < 0) {
        std::cerr << "[ERROR] dt must be positive" << std::endl;
        return phi_next;
    }
    for (size_t i = 0; i < N_X; ++i) {
        if (i == 0) {
            phi_next[i] = phi[i];
        }
        else if (i == N_X - 1) {
            phi_next[N_X - 1] = phi[N_X - 1];
        }
        else {
        phi_next[i] = phi[i] + dt * ddphi(phi, i);
        }
    }
    return phi_next;
}

inline void save_phi(const std::vector<double>& phi, const std::string filename) {
    std::ofstream file(filename);

    for (size_t i = 0; i < phi.size(); ++i) {
        file << X_GRID[i] << ", " << phi[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

inline void flow_phi(const std::vector<double>& phi_init, double dt) {
    if (dt < 0) {
        std::cerr << "[ERROR] dt must be positive" << std::endl;
        return;
    }
    std::vector<double> phi = phi_init;
    const double total_t = T_F - T_I;
    size_t N = static_cast<size_t>(std::ceil(total_t / std::abs(dt))) + 1;
    for (size_t i = 0; i < N; ++i) {
        double t = T_START + i*(T_END - T_START)/(N - 1);
        if (i % 1000 == 0) {
            save_phi(phi, "results/HE/phi_" + std::to_string(std::abs(t)) + ".txt");
        }
        if (i + 1 < N) {
            phi = integrate_phi(phi, dt);
        }
    }
}

    
} // namespace HEATEQ
