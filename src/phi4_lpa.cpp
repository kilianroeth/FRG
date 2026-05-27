#include "phi4_lpa.hpp"


// Finite difference derivatives -------

double dV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[1] - V[0])/D_RHO; }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - V[N_RHO-2])/D_RHO; }
    return (V[i+1] - V[i-1])/(2*D_RHO);
}

double ddV(const std::vector<double>& V, size_t i) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(D_RHO*D_RHO); }
    if (i == N_RHO - 1) { return (V[N_RHO-1] - 2.0*V[N_RHO-2] + V[N_RHO-3])/(D_RHO*D_RHO); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(D_RHO*D_RHO);
}

// classical potential -----------------
std::vector<double> V_classical() {
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

std::vector<double> RHS(const std::vector<double>& V, double k) {
    std::vector<double> RHS_vals(N_RHO);

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
        RHS_vals[i] = prefactor / denom;
    }

    return RHS_vals;
}

// Integrate RG time step --------------

std::vector<double> step(const std::vector<double>& V, double k, double dt) {
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

void save_V(const std::vector<double>& V, const std::string filename) {
    std::ofstream file(filename);

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << i*D_RHO << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

// Integrate complete RG flow ----------

void integrate_flow(const std::vector<double>& V_init, double dt) {
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
            V = step(V, k, dt_t);
        }
    }
}

