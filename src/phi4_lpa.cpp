#include "phi4_lpa.hpp"


// Finite difference derivatives -------

double dV(const std::vector<double>& V, size_t i, const Params& p) {
    if (i == 0) { return (V[1] - V[0])/p.drho(); }
    if (i == p.n_rho - 1) { return (V[p.n_rho-1] - V[p.n_rho-2])/p.drho(); }
    return (V[i+1] - V[i-1])/(2*p.drho());
}

double ddV(const std::vector<double>& V, size_t i, const Params& p) {
    if (i == 0) { return (V[2] - 2.0*V[1] + V[0])/(p.drho()*p.drho()); }
    if (i == p.n_rho - 1) { return (V[p.n_rho-1] - 2.0*V[p.n_rho-2] + V[p.n_rho-3])/(p.drho()*p.drho()); }
    return (V[i+1] - 2.0*V[i] + V[i-1])/(p.drho()*p.drho());
}

// classical potential -----------------
std::vector<double> V_classical(const Params& p) {
    std::vector<double> V;
    V.reserve(p.n_rho);
    double rho;

    for (size_t i = 0; i < p.n_rho; ++i) {
        rho = i * p.drho();
        V.push_back(p.m2*rho + p.lambda/6.0*rho*rho);
    }

    return V;
}

double V_min_classical(const Params& p) {
    return -6.0*p.m2/p.lambda;
}

// Compute RHS -------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p) {
    std::vector<double> RHS_vals(p.n_rho);

    double prefactor = k*k*k/(6.0*M_PI*M_PI);
    // std::cout << "k = " << k << std::endl;

    for (size_t i = 0; i < p.n_rho; ++i) {
        const double rho = i * p.drho();
        double denom = (1.0 + (dV(V, i, p) + 2.0*rho*ddV(V, i, p))/(k*k));
        if (k == 0) {
            denom = 1.0;
        }
        if (std::abs(denom) < 1e-5) {
            std::cerr << "[WARNING] |denom| = " << denom << "rho = " << rho << std::endl;
            denom = 1e-5;
        }
        // if (denom < 0) {
        //     std::cout << "[WARNING] Γ(2) < 0 at k = " << k << ",  = " << rho << "\n";
        // }
        RHS_vals[i] = prefactor / denom;
    }

    return RHS_vals;
}

// Integrate RG time step --------------

std::vector<double> step(const std::vector<double>& V, double k, double dt, const Params& p) {
    std::vector<double> V_next(p.n_rho);
    if (dt >= 0) {
        std::cerr << "[ERROR] dt must be negative" << std::endl;
        return V_next;
    }

    std::vector<double> RHS_vals = RHS(V, k, p);
    // std::cout << "RHS_0 = " << RHS_vals[0] << "\nRHS_15 = " << RHS_vals[15] <<  std::endl;
    for (size_t i = 0; i < p.n_rho; ++i) {
        V_next[i] = V[i] + dt * RHS_vals[i];
    }
    return V_next;
}

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string filename, const Params& p) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; }

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << p.rho_at(i) << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

void save_all(const std::vector<std::vector<double>>& snapshots, const std::vector<double>& k_values, const Params& p, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; }

    // metadata
    file << "# Wetterich LPA flow, phi^4, d=3\n";
    file << "# m2 = " << p.m2 << ", lambda = " << p.lambda << "\n";
    file << "# rho_max = " << p.rho_max << ", n_rho = " << p.n_rho << "\n";

    // header
    file << "rho";
    for (double k : k_values) {
        file << ", k = " << std::fixed << std::setprecision(6) << k;
    }
    file << "\n";

    file << std::scientific << std::setprecision(10);

    // data rows
    for (size_t i = 0; i < p.n_rho; ++i) {
        file << p.rho_at(i);
        for (const auto& V : snapshots) {
            file << ", " << V[i];
        }
        file << "\n";
    }

    std::cout << "Saved " << snapshots.size() << " snapshots -> " << filename << "\n";

}

// Integrate complete RG flow ----------

void integrate_flow(const std::vector<double>& V_init, double dt, const Params& p, const std::string& filename, int n_snapshots) {
    if (dt >= 0) {
        std::cerr << "[ERROR] dt must be negative" << std::endl;
        return;
    }

    const double total_t = p.t_start - p.t_end;
    size_t N = static_cast<size_t>(std::ceil(total_t / std::abs(dt))) + 1;
    double dt_t = -total_t / (N - 1);

    std::vector<size_t> snap_indices;
    for (int s = 0; s < n_snapshots; ++s) {
        double frac = static_cast<double>(s) / (n_snapshots - 1);
        size_t idx = static_cast<size_t>(std::round(frac * (N - 1)));
        snap_indices.push_back(idx);
    }

    std::vector<std::vector<double>> snapshots;
    std::vector<double> k_values;

    std::vector<double> V = V_init;
    for (size_t i = 0; i < N; ++i) {
        double t = p.t_start + i*(p.t_end - p.t_start)/(N - 1);
        double k = exp(t);

        if (!snap_indices.empty() && i == snap_indices.front()) {
            snapshots.push_back(V);
            k_values.push_back(k);
            snap_indices.erase(snap_indices.begin());
        }

        if (i + 1 < N) {
            V = step(V, k, dt_t, p);
        }
    }
    save_all(snapshots, k_values, p, filename);
}

