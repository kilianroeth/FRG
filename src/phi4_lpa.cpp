#include "phi4_lpa.hpp"


// classical potential -----------------
std::vector<double> V_classical(const Params& p) {
    std::vector<double> V(p.grid.n_rho());
    double rho;

    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        rho = i * p.grid.d_rho();
        V[i] = p.m2*rho + p.lambda/6.0*rho*rho;
    }

    return V;
}

double V_min_classical(const Params& p) {
    return -3.0*p.m2/p.lambda;
}


std::vector<double> u_classical(const Params& p) {
    std::vector<double> u(p.grid.n_rho());
    double rho;
    double k = exp(p.t_start);

    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        rho = std::pow(k,p.d - 2) * i * p.grid.d_rho();
        u[i] = std::pow(k, -p.d) * (p.m2*rho + p.lambda/6.0*rho*rho);
    }

    return u;
}

double u_min_classical(const Params& p) {
    double k = exp(p.t_start);
    return -std::pow(k, -p.d) * 3.0*p.m2/p.lambda;
}


// Compute RHS -------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p) {
    std::vector<double> RHS_vals(p.grid.n_rho());

    double prefactor = k*k*k/(6.0*M_PI*M_PI);
    // std::cout << "k = " << k << std::endl;

    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        const double rho = i * p.grid.d_rho();
        double denom = (1.0 + (p.grid.d1(V, i) + 2.0*rho*p.grid.d2(V, i))/(k*k));
        if (k == 0) {
            denom = 1.0;
        }
        if (!std::isfinite(denom) || std::abs(denom) < 1e-12) {
            #pragma omp critical
            std::cerr << "[WARNING] |denom| = " << denom << ", rho = " << rho << std::endl;
            denom = 1e-12;
        }
        // if (denom < 0) {
        //     std::cout << "[WARNING] Γ(2) < 0 at k = " << k << ",  = " << rho << "\n";
        // }
        RHS_vals[i] = prefactor / denom;
    }
    return RHS_vals;
}

// dimensionless quantities
// ̄ρ = k^2-d ρ
// u = k^-d V_k(k^d-2 ̄ρ)
std::vector<double> RHS_dimless(const std::vector<double>& u, const Params& p) {
    std::vector<double> RHS_vals(p.grid.n_rho());
    std::vector<double> RHS_remainder(p.grid.n_rho());
    std::vector<double> goldstone_propagator(p.grid.n_rho());
    std::vector<double> massive_propagator(p.grid.n_rho());

    double prefactor = Ω(p.d)/(p.d * std::pow(2*M_PI, p.d));

    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        double du = p.grid.d1(u, i);
        double ddu = p.grid.d2(u, i);

        const double rho = i * p.grid.d_rho();
        // LHS remainings
        RHS_remainder[i] = p.d * u[i] + (2 - p.d) * rho * du;
        
        // goldstone modes
        double goldstone_denom = 1 + du;
        if (!std::isfinite(goldstone_denom) || std::abs(goldstone_denom) < 1e-12) {
            #pragma omp critical 
            std::cerr << "[WARNING] |goldstone denom| = " << goldstone_denom << ", rho = " << rho << std::endl;
            goldstone_denom = 1e-12;
        }
        goldstone_propagator[i] = (p.N - 1)/goldstone_denom;

        // massive modes
        double massive_denom = 1 + du + 2 * rho * ddu;
        if (!std::isfinite(massive_denom) || std::abs(massive_denom) < 1e-12) {
            #pragma omp critical 
            std::cerr << "[WARNING] |massive denom| = " << massive_denom << ", rho = " << rho << std::endl;
            massive_denom = 1e-12;
        }
        massive_propagator[i] = 1. / massive_denom;


        RHS_vals[i] = -RHS_remainder[i] + prefactor * (goldstone_propagator[i] + massive_propagator[i]);
    }
    return RHS_vals;
}



// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Params& p) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; }

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << p.grid.rho_vals(i) << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

void save_all(const std::vector<std::vector<double>>& snapshots, const std::vector<std::vector<double>>& rhs_snapshots, const std::vector<double>& k_values, const Params& p, const std::string& filename) {
    if (filename.empty()) { return; }
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; return; }

    // metadata
    file << "# Wetterich LPA flow, phi^4, d=3\n";
    file << "# m2 = " << p.m2 << ", lambda = " << p.lambda << "\n";
    file << "# rho_max = " << p.grid.rho_max() << ", n_rho = " << p.grid.n_rho() << "\n";

    // V block
    file << "# block: V\n";
    file << "rho";
    for (double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    file << std::scientific << std::setprecision(10);
    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        file << i * p.grid.d_rho();
        for (const auto& V : snapshots)
            file << ", " << V[i];
        file << "\n";
    }

    // RHS
    file << "# block: RHS\n";
    file << "rho";
    for (double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    for (size_t i = 0; i < p.grid.n_rho(); ++i) {
        file << i * p.grid.d_rho();
        for (const auto& R : rhs_snapshots)
            file << ", " << R[i];
        file << "\n";
    }

    file.close();

    std::cout << "Saved " << snapshots.size() << " snapshots -> " << filename << "\n";
}

void save_dt_hist(const std::vector<double>& dt_values, const std::vector<double>& k_values, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open" << filename << "\n"; return; }
    
    // metadata
    file << "# time steps dt for adaptive RK4 time stepper\n";
    file << "# k, dt\n";
    // data
    for (size_t i = 0; i < dt_values.size(); ++i) {
        file << k_values[i] << ", " << dt_values[i] << "\n";
    }

    file.close();
}


// Integrate complete RG flow ----------

void integrate_flow(const std::vector<double>& V_init, double dt, const Params& p, const std::string& filename, int n_snapshots) {
    if (dt >= 0) {
        std::cerr << "[ERROR] dt must be negative" << std::endl;
        return;
    }
    std::cout << "Solving flow equation...\n";

    const RHSfunc rhs = [&p](const std::vector<double>& state, double t) {
        return RHS(state, std::exp(t), p);
    };

    const double total_t = p.t_start - p.t_end;
    size_t N = static_cast<size_t>(std::ceil(total_t / std::abs(dt))) + 1;
    double dt_t = -total_t / (N - 1);

    std::vector<size_t> snap_indices;
    for (int s = 0; s < n_snapshots; ++s) {
        double frac = static_cast<double>(s) / (n_snapshots - 1);
        size_t idx = static_cast<size_t>(std::round(frac * (N - 1)));
        snap_indices.push_back(idx);
    }

    std::vector<std::vector<double>> snapshots, rhs_snapshots;
    std::vector<double> k_values;

    std::vector<double> V = V_init;
    size_t next_snap = 0;
    for (size_t i = 0; i < N; ++i) {
        double t = p.t_start + i*(p.t_end - p.t_start)/(N - 1);
        double k = exp(t);
        progressBar(i + 1, N);

        if (next_snap < snap_indices.size() && i == snap_indices[next_snap]) {
            snapshots.push_back(V);
            rhs_snapshots.push_back(RHS(V, k, p));
            k_values.push_back(k);
            ++next_snap;
        }

        if (i + 1 < N) {
            V = step_euler(V, t, dt_t, rhs);
        }
    }
    save_all(snapshots, rhs_snapshots, k_values, p, filename);
}


// Adaptive integrator using RK4 + step-doubling
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p, const StepperConfig& cfg, const std::string& filename, int n_snapshots) {
    if (dt_init >= 0) {
        std::cerr << "[ERROR] dt_init must be negative" << std::endl;
        return;
    }

    std::cout << "Solving flow equation with adaptive time step...\n";

    const RHSfunc rhs = [&p](const std::vector<double>& state, double t) {
        return RHS(state, std::exp(t), p);
    };

    std::vector<double> snap_targets(n_snapshots);
    for (int s = 0; s < n_snapshots; ++s) {
        double fraction = static_cast<double>(s) / (n_snapshots - 1);
        snap_targets[s] = p.t_start + fraction * (p.t_end - p.t_start);
    }

    std::vector<std::pair<double, std::vector<double>>> snapshot_pairs;
    std::vector<std::pair<double, double>> step_history;

    integrate_adaptive(V_init, p.t_start, p.t_end, dt_init, rhs, cfg, snap_targets, &snapshot_pairs, &step_history);

    std::vector<std::vector<double>> snapshots;
    std::vector<std::vector<double>> rhs_snapshots;
    std::vector<double> k_values;
    snapshots.reserve(snapshot_pairs.size());
    rhs_snapshots.reserve(snapshot_pairs.size());
    k_values.reserve(snapshot_pairs.size());

    for (const auto& [t, V] : snapshot_pairs) {
        snapshots.push_back(V);
        rhs_snapshots.push_back(RHS(V, std::exp(t), p));
        k_values.push_back(std::exp(t));
    }

    std::vector<double> dt_values;
    std::vector<double> dt_k_values;
    dt_values.reserve(step_history.size());
    dt_k_values.reserve(step_history.size());
    for (const auto& [t, dt] : step_history) {
        dt_values.push_back(dt);
        dt_k_values.push_back(std::exp(t));
    }

    save_all(snapshots, rhs_snapshots, k_values, p, filename);
    if (!filename.empty()) {
        save_dt_hist(dt_values, dt_k_values, "results/dt_values.txt");
    }
}

