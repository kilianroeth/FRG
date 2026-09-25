#include "phi4_lpa.hpp"

namespace phi4 {

// classical potential -----------------
std::vector<double> V_classical(const Params& p, const Grid& grid) {
    std::vector<double> V(grid.n_rho());
    double rho;
    for(size_t i = 0; i < grid.n_rho(); ++i) {
        rho = i * grid.d_rho();
        V[i] = p.m2 * rho + p.lambda / 6.0 * rho * rho;
    }
    return V;
}

double V_min_classical(const Params& p) {
    return -3.0 * p.m2 / p.lambda;
}

std::vector<double> u_classical(const Params& p, const Grid& grid) {
    std::vector<double> u(grid.n_rho());
    double rho;
    double k = exp(p.t_start);

    for(size_t i = 0; i < grid.n_rho(); ++i) {
        rho = std::pow(k, p.d - 2) * i * grid.d_rho();
        u[i] = std::pow(k, -p.d) * (p.m2 * rho + p.lambda / 6.0 * rho * rho);
    }

    return u;
}

double u_min_classical(const Params& p) {
    double k = exp(p.t_start);
    return -std::pow(k, -p.d) * 3.0 * p.m2 / p.lambda;
}

// Compute RHS -------------------------

// computes RHS of Wetterich equation
void RHS(const std::vector<double>& V, double t, std::vector<double>& out, const Params& p,
         const Grid& grid) {
    std::vector<double> RHS_vals(grid.n_rho());

    const double k = std::exp(t);
    const double k2 = k * k;
    const size_t N_grid = grid.n_rho();
    const double prefactor = Ω(p.d) / std::pow(2 * M_PI, p.d) * std::pow(k, p.d + 2.) / p.d;

    bool warning_triggered = false;

#pragma omp parallel for schedule(static)
    for(size_t i = 0; i < N_grid; ++i) {
        const double rho = grid.rho_vals(i);
        const double d1 = grid.d1(V, i);
        const double d2 = grid.d2(V, i);

        // goldstone propagator
        double denom_goldstone = k2 + d1;
        if(!std::isfinite(denom_goldstone) || std::abs(denom_goldstone) < 1e-13) {
            warning_triggered = true;
            denom_goldstone = 1e-13;
        }

        // massive propagator
        double denom_massive = k2 + d1 + 2.0 * rho * d2;
        if(!std::isfinite(denom_massive) || std::abs(denom_massive) < 1e-13) {
            warning_triggered = true;
            denom_massive = 1e-13;
        }

        out[i] = prefactor * ((p.N - 1.) / denom_goldstone + 1. / denom_massive);

        if(warning_triggered && p.warning_level >= 1) {
            std::cerr << "[WARNING] Small/non-finite propagator denominator detected at t = " << t
                      << "\n";
        }
    }
}

// dimensionless quantities
// ̄ρ = k^2-d ρ
// u = k^-d V_k(k^d-2 ̄ρ)
// computes RHS of Wetterich equation minus all terms of the LHS that is not the RG-time derivative
void RHS_dimless(const std::vector<double>& u, const Params& p, const Grid& grid) {
    std::vector<double> RHS_vals(grid.n_rho());
    std::vector<double> LHS_remainder(grid.n_rho());
    std::vector<double> goldstone_propagator(grid.n_rho());
    std::vector<double> massive_propagator(grid.n_rho());

    const double d = static_cast<double>(p.d);
    const double N = static_cast<double>(p.N);

    double prefactor = Ω(d) / (d * std::pow(2 * M_PI, d));

    for(size_t i = 0; i < grid.n_rho(); ++i) {
        double du = grid.d1(u, i);
        double ddu = grid.d2(u, i);

        const double rho = i * grid.d_rho();
        // LHS remainings
        LHS_remainder[i] = d * u[i] + (2 - d) * rho * du;

        // goldstone modes
        double goldstone_denom = 1 + du;
        if(!std::isfinite(goldstone_denom) || std::abs(goldstone_denom) < 1e-12) {
            if(p.warning_level >= 1) {
                std::cerr << "[WARNING] |goldstone denom| = " << abs(goldstone_denom)
                          << ", rho = " << rho << std::endl;
            }
            goldstone_denom = 1e-12;
        }
        goldstone_propagator[i] = (N > 1) ? (N - 1) / goldstone_denom : 0.0;

        // massive modes
        double massive_denom = 1 + du + 2 * rho * ddu;
        if(!std::isfinite(massive_denom) || std::abs(massive_denom) < 1e-12) {
            if(p.warning_level >= 1) {
                std::cerr << "[WARNING] |massive denom| = " << abs(massive_denom)
                          << ", rho = " << rho << std::endl;
            }
            massive_denom = 1e-12;
        }
        massive_propagator[i] = 1. / massive_denom;

        RHS_vals[i] =
            -LHS_remainder[i] + prefactor * (goldstone_propagator[i] + massive_propagator[i]);
    }
}

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Grid& grid) {
    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
    }

    file << "ρ = 1/2 φ², V(ρ)\n";
    for(size_t i = 0; i < V.size(); ++i) {
        file << grid.rho_vals(i) << ", " << V[i] << "\n";
    }
    file << std::endl;

    file.close();
}

void save_all(const std::vector<std::vector<double>>& snapshots,
              const std::vector<std::vector<double>>& rhs_snapshots,
              const std::vector<double>& k_values, const Params& p, const Grid& grid,
              const std::string& filename) {
    if(filename.empty()) {
        return;
    }
    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
        return;
    }

    // metadata
    file << "# Wetterich LPA flow, phi^4, d=3, N=";
    file << p.N << "\n";
    file << "# m2 = " << p.m2 << ", lambda = " << p.lambda << "\n";
    file << "# rho_max = " << grid.rho_max() << ", n_rho = " << grid.n_rho() << "\n";

    // V block
    file << "# block: V\n";
    file << "rho";
    for(double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    file << std::scientific << std::setprecision(10);
    for(size_t i = 0; i < grid.n_rho(); ++i) {
        file << i * grid.d_rho();
        for(const auto& V : snapshots)
            file << ", " << V[i];
        file << "\n";
    }

    // RHS
    file << "# block: RHS\n";
    file << "rho";
    for(double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    for(size_t i = 0; i < grid.n_rho(); ++i) {
        file << i * grid.d_rho();
        for(const auto& R : rhs_snapshots)
            file << ", " << R[i];
        file << "\n";
    }

    file.close();

    std::cout << "Saved " << snapshots.size() << " snapshots -> " << filename << "\n";
}

void save_dt_hist(const std::vector<double>& dt_values, const std::vector<double>& k_values,
                  const std::string& filename) {
    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
        return;
    }

    // metadata
    file << "# time steps dt for adaptive RK4 time stepper\n";
    file << "# k, dt\n";
    // data
    for(size_t i = 0; i < dt_values.size(); ++i) {
        file << k_values[i] << ", " << dt_values[i] << "\n";
    }

    file.close();
}

// Integrate complete RG flow ----------

// Adaptive integrator using RK4 + step-doubling
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p,
                             const Grid& grid, const StepperConfig& cfg,
                             const std::string& filename, int n_snapshots) {
    if(dt_init >= 0) {
        std::cerr << "[ERROR] dt_init must be negative" << std::endl;
        return;
    }

    std::cout << "Solving flow equation with adaptive time step...\n";

    const RHSfunc rhs_func = [&p, &grid](const std::vector<double>& state, double t,
                                         std::vector<double>& out) { RHS(state, t, out, p, grid); };

    std::vector<double> snap_targets(n_snapshots);
    for(int s = 0; s < n_snapshots; ++s) {
        double fraction = static_cast<double>(s) / (n_snapshots - 1);
        snap_targets[s] = p.t_start + fraction * (p.t_end - p.t_start);
    }

    std::vector<std::pair<double, std::vector<double>>> snapshot_pairs;
    std::vector<std::pair<double, double>> step_history;

    integrate_adaptive(V_init, p.t_start, p.t_end, dt_init, rhs_func, cfg, snap_targets,
                       &snapshot_pairs, &step_history);

    std::vector<std::vector<double>> snapshots;
    std::vector<std::vector<double>> rhs_snapshots;
    std::vector<double> k_values;
    snapshots.reserve(snapshot_pairs.size());
    rhs_snapshots.reserve(snapshot_pairs.size());
    k_values.reserve(snapshot_pairs.size());

    std::vector<double> rhs_buf(grid.n_rho());

    for(const auto& [t, V] : snapshot_pairs) {
        snapshots.push_back(V);
        RHS(V, t, rhs_buf, p, grid);
        rhs_snapshots.push_back(rhs_buf);
        k_values.push_back(std::exp(t));
    }

    std::vector<double> dt_values;
    std::vector<double> dt_k_values;
    dt_values.reserve(step_history.size());
    dt_k_values.reserve(step_history.size());
    for(const auto& [t, dt] : step_history) {
        dt_values.push_back(dt);
        dt_k_values.push_back(std::exp(t));
    }

    save_all(snapshots, rhs_snapshots, k_values, p, grid, filename);
    if(!filename.empty()) {
        save_dt_hist(
            dt_values, dt_k_values,
            "results/phi4/dt_values.txt"); // TODO use the first part of string of filename of
                                           // integrate_flow_adaptive to have a dynamical directory
    }
}

} // namespace phi4
