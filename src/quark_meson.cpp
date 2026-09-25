#include "quark_meson.hpp"

namespace QM {

// Compute RHS ------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p) {
    std::vector<double> RHS_vals(p.grid.n_rho());

    phi4::Params phi4_p;
    phi4_p.d = p.d;
    phi4_p.grid = p.grid;
    phi4_p.lambda = p.lambda;
    phi4_p.m2 = p.m2;
    phi4_p.N = p.N;
    phi4_p.t_end = p.t_end;
    phi4_p.t_start = p.t_start;
    phi4_p.warning_level = p.warning_level;

    RHS_vals = phi4::RHS(V, k, phi4_p);
    for(size_t i = 0; i < p.grid.n_rho(); ++i) {
        const double rho = p.grid.rho_vals(i);
        RHS_vals[i] -= 4 * p.Nc * Ω(p.d) / pow(2 * M_PI, p.d) * pow(k, p.d + 2) / p.d * (1) /
                       (k * k + p.h * p.h * rho);
    }

    return RHS_vals;
}

// save current potential --------------
void save_V(const std::vector<double>& V, const std::string& filename, const Params& p) {
    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
    }

    file << "ρ = 1/2 φ², V(ρ)\n";
    for(size_t i = 0; i < V.size(); ++i) {
        file << p.grid.rho_vals(i) << ", " << V[i] << "\n";
    }
    file << std::endl;

    file.close();
}

void save_all(const std::vector<std::vector<double>>& snapshots,
              const std::vector<std::vector<double>>& rhs_snapshots,
              const std::vector<double>& k_values, const Params& p, const std::string& filename) {
    if(filename.empty()) {
        return;
    }
    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
        return;
    }

    // metadata
    file << "# Wetterich QM LPA flow, phi^4, d=3, N=";
    file << p.N << ", Nc = " << p.Nc << "\n";
    file << "# m2 = " << p.m2 << ", lambda = " << p.lambda << "\n";
    file << "# rho_max = " << p.grid.rho_max() << ", n_rho = " << p.grid.n_rho() << "\n";

    // V block
    file << "# block: V\n";
    file << "rho";
    for(double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    file << std::scientific << std::setprecision(10);
    for(size_t i = 0; i < p.grid.n_rho(); ++i) {
        file << i * p.grid.d_rho();
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
    for(size_t i = 0; i < p.grid.n_rho(); ++i) {
        file << i * p.grid.d_rho();
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

// Adaptive integrator (RK4 with step-doubling error estimate)
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p,
                             const StepperConfig& cfg, const std::string& filename,
                             int n_snapshots) {
    if(dt_init >= 0) {
        std::cerr << "[ERROR] dt_init must be negative" << std::endl;
        return;
    }

    std::cout << "Solving flow equation with adaptive time step...\n";

    const RHSfunc rhs = [&p](const std::vector<double>& state, double t) {
        return RHS(state, std::exp(t), p);
    };

    std::vector<double> snap_targets(n_snapshots);
    for(int s = 0; s < n_snapshots; ++s) {
        double fraction = static_cast<double>(s) / (n_snapshots - 1);
        snap_targets[s] = p.t_start + fraction * (p.t_end - p.t_start);
    }

    std::vector<std::pair<double, std::vector<double>>> snapshot_pairs;
    std::vector<std::pair<double, double>> step_history;

    integrate_adaptive(V_init, p.t_start, p.t_end, dt_init, rhs, cfg, snap_targets, &snapshot_pairs,
                       &step_history);

    std::vector<std::vector<double>> snapshots;
    std::vector<std::vector<double>> rhs_snapshots;
    std::vector<double> k_values;
    snapshots.reserve(snapshot_pairs.size());
    rhs_snapshots.reserve(snapshot_pairs.size());
    k_values.reserve(snapshot_pairs.size());

    for(const auto& [t, V] : snapshot_pairs) {
        snapshots.push_back(V);
        rhs_snapshots.push_back(RHS(V, std::exp(t), p));
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

    save_all(snapshots, rhs_snapshots, k_values, p, filename);
    if(!filename.empty()) {
        save_dt_hist(
            dt_values, dt_k_values,
            "results/phi4/dt_values.txt"); // TODO use the first part of string of filename of
                                           // integrate_flow_adaptive to have a dynamical directory
    }
}

} // namespace QM
