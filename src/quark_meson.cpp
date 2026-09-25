#include "quark_meson.hpp"

namespace QM {

// classical potential -----------------
std::vector<double> V_classical(const Params& p) {
    std::vector<double> V(p.grid.n_rho());
    double rho;

    for(size_t i = 0; i < p.grid.n_rho(); ++i) {
        rho = i * p.grid.d_rho();
        V[i] = p.m2 * rho + p.lambda / 6.0 * rho * rho;
    }

    return V;
}

double V_min_classical(const Params& p) {
    return -3.0 * p.m2 / p.lambda;
}

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

// Compute observables of the QM model
Observables compute_observables(std::vector<double>& V, const Params& p) {
    Observables obs;
    const size_t N = V.size();
    if(N < 3) {
        return obs;
    }

    const double drho = p.grid.d_rho();

    // 1. Find minimum index
    auto min_it = std::min_element(V.begin(), V.end());
    size_t idx = std::distance(V.begin(), min_it);

    if(idx == 0 || idx >= N - 1) {
        obs.rho0 = 0.0;

        double V_prime_0 = (V[1] - V[0]) / drho;
        obs.m2_pi = V_prime_0;
        obs.m2_sigma = 0.;
        return obs;
    }

    // 2. Parabolic Interpolation around local minimum
    double y1 = V[idx - 1];
    double y2 = V[idx];
    double y3 = V[idx + 1];

    double denom = y3 - 2.0 * y2 + y1;
    double delta = (denom > 1e-12) ? -0.5 * (y3 - y1) / denom : 0.;
    obs.rho0 = (static_cast<double>(idx) + delta) * drho;

    // 3. Compute derivatives
    double V_prime_idx = (y3 - y1) / (2. * drho);
    double V_double_prime = denom / (drho * drho);
    double V_prime_at_rho0 = V_prime_idx + delta * drho * V_double_prime;

    // 4. Compute masses
    obs.m2_pi = V_prime_at_rho0;
    obs.m2_sigma = V_prime_at_rho0 + 2.0 * obs.rho0 * V_double_prime;

    return obs;
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

void sweep_params(const std::vector<double>& m2, const std::vector<double>& lambda,
                  const std::vector<double>& h, const Params& p, const StepperConfig& cfg,
                  const std::string& filename) {
    std::cout << "============== Quark Meson Model ==============\n";
    std::cout << "Sweep UV params..." << "\n";
    const int n_m2 = static_cast<int>(m2.size());
    const int n_lambda = static_cast<int>(lambda.size());
    const int n_h = static_cast<int>(h.size());
    int number_of_sweeps = n_m2 * n_lambda * n_h;

    std::cout << "Number of sweeps = " << number_of_sweeps << "\n";
    std::cout << "Output file: " << filename << std::endl;

    std::ofstream file(filename);
    if(!file) {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
        return;
    }

    // metadata
    file << "# Wetterich QM LPA flow, d=3, N=";
    file << p.N << ", Nc = " << p.Nc << "\n";
    file << "# rho_max = " << p.grid.rho_max() << ", n_rho = " << p.grid.n_rho() << "\n";
    file << "# ------------------------------\n";
    file << "# m2, lambda, h, rho0, m2_sigma, m2_pi \n";

    int sweeps_done = 0;
    double total_comp_time = 0.;

    std::vector<double> rho0_vals(number_of_sweeps);
    std::vector<double> m2_sigma_vals(number_of_sweeps);
    std::vector<double> m2_pi_vals(number_of_sweeps);

#pragma omp parallel for collapse(3) schedule(dynamic)
    for(int i_m2 = 0; i_m2 < n_m2; ++i_m2) {
        for(int i_lambda = 0; i_lambda < n_lambda; ++i_lambda) {
            for(int i_h = 0; i_h < n_h; ++i_h) {
                const int idx = (i_m2 * n_lambda + i_lambda) * n_h + i_h;

                Params p_sweep = p;
                p_sweep.m2 = m2[i_m2];
                p_sweep.lambda = lambda[i_lambda];
                p_sweep.h = h[i_h];

                double dt_init = -0.0001;
                const std::vector<double>& V_init = V_classical(p_sweep);

                const RHSfunc rhs = [&p_sweep](const std::vector<double>& state, double t) {
                    return RHS(state, std::exp(t), p_sweep);
                };
                // std::cout << "Solving flow equation...\n";
                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                std::vector<double> V = integrate_adaptive(V_init, p_sweep.t_start, p_sweep.t_end,
                                                           dt_init, rhs, cfg, {}, nullptr, nullptr);
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                double duration = std::chrono::duration<double>(end - begin).count();

                // Compute observables
                Observables obs = compute_observables(V, p_sweep);
                rho0_vals[idx] = obs.rho0;
                m2_sigma_vals[idx] = obs.m2_sigma;
                m2_pi_vals[idx] = obs.m2_pi;

                int done = ++sweeps_done;
#pragma omp critical
                {
                    total_comp_time += duration;
                    double avg_time = total_comp_time / done;
                    int remaining_sweeps = number_of_sweeps - done;
                    double remaining_time = avg_time * remaining_sweeps;
                    std::cout << "---------- [" << done << "/" << number_of_sweeps
                              << "] ----------\n";
                    std::cout << "m2 = " << p_sweep.m2 << ", lambda = " << p_sweep.lambda
                              << ", h = " << p_sweep.h << "\n";
                    std::cout << "ρ0 = " << obs.rho0 << "\n";
                    std::cout << "Comptime          = " << duration << "[s]\n";
                    std::cout << "avg. comptime     = " << avg_time << "[s]\n";
                    std::cout << "Remaining time    = " << remaining_time / 60 << "[min]\n";
                };
            }
        }
    }

    // save data
    int i = 0;
    for(double m2_val : m2) {
        for(double lambda_val : lambda) {
            for(double h_val : h) {
                file << m2_val << ", " << lambda_val << ", " << h_val << ", " << rho0_vals[i]
                     << ", " << m2_sigma_vals[i] << ", " << m2_pi_vals[i] << "\n";
                i += 1;
            }
        }
    }

    file.close();
}

} // namespace QM
