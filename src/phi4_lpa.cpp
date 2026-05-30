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
    return -3.0*p.m2/p.lambda;
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
        if (std::abs(denom) < 1e-12) {
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

// Integrate RG time step --------------

// Simple forward euler
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

// Classical RK4 step
std::vector<double> step_rk4(const std::vector<double>& V, double k, double dt, const Params& p) {
    std::vector<double> V_next(p.n_rho);

    // first step
    std::vector<double> K1 = RHS(V, k, p);

    std::vector<double> K2, K3, K4;
    std::vector<double> V_temp = V;
    double k_half = k*exp(dt/2.0);
    double k_next = k*exp(dt);

    // second step at dt/2
    for (size_t i = 0; i < p.n_rho; ++i) {
        V_temp[i] = V[i] + dt/2.0 * K1[i];
    }
    K2 = RHS(V_temp, k_half, p);

    // third step at dt/2
    for (size_t i = 0; i < p.n_rho; ++i) {
        V_temp[i] = V[i] + dt/2.0 * K2[i];
    }
    K3 = RHS(V_temp, k_half, p);

    // fourth step at dt
    for (size_t i = 0; i < p.n_rho; ++i) {
        V_temp[i] = V[i] + dt * K3[i];
    }
    K4 = RHS(V_temp, k_next, p);

    for (size_t i = 0; i < p.n_rho; ++i) {
        V_next[i] = V[i] + dt/6.0 * (K1[i] + 2.0 * K2[i] + 2.0 * K3[i] + K4[i]);
    }

    return V_next;
}

// save current potential --------------

void save_V(const std::vector<double>& V, const std::string& filename, const Params& p) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; }

    file << "ρ = 1/2 φ², V(ρ)\n";
    for (size_t i = 0; i < V.size(); ++i) {
        file << p.rho_at(i) << ", " << V[i] <<"\n";
    }
    file << std::endl;

    file.close();
}

void save_all(const std::vector<std::vector<double>>& snapshots, const std::vector<std::vector<double>>& rhs_snapshots, const std::vector<double>& k_values, const Params& p, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open " << filename << "\n"; return; }

    // metadata
    file << "# Wetterich LPA flow, phi^4, d=3\n";
    file << "# m2 = " << p.m2 << ", lambda = " << p.lambda << "\n";
    file << "# rho_max = " << p.rho_max << ", n_rho = " << p.n_rho << "\n";

    // V block
    file << "# block: V\n";
    file << "rho";
    for (double k : k_values)
        file << ", k=" << std::fixed << std::setprecision(6) << k;
    file << "\n";
    file << std::scientific << std::setprecision(10);
    for (size_t i = 0; i < p.n_rho; ++i) {
        file << i * p.drho();
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
    for (size_t i = 0; i < p.n_rho; ++i) {
        file << i * p.drho();
        for (const auto& R : rhs_snapshots)
            file << ", " << R[i];
        file << "\n";
    }

    std::cout << "Saved " << snapshots.size() << " snapshots -> " << filename << "\n";
}

void save_dt_hist(const std::vector<double>& dt_values, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open" << filename << "\n"; return; }
    
    // metadata
    file << "#time steps dt for adaptive RK4 time stepper\n";
    // data
    for (double dt : dt_values) {
        file << dt << "\n";
    }
}

void save_m2_flow(const std::vector<double>& m2_values, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) { std::cerr << "[ERROR] Cannot open" << filename << "\n"; return; }

    // metadata
    file << "# flow of m2\n";
    // data
    for (size_t i = 0; i < m2_values.size(); ++i) {
        file << m2_values[i] << "\n";
    }
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
            V = step(V, k, dt_t, p);
        }
    }
    save_all(snapshots, rhs_snapshots, k_values, p, filename);
}


// Adaptive integrator using RK4 + step-doubling
void integrate_flow_adaptive(const std::vector<double>& V_init, double dt_init, const Params& p, const std::string& filename, int n_snapshots, double absolute_tolerance, double relative_tolerance) {
    if (dt_init >= 0) {
        std::cerr << "[ERROR] dt_init must be negative" << std::endl;
        return;
    }

    std::vector<std::vector<double>> snapshots, rhs_snapshots, m2_snapshots;
    // target times for snapshots
    std::vector<double> snap_targets(n_snapshots);
    for (int s = 0; s < n_snapshots; ++s) {
        double fraction = static_cast<double>(s) / (n_snapshots - 1);
        snap_targets[s] = p.t_start + fraction * (p.t_end - p.t_start);
    }
    int next_snap = 0;
    std::vector<double> k_values;
    std::vector<double> dt_values;
    
    std::vector<double> V = V_init;

    double t = p.t_start;
    double dt = dt_init;
    double sign = (dt > 0.0) ? 1.0 : -1.0;

    snapshots.push_back(V);
    rhs_snapshots.push_back(RHS(V, exp(t), p));
    k_values.push_back(exp(t));
    ++next_snap;

    // acceptance statistics
    int n_accepted = 0, n_rejected = 0;

    while (sign * t < sign * p.t_end) {
        // avoid overshooting
        if (sign * (t + dt) > sign * p.t_end) {
            dt = p.t_end - t;
        }

        // full step
        double k = exp(t);
        std::vector<double> V_full = step_rk4(V, k, dt, p);

        // two half steps
        double k_half = k * exp(dt/2.0);
        std::vector<double> V_mid = step_rk4(V, k, dt/2.0, p);
        std::vector<double> V_half = step_rk4(V_mid, k_half, dt/2.0, p);

        // estimate error
        double error = compute_error(V_full, V_half, absolute_tolerance, relative_tolerance);

        // check if error is small enough
        if (error < 1.0) {
            V = V_half;
            t += dt;
            dt_values.push_back(dt);
            ++n_accepted;

            // progress bar
            double progress = std::abs(t - p.t_start) / std::abs(p.t_end - p.t_start);
            progressBar(static_cast<size_t>(progress * 1000), 1000);

            while (next_snap < n_snapshots && sign * t >= sign * snap_targets[next_snap]) {
                snapshots.push_back(V);
                rhs_snapshots.push_back(RHS(V, exp(t), p));
                k_values.push_back(exp(t));
                ++next_snap;
            }
        }
        else {
            ++n_rejected;
        }

        // update time step dt, safety factor 0.9
        double factor = (error > 0) ? std::clamp(0.9 * std::pow(1.0/error, 1.0/4.0), 0.1, 5.0) : 5.0;
        dt *= factor;

        if (std::abs(dt) < 1e-15) {
            std::cerr << "[ERROR] dt too small, aboritng\n";
            break;
        }
    }

    std::cout << "Accepted: " << n_accepted << ", Rejected: " << n_rejected << "\n";
    save_all(snapshots, rhs_snapshots, k_values, p, filename);
    save_dt_hist(dt_values, "results/dt_values.txt");
}

double compute_error(const std::vector<double>& V1, const std::vector<double>& V2, double absolute_tolerance, double relative_tolerance) {
    if (V1.size() != V2.size()) {
        std::cerr << "Arrays don't have the samve size. V1.size() = " << V1.size() << ", V2.size() = " << V2.size() << "\n";
        throw std::runtime_error("V1 and V2 size mismach");
    }
    
    // sum of squared errors of each array entry
    double squared_errors = 0.0;

    // iterate over all entries
    for (size_t i = 0; i < V1.size(); ++i) {
        // scale = a_tol + r_rol * max(V1,V2)
        double scale = absolute_tolerance + relative_tolerance * std::max(std::abs(V1[i]), std::abs(V2[i]));
        double error = std::abs(V1[i] - V2[i]) / scale;
        // if error < 1, we are within the tolerances
        squared_errors += error * error;
    }

    // return RMS of sum of array entries
    return sqrt(squared_errors / V1.size());
}

size_t find_min(const std::vector<double>& V) {
    size_t i_min = 0;
    for (size_t i = 0; i < V.size(); ++i) {
        if(V[i] < V[i_min]) {
            i_min = i;
        }
    }
    return i_min;
}

double m2(const std::vector<double>& V, size_t i, const Params& p) {
    double m2;
    if (i == 0) {
        return dV(V, i, p);
    }
    else {
        return 2.0 * p.rho_at(i) * ddV(V, i, p);
    }
    return m2;
}

