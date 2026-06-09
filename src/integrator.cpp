#include "integrator.hpp"

// error computation --------------------------------

double compute_error(const std::vector<double>& u1, const std::vector<double>& u2, double absolute_tolerance, double relative_tolerance) {
    if (u1.size() != u2.size()) {
        std::cerr << "Arrays don't have the samve size. V1.size() = " << u1.size() << ", V2.size() = " << u2.size() << "\n";
        throw std::runtime_error("V1 and V2 size mismach");
    }
    
    // sum of squared errors of each array entry
    double squared_errors = 0.0;

    // iterate over all entries
    for (size_t i = 0; i < u1.size(); ++i) {
        // scale = a_tol + r_rol * max(V1,V2)
        double scale = absolute_tolerance + relative_tolerance * std::max(std::abs(u1[i]), std::abs(u2[i]));
        double error = std::abs(u1[i] - u2[i]) / scale;
        // if error < 1, we are within the tolerances
        squared_errors += error * error;
    }

    // return RMS of sum of array entries
    return sqrt(squared_errors / u1.size());
}

// time stepper -------------------------------------

std::vector<double> step_euler(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs) {
    auto f = rhs(state, t);
    size_t N = state.size();
    std::vector<double> next(N);
    for (size_t i = 0; i < N; ++i) {
        next[i] = state[i] + dt * f[i];
    }
    return next;
}


std::vector<double> step_rk4(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs) {
    size_t N = state.size();
    std::vector<double> tmp(N), next(N);

    auto k1 = rhs(state, t);

    for (size_t i = 0; i < N; ++i) { tmp[i] = state[i] + dt/2. * k1[i]; }
    auto k2 = rhs(tmp, t + dt/2.);

    for (size_t i = 0; i < N; ++i) { tmp[i] = state[i] + dt/2. * k2[i]; }
    auto k3 = rhs(tmp, t + dt/2.);

    for (size_t i = 0; i < N; ++i) { tmp[i] = state[i] + dt    * k3[i]; }
    auto k4 = rhs(tmp, t + dt);

    for (size_t i = 0; i < N; ++i) {
        next[i] = 1./6. * (k1[i] + 2.*k2[2] + 2.*k3[i] + k4[i]);
    }

    return next;
}

// adaptive flow integrator -----------------------------

std::vector<double> integrate_adaptive(
    const std::vector<double>& state_init,
    double t_start, double t_end,
    double dt_init,
    const RHSfunc& rhs,
    const StepperConfig& cfg,
    const std::vector<double>& snap_targets,
    std::vector<std::pair<double, std::vector<double>>>* snapshots) {
        double sign = (t_end > t_start) ? 1.0 : -1.0;
        double dt   = dt_init;
        double t    = t_start;
        auto V      = state_init;

        size_t next_snap = 0;
        int n_accepted = 0, n_rejected = 0;

        while (sign * t < sign * t_end) {
            // avoid overshooting
            if (sign * (t + dt) > sign * t_end) {
                dt = t_end - t;
            }

            // full step and two half steps to compare dt and dt/2 deviation
            auto V_full = step_rk4(V, t, dt, rhs);
            auto V_mid  = step_rk4(V, t, dt/2., rhs);
            auto V_half = step_rk4(V_mid, t + dt/2., dt/2., rhs);

            double error = compute_error(V_full, V_half, cfg.abs_tol, cfg.rel_tol);

            if (error < 1.0) {
                V = V_half;
                t += dt;
                ++n_accepted;

                // collect snaphost
                if (snapshots) {
                    while (next_snap < snap_targets.size() && snap_targets[next_snap]) {
                        snapshots->push_back({t,V});
                        ++next_snap;
                    }
                }

                if(cfg.show_progress) {
                    double progress = std::abs(t - t_start) / std::abs(t_end - t_start);
                    progressBar(static_cast<size_t>(progress * 1000), 1000);
                }
            }
            else {
                ++n_rejected;
            }

            // rescale dt - rk4 is a fourth order methods, so exponents is 1/4
            double factor = (error > 0)
                ? std::clamp(cfg.safety * std::pow(1.0/error, 0.25), cfg.factor_min, cfg.factor_max)
                :cfg.factor_max;
            dt *= factor;

            if (std::abs(dt) < cfg.dt_min) {
                std::cerr << "[ERROR] dt below minimum, aboritng\n";
                break;
            }
        }

        if (cfg.show_progress) {
            std::cout << "\nAccepted: " << n_accepted << ", Rejected: " << n_rejected << "\n";
        }

        return V;
    }
