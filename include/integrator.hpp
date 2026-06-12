#pragma once

#include <functional>
#include <vector>
#include <cmath>

#include "utils.hpp"

// General flow equation structurs: ∂_t f(state,t) = RHS(state,t)
// with t = ln(k/Λ) the RG time

using RHSfunc = std::function<std::vector<double>(const std::vector<double>&, double)>;

struct StepperConfig {
    double abs_tol      = 1e-8;     // absolute tolerance for adaptive time step error acceptance
    double rel_tol      = 1e-8;     // relative tolerance for adaptive time step error acceptance
    double dt_min       = 1e-15;    // minimal allowed time step
    double safety       = 0.9;      // safety factor for updating adaptive time step dt -> dt * 0.9 * scaling factor
    double factor_max   = 5.0;      // upper bound for scaling factor
    double factor_min   = 0.1;      // lower bound for scaling factor
    bool show_progress  = true;     // print progress
};


// difference between two arrays u1 and u2, rms of error of entries, error > 1 if difference > tolerances, error < 1 if difference < tolerances
double compute_error(const std::vector<double>& u1, const std::vector<double>& u2, double absolute_tolerance, double relative_tolerance);

// time stepper
std::vector<double> step_euler(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs);
std::vector<double> step_rk4(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs);

// adaptive flow integrator
std::vector<double> integrate_adaptive(
    const std::vector<double>& state_init,
    double t_start, double t_end, 
    double dt_init,
    const RHSfunc& rhs,
    const StepperConfig& cfg,
    const std::vector<double>& snap_targets = {},
    std::vector<std::pair<double, std::vector<double>>>* snapshots = nullptr,
    std::vector<std::pair<double, double>>* step_history = nullptr);