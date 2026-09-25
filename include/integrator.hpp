#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "utils.hpp"

// General flow equation structurs: ∂_t f(state,t) = RHS(state,t)
// with t = ln(k/Λ) the RG time

using RHSfunc =
    std::function<void(const std::vector<double>& state, double t, std::vector<double>& out)>;

struct StepperConfig {
    double abs_tol = 1e-8; // absolute tolerance for adaptive time step error acceptance
    double rel_tol = 1e-8; // relative tolerance for adaptive time step error acceptance
    double dt_min = 1e-15; // minimal allowed time step
    double safety =
        0.9; // safety factor for updating adaptive time step dt -> dt * 0.9 * scaling factor
    double factor_max = 5.0;   // upper bound for scaling factor
    double factor_min = 0.1;   // lower bound for scaling factor
    bool show_progress = true; // print progress
};

// Scratch space for RK4 intermediate steps
struct RK4Scratch {
    std::vector<double> tmp;
    std::vector<double> k1;
    std::vector<double> k2;
    std::vector<double> k3;
    std::vector<double> k4;

    explicit RK4Scratch(size_t size) : tmp(size), k1(size), k2(size), k3(size), k4(size) {}
};

// Scratch space for adaptive step ocmparison
struct AdaptiveScratch {
    std::vector<double> V_full;
    std::vector<double> V_mid;
    std::vector<double> V_half;
    RK4Scratch rk4_scratch;

    explicit AdaptiveScratch(size_t size) : V_full(size), V_mid(size), rk4_scratch(size) {}
};

// difference between two arrays u1 and u2, rms of error of entries, error > 1 if difference >
// tolerances, error < 1 if difference < tolerances
double compute_error(const std::vector<double>& u1, const std::vector<double>& u2,
                     double absolute_tolerance, double relative_tolerance);

// time stepper
void step_euler(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs,
                std::vector<double>& out, std::vector<double>& scratch_k);
void step_rk4(const std::vector<double>& state, double t, double dt, const RHSfunc& rhs,
              std::vector<double>& out, RK4Scratch& scratch);

// adaptive flow integrator
std::vector<double>
integrate_adaptive(const std::vector<double>& state_init, double t_start, double t_end,
                   double dt_init, const RHSfunc& rhs, const StepperConfig& cfg,
                   const std::vector<double>& snap_targets = {},
                   std::vector<std::pair<double, std::vector<double>>>* snapshots = nullptr,
                   std::vector<std::pair<double, double>>* step_history = nullptr);