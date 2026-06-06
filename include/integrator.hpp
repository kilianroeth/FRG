#pragma once

#include <functional>
#include <vector>

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