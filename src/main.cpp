#include <chrono>
#include <iostream>

#include "phi4_lpa.hpp"
#include "quark_meson.hpp"
#include "utils.hpp"

int main() {

    ///////////////////////////////
    // QM model
    ///////////////////////////////

    Grid grid;
    grid.set_rho_vals(linspace(0., 0.75, 500));

    QM::Params p;
    p.m2 = -0.08;
    p.d = 3;
    p.N = 4;
    p.Nc = 3;
    p.lambda = 1.0;
    p.h = 5;
    p.t_start = 0.0;
    p.t_end = -10.0;
    p.warning_level = 0;

    StepperConfig cfg;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-10;
    cfg.show_progress = false;

    std::vector<double> V_init = QM::V_classical(p, grid);
    std::cout << "Classical minimum: V_min = " << QM::V_min_classical(p) << std::endl;

    std::vector<double> dV_vals(V_init.size()), ddV_vals(V_init.size());
    for(size_t i = 0; i < V_init.size(); ++i) {
        dV_vals[i] = grid.d1(V_init, i);
        ddV_vals[i] = grid.d2(V_init, i);
    }

    // Allocate bufer and call QM::RHS in-place
    std::vector<double> RHS_vals(grid.n_rho());
    QM::RHS(V_init, 1, RHS_vals, p, grid);
    QM::save_V(V_init, "results/QM/V_classical.txt", grid);
    QM::save_V(dV_vals, "results/QM/V_classical_prime.txt", grid);
    QM::save_V(ddV_vals, "results/QM/V_classical_doubleprime.txt", grid);
    QM::save_V(RHS_vals, "results/QM/RHS.txt", grid);

    // QM::integrate_flow_adaptive(V_init, -0.0001, p, cfg, "results/QM/flow_adaptive.csv", 100);
    QM::sweep_params(linspace(-0.1, 0.75, 5), linspace(1.0, 1.0, 1), linspace(0., 7.5, 5), p, grid,
                     cfg, "results/QM/UV_params.csv");

    return 0;
}
