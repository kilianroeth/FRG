#include <iostream>
#include <chrono>

#include "phi4_lpa.hpp"
#include "utils.hpp"

int main() {

    Grid grid;
    grid.set_rho_vals(linspace(0.,0.5,500));

    Params p;
    p.grid = grid;
    p.m2 = -0.02;
    p.d = 3;
    p.N = 4;
    p.lambda = 1.0;
    p.t_start = 0.0;
    p.t_end = -10.0;

    StepperConfig cfg;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-10;
    cfg.show_progress = true;

    
    std::vector<double> V_init = V_classical(p);
    std::cout << "Classical minimum: V_min = " << V_min_classical(p) << std::endl;
    std::vector<double> dV_vals(V_init.size()), ddV_vals(V_init.size());
    for (size_t i = 0; i < V_init.size(); ++i) {
        dV_vals[i] = grid.d1(V_init,i);
        ddV_vals[i] = grid.d2(V_init,i);
    }
    std::vector<double> RHS_vals = RHS(V_init,1,p);
    save_V(V_init, "results/V_classical.txt",p);
    save_V(dV_vals, "results/V_classical_prime.txt",p);
    save_V(ddV_vals, "results/V_classical_doubleprime.txt",p);
    save_V(RHS_vals, "results/RHS.txt",p);

    integrate_flow_adaptive(V_init, -0.0001, p, cfg, "results/flow_adaptive.csv", 100);

    return 0;
}
