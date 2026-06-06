#include <iostream>
#include <chrono>

#include "phi4_lpa.hpp"
#include "utils.hpp"

int main() {

    Params p;
    p.n_rho = 10000;
    p.rho_max = 0.5;
    p.m2 = -0.02;
    p.lambda = 1.0;
    p.t_start = 0.0;
    p.t_end = -10.0;

    
    std::vector<double> V = V_classical(p);
    std::cout << "Classical minimum: V_min = " << V_min_classical(p) << std::endl;
    std::vector<double> dV_vals(V.size()), ddV_vals(V.size());
    for (size_t i = 0; i < V.size(); ++i) {
        dV_vals[i] = dV(V,i,p);
        ddV_vals[i] = ddV(V,i,p);
    }
    std::vector<double> RHS_vals = RHS(V,1,p);
    save_V(V, "results/V_classical.txt",p);
    save_V(dV_vals, "results/V_classical_prime.txt",p);
    save_V(ddV_vals, "results/V_classical_doubleprime.txt",p);
    save_V(RHS_vals, "results/RHS.txt",p);

    std::vector<double> m2_values = linspace(0.05,-0.05,20);
    // sweep_m2(p, m2_values, "results/m2_sweep.txt", -0.0001, 1e-10, 1e-10);

    integrate_flow_adaptive(V, -0.0001, p, "results/flow_adaptive.csv", 100, 1e-12, 1e-12);

    return 0;
}
