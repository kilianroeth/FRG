#include <iostream>

#include "phi4_lpa.hpp"

int main() {
    
    std::vector<double> V = V_classical();
    std::vector<double> dV_vals(V.size()), ddV_vals(V.size());
    for (size_t i = 0; i < V.size(); ++i) {
        dV_vals[i] = dV(V,i);
        ddV_vals[i] = ddV(V,i);
    }
    std::vector<double> RHS_vals = RHS(V,1);
    save_V(V, "results/V_classical.txt");
    save_V(dV_vals, "results/V_classical_prime.txt");
    save_V(ddV_vals, "results/V_classical_doubleprime.txt");
    save_V(RHS_vals, "results/RHS.txt");
    integrate_flow(V, -0.001);

    return 0;
}
