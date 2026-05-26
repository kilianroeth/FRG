#include <iostream>

#include "phi4_lpa.cpp"

int main() {
    

    std::vector<double> V = V_classical();
    save_V(V, "results/V_classical.txt");
    integrate_flow(V, -0.001);


    return 0;
}
