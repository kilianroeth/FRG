#include <iostream>

#include "phi4_lpa.cpp"

int main() {
    

    std::vector<double> V = V_classical();

    std::ofstream out("vector.txt");
    for (size_t i = 0; i < V.size(); ++i) {
        out << i/50. << ", " << V[i] << "\n";
    }

    std::cout << "Hello FRG" << std::endl;

    return 0;
}
