#include "quark_meson.hpp"

namespace QM {

// Compute RHS ------------------------

std::vector<double> RHS(const std::vector<double>& V, double k, const Params& p) {
    std::vector<double> RHS_vals(p.grid.n_rho());
    return RHS_vals;
}

} // namespace QM
