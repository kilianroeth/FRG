#include "grid.hpp"

double Grid::d1(const std::vector<double>& u, size_t i) {

}

double Grid::d2(const std::vector<double>& u, size_t i) {

}



std::vector<double> Grid::d1(const std::vector<double>& u) {

}

std::vector<double> Grid::d2(const std::vector<double>& u) {

}

// setters
// --------------------

void Grid::set_rho_vals(std::vector<double> rho_vals) {
    m_rho_vals = rho_vals;
}

void Grid::set_rho_min(double rho_min) {
    m_rho_max = rho_min;
}

void Grid::set_rho_max(double rho_max) {
    m_rho_max = rho_max;
}

void Grid::set_rho_uniform(bool rho_uniform) {
    m_rho_uniform = rho_uniform;
}


void Grid::set_d_rho(double d_rho) {
    m_d_rho = d_rho;
}

void Grid::set_n_rho(double n_rho) {
    m_n_rho = n_rho;
}

void Grid::set_n10_rho(double n10_rho) {
    m_n10_rho = n10_rho;
}
