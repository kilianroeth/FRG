#include "grid.hpp"

double Grid::d1(const std::vector<double>& u, size_t i) {
    double d1;
    return d1;
}

double Grid::d2(const std::vector<double>& u, size_t i) {
    double d2;
    return d2;
}



std::vector<double> Grid::d1(const std::vector<double>& u) {
    std::vector<double> d1;
    return d1;
}

std::vector<double> Grid::d2(const std::vector<double>& u) {
    std::vector<double> d2;
    return d2;
}

void Grid::update_grid() {
    if (m_rho_uniform) {
        m_rho_vals = linspace(m_rho_min, m_rho_max, m_n_rho);
    }
    else {
        m_rho_vals = logspace(m_rho_min_pot, m_rho_max_pot, m_n10_rho);
    }

}

// setters
// --------------------

void Grid::set_rho_vals(std::vector<double> rho_vals) {
    m_rho_vals = rho_vals;
}

void Grid::set_rho_min(double rho_min) {
    m_rho_max = rho_min;
    update_grid();
}

void Grid::set_rho_max(double rho_max) {
    m_rho_max = rho_max;
    update_grid();
}

void Grid::set_rho_uniform(bool rho_uniform) {
    m_rho_uniform = rho_uniform;
}


void Grid::set_d_rho(double d_rho) {
    m_d_rho = d_rho;
    m_n_rho = (m_rho_max - m_rho_min) / m_d_rho;
    update_grid();
}

void Grid::set_n_rho(double n_rho) {
    m_n_rho = n_rho;
    update_grid();
}

void Grid::set_n10_rho(double n10_rho) {
    m_n10_rho = n10_rho;
}

void Grid::set_rho_min_pot(double rho_min_pot) {
    m_rho_min_pot = rho_min_pot;
    update_grid();
}

void Grid::set_rho_max_pot(double rho_max_pot) {
    m_rho_max_pot = rho_max_pot;
    update_grid();
}