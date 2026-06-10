#include "grid.hpp"

// dervatives -------------------------------------


double Grid::d1(const std::vector<double>& u, size_t i) {
    double d1;
    if (m_rho_uniform) {
        // interior: 4th order finite difference
        if (i >= 2 && i <= m_n_rho - 3)  {
            return (-u[i+2] + 8.0*u[i+1] - 8.0*u[i-1] + u[i-2])/(12.0*m_d_rho);
        }
        // one point from boundary: 2nd order finite difference
        if (i == 1 || i == m_n_rho - 2) {
            return (u[i+1] - u[i-1])/(2*m_d_rho);
        }
        // boundary: 1st order one sided
        if (i == 0) { return (u[1] - u[0])/m_d_rho; }
        if (i == m_n_rho - 1) { return (u[m_n_rho-1] - u[m_n_rho-2])/m_d_rho; }
    }
    // for log grid we have a non uniform spacing in rho, but a uniform spacing in t = log(rho), du / drho = du / dt * dt / drho = du / dt 1 / rho
    else {
        double dt = std::log(10.0) / m_n10_rho;
        double rho = m_rho_vals[i];
        double du_dt;
        // interior: 4th order finite difference
        if (i >= 2 && i <= m_n_rho - 3)  {
            du_dt =  (-u[i+2] + 8.0*u[i+1] - 8.0*u[i-1] + u[i-2])/(12.0*dt);
        }
        // one point from boundary: 2nd order finite difference
        if (i == 1 || i == m_n_rho - 2) {
            du_dt = (u[i+1] - u[i-1])/(2*dt);
        }
        // boundary: 1st order one sided
        if (i == 0) { du_dt =  (u[1] - u[0])/dt; }
        if (i == m_n_rho - 1) { du_dt = (u[m_n_rho-1] - u[m_n_rho-2])/dt; }

        return du_dt / rho;
    }
    return d1;
}

double Grid::d2(const std::vector<double>& u, size_t i) {
    double d2;
    if (m_rho_uniform) {
        // interions: 4th order finite difference
        if (i >= 2 && i <= m_n_rho - 3) {
            return (-u[i+2] + 16.0*u[i+1] - 30.0*u[i] + 16.0*u[i-1] - u[i-2])/(12.0*m_d_rho*m_d_rho);
        }
        // one point from boundary: 2nd order finite difference
        if (i == 1 || i == m_n_rho - 2) {
            return (u[i+1] - 2.0*u[i] + u[i-1])/(m_d_rho*m_d_rho);
        }
        // boundary: 1st order one sided
        if (i == 0) { return (u[2] - 2.0*u[1] + u[0])/(m_d_rho*m_d_rho); }
        if (i == m_n_rho - 1) { return (u[m_n_rho-1] - 2.0*u[m_n_rho-2] + u[m_n_rho-3])/(m_d_rho*m_d_rho); }
    }
    else {
        double dt = std::log(10.0) / m_n10_rho;
        double rho = m_rho_vals[i];
        double du_dt;
        // interions: 4th order finite difference
        if (i >= 2 && i <= m_n_rho - 3) {
            du_dt = (-u[i+2] + 16.0*u[i+1] - 30.0*u[i] + 16.0*u[i-1] - u[i-2])/(12.0*dt*dt);
        }
        // one point from boundary: 2nd order finite difference
        if (i == 1 || i == m_n_rho - 2) {
            du_dt = (u[i+1] - 2.0*u[i] + u[i-1])/(dt*dt);
        }
        // boundary: 1st order one sided
        if (i == 0) { du_dt =  (u[2] - 2.0*u[1] + u[0])/(dt*dt); }
        if (i == m_n_rho - 1) { du_dt = (u[m_n_rho-1] - 2.0*u[m_n_rho-2] + u[m_n_rho-3])/(dt*dt); }

        return du_dt / rho;

    }
    return d2;
}



std::vector<double> Grid::d1(const std::vector<double>& u) {
    std::vector<double> result(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
        result[i] = d1(u,i);
    }
    return result;
}

std::vector<double> Grid::d2(const std::vector<double>& u) {
    std::vector<double> result(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
        result[i] = d2(u,i);
    }
    return result;
}

// grid updater -------------------------------------

void Grid::update_grid() {
    if (m_rho_uniform) {
        m_rho_vals = linspace(m_rho_min, m_rho_max, m_n_rho);
    }
    else {
        m_rho_vals = logspace(m_rho_min_pot, m_rho_max_pot, m_n10_rho);
    }

}

// setters --------------------

void Grid::set_rho_vals(std::vector<double> rho_vals) {
    m_rho_vals = rho_vals;
}

void Grid::set_rho_min(double rho_min) {
    m_rho_min = rho_min;
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