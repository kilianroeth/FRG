#pragma once

#include <vector>
#include <cmath>

#include "utils.hpp"

class Grid {
    private:
        std::vector<double> m_rho_vals;     // stores the rho grid
        double  m_rho_min = 0.0;            // smallest rho value
        double  m_rho_max = 0.0;            // largest rho value
        bool    m_rho_uniform  = true;      // linear grid / log grid
        double  m_d_rho = 1.0;              // rho grid spacing for linear grid
        size_t  m_n_rho = 0;                // number of rho points for the linear grid
        double  m_rho_min_pot = 0.0;        // log grid: smallest 10^pot
        double  m_rho_max_pot = 0.0;        // log grid largest 10^pot
        double  m_n10_rho = 0.0;            // number of rho points per decade for the log grid

        void update_grid();

    public:
        double d1(const std::vector<double>& u, size_t i) const;
        double d2(const std::vector<double>& u , size_t i) const;

        std::vector<double> d1(const std::vector<double>& u) const;
        std::vector<double> d2(const std::vector<double>& u) const;

        // setters ------------------------------
        void set_rho_vals(std::vector<double> rho_vals);
        void set_rho_min(double rho_min);
        void set_rho_max(double rho_max);
        void set_rho_uniform(bool rho_uniform);
        void set_d_rho(double d_rho);
        void set_n_rho(size_t n_rho);
        void set_n10_rho(double n10_rho);
        void set_rho_min_pot(double rho_min_pot);
        void set_rho_max_pot(double rho_max_pot);
        
        // getters ------------------------------
        std::vector<double> rho_vals() const { return m_rho_vals; }
        double rho_min() const { return m_rho_min; }
        double rho_max() const { return m_rho_max; }
        bool is_uniform() const { return m_rho_uniform; }
        double d_rho() const { return m_d_rho; }
        size_t n_rho() const { return m_n_rho; }
        double n10_rho() const { return m_n10_rho; }
        double rho_min_pot() const { return m_rho_min_pot; }
        double rho_max_pot() const { return m_rho_max_pot; }
        double rho_vals(size_t i) const { return m_rho_vals[i];}

};