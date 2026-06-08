#pragma once

#include <vector>

class Grid {
    private:
        std::vector<double> m_rho_vals;     // stores the rho grid
        double  m_rho_min;                  // smallest rho value
        double  m_rho_max;                  // largest rho value
        bool    m_rho_uniform;              // linear grid / log grid
        double  m_d_rho;                    // rho grid spacing for linear grid
        double  m_n_rho;                    // number of rho points for the linear grid
        double  m_n10_rho;                  // number of rho points per decade for the log grid

    public:
        double d1(const std::vector<double>& u, size_t i);
        double d2(const std::vector<double>& u , size_t i);

        std::vector<double> d1(const std::vector<double>& u);
        std::vector<double> d2(const std::vector<double>& u);

        // setters
        // ------------------------------
        void set_rho_vals(std::vector<double> rho_vals);
        void set_rho_min(double rho_min);
        void set_rho_max(double rho_max);
        void set_d_rho(double d_rho);
        void set_n_rho(double n_rho);
        void set_n10_rho(double n10_rho);
        
};