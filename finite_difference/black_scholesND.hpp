#ifndef BLACKSCHOLESND_H
#define BLACKSCHOLESND_H

#include "random_process.hpp"
#include "../utils/matrix.hpp"

class BlackScholesND : public RandomProcess {
    protected:
        Matrix<double> m_initial_spots_matrix;
        Matrix<double> m_rates_matrix;
        Matrix<double> m_vols_matrix;
        Matrix<double> m_corrs_matrix;
        Matrix<double> *m_spot_at_maturity;
        int m_dimension;

    public:
        BlackScholesND(Normal* gen, const Matrix<double> initial_spots, const Matrix<double> rates, const Matrix<double> vols, const Matrix<double> corrs);
        void simulate(double start_time, double end_time, size_t nb_steps);
        Matrix<double> spot_at_inception() const {return m_initial_spots_matrix;};
        Matrix<double> spot_at_maturity() const {return *m_spot_at_maturity;};
        ~BlackScholesND();

};

#endif
