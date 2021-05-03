#ifndef STATIC_CONTROL_VARIATE_H
#define STATIC_CONTROL_VARIATE_H

#include "black_scholesND.hpp"
#include "pde_grid_2d.hpp"
#include "../utils/matrix.hpp"
#include <math.h>

class StaticControlVariate : public BlackScholesND {
    protected:
        Matrix<double> m_weights;

    public:
        StaticControlVariate(Normal* gen,
                             Matrix<double> initial_spots,
                             Matrix<double> rates,
                             Matrix<double> vols,
                             Matrix<double> corrs,
                             Matrix<double> weights);
        ~StaticControlVariate();

        double joint_simulation(double start_time, double end_time, size_t nb_steps, double strike);
};

#endif
