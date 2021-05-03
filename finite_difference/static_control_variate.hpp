#ifndef STATIC_CONTROL_VARIATE_H
#define STATIC_CONTROL_VARIATE_H

#include "MonteCarlo.hpp"
#include "pde_grid_2d.hpp"
#include "../utils/matrix.hpp"
#include "rnr1_function.hpp"
#include <math.h>

class StaticControlVariate : public MonteCarlo {
    protected:
        double m_exp_htp;
        double m_start_time;
        double m_end_time;
        size_t m_nb_steps;
        double m_r;
    
    public:
        StaticControlVariate(BasketPayoff* payoff,
                             BlackScholesND* target_process,
                             double start,
                             double end,
                             size_t nb_steps);
        ~StaticControlVariate();

        void compute_exp_htp();
        double joint_simulation();
};

#endif
