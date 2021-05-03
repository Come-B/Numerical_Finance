#ifndef MONTECARLO_H
#define MONTECARLO_H
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "rnr1_function.hpp"
#include "black_scholesND.hpp"

class MonteCarlo{
    public:
        MonteCarlo(RNR1Function* payoff, BlackScholesND* target_process);
        double compute_payoff_fixed_number(double start_time, double end_time, size_t nb_steps, int nb_simulation);
        double compute_payoff_conf(double start_time, double end_time, size_t nb_steps, double eps, double conf, int simu_step_size=1000);
        double static inv_cdf(const double& u);
        virtual ~MonteCarlo();

    protected:
        RNR1Function* m_payoff_function;
        BlackScholesND* m_target_process;
        const int MAX_ITER = 50000;
    private:
};

#endif // MONTECARLO_H
