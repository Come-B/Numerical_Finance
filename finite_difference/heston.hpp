#ifndef HESTON_H
#define HESTON_H

#include "random_process.hpp"

class Heston : public RandomProcess {
    protected:
        double spot;
        double init_variance;
        double mu;
        double theta;
        double kappa;
        double sigma;
        double rho;
    
    public:
        Heston(Normal* gen, double spot, double init_variance, double mu, double theta, double kappa, double sigma, double rho);
        void simulate(double strat_time, double end_time, size_t nb_steps);
        ~Heston();
};

#endif
