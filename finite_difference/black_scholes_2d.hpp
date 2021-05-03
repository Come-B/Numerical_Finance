#ifndef BLACK_SCHOLES_2D_H
#define BLACK_SCHOLES_2D_H

#include "random_process.hpp"

class BlackScholes2D : public RandomProcess {
    protected:
        double spot1, spot2;
        double rate1, rate2;
        double vol1, vol2;
        double rho;
    
    public:
        BlackScholes2D(Normal* gen, double spot1, double spot2, double rate1, double rate2, double vol1, double vol2, double rho);
        ~BlackScholes2D();
};

class BSMilstein2D : public BlackScholes2D {
    public:
        BSMilstein2D(Normal* gen, double spot1, double spot2, double rate1, double rate2, double vol1, double vol2, double rho);
        void simulate(double start_time, double end_time, size_t nb_steps);
        ~BSMilstein2D();
};

#endif
