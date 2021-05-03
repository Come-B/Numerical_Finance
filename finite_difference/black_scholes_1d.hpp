#ifndef BLACK_SCHOLES_1D_H
#define BLACK_SCHOLES_1D_H

#include "random_process.hpp"

class BlackScholes1D : public RandomProcess {
    protected:
        double spot;
        double rate;
        double vol;
    
    public:
        BlackScholes1D(Normal* gen, double spot, double rate, double vol);
        ~BlackScholes1D();
};

class BSEuler1D : public BlackScholes1D {
    public:
        BSEuler1D(Normal* gen, double spot, double rate, double vol);
        void simulate(double start_time, double end_time, size_t nb_steps);
        ~BSEuler1D();
};

class BSMilstein1D : public BlackScholes1D {
    public:
        BSMilstein1D(Normal* gen, double spot, double rate, double vol);
        void simulate(double start_time, double end_time, size_t nb_steps);
        ~BSMilstein1D();
};

#endif
