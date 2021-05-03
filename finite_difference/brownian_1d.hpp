#ifndef BROWNIAN_1D_H
#define BROWNIAN_1D_H

#include "random_process.hpp"

class Brownian1D : public RandomProcess {
    public:
        Brownian1D(Normal* gen);
        void simulate(double start_time, double end_time, size_t nb_steps);
        ~Brownian1D();
};

#endif
