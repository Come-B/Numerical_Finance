#ifndef BROWNIAN_ND_H
#define BROWNIAN_ND_H

#include "random_process.hpp"
#include "../utils/matrix.hpp"

class BrownianND : public RandomProcess {
    protected:
        Matrix<double>* correlation_matrix;

    public:
        BrownianND(Normal* gen, int dim, Matrix<double>* cor);
        void simulate(double start_time, double end_time, size_t nb_steps);
        ~BrownianND();
};

#endif
