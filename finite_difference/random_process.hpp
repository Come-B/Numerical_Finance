#ifndef RANDOM_PROCESS_H
#define RANDOM_PROCESS_H

#include "single_path.hpp"
#include "../random_generator/continuous_generator.hpp"
#include <cmath>

class RandomProcess {
    protected:
        RandomGenerator* generator;
        vector<SinglePath*> paths;
        int dimension;

    public:
        RandomProcess(RandomGenerator* gen, int dim = 1);
        virtual void simulate(double start_time, double end_time, size_t nb_steps) = 0;
        SinglePath* get_path(int index);
        ~RandomProcess();
};

#endif
