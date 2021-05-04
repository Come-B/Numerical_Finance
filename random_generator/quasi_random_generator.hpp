#ifndef QUASI_RANDOM_GENERATOR_H
#define QUASI_RANDOM_GENERATOR_H

#include "uniform_generator.hpp"
#include "../utils/matrix.hpp"
#include <vector>

class QuasiRandomGenerator : public UniformGenerator {
    public:
        QuasiRandomGenerator(double target_mean, double target_variance);
        QuasiRandomGenerator();
        ~QuasiRandomGenerator();
};

class VDCHSequence : public QuasiRandomGenerator {
    protected:
        int p = 2;
        myLong n = 1;

    public:
        VDCHSequence(int p);
        VDCHSequence();
        ~VDCHSequence();

        virtual double generate();
};

class HaltonSequence {
    protected:
        int m_dimension;
        std::vector<VDCHSequence*> m_sequences;

    public:
        HaltonSequence(int dim = 1);
        double generate_at_dim(int dim);
        UniformGenerator* generator_at(int dim);
        ~HaltonSequence();

        static const int primes[15];
};

#endif
