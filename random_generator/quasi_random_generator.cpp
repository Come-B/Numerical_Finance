#include "quasi_random_generator.hpp"
#include <iostream>

QuasiRandomGenerator::QuasiRandomGenerator() {};

QuasiRandomGenerator::~QuasiRandomGenerator() {};

VDCHSequence::VDCHSequence(int p) : p(p) {};

VDCHSequence::VDCHSequence() {};

double VDCHSequence::generate() {
    double q = 0.;
    double pk_1 = (double) 1 / p;
    myLong n_gen = n;

    while (n_gen > 0) {
      q += (n_gen % p) * pk_1;
      n_gen /= p;
      pk_1 /= p;
    }
    n++;
    return q;
};

VDCHSequence::~VDCHSequence() {};

const int HaltonSequence::primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
HaltonSequence::HaltonSequence(int dim): m_dimension(dim){
    for(int i=0;i<m_dimension;i++)
        m_sequences.push_back(new VDCHSequence(primes[i]));
}

double HaltonSequence::generate_at_dim(int dim){
    return m_sequences.at(dim)->generate();
}

UniformGenerator* HaltonSequence::generator_at(int dim){
    return m_sequences.at(dim);
}

HaltonSequence::~HaltonSequence() {};
