#include "quasi_random_generator.hpp"

QuasiRandomGenerator::QuasiRandomGenerator() {};

QuasiRandomGenerator::~QuasiRandomGenerator() {};

VDCHSequence::VDCHSequence(int p, myLong n) : p(p), n(n) {};

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
