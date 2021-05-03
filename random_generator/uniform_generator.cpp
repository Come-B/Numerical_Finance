#include "uniform_generator.hpp"

UniformGenerator::UniformGenerator() : RandomGenerator(0.5, 1 / 12.) {};

UniformGenerator::~UniformGenerator() {};

PseudoGenerator::PseudoGenerator(myLong input_seed) : seed(input_seed), current(input_seed) {};

PseudoGenerator::PseudoGenerator() {};

PseudoGenerator::~PseudoGenerator() {};

LinearCongruential::LinearCongruential(myLong seed, myLong multiplier,
                                       myLong increment, myLong modulus) :
                                       PseudoGenerator(seed), multiplier(multiplier),
                                       increment(increment), modulus(modulus) {};

LinearCongruential::LinearCongruential(myLong seed) : PseudoGenerator(seed) {};

LinearCongruential::LinearCongruential() {};

LinearCongruential::~LinearCongruential() {};

double LinearCongruential::generate() {
    double result = (double) current / modulus;
    current = (current * multiplier + increment) % modulus;
    return result;
};

myLong LinearCongruential::getModulus() {
    return modulus;
};

myLong LinearCongruential::getCurrent() {
    return current;
};

EcuyerCombined::EcuyerCombined(LinearCongruential lcm1, LinearCongruential lcm2) : lcm1(lcm1), lcm2(lcm2) {};

EcuyerCombined::EcuyerCombined(myLong seed) : seed(seed) {};

EcuyerCombined::EcuyerCombined() {};

EcuyerCombined::~EcuyerCombined() {};

double EcuyerCombined::generate() {
    lcm1.generate();
    lcm2.generate();
    myLong lcm1_current = lcm1.getCurrent();
    myLong lcm2_current = lcm2.getCurrent();
    myLong lcm1_modulus = lcm1.getModulus();
    myLong result = (lcm1_current - lcm2_current) % (lcm1_modulus - 1);
    if (result == 0)
        return (double) (lcm1_modulus - 1) / lcm1_modulus;
    return (double) result / lcm1_modulus;
};
