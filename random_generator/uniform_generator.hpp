#ifndef UNIFORM_GENERATOR_H
#define UNIFORM_GENERATOR_H

#include "random_generator.hpp"

class UniformGenerator : public RandomGenerator {
    public:
        UniformGenerator();
        ~UniformGenerator();
};

class PseudoGenerator : public UniformGenerator {
    protected:
        myLong seed = 272727272;
        myLong current = seed;

    public:
        PseudoGenerator(myLong input_seed);
        PseudoGenerator();
        ~PseudoGenerator();
};

class LinearCongruential : public PseudoGenerator {
    private:
        myLong multiplier = 40014;
        myLong increment = 0;
        myLong modulus = 2147483563;

    public:
        LinearCongruential(myLong seed,
                           myLong multiplier,
                           myLong increment,
                           myLong modulus);
        LinearCongruential(myLong seed);
        LinearCongruential();
        ~LinearCongruential();
        
        virtual double generate();
        myLong getModulus();
        myLong getCurrent();
};

class EcuyerCombined : public PseudoGenerator {
    private:
        myLong seed = 272727272;
        LinearCongruential lcm1 = LinearCongruential(seed);
        LinearCongruential lcm2 = LinearCongruential(seed, 40692, 0, 2147483399);

    public:
        EcuyerCombined(LinearCongruential lcm1, LinearCongruential lcm2);
        EcuyerCombined(myLong seed);
        EcuyerCombined();
        ~EcuyerCombined();

        virtual double generate();
};

#endif
