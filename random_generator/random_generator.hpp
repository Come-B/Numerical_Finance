#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#define _USE_MATH_DEFINES

typedef unsigned long long myLong;

class RandomGenerator {
    protected:
        double target_mean;
        double target_variance;

    public:
        RandomGenerator(double target_mean, double target_variance);
        RandomGenerator();
        ~RandomGenerator();

        virtual double generate() = 0;

        double mean(myLong nb_sim);
        double get_target_mean() {return target_mean;};
        bool test_mean(myLong nb_sim, double tol, double* res_ptr = 0);

        double variance(myLong nb_sim);
        double get_target_variance() {return target_variance;};
        bool test_variance(myLong nb_sim, double tol, double* res_ptr = 0);
};

#endif
