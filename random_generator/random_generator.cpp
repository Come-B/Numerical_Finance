#include "random_generator.hpp"
#include <iostream>
#include <vector>
#include <cmath>

RandomGenerator::RandomGenerator(double target_mean, double target_variance) : target_mean(target_mean), target_variance(target_variance) {};

RandomGenerator::RandomGenerator() : target_mean(0.), target_variance(1.) {};

RandomGenerator::~RandomGenerator() {};

double RandomGenerator::mean(myLong nb_sim){
    double sum = 0.;
    for (myLong i = 0; i < nb_sim; ++i)
        sum += generate();
    return sum / nb_sim;
};

bool RandomGenerator::test_mean(myLong nb_sim, double tol, double* res_ptr){
    double experimental_mean = mean(nb_sim);
    if(res_ptr)
        *res_ptr = experimental_mean;
    return (std::abs(experimental_mean - target_mean) <= tol);
};

double RandomGenerator::variance(myLong nb_sim){
    double mean = 0.;
    std::vector<double> generated_numbers;
    for (myLong i = 0; i < nb_sim; ++i){
        double generated = generate();
        generated_numbers.push_back(generated);
        mean += generated;
    }
    mean = mean / nb_sim;

    double sum_dist = 0.;
    for (myLong i = 0; i < nb_sim; ++i){
        double dista = generated_numbers[i] - mean;
        sum_dist += (dista*dista);
    }
    return sum_dist / nb_sim;
};

bool RandomGenerator::test_variance(myLong nb_sim, double tol, double* res_ptr){
    double experimental_var = variance(nb_sim);
    if(res_ptr)
        *res_ptr = experimental_var;
    return (std::abs(experimental_var - target_variance) <= tol);
};
