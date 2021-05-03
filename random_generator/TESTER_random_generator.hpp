#ifndef TESTER_RANDOM_GENERATOR_C
#define TESTER_RANDOM_GENERATOR_C

#include "uniform_generator.hpp"
#include "discrete_generator.hpp"
#include "continuous_generator.hpp"
#include "quasi_random_generator.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <type_traits>

template<typename Base, typename T>
inline bool instanceof(const T*) {
   return std::is_base_of<Base, T>::value;
};

void compute_generation(RandomGenerator* generator, std::string generator_name);

void test_random_generator();

#endif
