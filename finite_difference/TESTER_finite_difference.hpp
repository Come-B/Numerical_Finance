#ifndef TESTER_FINITE_DIFFERENCE_C
#define TESTER_FINITE_DIFFERENCE_C

#include "../random_generator/quasi_random_generator.hpp"
#include "single_path.hpp"
#include "random_process.hpp"
#include "brownian_1d.hpp"
#include "brownian_nd.hpp"
#include "black_scholes_1d.hpp"
#include "black_scholes_2d.hpp"
#include "black_scholesND.hpp"
#include "heston.hpp"
#include "r2r1_function.hpp"
#include "bs.hpp"
#include "vanilla_boundaries.hpp"
#include "vanilla_terminal_condition.hpp"
#include "null_function.hpp"
#include "pde_grid_2d.hpp"
#include "MonteCarlo.hpp"
#include "static_control_variate.hpp"
#include "../utils/matrix.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>

using namespace std;

void display_values(SinglePath* path, string name, bool five_first_five_last);

void test_finite_difference();

#endif
