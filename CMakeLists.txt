cmake_minimum_required(VERSION 3.1)
project(NumericalFinance)

include(CheckCXXCompilerFlag)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion -Wsign-conversion")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wunused-parameter -Wextra -Wreorder -Wconversion")
    CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

    if (HAS_CPP14_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    else()
        message(FATAL_ERROR "Unsupported compiler -- NumericalFinance requires C++14 support!")
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif()

set(NUMERICAL_FINANCE_SOURCES
    "random_generator/random_generator.hpp"
    "random_generator/random_generator.cpp"
    "random_generator/uniform_generator.hpp"
    "random_generator/uniform_generator.cpp"
    "random_generator/discrete_generator.hpp"
    "random_generator/discrete_generator.cpp"
    "random_generator/continuous_generator.hpp"
    "random_generator/continuous_generator.cpp"
    "random_generator/quasi_random_generator.hpp"
    "random_generator/quasi_random_generator.cpp"
    "random_generator/TESTER_random_generator.hpp"
    "random_generator/TESTER_random_generator.cpp"

    "finite_difference/single_path.hpp"
    "finite_difference/single_path.cpp"
    "finite_difference/random_process.hpp"
    "finite_difference/random_process.cpp"
    "finite_difference/r1r1_function.hpp"
    "finite_difference/r1r1_function.cpp"
    "finite_difference/r2r1_function.hpp"
    "finite_difference/r2r1_function.cpp"
    "finite_difference/vanilla_terminal_condition.hpp"
    "finite_difference/vanilla_terminal_condition.cpp"
    "finite_difference/vanilla_boundaries.hpp"
    "finite_difference/vanilla_boundaries.cpp"
    "finite_difference/null_function.hpp"
    "finite_difference/null_function.cpp"
    "finite_difference/black_scholes_1d.hpp"
    "finite_difference/black_scholes_1d.cpp"
    "finite_difference/black_scholes_2d.hpp"
    "finite_difference/black_scholes_2d.cpp"
    "finite_difference/black_scholesND.hpp"
    "finite_difference/black_scholesND.cpp"
    "finite_difference/brownian_1d.hpp"
    "finite_difference/brownian_1d.cpp"
    "finite_difference/brownian_nd.hpp"
    "finite_difference/brownian_nd.cpp"
    "finite_difference/bs.hpp"
    "finite_difference/bs.cpp"
    "finite_difference/heston.hpp"
    "finite_difference/heston.cpp"
    "finite_difference/pde_grid_2d.hpp"
    "finite_difference/pde_grid_2d.cpp"
    "finite_difference/MonteCarlo.hpp"
    "finite_difference/MonteCarlo.cpp"
    "finite_difference/rnr1_function.hpp"
    "finite_difference/rnr1_function.cpp"
    "finite_difference/TESTER_finite_difference.hpp"
    "finite_difference/TESTER_finite_difference.cpp"

    "utils/matrix.hpp"
    "utils/matrix.cpp"
    
    main.cpp
)

set(NUMERICAL_FINANCE_TARGET NumericalFinance)
add_executable(${NUMERICAL_FINANCE_TARGET} ${NUMERICAL_FINANCE_SOURCES})

