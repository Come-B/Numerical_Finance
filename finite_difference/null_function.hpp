#ifndef NULL_FUNCTION_H
#define NULL_FUNCTION_H

#include "r1r1_function.hpp"
#include "r2r1_function.hpp"

class NullFunction1D : public R1R1Function {
    public:
        NullFunction1D();
        double operator()(double t);
        ~NullFunction1D();
};

class NullFunction2D : public R2R1Function {
    public:
        NullFunction2D();
        double operator()(double x, double t);
        ~NullFunction2D();
};

#endif
