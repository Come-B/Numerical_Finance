#ifndef BS_H
#define BS_H

#include "r2r1_function.hpp"

class BSActualization : public R2R1Function {
    private:
        double rate;
    
    public:
        BSActualization(double rate);
        double operator()(double x, double t);
        ~BSActualization();
};

class BSVariance : public R2R1Function {
    private:
        double sigma;
    
    public:
        BSVariance(double sigma);
        double operator()(double x, double t);
        ~BSVariance();
};

class BSTrend : public R2R1Function {
    private:
        double rate;
    
    public:
        BSTrend(double rate);
        double operator()(double x, double t);
        ~BSTrend();
};

#endif
