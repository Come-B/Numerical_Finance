#ifndef VANILLA_BOUNDARIES_H
#define VANILLA_BOUNDARIES_H

#include "r1r1_function.hpp"

class CallTopBoundary : public R1R1Function {
    private:
        double s_max;
        double strike;
    
    public:
        CallTopBoundary(double s_max, double strike);
        double operator()(double t);
        ~CallTopBoundary();
};

class PutTopBoundary : public R1R1Function {
    private:
        double s_max;
        double strike;
    
    public:
        PutTopBoundary(double s_max, double strike);
        double operator()(double t);
        ~PutTopBoundary();
};

class CallBottomBoundary : public R1R1Function {
    private:
        double s_min;
        double strike;
    
    public:
        CallBottomBoundary(double s_min, double strike);
        double operator()(double t);
        ~CallBottomBoundary();
};

class PutBottomBoundary : public R1R1Function {
    private:
        double s_min;
        double strike;
    
    public:
        PutBottomBoundary(double s_min, double strike);
        double operator()(double t);
        ~PutBottomBoundary();
};

#endif
