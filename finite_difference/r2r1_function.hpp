#ifndef R2R1_FUNCTION_H
#define R2R1_FUNCTION_H

using namespace std;

class R2R1Function {
    public:
        R2R1Function();
        virtual double operator()(double x, double t) = 0;
        ~R2R1Function();
};

#endif