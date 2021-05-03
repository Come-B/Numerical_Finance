#ifndef R1R1_FUNCTION_H
#define R1R1_FUNCTION_H

using namespace std;

class R1R1Function {
    public:
        R1R1Function();
        virtual double operator()(double x) = 0;
        ~R1R1Function();
};

#endif