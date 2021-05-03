#ifndef VANILLA_TERMINAL_CONDITION_H
#define VANILLA_TERMINAL_CONDITION_H

#include "r1r1_function.hpp"

class VanillaTerminalCondition : public R1R1Function {
    protected:
        double strike;
    
    public:
        VanillaTerminalCondition(double strike);
        ~VanillaTerminalCondition();
};

class CallTerminalCondition : public VanillaTerminalCondition {
    public:
        CallTerminalCondition(double strike);
        double operator()(double x);
        ~CallTerminalCondition();
};

class PutTerminalCondition : public VanillaTerminalCondition {
    public:
        PutTerminalCondition(double strike);
        double operator()(double x);
        ~PutTerminalCondition();
};

#endif
