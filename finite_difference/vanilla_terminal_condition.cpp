#include "vanilla_terminal_condition.hpp"
#include <algorithm>

using namespace std;

VanillaTerminalCondition::VanillaTerminalCondition(double strike) : R1R1Function(), strike(strike) {};

VanillaTerminalCondition::~VanillaTerminalCondition() {};

CallTerminalCondition::CallTerminalCondition(double strike) : VanillaTerminalCondition(strike) {};

double CallTerminalCondition::operator()(double x) {
    return max(x - strike, 0.);
}

CallTerminalCondition::~CallTerminalCondition() {};

PutTerminalCondition::PutTerminalCondition(double strike) : VanillaTerminalCondition(strike) {};

double PutTerminalCondition::operator()(double x) {
    return max(strike - x, 0.);
}

PutTerminalCondition::~PutTerminalCondition() {};
