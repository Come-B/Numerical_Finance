#include "vanilla_boundaries.hpp"
#include <algorithm>

using namespace std;

CallTopBoundary::CallTopBoundary(double s_max, double strike) : R1R1Function(), s_max(s_max), strike(strike) {};

double CallTopBoundary::operator()(double t) {
    return max(s_max - strike, 0.);
};

CallTopBoundary::~CallTopBoundary() {};

PutTopBoundary::PutTopBoundary(double s_max, double strike) : R1R1Function(), s_max(s_max), strike(strike) {};

double PutTopBoundary::operator()(double t) {
    return max(strike - s_max, 0.);
};

PutTopBoundary::~PutTopBoundary() {};

CallBottomBoundary::CallBottomBoundary(double s_min, double strike) : R1R1Function(), s_min(s_min), strike(strike) {};

double CallBottomBoundary::operator()(double t) {
    return max(s_min - strike, 0.);
};

CallBottomBoundary::~CallBottomBoundary() {};

PutBottomBoundary::PutBottomBoundary(double s_min, double strike) : R1R1Function(), s_min(s_min), strike(strike) {};

double PutBottomBoundary::operator()(double t) {
    return max(strike - s_min, 0.);
};

PutBottomBoundary::~PutBottomBoundary() {};
