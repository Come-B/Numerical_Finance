#include "bs.hpp"

BSActualization::BSActualization(double rate) : rate(rate) {};

double BSActualization::operator()(double x, double t) {
    return rate;
};

BSActualization::~BSActualization() {};

BSVariance::BSVariance(double sigma) : sigma(sigma) {};

double BSVariance::operator()(double x, double t) {
    return sigma * sigma * x * x;
};

BSVariance::~BSVariance() {};

BSTrend::BSTrend(double rate) : rate(rate) {};

double BSTrend::operator()(double x, double t) {
    return rate * x;
};

BSTrend::~BSTrend() {};
