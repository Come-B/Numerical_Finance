#include "null_function.hpp"

NullFunction1D::NullFunction1D() {};

double NullFunction1D::operator()(double t) {
    return 0.;
}

NullFunction1D::~NullFunction1D() {};

NullFunction2D::NullFunction2D() {};

double NullFunction2D::operator()(double x, double t) {
    return 0.;
}

NullFunction2D::~NullFunction2D() {};
