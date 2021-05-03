#include "rnr1_function.hpp"
#include <iostream>

RNR1Function::RNR1Function() {}

RNR1Function::~RNR1Function() {}

BasketPayoff::BasketPayoff(Matrix<double> coeffs, double strike){
    m_coeffs = coeffs;
    m_strike = strike;
}

double BasketPayoff::operator()(Matrix<double> X){
    return (m_coeffs.dot(X)).elem_at(0,0)-m_strike;
}

BasketPayoff::~BasketPayoff() {}
