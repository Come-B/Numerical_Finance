#ifndef RNR1FUNCTION_H
#define RNR1FUNCTION_H
#include "../utils/matrix.hpp"

class RNR1Function {
    public:
        RNR1Function();
        virtual double operator()(Matrix<double> X) = 0;
        ~RNR1Function();
};

class BasketPayoff: public RNR1Function {
    protected:
        Matrix<double> m_coeffs;
        double m_strike;
    public:
        BasketPayoff(Matrix<double> coeffs, double strike=0.);
        double operator()(Matrix<double> X);
        double get_strike() {return m_strike;};
        Matrix<double> get_coeffs() {return m_coeffs;};
        ~BasketPayoff();
};

#endif
