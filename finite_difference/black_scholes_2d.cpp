#include "black_scholes_2d.hpp"
#include "random_process.hpp"

BlackScholes2D::BlackScholes2D(Normal* gen, double spot1, double spot2,
                                            double rate1, double rate2,
                                            double vol1, double vol2, double rho) :
                                            RandomProcess(gen, 2), spot1(spot1), spot2(spot2),
                                            rate1(rate1), rate2(rate2), vol1(vol1), vol2(vol2), rho(rho) {};

BlackScholes2D::~BlackScholes2D() {};

BSMilstein2D::BSMilstein2D(Normal* gen, double spot1, double spot2, double rate1, double rate2, double vol1, double vol2, double rho) : BlackScholes2D(gen, spot1, spot2, rate1, rate2, vol1, vol2, rho) {};

void BSMilstein2D::simulate(double start_time, double end_time, size_t nb_steps) {
    // To implement
};

BSMilstein2D::~BSMilstein2D() {};
