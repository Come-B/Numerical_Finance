#include "black_scholes_1d.hpp"
#include "random_process.hpp"

BlackScholes1D::BlackScholes1D(Normal* gen, double spot, double rate, double vol) : RandomProcess(gen), spot(spot), rate(rate), vol(vol) {};

BlackScholes1D::~BlackScholes1D() {};

BSEuler1D::BSEuler1D(Normal* gen, double spot, double rate, double vol) : BlackScholes1D(gen, spot, rate, vol) {}

void BSEuler1D::simulate(double start_time, double end_time, size_t nb_steps) {
    SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
    path->add_value(spot);
    double dt = path->get_time_step();
    double last_inserted = spot;
    for (size_t i = 0; i < nb_steps; ++i) {
        double next_value = last_inserted * (1 + rate*dt + vol*generator->generate()*std::sqrt(dt));
        path->add_value(next_value);
        last_inserted = next_value;
    }
    paths.push_back(path);
}

BSEuler1D::~BSEuler1D() {};

BSMilstein1D::BSMilstein1D(Normal* gen, double spot, double rate, double vol) : BlackScholes1D(gen, spot, rate, vol) {}

void BSMilstein1D::simulate(double start_time, double end_time, size_t nb_steps) {
    SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
    path->add_value(spot);
    double dt = path->get_time_step();
    double last_inserted = spot;
    double generated;
    for (size_t i = 0; i < nb_steps; ++i) {
        generated = generator->generate();
        double next_value = last_inserted * (1 + rate*dt + vol*generated*std::sqrt(dt));
        path->add_value(next_value);
        last_inserted = next_value;
    }
    paths.push_back(path);
}

BSMilstein1D::~BSMilstein1D() {};
