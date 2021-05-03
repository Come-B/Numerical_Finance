#include "brownian_1d.hpp"
#include "random_process.hpp"

Brownian1D::Brownian1D(Normal* gen) : RandomProcess(gen) {};

void Brownian1D::simulate(double start_time, double end_time, size_t nb_steps) {
    SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
    path->add_value(0.);
    double dt = path->get_time_step();
    double last_inserted = 0.;
    double next_value;
    for (size_t i = 0; i < nb_steps; i++) {
        next_value = last_inserted + generator->generate() * std::sqrt(dt);
        path->add_value(next_value);
        last_inserted = next_value;
    }
    paths.push_back(path);
}

Brownian1D::~Brownian1D() {};
