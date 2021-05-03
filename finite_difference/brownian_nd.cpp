#include "brownian_nd.hpp"
#include "random_process.hpp"

BrownianND::BrownianND(Normal* gen, int dim, Matrix<double>* cor) : RandomProcess(gen, dim), correlation_matrix(cor) {};

void BrownianND::simulate(double start_time, double end_time, size_t nb_steps) {
    for (int i = 0; i < dimension; i++)
        paths.push_back(new SinglePath(start_time, end_time, nb_steps));

    Matrix<double> cholesky = correlation_matrix->cholesky();
    Matrix<double> draws(dimension, nb_steps);
    for(int i = 0; i < dimension; i++)
        for(int j = 0; j < nb_steps; j++)
            draws.set_elem_at(i, j, generator->generate());

    Matrix<double> prod(dimension, nb_steps);
    prod = cholesky.dot(draws);
    prod = prod * std::sqrt(paths.back()->get_time_step());
    double next_value;
    for(int i = 0; i < dimension; i++) {
        double last_inserted = 0.;
        for(int j = 0; j < nb_steps; j++) {
            next_value = last_inserted + prod.elem_at(i, j);
            paths[i]->add_value(next_value);
            last_inserted = next_value;
        }
    }
}

BrownianND::~BrownianND() {};
