#include "black_scholesND.hpp"
#include <iostream>

BlackScholesND::BlackScholesND(Normal* gen, const Matrix<double> initial_spots, const Matrix<double> rates, const Matrix<double> vols, const Matrix<double> corrs) :
    RandomProcess(gen, initial_spots.get_num_rows()),
    m_dimension(initial_spots.get_num_rows()){
    m_initial_spots_matrix = initial_spots;
    m_rates_matrix = rates;
    m_vols_matrix = vols;
    m_corrs_matrix = corrs;
    m_spot_at_maturity = new Matrix<double>(initial_spots.get_num_rows(),1);

}

void BlackScholesND::simulate(double start_time, double end_time, size_t nb_steps) {
    paths.clear();
    for (int i = 0; i < m_dimension; i++)
        paths.push_back(new SinglePath(start_time, end_time, nb_steps));

    double dt = paths.back()->get_time_step();
    Matrix<double> cholesky = m_corrs_matrix.cholesky();
    Matrix<double> random_law(dimension, nb_steps);
    for(int i = 0; i < dimension; i++)
        for(size_t j = 0; j < nb_steps; j++)
            random_law.set_elem_at(i, j, generator->generate());

    Matrix<double> current_spot(m_dimension,1);
    for(int i = 0; i < m_dimension; i++)
        current_spot.set_elem_at(i,0,m_initial_spots_matrix.elem_at(i,0));

    for(size_t j = 0; j < nb_steps; j++) {
        for(int i = 0; i < m_dimension; i++)
           paths.at(i)->add_value(current_spot.elem_at(i,0));
        current_spot  = current_spot*(m_rates_matrix*dt + m_vols_matrix*(cholesky.dot(random_law.extract_column(j)))*std::sqrt(dt) + 1.);
    }
    for(int i = 0; i < m_dimension; i++)
        m_spot_at_maturity->set_elem_at(i,0,paths.at(i)->get_value(end_time));
}

BlackScholesND::~BlackScholesND(){

}
