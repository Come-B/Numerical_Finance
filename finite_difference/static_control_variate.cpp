#include "static_control_variate.hpp"
#include "pde_grid_2d.hpp"
#include <cmath>

StaticControlVariate::StaticControlVariate(Normal* gen,
                                           Matrix<double> initial_spots,
                                           Matrix<double> rates,
                                           Matrix<double> vols,
                                           Matrix<double> corrs,
                                           Matrix<double> weights)
                     : BlackScholesND(gen,
                                      initial_spots,
                                      rates,
                                      vols,
                                      corrs){
    m_weights = weights;
};

double StaticControlVariate::joint_simulation(double start_time, double end_time, size_t nb_steps, double strike) {
    simulate(start_time, end_time, nb_steps);

    double price_products = 1.;
    for (int i = 0; i < m_dimension; i++) {
//        printf("m_initial_spots_matrix->elem_at(i, 0): %6.6f\n", m_initial_spots_matrix.elem_at(i, 0));
//        printf("m_weights.elem_at(i, 0): %6.6f\n", m_weights.elem_at(i, 0));
        price_products *= std::pow(m_initial_spots_matrix.elem_at(i, 0), m_weights.elem_at(i, 0));
    }
//    printf("Products: %6.6f\n", price_products);

    double r = 0.01;
    Matrix<double> weitghs_t = m_weights.transpose();
    Matrix<double> m_vols_matrix_t = m_vols_matrix.transpose();
    double implied_rate = weitghs_t.dot(m_vols_matrix.dot(m_vols_matrix_t.dot(m_weights))).elem_at(0, 0) * (1 / 2.) + r;
    for (int i = 0; i < m_dimension; i++) {
        implied_rate += m_weights.elem_at(i, 0) * m_vols_matrix.elem_at(i, 0) * m_vols_matrix.elem_at(i, 0) / 2.;
    }
    printf("Rate: %6.6f\n", implied_rate);

    double implied_vol = std::sqrt(weitghs_t.dot(m_vols_matrix.dot(m_vols_matrix_t.dot(m_weights))).elem_at(0, 0));
    printf("Vol: %6.6f\n", implied_vol);

    double exp_htp = (std::log(price_products / strike) + implied_rate * end_time) / (implied_vol * std::sqrt(end_time));
    printf("exp_htp: %6.6f\n", exp_htp);

    double sum_prices = 0.;
    for (int i = 0; i < m_dimension; i++) {
        sum_prices += m_weights.elem_at(i, 0) * get_path(i)->get_value(end_time);
    }
    printf("sum_prices: %6.6f\n", sum_prices);
    double sum_log_prices = 0.;
    for (int i = 0; i < m_dimension; i++) {
        sum_log_prices += m_weights.elem_at(i, 0) * std::log(get_path(i)->get_value(end_time));
    }
    printf("sum_log_prices: %6.6f\n", sum_log_prices);
    double exp_sum_log_prices = std::exp(sum_log_prices);
    printf("exp_sum_log_prices: %6.6f\n", exp_sum_log_prices);

    double ht;
    if (sum_prices - strike >= 0) {
        ht = sum_prices - strike;
    } else {
        ht = 0.;
    }
    printf("ht: %6.6f\n", ht);

    double htp;
    if (exp_sum_log_prices - strike >= 0) {
        htp = exp_sum_log_prices - strike;
    } else {
        htp = 0.;
    }
    printf("htp: %6.6f\n", htp);

    return std::exp(-r * end_time) * (ht - htp + exp_htp);
};

StaticControlVariate::~StaticControlVariate() {};
