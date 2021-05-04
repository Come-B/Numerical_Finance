#include "static_control_variate.hpp"
#include "pde_grid_2d.hpp"
#include <cmath>
#include <iostream>

StaticControlVariate::StaticControlVariate(BasketPayoff* payoff,
                                           BlackScholesND* target_process,
                                           double start,
                                           double end,
                                           size_t nb_steps)
                     : MonteCarlo(payoff, target_process),
                       m_start_time(start),
                       m_end_time(end),
                       m_nb_steps(nb_steps) {
    m_r = 1.;
    for (int i = 0; i < m_target_process->get_dim(); i++) {
        m_r *= std::pow(m_target_process->get_rates().elem_at(i, 0), m_payoff_function->get_coeffs().elem_at(0,i));
    }
    compute_exp_htp();
};

void StaticControlVariate::compute_exp_htp() {
    double price_products = 1.;
    for (int i = 0; i < m_target_process->get_dim(); i++) {
        price_products *= std::pow(m_target_process->get_initial_spots().elem_at(i, 0), m_payoff_function->get_coeffs().elem_at(0, i));
    }
    printf("Products: %6.6f\n", price_products);

    Matrix<double> chol = m_target_process->get_corrs().cholesky();
    Matrix<double> chol_t = chol.transpose();
    Matrix<double> weitghs_t = m_payoff_function->get_coeffs().transpose();
    double var = m_payoff_function->get_coeffs().dot(chol.dot(chol_t.dot(weitghs_t))).elem_at(0,0);
    
    double implied_rate = var / 2. + m_r;
    for (int i = 0; i < m_target_process->get_dim(); i++) {
        implied_rate -= m_payoff_function->get_coeffs().elem_at(0, i) * chol.elem_at(i, 0) * chol.elem_at(i, 0) / 2.;
    }
    printf("Rate: %6.6f\n", implied_rate);

    double implied_vol = std::sqrt(var);
    printf("Vol: %6.6f\n", implied_vol);

    m_exp_htp = (std::log(price_products / m_payoff_function->get_strike()) + implied_rate * m_end_time) / (implied_vol * std::sqrt(m_end_time));
    printf("m_exp_htp: %6.6f\n", m_exp_htp);
};

double StaticControlVariate::joint_simulation() {
    m_target_process->simulate(m_start_time, m_end_time, m_nb_steps);

    double sum_prices = 0.;
    for (int i = 0; i < m_target_process->get_dim(); i++) {
        sum_prices += m_payoff_function->get_coeffs().elem_at(0, i) * m_target_process->get_path(i)->get_value(m_end_time);
    }
    printf("sum_prices: %6.6f\n", sum_prices);

    double sum_log_prices = 0.;
    for (int i = 0; i < m_target_process->get_dim(); i++) {
        sum_log_prices += m_payoff_function->get_coeffs().elem_at(0, i) * std::log(m_target_process->get_path(i)->get_value(m_end_time));
    }
    printf("sum_log_prices: %6.6f\n", sum_log_prices);
    double exp_sum_log_prices = std::exp(sum_log_prices);
    printf("exp_sum_log_prices: %6.6f\n", exp_sum_log_prices);
    
    double ht;
    if (sum_prices - m_payoff_function->get_strike() >= 0) {
        ht = sum_prices - m_payoff_function->get_strike();
    } else {
        ht = 0.;
    }
    printf("ht: %6.6f\n", ht);

    double htp;
    if (exp_sum_log_prices - m_payoff_function->get_strike() >= 0) {
        htp = exp_sum_log_prices - m_payoff_function->get_strike();
    } else {
        htp = 0.;
    }
    printf("htp: %6.6f\n", htp);

    return std::exp(-m_r * m_end_time) * (ht - htp + m_exp_htp);
};

StaticControlVariate::~StaticControlVariate() {};
