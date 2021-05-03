#include "pde_grid_2d.hpp"
#include <stdio.h>
#include <stdexcept>

PDEGrid2D::PDEGrid2D(double maturity,
                     double min_underlying_value,
                     double max_underlying_value,
                     int nbr_time_steps,
                     double underlying_step,
                     R2R1Function* variance_function,
                     R2R1Function* trend_function,
                     R2R1Function* actualization_function,
                     R2R1Function* source_term_function,
                     R1R1Function* top_boundary_function,
                     R1R1Function* bottom_boundary_function,
                     R1R1Function* right_boundary_function) :
                     T(maturity),
                     x_min(min_underlying_value),
                     x_max(max_underlying_value),
                     h0(maturity / nbr_time_steps),
                     h1(underlying_step),
                     nodes_width(nbr_time_steps + 1),
                     a(variance_function),
                     b(trend_function),
                     r(actualization_function),
                     f(source_term_function),
                     top_boundary_function(top_boundary_function),
                     bottom_boundary_function(bottom_boundary_function),
                     right_boundary_function(right_boundary_function) {
    nodes_height = (max_underlying_value - min_underlying_value) / h1 + 1;
    double norm_a = (*a)(x_max, T);
    if (h0 > (h1 * h1) / norm_a) {
        throw std::invalid_argument("the scheme cannot be stable");
    }
    nodes = vector< vector<double> > (nodes_width);
    for (size_t i = 0; i < nodes_width; ++i) {
        nodes[i] = vector<double> (nodes_height);
    }
};

void PDEGrid2D::fill_nodes() {
    fill_right_boundary();
    fill_top_and_bottom_boundaries();
};

double PDEGrid2D::get_time_zero_value(double spot) {
    return nodes[0][(int) (spot - x_min) / h1];
};

double PDEGrid2D::get_value(double spot, double time) {
    int spot_pos = (int) (spot - x_min) / h1;
    int time_pos = (int) ((double) time / h0);
    return nodes[time_pos][spot_pos];
};

void PDEGrid2D::fill_top_and_bottom_boundaries() {
    for (size_t i = 0; i < nodes_width; ++i) {
        nodes[i][nodes_height - 1] = (*top_boundary_function)(i * h0);
        nodes[i][0] = (*bottom_boundary_function)(i * h0);
    }
};

void PDEGrid2D::fill_right_boundary() {
    for (size_t j = 0; j < nodes_height; ++j) {
        nodes[nodes_width - 1][j] = (*right_boundary_function)(x_min + j * h1);
    }
};

void PDEGrid2D::display_nodes() {
    printf("\n The nodes matrix is :\n");
    printf("[");
    for (size_t j = 0; j < nodes_height; ++j) {
        if (j == 0) {
            printf("[");
        } else {
            printf(" [");
        }
        for (size_t i = 0; i < nodes_width; ++i) {
            if (i == nodes_width - 1) {
                printf("% 07.3f", nodes[i][j]);
            } else {
                printf("% 07.3f, ", nodes[i][j]);
            }
        }
        if (j == nodes_height - 1) {
            printf("]");
        } else {
            printf("],\n");
        }
    }
    printf("]\n\n");
}

PDEGrid2D::~PDEGrid2D() {};

PDEGridExplicit::PDEGridExplicit(double maturity,
                                 double minimum_underlying_value,
                                 double maximum_underlying_value,
                                 int nbr_time_steps,
                                 double underlying_step,
                                 R2R1Function* variance_function,
                                 R2R1Function* trend_function,
                                 R2R1Function* actualization_function,
                                 R2R1Function* source_term_function,
                                 R1R1Function* top_boundary_function,
                                 R1R1Function* bottom_boundary_function,
                                 R1R1Function* right_boundary_function) :
                                 PDEGrid2D(maturity,
                                           minimum_underlying_value,
                                           maximum_underlying_value,
                                           nbr_time_steps,
                                           underlying_step,
                                           variance_function,
                                           trend_function,
                                           actualization_function,
                                           source_term_function,
                                           top_boundary_function,
                                           bottom_boundary_function,
                                           right_boundary_function) {};

void PDEGridExplicit::fill_nodes() {
    PDEGrid2D::fill_nodes();
    double x;
    double t;
    double ajk_h0_to_h12;
    double bjk_h0_to_h1;
    for (size_t k = nodes_width - 1; k > 0; --k) {
        for (size_t j = 1; j < nodes_height - 1; ++j) {
            x = x_min + j * h1;
            t = k * h0;
            ajk_h0_to_h12 = h0 * (*a)(x, t) / (h1 * h1);
            bjk_h0_to_h1 = h0 * (*b)(x, t) / h1;
            nodes[k - 1][j] = nodes[k][j] * (1 - ajk_h0_to_h12 - bjk_h0_to_h1 - h0 * (*r)(x, t))
                            + nodes[k][j + 1] * (bjk_h0_to_h1 + 0.5 * ajk_h0_to_h12)
                            + nodes[k][j - 1] * (0.5 * ajk_h0_to_h12)
                            + h0 * (*f)(x, t);
        }
    }
};

PDEGridExplicit::~PDEGridExplicit() {};

PDEGridImplicit::PDEGridImplicit(double maturity,
                                 double minimum_underlying_value,
                                 double maximum_underlying_value,
                                 int nbr_time_steps,
                                 double underlying_step,
                                 R2R1Function* variance_function,
                                 R2R1Function* trend_function,
                                 R2R1Function* actualization_function,
                                 R2R1Function* source_term_function,
                                 R1R1Function* top_boundary_function,
                                 R1R1Function* bottom_boundary_function,
                                 R1R1Function* right_boundary_function) :
                                 PDEGrid2D(maturity,
                                           minimum_underlying_value,
                                           maximum_underlying_value,
                                           nbr_time_steps,
                                           underlying_step,
                                           variance_function,
                                           trend_function,
                                           actualization_function,
                                           source_term_function,
                                           top_boundary_function,
                                           bottom_boundary_function,
                                           right_boundary_function) {};

void PDEGridImplicit::fill_nodes() {
    PDEGrid2D::fill_nodes();
    // To implement
};

PDEGridImplicit::~PDEGridImplicit() {};

PDEGridTheta::PDEGridTheta(double maturity,
                           double minimum_underlying_value,
                           double maximum_underlying_value,
                           int nbr_time_steps,
                           double underlying_step,
                           R2R1Function* variance_function,
                           R2R1Function* trend_function,
                           R2R1Function* actualization_function,
                           R2R1Function* source_term_function,
                           R1R1Function* top_boundary_function,
                           R1R1Function* bottom_boundary_function,
                           R1R1Function* right_boundary_function,
                           double theta) :
                           PDEGrid2D(maturity,
                                     minimum_underlying_value,
                                     maximum_underlying_value,
                                     nbr_time_steps,
                                     underlying_step,
                                     variance_function,
                                     trend_function,
                                     actualization_function,
                                     source_term_function,
                                     top_boundary_function,
                                     bottom_boundary_function,
                                     right_boundary_function) {};

void PDEGridTheta::fill_nodes() {
    PDEGrid2D::fill_nodes();
    // To implement
};

PDEGridTheta::~PDEGridTheta() {};

