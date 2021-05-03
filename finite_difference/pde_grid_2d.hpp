#ifndef PDE_GRID_2D_H
#define PDE_GRID_2D_H

#include "r2r1_function.hpp"
#include "r1r1_function.hpp"
#include <vector>

using namespace std;

class PDEGrid2D {
    protected:
        double T;
        double x_min;
        double x_max;
        double h0;
        double h1;
        R2R1Function* a;
        R2R1Function* b;
        R2R1Function* r;
        R2R1Function* f;
        R1R1Function* top_boundary_function;
        R1R1Function* bottom_boundary_function;
        R1R1Function* right_boundary_function;
        vector< vector<double> > nodes;
        size_t nodes_height;
        size_t nodes_width;
    
    public:
        PDEGrid2D(double maturity,
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
                  R1R1Function* right_boundary_function);
        virtual void fill_nodes();
        double get_time_zero_value(double spot);
        double get_value(double spot, double time);
        void display_nodes();
        ~PDEGrid2D();

    private:
        void fill_top_and_bottom_boundaries();
        void fill_right_boundary();
};

class PDEGridExplicit : public PDEGrid2D {
    public:
        PDEGridExplicit(double maturity,
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
                        R1R1Function* right_boundary_function);
        void fill_nodes();
        ~PDEGridExplicit();
};

class PDEGridImplicit : public PDEGrid2D {
    public:
        PDEGridImplicit(double maturity,
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
                        R1R1Function* right_boundary_function);
        void fill_nodes();
        ~PDEGridImplicit();
};

class PDEGridTheta : public PDEGrid2D {
    public:
        PDEGridTheta(double maturity,
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
                     double theta);
        void fill_nodes();
        ~PDEGridTheta();
};

#endif
