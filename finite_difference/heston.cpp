#include "heston.hpp"
#include "random_process.hpp"

Heston::Heston(Normal* gen, double spot, double init_variance,
                                     double mu, double theta, double kappa,
                                     double sigma, double rho) :
                                     RandomProcess(gen, 1), spot(spot),
                                     init_variance(init_variance), mu(mu),
                                     theta(theta), kappa(kappa), sigma(sigma), rho(rho) {};

void Heston::simulate(double start_time, double end_time, size_t nb_steps) {
    SinglePath* spot_path = new SinglePath(start_time, end_time, nb_steps);
    SinglePath* var_path = new SinglePath(start_time, end_time, nb_steps);
    spot_path->add_value(spot);
    var_path->add_value(init_variance);
    double dt = spot_path->get_time_step();
    double last_spot_inserted = spot;
    double last_var_inserted = init_variance;
    double z1, z2, generated_x, generated_y, next_spot_value, next_var_value;
    for (size_t i = 0; i < nb_steps; ++i) {
        z1 = generator->generate();
        z2 = generator->generate();
        generated_x = z1;
        generated_y = rho * z1 + sqrt(1 - pow(rho, 2)) * z2;
        next_spot_value = last_spot_inserted * (1 + mu * dt + sqrt(last_var_inserted) * generated_x * sqrt(dt));
        next_var_value = last_var_inserted + kappa * (theta - last_var_inserted) * dt + sigma * sqrt(last_var_inserted) * generated_y * sqrt(dt);
        spot_path->add_value(next_spot_value);
        var_path->add_value(next_var_value);
        last_spot_inserted = next_spot_value;
        last_var_inserted = next_var_value;
    }
    paths.push_back(spot_path);
    paths.push_back(var_path);
};

Heston::~Heston() {};