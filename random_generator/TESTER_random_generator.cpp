#include "TESTER_random_generator.hpp"

void compute_generation(RandomGenerator* generator, string generator_name) {
    myLong n_sim = 100000;
    double tol = 0.03;
    std::cout << generator_name << " :" << std::endl;
    std::cout << generator->generate() << std::endl;
    std::cout << generator->generate() << std::endl;
    std::cout << generator->generate() << std::endl;

    double* var_ptr = new double;
    double* mean_ptr = new double;
    string correct_mean = generator->test_mean(n_sim, tol, mean_ptr) ? "________ OK" : "________ ERROR";
    printf("Mean : %.4f vs theo. mean: %.4f %s\n", *mean_ptr,generator->get_target_mean(), correct_mean.c_str());
    string correct_variance = generator->test_variance(n_sim, tol, var_ptr) ? "____ OK" : "____ ERROR";
    printf("Variance : %.4f vs theo. var: %.4f %s\n", *var_ptr,generator->get_target_variance(), correct_variance.c_str());
    std::cout << std::endl;
};

void test_random_generator() {
    RandomGenerator* lcm = new LinearCongruential();
    compute_generation(lcm, "Linear congruential");

    RandomGenerator* clcg = new EcuyerCombined();
    compute_generation(clcg, "Ecuyer combined");

    RandomGenerator* head_tail = new HeadTail();
    compute_generation(head_tail, "Head or tail");

    RandomGenerator* bernoulli = new Bernoulli(0.8);
    compute_generation(bernoulli, "Bernoulli");

    RandomGenerator* binomial = new Binomial(100, 0.7);
    compute_generation(binomial, "Binomial");

    /*
    vector<double> probas {0.2, 0.5, 0.3};
    RandomGenerator* finite_set = new FiniteSet(probas);
    compute_generation(finite_set, "Finite set");
    */

    RandomGenerator* poisson1 = new Poisson(2.);
    compute_generation(poisson1, "Poisson (1st algo)");

    RandomGenerator* poisson2 = new Poisson(2., false);
    compute_generation(poisson2, "Poisson (2nd algo)");

    RandomGenerator* exponential1 = new Exponential(2.);
    compute_generation(exponential1, "Exponential (inverse function)");

    RandomGenerator* exponential2 = new Exponential(2., false);
    compute_generation(exponential2, "Exponential (rejection sampling)");

    RandomGenerator* normal1 = new Normal(1., 1.);
    compute_generation(normal1, "Normal (Box-Müller algorithm)");

    RandomGenerator* normal2 = new Normal(1., 1., 2);
    compute_generation(normal2, "Normal (Central Limit Theorem)");

    RandomGenerator* normal3 = new Normal(1., 1., 3);
    compute_generation(normal3, "Normal (rejection sampling)");

    RandomGenerator* normal7 = new Normal(1., 1., 4);
    compute_generation(normal7, "Normal (Inverse distribution)");

    RandomGenerator* normal9 = new Normal(1., 1., 5);
    compute_generation(normal9, "Normal (Marsaglia polar)");

    RandomGenerator* vdch_sequence = new VDCHSequence(5);
    compute_generation(vdch_sequence, "Van Der Corput and Halton sequence 5");

    RandomGenerator* normal4 = new Normal(1., 1.,1,new VDCHSequence(5));
    compute_generation(normal4, "Normal (Box-Müller algorithm) with VDSH");

    RandomGenerator* normal5 = new Normal(1., 1.,2,new VDCHSequence(5));
    compute_generation(normal5, "Normal (Central Limit Theorem) with VDSH");

    RandomGenerator* normal6 = new Normal(1., 1.,3,new VDCHSequence(5));
    compute_generation(normal6, "Normal (rejection sampling) with VDSH");

    RandomGenerator* normal8 = new Normal(1., 1.,4,new VDCHSequence(5));
    compute_generation(normal8, "Normal (Inverse distribution) with VDSH");

    RandomGenerator* normal10 = new Normal(1., 1.,5,new VDCHSequence(5));
    compute_generation(normal10, "Normal (Marsaglia polar) with VDSH");
};
