#ifndef CONTINUOUS_GENERATOR_H
#define CONTINUOUS_GENERATOR_H

#include "random_generator.hpp"
#include "uniform_generator.hpp"
#include "quasi_random_generator.hpp"
#include <vector>

using namespace std;

class ContinuousGenerator : public RandomGenerator {
    public:
        ContinuousGenerator(double target_mean, double target_variance);
        ContinuousGenerator();
        ~ContinuousGenerator();
};

class Exponential : public ContinuousGenerator {
    protected:
        double m_lambda;
        bool m_expo_algo = true;
        UniformGenerator* m_uniform = new EcuyerCombined();
        const myLong MAX_ITER = 10000;

    public:
        Exponential(double lambda);
        Exponential(double lambda, bool expo_algo);
        Exponential(double lambda, bool expo_algo, UniformGenerator* uniform1);
        Exponential();
        double pdf(double x);
        ~Exponential();

        virtual double generate();
};

class Normal : public ContinuousGenerator {
    protected:
        double m_mu = 0.;
        double m_sigma = 1.;
        int m_normal_algo = 1;
        UniformGenerator* m_uniform = new EcuyerCombined();
        Exponential* m_exponential = new Exponential(1,true,m_uniform);

        bool m_has_spare_value = false;
        double m_spare_value;
        const myLong MAX_ITER = 100000;

        double pdf(double x);

    public:
        Normal(double mu, double sigma);
        Normal(double mu, double sigma, int normal_algo);
        Normal(double mu, double sigma, int normal_algo, UniformGenerator* uniform);
        Normal();
        double get_mu() {return m_mu;};
        double get_sigma() {return m_sigma;};
        ~Normal();

        virtual double generate();
};

class NormalND {
    protected:
        double m_mu = 0.;
        double m_sigma = 1.;
        int m_dimension= 1;
        HaltonSequence* m_uniform;

    public:
        NormalND(double mu, double sigma, HaltonSequence* uniform);
        double get_mu() {return m_mu;};
        double get_sigma() {return m_sigma;};
        ~NormalND();

        double generate_at_dim(int dim);
        RandomGenerator* generator_at(int dim);
};

#endif
