#ifndef DISCRETE_GENERATOR_H
#define DISCRETE_GENERATOR_H

#include "random_generator.hpp"
#include "uniform_generator.hpp"
#include "continuous_generator.hpp"
#include <vector>

using namespace std;

class DiscreteGenerator : public RandomGenerator {
    public:
        DiscreteGenerator(double target_mean, double target_variance);
        DiscreteGenerator();
        ~DiscreteGenerator();
};

class Bernoulli : public DiscreteGenerator {
    private:
        double m_p;
        UniformGenerator* m_uniform = new EcuyerCombined();

    public:
        Bernoulli(double p);
        Bernoulli(double p, UniformGenerator* uniform);
        ~Bernoulli();

        virtual double generate();
};

class HeadTail : public Bernoulli {
    public:
        HeadTail(UniformGenerator* uniform);
        HeadTail();
        ~HeadTail();
};

class Binomial : public DiscreteGenerator {
    private:
        myLong m_n;
        double m_p;
        UniformGenerator* m_uniform = new EcuyerCombined();
        Bernoulli m_bernoulli = Bernoulli(m_p, m_uniform);

    public:
        Binomial(myLong n, double p);
        Binomial(myLong n, double p, UniformGenerator* uniform);
        ~Binomial();

        virtual double generate();
};

class FiniteSet : public DiscreteGenerator {
    protected:
        vector<double> m_probabilities;
        UniformGenerator* m_uniform = new EcuyerCombined();
        void adjust_moments();

    public:
        FiniteSet(const vector<double>& probabilities);
        FiniteSet(const vector<double>& probabilities, UniformGenerator* uniform);
        ~FiniteSet();

        virtual double generate();
};

class Poisson : public DiscreteGenerator {
    protected:
        double m_lambda;
        bool m_poisson_algo = true;
        UniformGenerator* m_uniform = new EcuyerCombined();
        Exponential m_exponential = Exponential(m_lambda, m_uniform);

    public:
        Poisson(double lambda);
        Poisson(double lambda, bool poisson_algo);
        Poisson(double lambda, bool poisson_algo, UniformGenerator* uniform);
        Poisson();
        ~Poisson();

        virtual double generate();
};

#endif
