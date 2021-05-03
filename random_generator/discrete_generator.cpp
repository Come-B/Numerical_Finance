#include "discrete_generator.hpp"
#include <stdexcept>
#include <cmath>

DiscreteGenerator::DiscreteGenerator(double target_mean, double target_variance) : RandomGenerator(target_mean, target_variance) {};

DiscreteGenerator::DiscreteGenerator() {};

DiscreteGenerator::~DiscreteGenerator() {};

Bernoulli::Bernoulli(double p) : DiscreteGenerator(p, p * (1 - p)), m_p(p) {};

Bernoulli::Bernoulli(double p, UniformGenerator* uniform) : DiscreteGenerator(p, p * (1 - p)), m_p(p), m_uniform(uniform) {};

Bernoulli::~Bernoulli() {};

double Bernoulli::generate() {
    if (m_uniform->generate() <= m_p)
        return 1.;
    return 0.;
};

HeadTail::HeadTail(UniformGenerator* uniform) : Bernoulli(0.5, uniform){};

HeadTail::HeadTail() : Bernoulli(0.5) {};

HeadTail::~HeadTail() {};

Binomial::Binomial(myLong n, double p) : DiscreteGenerator(n * p, n * p * (1 - p)), m_n(n), m_p(p) {};

Binomial::Binomial(myLong n, double p, UniformGenerator* uniform) : DiscreteGenerator(n * p, n * p * (1 - p)), m_n(n), m_p(p), m_uniform(uniform) {};

Binomial::~Binomial() {};

double Binomial::generate() {
    double sum = 0.;
    for (myLong i = 0; i < m_n; ++i)
        sum += m_bernoulli.generate();
    return sum;
};

void check_probabilities(const vector<double>& probabilities) {
    double s = 0.;
    double pi;
    for (size_t i = 0; i < probabilities.size(); i++) {
        pi = probabilities[i];
        if (pi < 0 || pi > 1) throw std::invalid_argument("a probability given is not between 0 and 1");
        s += pi;
    }
    if (s != 1.) throw std::invalid_argument("the sum of probabilities is not 1");
};

FiniteSet::FiniteSet(const vector<double>& input_probabilities) : m_probabilities(input_probabilities) {
    check_probabilities(m_probabilities);
    adjust_moments();
};

FiniteSet::FiniteSet(const vector<double>& input_probabilities, UniformGenerator* uniform) : m_probabilities(input_probabilities), m_uniform(uniform) {
    check_probabilities(m_probabilities);
    adjust_moments();
};

void FiniteSet::adjust_moments(){
    double mean = 0.;
    for (size_t i = 0; i < m_probabilities.size(); i++)
        mean += (i + 1.)*m_probabilities[i];
    target_mean = mean;

    double var = 0.;
    for (size_t i = 0; i < m_probabilities.size(); i++){
        double temp = ((i + 1.) - target_mean);
        var += temp*temp*m_probabilities[i];
    }
    target_variance = var;
}

FiniteSet::~FiniteSet() {};

double FiniteSet::generate() {
    double u = m_uniform->generate();
    double Pk = 0.;
    double Pk_1 = 0.;
    for (myLong i = 0; i < m_probabilities.size(); ++i) {
        Pk += m_probabilities[i];
        if (Pk_1 <= u && u < Pk)
            return i+1.;
        Pk_1 = Pk;
    }
   throw std::range_error("didn't encounter target probability");
};

Poisson::Poisson(double lambda) : DiscreteGenerator(lambda, lambda), m_lambda(lambda) {};

Poisson::Poisson(double lambda, bool poisson_algo) : DiscreteGenerator(lambda, lambda), m_lambda(lambda), m_poisson_algo(poisson_algo) {};

Poisson::Poisson(double lambda, bool poisson_algo, UniformGenerator* uniform) : DiscreteGenerator(lambda, lambda), m_lambda(lambda), m_poisson_algo(poisson_algo), m_uniform(uniform) {};

Poisson::Poisson() {};

Poisson::~Poisson() {};

double Poisson::generate() {
    if (m_poisson_algo) {
        double u = m_uniform->generate();
        double Pk = exp(-m_lambda);
        double Pk_1 = 0.;
        double pk = Pk;
        double k = 0.;
        while (u > Pk and u != Pk_1) {
            Pk_1 = Pk;
            pk = (m_lambda / (k + 1)) * pk;
            Pk += pk;
            k += 1.;
        };
        return k;
    } else {
        double s = 0.;
        for(myLong i = 0; i < 100000; i++){
            s += m_exponential.generate();
            if(s>1.)
                return (double) i;
        }
        throw std::range_error("Poisson simulation got out of range (max tried: 100000)");
    }
};
