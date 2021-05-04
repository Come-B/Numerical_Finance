#include "continuous_generator.hpp"
#include <cmath>
#include <stdexcept>
#include "MonteCarlo.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_E
    #define M_E 2.71828182845904523536
#endif

ContinuousGenerator::ContinuousGenerator(double target_mean, double target_variance) : RandomGenerator(target_mean, target_variance) {};

ContinuousGenerator::ContinuousGenerator() {};

ContinuousGenerator::~ContinuousGenerator() {};

Exponential::Exponential(double lambda) : ContinuousGenerator(1 / lambda, 1 / (lambda * lambda)), m_lambda(lambda) {};

Exponential::Exponential(double lambda, bool expo_algo) : ContinuousGenerator(1 / lambda, 1 / (lambda * lambda)), m_lambda(lambda), m_expo_algo(expo_algo) {};

Exponential::Exponential(double lambda, bool expo_algo, UniformGenerator* uniform1) : ContinuousGenerator(1 / lambda, 1 / (lambda * lambda)), m_lambda(lambda), m_expo_algo(expo_algo), m_uniform(uniform1) {};

double Exponential::pdf(double x){
    return m_lambda*exp(-m_lambda*x);
}

double Exponential::generate() {
    if (m_expo_algo)
        return -log(m_uniform->generate()) / m_lambda;
    else {
        double x_max = -log(0.000001)/m_lambda;
        for(myLong it=0; it<MAX_ITER; it++){
            double X = x_max*m_uniform->generate();
            double Y = m_lambda*m_uniform->generate();
            if(Y<=pdf(X))
                return X;
        }
        throw std::overflow_error("no y have been lower than the density function of x");
    };
    return -42.;
};

Exponential::~Exponential() {};

Normal::Normal(double mu, double sigma) : ContinuousGenerator(mu, sigma), m_mu(mu), m_sigma(sigma) {};

Normal::Normal(double mu, double sigma, int normal_algo) : ContinuousGenerator(mu, sigma), m_mu(mu), m_sigma(sigma), m_normal_algo(normal_algo) {};

Normal::Normal(double mu, double sigma, int normal_algo, UniformGenerator* uniform) : ContinuousGenerator(mu, sigma), m_mu(mu), m_sigma(sigma), m_normal_algo(normal_algo), m_uniform(uniform), m_exponential(new Exponential(1,true,uniform)) {};

Normal::Normal() {};

double Normal::pdf(double x){
    return std::exp(-x*x/2.)/(2.*M_PI);
}

double Normal::generate() {
    if(m_has_spare_value){
        m_has_spare_value = false;
        return m_spare_value;
    }
    if (m_normal_algo == 2) { // Central Limit Theorem
        double s = 0.;
        for (int i = 0; i < 12; ++i)
            s += m_uniform->generate();
        return m_mu + m_sigma * (s - 6.);
    } else if (m_normal_algo == 3) { // Rejection sampling method
        for(myLong it=0; it<MAX_ITER; it++){
            double X = (2*(m_uniform->generate() < 0.5)-1)*m_exponential->generate();
            double Y = std::sqrt(2*M_E/M_PI)*m_exponential->pdf(std::abs(X))*m_uniform->generate()/2.;
            if(Y<=pdf(X))
                return m_mu + m_sigma * X;
        }
        throw std::overflow_error("no y have been lower than the density function of x");
    } else if (m_normal_algo == 4) { //Inverse distribution
        return m_mu + m_sigma * MonteCarlo::inv_cdf(m_uniform->generate());
    } else if (m_normal_algo == 5) { //Marsaglia polar
         double u,v,s;
        do{
            u = m_uniform->generate()*2. - 1.;
            v = m_uniform->generate()*2. - 1.;
            s = u*u + v*v;
        }while(s>1);
        m_spare_value = m_mu + m_sigma * v*std::sqrt(-2.*std::log(s)/s);
        m_has_spare_value = true;
        return m_mu + m_sigma * u*std::sqrt(-2.*std::log(s)/s);
    } else { // If normal_algo is anything else, we use Box-Müller algorithm
        double u1 = m_uniform->generate();
        double u2 = m_uniform->generate();
        double r = sqrt(-2 * log(u1));
        double theta = 2 * M_PI * u2;
        m_spare_value = m_mu + m_sigma * r * sin(theta);
        m_has_spare_value = true;
        return m_mu + m_sigma * r * cos(theta);
    };
};

Normal::~Normal() {};

NormalND::NormalND(double mu, double sigma, HaltonSequence* uniform) : m_uniform(uniform) {};

double NormalND::generate_at_dim(int dim){
    return m_mu + m_sigma * MonteCarlo::inv_cdf(m_uniform->generate_at_dim(dim));
};

RandomGenerator* NormalND::generator_at(int dim){
    return new Normal(m_mu,m_sigma,4,m_uniform->generator_at(dim));
}

NormalND::~NormalND() {};
