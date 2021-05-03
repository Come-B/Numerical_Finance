#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(RNR1Function* payoff, BlackScholesND* target_process):
    m_payoff_function(payoff),
    m_target_process(target_process){

}

double MonteCarlo::compute_payoff_fixed_number(double start_time, double end_time, size_t nb_steps, int nb_simulation){
    vector<double> res;
    double mean=0.;
    for(int i=0;i<nb_simulation;i++){
        m_target_process->simulate(start_time,end_time,nb_steps);
        res.push_back((*m_payoff_function)(m_target_process->spot_at_maturity()));
        mean+=res.back();
        //printf("Final spots:\n%s\n",m_target_process->spot_at_maturity().to_string().c_str());
    }
    mean /= nb_simulation;
    double var=0.;
    for(int i=0;i<nb_simulation;i++)
        var+=(res.at(i)-mean)*(res.at(i)-mean);

    var /= (nb_simulation-1);
    printf("Finished Monte Carlo with %i simulations: average=%4.2f variance=%4.2f\n", nb_simulation, mean, var);

    double normal_quantile_95 = inv_cdf(0.975);
    Matrix<double> confidence_interval(1,2,new double[2]{-1,1});
    confidence_interval=confidence_interval*(normal_quantile_95*std::sqrt(var)/std::sqrt(nb_simulation));
    confidence_interval=confidence_interval+mean;
    printf("... 95%% confidence interval: %s\n", confidence_interval.to_string().c_str());
    return mean;
}

double MonteCarlo::inv_cdf(const double& u){
     static double a[4]={ 2.50662823884,
                         -18.61500062529,
                         41.39119773534,
                         -25.44106049637};
     static double b[4]={-8.47351093090,
                         23.08336743743,
                         -21.06224101826,
                         3.13082909833};
     static double c[9]={0.3374754822726147,
                         0.9761690190917186,
                         0.1607979714918209,
                         0.0276438810333863,
                         0.0038405729373609,
                         0.0003951896511919,
                         0.0000321767881768,
                         0.0000002888167364,
                         0.0000003960315187};
     double x=u-0.5;
     double r;
     if(fabs(x)<0.42){  // Beasley-Springer
        double y=x*x;
        r=x*(((a[3]*y+a[2])*y+a[1])*y+a[0])/((((b[3]*y+b[2])*y+b[1])*y+b[0])*y+1.0);
     }else{ // Moro
        r=u;
        if (x>0.0)
            r=1.0-u;
        r=log(-log(r));
        r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
        if (x<0.0)
            r=-r;
    }
    return r;
}

double MonteCarlo::compute_payoff_conf(double start_time, double end_time, size_t nb_steps, double eps, double conf, int simu_step_size){
    vector<double> res;
    double running_tot=0.;
    double running_mean=0.;
    double running_dist=0.;
    double running_var=0.;
    double normal_quantile = inv_cdf((1+conf)/2);
    double boundary = normal_quantile*normal_quantile/(eps*eps);
    int iteration=1;
    do{
        m_target_process->simulate(start_time,end_time,nb_steps);
        res.push_back((*m_payoff_function)(m_target_process->spot_at_maturity()));
        running_tot+=res.back();
        running_mean = running_tot/iteration;
        if(iteration%simu_step_size==0){
            running_dist=0.;
        for(int i=0;i<res.size();i++)
            running_dist+=(res.at(i)-running_mean)*(res.at(i)-running_mean);
        running_var = running_dist/(iteration-1);
        }
    }while((iteration<simu_step_size || iteration<boundary*running_var) && iteration++<MAX_ITER);

    printf("Finished Monte Carlo for error=%4.2f and confidence=%4.2f with %i simulations: average=%4.2f variance=%4.2f\n",eps,conf,iteration, running_mean, running_var);
    Matrix<double> confidence_interval(1,2,new double[2]{-1,1});
    confidence_interval=confidence_interval*(normal_quantile*std::sqrt(running_var)/std::sqrt(iteration));
    confidence_interval=confidence_interval+running_mean;
    printf("... %i%% confidence interval: %s\n", std::lround(100*conf),confidence_interval.to_string().c_str());
    return running_mean;
}

MonteCarlo::~MonteCarlo(){
}
