#ifndef SINGLE_PATH_H
#define SINGLE_PATH_H

#include <vector>

using namespace std;

class SinglePath {
    protected:
        vector<double> values;
        vector<double> times;
        double start_time;
        double end_time;
        size_t nb_steps;
        double time_step;

    public:
        SinglePath(double start, double end, size_t nb_steps);
        void add_value(double val);
        double get_value(double time);
        vector<double>& get_values() {return values;};
        vector<double>& get_times() {return times;};
        double get_time_step() {return time_step;};
        ~SinglePath();
};

#endif
