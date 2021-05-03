#include "single_path.hpp"
#include <vector>
#include <stdexcept>

SinglePath::SinglePath(double start, double end, size_t nb_steps) :
    start_time(start),
    end_time(end),
    nb_steps(nb_steps) {
        if (nb_steps == 0)
            throw invalid_argument("the number of steps is null");
        time_step = (end_time - start_time) / nb_steps;
}

void SinglePath::add_value(double val) {
    values.push_back(val);
    if (times.size() == 0) {
        times.push_back(start_time);
    } else {
        double previous_time = times.back();
        times.push_back(previous_time + time_step);
    }
}

double SinglePath::get_value(double time) {
    if (time <= start_time) {
        return values[0];
    } else if (time >= end_time) {
        return values.back();
    } else {
        for (size_t i = 1; i < times.size(); ++i) {
            if (times[i] == time) {
                return values[i];
            } else if (times[i] > time) {
                double u_time = times[i];
                double l_time = times[i - 1];
                double u_value = values[i];
                double l_value = values[i - 1];
                return l_value * (u_time - time) / (u_time - l_time) + u_value * (time - l_time) / (u_time - l_time);
            }
        }
    }
}

SinglePath::~SinglePath() {
    values.clear();
    times.clear();
}
