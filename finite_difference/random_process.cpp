#include "random_process.hpp"

RandomProcess::RandomProcess(RandomGenerator* gen, int dim) : generator(gen), dimension(dim) {};

SinglePath* RandomProcess::get_path(int index) {
    return paths.at(index);
};

RandomProcess::~RandomProcess() {};
