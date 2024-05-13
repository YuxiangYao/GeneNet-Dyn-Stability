#ifndef RANDOM_SETTING_H_INCLUDED
#define RANDOM_SETTING_H_INCLUDED
#include <random>
std::mt19937 mt(12345);
std::uniform_int_distribution<int> r_binary(0,1);
std::uniform_real_distribution<double> runif(0,1);
#endif // RANDOM_SETTING_H_INCLUDED
