#ifndef RADIATOR_H_INCLUDED
#define RADIATOR_H_INCLUDED
#include<vector>

#define C_NUMBER_LAYERS 11
#define C_NUMBER_ACTIVE_LAYERS 2

double radiator (const std::vector<double>& x);
std::vector<double> var2struct(const std::vector<double> &x);


#endif // RADIATOR_H_INCLUDED
