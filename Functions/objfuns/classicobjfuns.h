/*
 *  GOProblems.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef CLASSICOBJFUNS_H
#define CLASSICOBJFUNS_H

#include <vector>

//NOTE: the functions here have passing by reference + const as they are called a lot of time during execution and thus
//it is worth trying to save time by avoiding to make a copy of the variable passed

double testfunction (const std::vector<double>& x);
double rastrigin (const std::vector<double>& x);
double schwefel (const std::vector<double>& x);
double ackley (const std::vector<double>& x);
double rosenbrock (const std::vector<double>& x);
double lennardjones (const std::vector<double>& x);
double levy (const std::vector<double>& x);

#endif
