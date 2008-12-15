/*
 *  GOProblems.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/17/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef TRAJOBJFUNS_H
#define TRAJOBJFUNS_H
#include<vector>

//NOTE: the functions here have passing by reference + const as they are called a lot of time during execution and thus
//it is worth trying to save time by avoiding to make a copy of the variable passed

double messenger (const std::vector<double>& x);
double rosetta (const std::vector<double>& x);
double gtoc1 (const std::vector<double>& x);
double cassini1 (const std::vector<double>& x);
double cassini2 (const std::vector<double>& x);
double sagas (const std::vector<double>& x);


#endif
