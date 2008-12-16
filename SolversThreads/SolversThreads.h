/*
 *  SolverThreads.h
 *  PaGMO
 *
 *  Created by Dario Izzo on 9/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef SOLVERSTHREADS_H
#define SOLVERSTHREADS_H

#include <pthread.h>

#include "GOproblem.h"
#include "population.h"

//Here we define the parameters needed to instanciate a thread. These contain
//datas that are algorithm specific, but also data that are needed for all aglorithms (LB,UB,objfun,mutex etc.)

struct threadParam{
	//Thread unique ID
	unsigned int threadID;

	//Solvers Data
	int NP;
	int generations;
	//DE
	int strategy;
	double F,CR;
	//PSO
	double omega,eta1,eta2,vcoeff;
	int nswarms;
	//GA
	double M,CRsga;
	int insert_best;
	//SA-AN
	double Ts,Tf;
	//pointers giving access to global resources
	GOProblem* problem;

	bool *isActive;
	pthread_mutex_t *TPmutex;
	pthread_cond_t *exit;
	Population *Ptr_pop;
	std::ofstream *Ptr_log;
	unsigned long randomSeed;
};

//Here we define the protoypes for each type of thread we may want to open
void *DEthread(void *data);
void *PSOthread(void *data);
void *MPSOthread(void *data);
void *SGAthread(void *data);
void *ASAthread(void *data);
#endif
