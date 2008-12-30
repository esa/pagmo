/*
 *  SolverThreads.cpp
 *  PaGMO
 *
 *  Created by Dario Izzo on 9/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <fstream>
#include <vector>

#include "ASA.h"
#include "DE.h"
#include "MPSO.h"
#include "PSO.h"
#include "SGA.h"
#include "SolversThreads.h"
#include "population.h"
#include "rng.h"

using namespace std;

// Shortcut definition for the lock type.
typedef boost::unique_lock<boost::mutex> lock_type;

//******************************************************************************************
//DE thread type
//******************************************************************************************

void *DEthread(void *data)
{
   threadParam *PtrTP = (threadParam *)data;
   Population deme;
   double oldfitness;
   vector<size_t> picks;
   GOProblem* problem;
   const vector<double> &LB = PtrTP->problem->getLB();
   DEalgorithm DE;
   rng_uint32 rng;
   rng_double drng;


	clock_t start,end;
	double dif;

    {
        lock_type lock(*PtrTP->TPmutex);
        rng.seed(PtrTP->randomSeed);
        drng.seed(PtrTP->randomSeed);
		deme = PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		problem = PtrTP->problem;
		DE.initDE(PtrTP->generations,LB.size(),PtrTP->F,PtrTP->CR,PtrTP->strategy, rng());
    }

    oldfitness = deme.extractBestIndividual().getFitness();

    start=clock();
   deme = DE.evolve(deme, *problem);
   end=clock();
   dif = (double)(end-start) / (double)CLOCKS_PER_SEC;

   {
        lock_type lock(*PtrTP->TPmutex);
        //insert deme in main population
		PtrTP->Ptr_pop->insertDeme(deme,picks);
		//log in cout
		cout << "Thread ID: " <<  PtrTP->threadID << endl ;
		cout << "\t\t\tDE:\t\t\t F "<< PtrTP->F  <<  "\t\tCR: " << PtrTP->CR << "\t\tStrategy: " << PtrTP->strategy << "\t\tGenerations: " << PtrTP->generations << endl;
		cout << "\t\t\tInitial fitness: " << oldfitness << endl;
		cout << "\t\t\tFinal fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
		cout << "\t\t\tSeconds elapsed: " << dif << endl;
		cout << "\t\t\tShutting down" << endl;
		//log in logfile
		*(PtrTP->Ptr_log) << PtrTP->generations * deme.size() << " " << (PtrTP->Ptr_pop->extractBestIndividual()).getFitness() << endl;
		//sinal exit
		--*(PtrTP->isActive);
		PtrTP->exit->notify_one();
   }
   return 0;
}





//******************************************************************************************
//MPSO thread type
//******************************************************************************************

void *MPSOthread(void *data)
{
   threadParam *PtrTP;
   PtrTP = (threadParam *)data;
   Population deme;
   double oldfitness;
   vector<size_t> picks;
   GOProblem* problem;
   const vector<double> &LB = PtrTP->problem->getLB();
   MPSOalgorithm MPSO;
   rng_uint32 rng;
   rng_double drng;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
        rng.seed(PtrTP->randomSeed);
        drng.seed(PtrTP->randomSeed);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		problem = PtrTP->problem;
		MPSO.initMPSO(PtrTP->generations,LB.size(),PtrTP->omega,PtrTP->eta1,PtrTP->eta2,PtrTP->vcoeff, PtrTP->nswarms, rng());
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = MPSO.evolve(deme, *problem);
   end=clock();
   dif = (double)(end-start) / (double)CLOCKS_PER_SEC;

    {
        lock_type lock(*PtrTP->TPmutex);
        //insert deme in main population
		PtrTP->Ptr_pop->insertDeme(deme,picks);
		//log in cout
		cout << "Thread ID: " <<  PtrTP->threadID << endl ;
		cout << "\t\t\tMPSO:\t\t omega "<< PtrTP->omega  <<  "\t\teta1: " << PtrTP->eta1 <<  "\t\teta2: " << PtrTP->eta2 << "\t\tVcoeff: " << PtrTP->vcoeff<< "\t\tNswarms" <<PtrTP->vcoeff << "\t\tGenerations: " << PtrTP->generations << endl;
		cout << "\t\t\tInitial fitness: " << oldfitness << endl;
		cout << "\t\t\tFinal fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
		cout << "\t\t\tSeconds elapsed: " << dif << endl;
		cout << "\t\t\tShutting down" << endl;
		//log in logfile
		*(PtrTP->Ptr_log) << PtrTP->generations * deme.size() << " " << (PtrTP->Ptr_pop->extractBestIndividual()).getFitness() << endl;
		//sinal exit
		--*(PtrTP->isActive);
        PtrTP->exit->notify_one();
    }
   return 0;
}

//******************************************************************************************
//PSO thread type
//******************************************************************************************

void *PSOthread(void *data)
{
   threadParam *PtrTP;
   PtrTP = (threadParam *)data;
   Population deme;
   double oldfitness;
   vector<size_t> picks;
   GOProblem* problem;
   const vector<double> &LB = PtrTP->problem->getLB();
   PSOalgorithm PSO;
   rng_uint32 rng;
   rng_double drng;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
        rng.seed(PtrTP->randomSeed);
        drng.seed(PtrTP->randomSeed);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		problem = PtrTP->problem;
		PSO.initPSO(PtrTP->generations,LB.size(),PtrTP->omega,PtrTP->eta1,PtrTP->eta2,PtrTP->vcoeff, rng());
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = PSO.evolve(deme, *problem);
   end=clock();
   dif = (double)(end-start) / (double)CLOCKS_PER_SEC;

   {
        lock_type lock(*PtrTP->TPmutex);
        //insert deme in main population
		PtrTP->Ptr_pop->insertDeme(deme,picks);
		//log in cout
		cout << "Thread ID: " <<  PtrTP->threadID << endl ;
		cout << "\t\t\tPSO:\t\t omega "<< PtrTP->omega  <<  "\t\teta1: " << PtrTP->eta1 <<  "\t\teta2: " << PtrTP->eta2 << "\t\tVcoeff: " << PtrTP->vcoeff << "\t\tGenerations: " << PtrTP->generations << endl;
		cout << "\t\t\tInitial fitness: " << oldfitness << endl;
		cout << "\t\t\tFinal fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
		cout << "\t\t\tSeconds elapsed: " << dif << endl;
		cout << "\t\t\tShutting down" << endl;
		//log in logfile
		*(PtrTP->Ptr_log) << PtrTP->generations * deme.size() << " " << (PtrTP->Ptr_pop->extractBestIndividual()).getFitness() << endl;
		//sinal exit
		--*(PtrTP->isActive);
		PtrTP->exit->notify_one();
   }
   return 0;
}



//******************************************************************************************
//SGA thread type
//******************************************************************************************

void *SGAthread(void *data)
{
   threadParam *PtrTP;
   PtrTP = (threadParam *)data;
   Population deme;
   double oldfitness;
   vector<size_t> picks;
   const vector<double> &LB = PtrTP->problem->getLB();
   SGAalgorithm SGA;
   rng_uint32 rng;
   rng_double drng;
   GOProblem* problem;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
        rng.seed(PtrTP->randomSeed);
        drng.seed(PtrTP->randomSeed);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		problem = PtrTP->problem;
		SGA.initSGA(PtrTP->generations,LB.size(),PtrTP->CRsga,PtrTP->M,PtrTP->insert_best, rng());
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = SGA.evolve(deme, *problem);
   end=clock();
   dif = (double)(end-start) / (double)CLOCKS_PER_SEC;

   {
        lock_type lock(*PtrTP->TPmutex);
        //insert deme in main population
		PtrTP->Ptr_pop->insertDeme(deme,picks);
		//log in cout
		cout << "Thread ID: " <<  PtrTP->threadID << endl ;
		cout << "\t\t\tSGA:\t\t CR "<< PtrTP->CRsga  <<  "\t\tM: " << PtrTP->M <<  "\t\tInsertBest: " << PtrTP->insert_best << "\t\tGenerations: " << PtrTP->generations << endl;
		cout << "\t\t\tInitial fitness: " << oldfitness << endl;
		cout << "\t\t\tFinal fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
		cout << "\t\t\tSeconds elapsed: " << dif << endl;
		cout << "\t\t\tShutting down" << endl;
		//log in logfile
		*(PtrTP->Ptr_log) << PtrTP->generations * deme.size() << " " << (PtrTP->Ptr_pop->extractBestIndividual()).getFitness() << endl;
		//sinal exit
		--*(PtrTP->isActive);
		PtrTP->exit->notify_one();
   }
   return 0;
}


//******************************************************************************************
//SA-AN thread type
//******************************************************************************************

void *ASAthread(void *data)
{
   threadParam *PtrTP;
   PtrTP = (threadParam *)data;
   Population deme;
   double oldfitness;
   vector<size_t> picks;
   const vector<double> &LB = PtrTP->problem->getLB();
   ASAalgorithm ASA;
   rng_uint32 rng;
   rng_double drng;
   GOProblem *problem = PtrTP->problem;

	clock_t start,end;
	double dif;

   {
	lock_type lock(*PtrTP->TPmutex);
	rng.seed(PtrTP->randomSeed);
	drng.seed(PtrTP->randomSeed);
	deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
	ASA.initASA(PtrTP->generations,LB.size(),PtrTP->Ts,PtrTP->Tf, rng());
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = ASA.evolve(deme[0], *problem);
   end=clock();
   dif = (double)(end-start) / (double)CLOCKS_PER_SEC;

   {
        lock_type lock(*PtrTP->TPmutex);
        //insert deme in main population
		PtrTP->Ptr_pop->insertDeme(deme,picks);
		//log in cout
		cout << "Thread ID: " <<  PtrTP->threadID << endl ;
		cout << "\t\t\tASA:\t\t\t Ts "<< PtrTP->Ts  <<  "\t\tTf: " << PtrTP->Tf << "\t\tFunction Evaluations: " << PtrTP->generations << endl;
		cout << "\t\t\tInitial fitness: " << oldfitness << endl;
		cout << "\t\t\tFinal fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
		cout << "\t\t\tSeconds elapsed: " << dif << endl;
		cout << "\t\t\tShutting down" << endl;
		//log in logfile
		*(PtrTP->Ptr_log) << PtrTP->generations * deme.size() << " " << (PtrTP->Ptr_pop->extractBestIndividual()).getFitness() << endl;
		//sinal exit
		--*PtrTP->isActive;
		PtrTP->exit->notify_one();
   }
   return 0;
}


