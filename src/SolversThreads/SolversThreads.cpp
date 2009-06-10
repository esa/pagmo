/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// 21/00/2008: Initial version by Dario Izzo.

#include <boost/scoped_ptr.hpp>
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
   double oldfitness;
   vector<size_t> picks;
   GOProblem* problem;
   Population deme(*PtrTP->problem,0);
   boost::scoped_ptr<DEalgorithm> DE;
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
		DE.reset(new DEalgorithm(PtrTP->generations, PtrTP->F,PtrTP->CR,PtrTP->strategy));
    }

    oldfitness = deme.extractBestIndividual().getFitness();

    start=clock();
   deme = DE->evolve(deme);
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
   double oldfitness;
   vector<size_t> picks;
   GOProblem* problem;
   Population deme(*PtrTP->problem,0);
   boost::scoped_ptr<MPSOalgorithm> MPSO;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		problem = PtrTP->problem;
		MPSO.reset(new MPSOalgorithm(PtrTP->generations,PtrTP->omega,PtrTP->eta1,PtrTP->eta2,PtrTP->vcoeff, PtrTP->nswarms));
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = MPSO->evolve(deme);
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
   double oldfitness;
   vector<size_t> picks;
   Population deme(*PtrTP->problem,0);
   boost::scoped_ptr<PSOalgorithm> PSO;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		PSO.reset(new PSOalgorithm(PtrTP->generations,PtrTP->omega,PtrTP->eta1,PtrTP->eta2,PtrTP->vcoeff));
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = PSO->evolve(deme);
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
   double oldfitness;
   vector<size_t> picks;
   Population deme(*PtrTP->problem,0);
   boost::scoped_ptr<SGAalgorithm> SGA;

	clock_t start,end;
	double dif;

   {
        lock_type lock(*PtrTP->TPmutex);
		deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
		SGA.reset(new SGAalgorithm (PtrTP->generations,PtrTP->CRsga,PtrTP->M,PtrTP->insert_best));
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = SGA->evolve(deme);
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
   double oldfitness;
   vector<size_t> picks;
   boost::scoped_ptr<ASAalgorithm> ASA_ptr;
   rng_uint32 rng;
   rng_double drng;
   GOProblem *problem = PtrTP->problem;
   Population deme(*problem,0);

	clock_t start,end;
	double dif;

   {
	lock_type lock(*PtrTP->TPmutex);
	rng.seed(PtrTP->randomSeed);
	drng.seed(PtrTP->randomSeed);
	deme=PtrTP->Ptr_pop->extractRandomDeme(PtrTP->NP,picks);
	ASA_ptr.reset(new ASAalgorithm(PtrTP->generations,PtrTP->Ts,PtrTP->Tf));
   }

   oldfitness = deme.extractBestIndividual().getFitness();

   start=clock();
   deme = ASA_ptr->evolve(deme);
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


