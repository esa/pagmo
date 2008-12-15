#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>	//for time()
#include "population.h"
#include "ASA.h"
#include "PSO.h"
#include "DE.h"
#include "SGA.h"
#include "MPSO.h"
#include "LOCAL.h"
#include "TrajectoryProblems.h"
#include "ClassicProblems.h"
#include <pthread.h>
#include "SolversThreads.h"
#include "PkRandom.h"

using namespace std;

int main(){

		int D =0;							//Problem dimension will be assigned later
		int choice=0;						//User choice
		cout.precision(9);

		//We prepare the pseudorandom sequence (TODO: check the randomnumbers of different threads are different)

		Pk::Random32 rng(time(0));

		//we set the problem
		messengerfullProb problem;
		//we extract its information into local variables
		const vector<double>& LB = problem.getLB();
		const vector<double>& UB = problem.getUB();
		D = problem.getDimension();

		//we declare populations and individuals
		Population demeDE,demePSO,demeASA,demeLOCA,demeSGA,pop;
		Individual x;
		vector <int> picksDE,picksPSO,picksASA,picksLOCAL,picksSGA;

		time_t start,end,start1,end1;
		double dif;

		//We create and open the logfile
		ofstream logfile("log.txt");



while (choice != -1) {

		//we choose the algorithm

		cout << "Choose: 1-ASA, 2-PSO, 3-MPSO, 4-DE, 5-SGA, 6-IslandModel(ring), Default-DiGMO: ";
		cin >> choice;

		switch (choice){
			case 1:
			{
			//Experiment Settings
			int NP = 1;					    //population size
			int trials = 300;				//number of trials
			int niterTot = 10000;			//generations per algorithm call
			int niterRange = 20;
			int niterTemp = 1;
			double T0 = 10;
			double Tf = T0/1000;
			double Tcoeff = 0.85;
			double StartStep =1;

			//stopping criteria
			int itermax = 120;					//Maximum number of iterations allowed (i.e. output printed on the screen)

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Instanciate the algorithm
			//Adaptive Simulated Annealing
			ASAalgorithm ASA;
			//ASA.initASA(niterTot,niterTemp,niterRange,LB.size(),T0,Tcoeff,StartStep, rng.next());
			ASA.initASA(niterTot,LB.size(),T0,Tf, rng.next());

			//Pruned bounds

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, rng);
				pop.evaluatePopulation(problem);
				int iter = 0;

				time(&start);
					while( iter < itermax){
						iter ++;
						//we print the best
						cout << "Initial fitness: " << pop.extractBestIndividual().getFitness() << endl;
						//we evolve it
						start1=clock();
						if (pop.extractBestIndividual().getFitness() < 5){
							ASA.initASA(niterTot,LB.size(),1,0.01, rng.next());
						}
						pop = ASA.evolve(pop[0],problem);
						end1=clock();
						dif = (double)(end1-start1) / (double)CLOCKS_PER_SEC;
						//we print the result
						cout << "Final fitness: " << pop.extractBestIndividual().getFitness() << endl;
						cout << "Worst fitness: " << pop.extractWorstIndividual().getFitness() << endl;
						cout << "Mean         : " << pop.evaluateMean() << endl;
						cout << "Std          : " << pop.evaluateStd() << endl;
						cout << "\t\tSeconds elapsed: " << dif << endl;
					}
				time(&end);
				results.push_back(pop.extractBestIndividual().getFitness());
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;
				for (int i=0;i<D;i++){
				logfile << pop.extractBestIndividual()[i] << " ";
				}
				logfile << exp(-pop.extractBestIndividual().getFitness()) <<endl;
			}

			//evaluate experiment results
			for (int i=0;i<trials;i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= trials;

			for (int i=0;i<trials;i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/trials);

			for (int i=0;i<trials;i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;

			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			}
			//Lines to test ASA
			break;

			case 2: //PSO sequential experiment
				{
			//Experiment Settings
			int NP = 20;					//population size
			int trials = 100;				//number of trials
			int gen = 500;					//generations per algorithm call
			double omega = 0.6;
			double eta1 = 2;
			double eta2 = 2;
			int vcoeff = 1;

			//stopping criteria
			int itermax = 120;					//Maximum number of iterations allowed (i.e. output printed on the screen)
			double stdmin = 1e-5;				//Standard deviation of the population that determines a stop

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Instanciate the algorithm
			PSOalgorithm PSO;
			PSO.initPSO(gen,LB.size(),omega,eta1,eta2,vcoeff, rng.next());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, rng);
				pop.evaluatePopulation(problem);
				int iter = 0;

				time(&start);
					while(iter < itermax){
						iter ++;
						//we print the best
						cout << "Initial fitness: " << pop.extractBestIndividual().getFitness() << endl;
						//we evolve it
						start1=clock();
						pop = PSO.evolve(pop,problem);
						end1=clock();
						dif = (double)(end1-start1) / (double)CLOCKS_PER_SEC;
						//we print the result
						cout << "Final fitness: " << pop.extractBestIndividual().getFitness() << endl;
						cout << "Worst fitness: " << pop.extractWorstIndividual().getFitness() << endl;
						cout << "Mean         : " << pop.evaluateMean() << endl;
						cout << "Std          : " << pop.evaluateStd() << endl;
						cout << "\t\tSeconds elapsed: " << dif << endl;
					}
				time(&end);
				results.push_back(pop.extractBestIndividual().getFitness());
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;
			}

			//evaluate experiment results
			for (int i=0;i<trials;i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= trials;

			for (int i=0;i<trials;i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/trials);

			for (int i=0;i<trials;i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;

			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			}
			break;

			case 3: //MPSO sequential experiment
				{
			//Experiment Settings
			int NP = 20;					//population size
			int trials = 100;				//number of trials
			int gen = 500;					//generations per algorithm call
			double omega = 0.65;
			double eta1 = 2.0;
			double eta2 = 2.0;
			int vcoeff = 1;
			int nswarms = 4;

			//stopping criteria
			int itermax = 100;					//Maximum number of iterations allowed (i.e. output printed on the screen)
			double stdmin = 1e-5;				//Standard deviation of the population that determines a stop

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Instanciate the algorithm
			MPSOalgorithm MPSO;
			MPSO.initMPSO(gen,LB.size(),omega,eta1,eta2,vcoeff,nswarms, rng.next());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, rng);
				pop.evaluatePopulation(problem);
				int iter = 0;

				time(&start);
					while(iter < itermax){
						iter ++;
						//we print the best
						cout << "Initial fitness: " << pop.extractBestIndividual().getFitness() << endl;
						//we evolve it
						start1=clock();
						pop = MPSO.evolve(pop,problem);
						end1=clock();
						dif = (double)(end1-start1) / (double)CLOCKS_PER_SEC;
						//we print the result
						cout << "Final fitness: " << pop.extractBestIndividual().getFitness() << endl;
						cout << "Worst fitness: " << pop.extractWorstIndividual().getFitness() << endl;
						cout << "Mean         : " << pop.evaluateMean() << endl;
						cout << "Std          : " << pop.evaluateStd() << endl;
						cout << "\t\tSeconds elapsed: " << dif << endl;
					}
				time(&end);
				results.push_back(pop.extractBestIndividual().getFitness());
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;
			}

			//evaluate experiment results
			for (int i=0;i<trials;i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= trials;

			for (int i=0;i<trials;i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/trials);

			for (int i=0;i<trials;i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;

			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			}
			break;



			case 4: //Sequential Differential Evolution Experiment
			{
			//Experiment Settings
			int NP = 20;					//population size for each island
			int trials = 100;				//number of trials
			int gen = 500;					//generations per algorithm call
			double F = 0.8;					//F in DE
			double CR = 0.8;				//CR in DE
			int strategy = 2;				//DE startegy

			//stopping criteria
			int itermax = 120;				//Maximum number of iterations allowed (i.e. output printed on the screen)

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Instanciate the algorithm
			DEalgorithm DE;
			DE.initDE(gen,LB.size(),F,CR,strategy, rng.next());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population


				pop.createRandomPopulation(LB,UB,NP, rng);
				pop.evaluatePopulation(problem);
				int iter = 0;

				time(&start);
					while( iter < itermax){
						iter ++;
						//we print the best
						cout << "Initial fitness: " << pop.extractBestIndividual().getFitness() << endl;
						//we evolve it
						start1=clock();
						pop = DE.evolve(pop,problem);
						end1=clock();
						dif = (double)(end1-start1) / (double)CLOCKS_PER_SEC;
						//we print the result
						cout << "Final fitness: " << pop.extractBestIndividual().getFitness() << endl;
						cout << "Worst fitness: " << pop.extractWorstIndividual().getFitness() << endl;
						cout << "Mean         : " << pop.evaluateMean() << endl;
						cout << "Std          : " << pop.evaluateStd() << endl;
						cout << "\t\tSeconds elapsed: " << dif << endl;
					}
				time(&end);
				//if (iter<itermax){
				//results.push_back(iter*NP*gen);
				results.push_back(pop.extractBestIndividual().getFitness());
				//}
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;
				for (int i=0;i<D;i++){
				logfile << pop.extractBestIndividual()[i] << " ";
				}
				logfile << exp(-pop.extractBestIndividual().getFitness()) <<endl;
			}

			//evaluate experiment results
			for (unsigned int i=0;i<results.size();i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= results.size();

			for (unsigned int i=0;i<results.size();i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/results.size());

			for (unsigned int i=0;i<results.size();i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;
			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			cout << "Success: " << results.size() << endl;
			}
			break;



			case 5: //SGA experiment
			{
			//Experiment Settings
			int NP = 20;					//population size
			int trials = 100;				//number of trials
			int gen = 500;					//generations per algorithm call
			double CR = 0.7;
			double M = 0.2;
			int insert_best = 1;

			//stopping criteria
			int itermax = 120;					//Maximum number of iterations allowed (i.e. output printed on the screen)

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Instanciate the algorithm
			//Simple Genetic Algorithm
			SGAalgorithm SGA;
			SGA.initSGA(gen,LB.size(),CR,M,insert_best, rng.next());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, rng);
				pop.evaluatePopulation(problem);
				int iter = 0;

				time(&start);
					while( iter < itermax){
						iter ++;
						//we print the best
						cout << "Initial fitness: " << pop.extractBestIndividual().getFitness() << endl;
						//we evolve it
						start1=clock();
						pop = SGA.evolve(pop,problem);
						end1=clock();
						dif = (double)(end1-start1) / (double)CLOCKS_PER_SEC;
						//we print the result
						cout << "Final fitness: " << pop.extractBestIndividual().getFitness() << endl;
						cout << "Worst fitness: " << pop.extractWorstIndividual().getFitness() << endl;
						cout << "Mean         : " << pop.evaluateMean() << endl;
						cout << "Std          : " << pop.evaluateStd() << endl;
						cout << "\t\tSeconds elapsed: " << dif << endl;
					}
				time(&end);
				results.push_back(pop.extractBestIndividual().getFitness());
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;
			}

			//evaluate experiment results
			for (int i=0;i<trials;i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= trials;

			for (int i=0;i<trials;i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/trials);

			for (int i=0;i<trials;i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;

			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			}
			//Lines to test SGA
			break;

			case 6: //Parallel asynchronous island model
			{

			//Experiment Settings
			int NP = 1;						//population size for each island
			int trials = 20;				//number of trials
			int gen = 500;					//generations per algorithm call
			double F = 0.8;					//F in DE
			double CR = 0.8;				//CR in DE
			int strategy = 2;				//DE startegy
			int islandsN  = 1;				//Number of Islands
			
			//stopping criteria
			int itermax = 120;				//Maximum number of iterations allowed (i.e. output printed on the screen)

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Main cycle creating threads
			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//We create a random population in each island (this is done by the main thread - i.e. no parallelization here)
				vector<Population> IslandPop(islandsN);
				vector<GOProblem*> parallelProblems(islandsN);
				for (int i=0;i<islandsN;i++){
					parallelProblems[i] = new tandemProb();
					
					IslandPop[i].createRandomPopulation(LB,UB,NP, rng);
					IslandPop[i].evaluatePopulation(*parallelProblems[i]);
				}

				//We instanciate the objects needed for pthreads allocating memory for the threads array
				pthread_t *threads;
				threads = new pthread_t [islandsN];
				pthread_mutex_t mutex;
				pthread_cond_t exitcond;
				pthread_mutex_init(&mutex, NULL);
				pthread_cond_init (&exitcond, NULL);

				//We allocate memory for the isActive array containing the status ofeach thread	and we initialise it to false (no threads opened)
				bool *isActive;
				isActive = new bool[islandsN];
				for (int i=0;i<islandsN;i++) isActive[i]=false;  //all threads are inactive

				//We allocate memory and initialise the data threadParam array containing the information to be passed to the different threads
				threadParam* data;
				data = new threadParam [islandsN];
				for (int i=0;i<islandsN;i++){
					data[i].threadID = i;

					data[i].problem = parallelProblems[i];
					data[i].Ptr_pop = &IslandPop[i];	//Each Island has a different population


					data[i].NP = NP;
					data[i].generations = gen;

					//Initialise DE specific data
					data[i].strategy = 1;
					data[i].F = F;
					data[i].CR = CR;

					//Initialise PSO/MPSO specific data
					data[i].omega = 0.65;
					data[i].eta1 = 2.0;
					data[i].eta2 = 2.0;
					data[i].vcoeff = 1;
					data[i].nswarms = 4;

					//Initialise SGA specific data
					data[i].CRsga = 0.95;
					data[i].M = 0.2;
					data[i].insert_best =1;

					//Initialise ASA specific data
					data[i].Ts = 1;
					data[i].Tf = data[i].Ts/1000;
					data[i].isActive = &isActive[i];
					data[i].TPmutex = &mutex;
					data[i].exit = &exitcond;

					data[i].Ptr_log = &logfile;
				}

				int iter=0;
				int rc;
				int IslandType = 0;
				double loweststd = 1000000;

				pthread_mutex_lock (&mutex);
				time(&start);
				while(iter < itermax){
					for (int i=0;i<islandsN;i++){
						if ( !*(data[i].isActive) ){
							//Create again the i-th thread to simulate an Island
							IslandType = 3;
							if (IslandType == 0){
								data[i].randomSeed = rng.next();
								cout << "\t\t\tPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, PSOthread, (void *)&data[i]);
							}
							else if (IslandType == 1){
								data[i].randomSeed = rng.next();
								cout << "\t\t\tMPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff  << "\t\tNswamrs: " << data[i].nswarms<< "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, MPSOthread, (void *)&data[i]);
							}
							else if (IslandType == 2){
								data[i].randomSeed = rng.next();
								cout << "\t\t\tSGA:\t\t CR: "<< data[i].CRsga  <<  "\t\tM: " << data[i].M <<  "\t\tInsertBest: " << data[i].insert_best << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, SGAthread, (void *)&data[i]);
							}
							else if (IslandType == 3){
								data[i].randomSeed = rng.next();
								data[i].generations=10000;
								data[i].NP = 1;
								cout << "\t\t\tASA:\t\t Ts: "<< data[i].Ts  <<  "\t\tTf: " << data[i].Tf << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, ASAthread, (void *)&data[i]);
							}
							else {
								data[i].randomSeed = rng.next();
								data[i].strategy = strategy;
								cout << "\t\t\tDE: \t\t F: "<< data[i].F  <<  "\t\tCR: " << data[i].CR << "\t\tStrategy: " << data[i].strategy << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, DEthread, (void *)&data[i]);
							}
							if (rc){
							//Problem creating the thread
								printf("ERROR; return code from pthread_create() is %d\n", rc);
								exit(-1);
							}
							else{
							//Thread Successfully Created
								iter += 1;
								*(data[i].isActive) = true;
								//thread is detached as it will never be joined
								//(on OSx Leopard this avoids that at a certain point
								//threads cannot be created anymore resulting in rc=35)
								pthread_detach(threads[i]);
							}

							//evaluate highest standard deviation across Islands
							//loweststd = IslandPop[0].evaluateStd();
							loweststd = IslandPop[0].extractBestIndividual().getFitness();
							for (int i2=1;i2<islandsN;i2++){
								if (loweststd < IslandPop[i2].extractBestIndividual().getFitness()) loweststd = IslandPop[i2].extractBestIndividual().getFitness();
							}

							//log islands best and std
							for (int i2=0;i2<islandsN;i2++){
							cout << "CPU#:" << i2 << ": " << IslandPop[i2].extractBestIndividual().getFitness() << " " << IslandPop[i2].evaluateStd() << endl;
							}
							cout << "HighestStd: " << loweststd << endl;
							cout << "Iter: " << iter+1 << endl;
							cout << "\nmain():" << endl <<  "\t\t\tcreating thread, ID " << i << endl;

							//ring topology migration
							if (Pk::nextDouble(rng) < 0.2){
								IslandPop[(i+1) % islandsN].substituteIndividual(IslandPop[i].extractBestIndividual(), rng.next() % data[i].NP);
							}

						}
					}
					pthread_cond_wait(&exitcond, &mutex);
				}
				pthread_mutex_unlock (&mutex);
				//The main cycle has finished: we wait for all threads to finish
				for (int i=0; i<islandsN;i++){
					while (*data[i].isActive); //infinite loop if a thread never ends
				}

				//we clean the thread objects and deallocate memory
				pthread_mutex_destroy(&mutex);
				pthread_cond_destroy(&exitcond);
				delete[] data;
				delete[] threads;
				delete[] isActive;
				time(&end);
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;

				double res=IslandPop[0].extractBestIndividual().getFitness();
				for (int i=1;i<islandsN;i++){
					if (res > IslandPop[i].extractBestIndividual().getFitness()) res = IslandPop[i].extractBestIndividual().getFitness();
				}
				//if (iter<itermax){
				results.push_back(res);
				//}
				
				for (int i = 0; i < islandsN; i++) {
					delete parallelProblems[i];
				}
			}

			//evaluate experiment results
			for (unsigned int i=0;i<results.size();i++){
				mean += results[i];
				if (results[i] > max) max = results[i];
				if (results[i] < min) min = results[i];
			}
			mean /= results.size();

			for (unsigned int i=0;i<results.size();i++){
				dev += (mean-results[i])*(mean-results[i]);
			}
			dev = sqrt(dev/results.size());

			for (unsigned int i=0;i<results.size();i++){
				cout << "\nTrial #" << i << ": " << results[i] << endl;
			}

			//print results
			cout << "\nMean: " << mean << endl;
			cout << "Std: " << dev << endl;
			cout << "Max: " << max << endl;
			cout << "Min: " << min << endl;
			cout << "Success: " << results.size() << endl;
			}
			break;

            case 9:
            {
                Individual x;
                x.createRandomIndividual(LB,UB,rng);
                cout << x << endl;
                x[0]=666;
                cout << x << endl;

                Population pop;
                pop.createRandomPopulation(LB,UB,5,rng);
                pop.evaluatePopulation(problem);
                cout << pop[3] << endl;
                pop[3][0]=666;
                cout << pop[3] << endl;
                cout << pop << endl;

            }
            break;

			default:
			{
			int NUM_THREADS=2;
			//We create a random population and evaluate its fitness (this is done by the main thread - i.e. no parallelization here)
			pop.createRandomPopulation(LB,UB,20*5, rng);
			pop.evaluatePopulation(problem);

			//We instanciate the objects needed for pthreads allocating memory for the threads array
			pthread_t *threads;
			threads = new pthread_t [NUM_THREADS];
			//pthread_attr_t attr;
			//pthread_attr_init(&attr);
			//pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			pthread_mutex_t mutex;
			pthread_cond_t exitcond;
			pthread_mutex_init(&mutex, NULL);
			pthread_cond_init (&exitcond, NULL);




			//We allocate memory for the isActive array containing the status ofeach thread	and we initialise it to false (no threads opened)
			bool *isActive;
			isActive = new bool[NUM_THREADS];		//all threads are inactive
			for (int i=0;i<NUM_THREADS;i++) {isActive[i]=false;}

			//We allocate memory and initialise the data threadParam array containing the information to be passed to the different threads
			threadParam *data;
			data = new threadParam [NUM_THREADS];
			for (int i=0;i<NUM_THREADS;i++){
				data[i].threadID = i;

				data[i].problem = &problem;
				data[i].Ptr_pop = &pop;

				data[i].generations = 500;

				//Initialise DE specific data
				data[i].strategy = 2;
				data[i].F = 0.8;
				data[i].CR = 0.8;

				//Initialise PSO specific data
				data[i].omega = 0.5;
				data[i].eta1 = 2.0;
				data[i].eta2 = 2.0;
				data[i].vcoeff = 1000;

				data[i].isActive = &isActive[i];
				data[i].TPmutex = &mutex;
				data[i].exit = &exitcond;

				data[i].Ptr_log = &logfile;
			}

			//Main cycle creating threads
			int globalCount=0;
			int count1=0,count2=0;
			int rc;
			int algorithmchoice = 0;

			pthread_mutex_lock (&mutex);
			time(&start);
			while (globalCount<120){
			//while(pop.extractBestIndividual().getFitness()!=0){
					for (int i=0;i<NUM_THREADS;i++){
					if ( !*(data[i].isActive) ){
							cout << "main():" << endl <<  "\t\t\tcreating thread, ID " << i << endl;
							algorithmchoice = 4;//r250()%10;
							if (algorithmchoice == 0){
								count1++;
								cout << "\t\t\tPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, PSOthread, (void *)&data[i]);
							}
							else {
								data[i].strategy = i+1;
								cout << "\t\t\tDE: \t\t F: "<< data[i].F  <<  "\t\tCR: " << data[i].CR << "\t\tStrategy: " << data[i].strategy << "\t\tGenerations: " << data[i].generations << endl;
								rc = pthread_create(&threads[i], NULL, DEthread, (void *)&data[i]);
								count2++;
							}

							if (rc){
								printf("ERROR; return code from pthread_create() is %d\n", rc);
								exit(-1);
							}
							else{
								globalCount += 1;
								*(data[i].isActive) = true;
							}
					}
				}
				pthread_cond_wait(&exitcond, &mutex);
			}
			pthread_mutex_unlock (&mutex);

			//The main cycle has finished: we wait for all threads to finish
			for (int i=0; i<NUM_THREADS;i++){
				pthread_join(threads[i],NULL);
			}

			//we clean the thread objects and deallocate memory
			pthread_mutex_destroy(&mutex);
			//pthread_attr_destroy(&attr);
			pthread_cond_destroy(&exitcond);
			delete[] data;
			delete[] threads;
			delete[] isActive;
			time(&end);
			dif = difftime(end,start);
			cout << "\nSeconds elapsed: " << dif << endl<<endl;
			cout << "\nCount1: " <<  count1 << endl;
		    cout << "\nCount2: " <<  count2 << endl;

			}
			break;
		}//end case

} //end while


			//The following lines are commented out as they implement exactly the DiGMO random strategy in a random
			//way. This creates a very unefficient cooperation between the algorithms as it does not have the implicit
			//learning of the best strategy. This is the mechanism that allow, in the distributed version of the algorithm,
			//the same individual to be evolved by different algorithms and only the most succesful one to be selected
			//
			//It seems that the reinsertion strategy (only improved individuals substitute their old versions) together with
			//the contemporary evolution of the same individuals is key to the efficiency of algorithm cooperation (mimicking co-evolution)
			//In a sequential version this can be  enforced by evolving with DE, PSO and SA and only at the end perform reinsertion
			//as in the uncommented lines above

			/*algorithmchoice = r250()%3;
			if ( algorithmchoice == 0 ){
				deme=pop.extractRandomDeme(20,picks);
				//we print the best
				cout << "Using DE" << endl;
				cout << "Initial fitness: " << deme.extractBestIndividual().getFitness() << endl;
				//DE.initDE(50,LB.size(),dr250(),dr250(),(int)dr250()*11+1);
				deme = DE.evolve(deme,&rosenbrock,LB,UB);
			}

			else if ( algorithmchoice == 1 ){
				deme=pop.extractRandomDeme(20,picks);
				//deme.resetVelocities(LB,UB);
			    //we print the best
				cout << "Using PSO" << endl;
				cout << "Initial fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
				//PSO.initPSO(500,LB.size(),(dr250()+1)/2,4.0*dr250(),4.0*dr250(),(int)dr250()*10000);
				deme = PSO.evolve(deme,&rosenbrock,LB,UB);
			}

			else if ( algorithmchoice == 2 ){
				deme=pop.extractRandomDeme(1,picks);
				//deme.resetVelocities(LB,UB);
				//we print the best
				cout << "Using ASA" << endl;
				cout << "Initial fitness: " << deme.extractBestIndividual().getFitness() << endl;
				//ASA.initASA(10000,1,8,LB.size(),dr250(),dr250(),1.0);
				deme = ASA.evolve(deme[0],&rosenbrock,LB,UB);
			}

			//we print the result
			cout << "Final fitness: " <<  deme.extractBestIndividual().getFitness() << endl;
			//we insert it in the main population
			pop.insertDeme(deme,picks);
			picks.clear();*/


		cout << "\nBest fitness found: " <<  pop.extractBestIndividual().getFitness() << endl;

        return 0;


    }

