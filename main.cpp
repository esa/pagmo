#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
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
#include "SolversThreads.h"
#include "rng.h"

using namespace std;

// Useful typedefs.
typedef boost::unique_lock<boost::mutex> lock_type;
typedef messengerfullProb problem_type;

int main(){

		int D =0;							//Problem dimension will be assigned later
		int choice=0;						//User choice
		cout.precision(9);

		//We prepare the pseudorandom sequence (TODO: check the randomnumbers of different threads are different)

		rng_uint32_type rng(time(0));
		rng_double_type drng(time(0));

		//we set the problem
		problem_type problem;
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

		cout << "Choose: 1-ASA, 2-PSO, 3-MPSO, 4-DE, 5-SGA, 6-IslandModel(ring), 7-DiGMO(simulator): ";
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
			//ASA.initASA(niterTot,niterTemp,niterRange,LB.size(),T0,Tcoeff,StartStep, rng());
			ASA.initASA(niterTot,LB.size(),T0,Tf, rng());

			//Pruned bounds

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, drng);
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
							ASA.initASA(niterTot,LB.size(),1,0.01, rng());
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
			PSO.initPSO(gen,LB.size(),omega,eta1,eta2,vcoeff, rng());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, drng);
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
			MPSO.initMPSO(gen,LB.size(),omega,eta1,eta2,vcoeff,nswarms, rng());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, drng);
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
			DE.initDE(gen,LB.size(),F,CR,strategy, rng());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population


				pop.createRandomPopulation(LB,UB,NP, drng);
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
			SGA.initSGA(gen,LB.size(),CR,M,insert_best, rng());

			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//we create a random population
				pop.createRandomPopulation(LB,UB,NP, drng);
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
			int islandsN  = 8;				//Number of Islands

			//stopping criteria
			int itermax = 120;				//Maximum number of iterations allowed (i.e. output printed on the screen)

			//Experiment Outputs
			double mean = 0, dev= 0, max = 0, min = 200000000;
			vector <double> results;

			//Double rng.
			rng_double_type drng(uint32_t(time(0)));

			//Main cycle creating threads
			for (int i=0;i<trials;i++){
				cout << "\nTrial number #" << i+1 << endl;
				//We create a random population in each island (this is done by the main thread - i.e. no parallelization here)
				vector<Population> IslandPop(islandsN);
				vector<GOProblem*> parallelProblems(islandsN);
				for (int i=0;i<islandsN;i++){
					parallelProblems[i] = new problem_type();

					IslandPop[i].createRandomPopulation(LB,UB,NP, drng);
					IslandPop[i].evaluatePopulation(*parallelProblems[i]);
				}

				//We instanciate the objects needed for threads
				boost::mutex mutex;
				boost::condition_variable exitcond;

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
				int IslandType = 0;
				double loweststd = 1000000;

                {
                    lock_type lock(mutex);
                    time(&start);
                    while(iter < itermax){
                        for (int i=0;i<islandsN;i++){
                            if ( !*(data[i].isActive) ){
                                //Create again the i-th thread to simulate an Island
                                IslandType = 3;
                                if (IslandType == 0){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(PSOthread,(void *)&data[i]).detach();
                                }
                                else if (IslandType == 1){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tMPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff  << "\t\tNswamrs: " << data[i].nswarms<< "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(MPSOthread, (void *)&data[i]).detach();
                                }
                                else if (IslandType == 2){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tSGA:\t\t CR: "<< data[i].CRsga  <<  "\t\tM: " << data[i].M <<  "\t\tInsertBest: " << data[i].insert_best << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(SGAthread,(void *)&data[i]).detach();
                                }
                                else if (IslandType == 3){
                                    data[i].randomSeed = rng();
                                    data[i].generations=10000;
                                    data[i].NP = 1;
                                    cout << "\t\t\tASA:\t\t Ts: "<< data[i].Ts  <<  "\t\tTf: " << data[i].Tf << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(ASAthread,(void *)&data[i]).detach();
                                }
                                else {
                                    data[i].randomSeed = rng();
                                    data[i].strategy = strategy;
                                    cout << "\t\t\tDE: \t\t F: "<< data[i].F  <<  "\t\tCR: " << data[i].CR << "\t\tStrategy: " << data[i].strategy << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(DEthread,(void *)&data[i]).detach();
                                }
                                //Thread Successfully Created
                                iter += 1;
                                *(data[i].isActive) = true;

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
                                if (drng() < 0.2){
                                    IslandPop[(i+1) % islandsN].substituteIndividual(IslandPop[i].extractBestIndividual(), rng() % data[i].NP);
                                }

                            }
                        }
                        exitcond.wait(lock);
                    }
                }
				//The main cycle has finished: we wait for all threads to finish
				for (int i=0; i<islandsN;i++){
				    //infinite loop if a thread never ends
				    while (__sync_bool_compare_and_swap(data[i].isActive,1,1));
				}

				//deallocate memory
				delete[] data;
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

            case 7:
           {

			//Experiment Settings
			int NP = 200;					//main population size
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
				Population IslandPop;
				vector<GOProblem*> parallelProblems(islandsN);
				for (int i=0;i<islandsN;i++) { parallelProblems[i] = new problem_type(); } //It is necessary to create a different object per CPU as to make sure to access different memory
																							//locations when the problem is mga or mga-1dsm for th r,v,and DV variables in the mgaproblem structures
				IslandPop.createRandomPopulation(LB,UB,NP, drng);
				IslandPop.evaluatePopulation(*parallelProblems[0]);							//all the problems are identical.... evaluation is done for the [0] one.


				//We instanciate the objects needed for threads
				boost::mutex mutex;
				boost::condition_variable exitcond;

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
					data[i].Ptr_pop = &IslandPop;			//Only one population in DiGMO


					data[i].NP = 20;
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
				int IslandType = 0;
				double loweststd = 1000000;

                {
                    lock_type lock(mutex);
                    time(&start);
                    while(iter < itermax){
                        for (int i=0;i<islandsN;i++){
                            if ( !*(data[i].isActive) ){
                                //Create again the i-th thread to simulate an Island
                                IslandType = 4;
                                if (IslandType == 0){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(PSOthread,(void *)&data[i]).detach();
                                }
                                else if (IslandType == 1){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tMPSO:\t\t omega: "<< data[i].omega  <<  "\t\teta1: " << data[i].eta1 <<  "\t\teta2: " << data[i].eta2 << "\t\tVcoeff: " << data[i].vcoeff  << "\t\tNswamrs: " << data[i].nswarms<< "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(MPSOthread,(void *)&data[i]).detach();
                                }
                                else if (IslandType == 2){
                                    data[i].randomSeed = rng();
                                    cout << "\t\t\tSGA:\t\t CR: "<< data[i].CRsga  <<  "\t\tM: " << data[i].M <<  "\t\tInsertBest: " << data[i].insert_best << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(SGAthread,(void *)&data[i]).detach();
                                }
                                else if (IslandType == 3){
                                    data[i].randomSeed = rng();
                                    data[i].generations=10000;
                                    data[i].NP = 1;
                                    cout << "\t\t\tASA:\t\t Ts: "<< data[i].Ts  <<  "\t\tTf: " << data[i].Tf << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(ASAthread,(void *)&data[i]).detach();
                                }
                                else {
                                    data[i].randomSeed = rng();
                                    data[i].strategy = strategy;
                                    cout << "\t\t\tDE: \t\t F: "<< data[i].F  <<  "\t\tCR: " << data[i].CR << "\t\tStrategy: " << data[i].strategy << "\t\tGenerations: " << data[i].generations << endl;
                                    boost::thread(DEthread,(void *)&data[i]).detach();
                                }

                                iter += 1;
                                *(data[i].isActive) = true;

                                //evaluate standard deviation in main population
                                //loweststd = IslandPop[0].evaluateStd();
                                loweststd = IslandPop.evaluateStd();

                                //log islands best and std
                                for (int i2=0;i2<islandsN;i2++){
                                cout << "CPU#:" << i2 << ": " << IslandPop.extractBestIndividual().getFitness() << " " << IslandPop.evaluateStd() << endl;
                                }
                                cout << "HighestStd: " << loweststd << endl;
                                cout << "Iter: " << iter+1 << endl;
                                cout << "\nmain():" << endl <<  "\t\t\tcreating thread, ID " << i << endl;

                                //ring topology migration


                            }
                        }
                        exitcond.wait(lock);
                    }
                }
				//The main cycle has finished: we wait for all threads to finish
				for (int i=0; i<islandsN;i++){
					while (*data[i].isActive); //infinite loop if a thread never ends
				}

				//deallocate memory
				delete[] data;
				delete[] isActive;
				time(&end);
				dif = difftime(end,start);
				cout << "\nSeconds elapsed: " << dif << endl<<endl;

				double res=IslandPop.extractBestIndividual().getFitness();
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

		}//end case

} //end while





		cout << "\nBest fitness found: " <<  pop.extractBestIndividual().getFitness() << endl;

        return 0;


    }

