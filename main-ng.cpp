#include "island.h"
#include "CS.h"
#include "ClassicProblems.h"
#include "TrajectoryProblems.h"
#include <iostream>


using namespace std;


int main(){
        CSalgorithm algo(0.001);
        messengerfullProb prob;
        island isl(prob,algo,20);
        cout << "Best: " << isl.best().getFitness() << endl;
        isl.evolve();
        isl.join();
        cout << "Best: " << isl.best().getFitness() << endl;
        return 0;
}
