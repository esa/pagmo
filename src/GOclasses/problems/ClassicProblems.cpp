/*
 *  GOProblem.cpp
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "ClassicProblems.h"
#include "GOproblem.h"
#include <math.h>
#include <vector>

//***********************************************************************************
//Classical problems
//***********************************************************************************

TestProb::TestProb(int dim):GOProblem(dim, std::vector<double>(dim, 0.0), std::vector<double>(dim, 1.0)) {
	//nothing else to be done
};

rastriginProb::rastriginProb(int dim):GOProblem(dim, std::vector<double>(dim, -5.12), std::vector<double>(dim, 5.12)) {
	//nothing else to be done
};

schwefelProb::schwefelProb(int dim):GOProblem(dim, std::vector<double>(dim, -500.0), std::vector<double>(dim, 500.0)) {
	//nothing else to be done
};

ackleyProb::ackleyProb(int dim):GOProblem(dim, std::vector<double>(dim, -15.0), std::vector<double>(dim, 30.0)) {
	//nothing else to be done
};

rosenbrockProb::rosenbrockProb(int dim):GOProblem(dim, std::vector<double>(dim, -5.0), std::vector<double>(dim, 10.0)) {
	//nothing else to be done
};

lennardjonesProb::lennardjonesProb(int atoms) {
		setDimension(3*atoms-6);
		//lennard jones bounds
		getLB().resize(getDimension());
		getUB().resize(getDimension());
		for (int i = 0; i < getDimension(); i++) {
			if ( (i != 0) && (i % 3) == 0 ) {
				getLB()[i] = 0.0;
				getUB()[i] = 6.0;
			}
			else {
				getLB()[i] = -3.0;
				getUB()[i] = 3.0;
			}
		}
};

levyProb::levyProb(int dim):GOProblem(dim, std::vector<double>(dim, -10.0), std::vector<double>(dim, 10.0)) {
	//nothing else to do
};
