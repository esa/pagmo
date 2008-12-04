/*
 *  LOCAL.cpp
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "LOCAL.h"
#include <vector>


void SIMPLELOCALalgorithm::initSIMPLELOCAL(double rangeInit, double reduxCoeffInit, double minrangeInit){
	range = rangeInit;
	reduxCoeff=reduxCoeffInit;
	minrange = minrangeInit;
}

Individual SIMPLELOCALalgorithm::evolve(Individual x0, double(*objfun)(const std::vector<double>&), std::vector<double> LB, std::vector<double> UB ){
	
	std::vector <double> x,newx;
	x = x0.getDecisionVector();
	
	double f = x0.getFitness();
	double newf;
	unsigned int D = x.size();
	bool flag = false;
	
	double newrange=range;
	
	while (newrange > minrange){
		flag = false;
		for (unsigned int i=0; i<D; i++){
			newx=x;
			
			
			newx[i] = x[i] + newrange * (UB[i]-LB[i]);
			//feasibility correction
			if (newx[i] > UB [i]) newx[i]=UB[i];
			
			newf = objfun(newx);
			if (newf < f) {
				f = newf;
				x = newx;
				flag=true;
				break; //accept
			}
		
			newx[i] = x[i] - newrange * (UB[i]-LB[i]);
			//feasibility correction
			if (newx[i] < LB [i]) newx[i]=LB[i];
			
			newf = objfun(newx);
			if (newf < f) {  //accept
				f = newf;
				x = newx;
				flag=true;
				break;
			}
		}
		if (!flag){
			newrange *= reduxCoeff;
		}
	} //end while
	
	Individual ret;
	ret.setFitness(f);
	ret.setDecisionVector(x);
	ret.setVelocity(x0.getVelocity());
	return ret;
}
	

