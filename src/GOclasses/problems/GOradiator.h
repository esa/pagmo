#ifndef GORADIATOR_H_INCLUDED
#define GORADIATOR_H_INCLUDED

#include <vector>
#include "rad_objfun.h"
#include "GOproblem.h"

//***********************************************************************************
//Nanostructured Radiator Problems....
//***********************************************************************************

class radiatorProb : public GOProblem {
public:
	radiatorProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return radiator(x); }
};	//end class radiatorProb

#endif // GORADIATOR_H_INCLUDED
