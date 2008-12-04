#include "GOradiator.h"
#include "rad_objfun.h"
#include <math.h>
#include <vector>

radiatorProb::radiatorProb(int tmp){ //dim needs to be an even number (IMPORTANT)
		int number_layers=C_NUMBER_LAYERS;
		int number_active_layers=C_NUMBER_ACTIVE_LAYERS;
		int dim = (number_layers-number_active_layers)*2+number_active_layers;
		setDimension(dim);
		//radiatorProb bounds
		std::vector <double> lb,ub;
		/*for ( int i=0; i<getDimension()/2; i++ ) {
				lb.push_back(0.0);
				ub.push_back(15.0);
		}
		for ( int i=getDimension()/2; i<getDimension(); i++ ) {
				lb.push_back(1e-9);
				ub.push_back(1e-6);
		}*/
		for ( int i=0; i<(number_layers-number_active_layers); i++ ) {
				lb.push_back(1.5);
				ub.push_back(15.0);
		}
		for ( int i=(number_layers-number_active_layers); i<dim; i++ ) {
				lb.push_back(1e-9);
				ub.push_back(1e-6);
		}




		setBounds(lb,ub);
};
