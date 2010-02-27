/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#include "GOradiator.h"
#include "rad_objfun.h"
#include <math.h>
#include <vector>

radiatorProb::radiatorProb(int tmp)  //dim needs to be an even number (IMPORTANT)
{
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
