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

#include "base_aco.h"

namespace pagmo { namespace problem {

base_aco::base_aco(int n, int ni, int nf, int nc, int nic, const double &c_tol):base(n,ni,nf,nc,nic,c_tol){};

//the default behaviour is to set to 1 all the values corresponding to values inside the bounds and 0 elsewhere
void base_aco::get_heuristic_information_matrix(std::vector<std::vector<std::vector<fitness_vector> > > &eta) {
	for(std::vector<std::vector<std::vector<fitness_vector> > >::size_type k = 0; k < eta.size(); ++k) {
		for(std::vector<std::vector<fitness_vector> >::size_type i=0; eta[0].size(); ++i) {
			for(std::vector<fitness_vector>::size_type  j = 0; j < eta[0][0].size(); ++j) {
				if(i<get_ub()[k] && j < get_ub()[k+1]) {
					eta[k][i][j][0] = 1;
				}
				else {
					eta[k][i][j][0] = 0;
				}
			}
		}
	}
}


bool base_aco::check_partial_feasibility(decision_vector x){
	x=x; //to avoid the  unused parameter ‘x’ warning by compiler
	return true;
}

}} //namespaces
