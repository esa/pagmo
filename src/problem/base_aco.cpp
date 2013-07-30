/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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
	
	/**
	 * @param[in] n integer dimension of the problem.
	 * @param[in] nc global constraints dimension
	 * @param[in] nic inequality constraints dimension
	 */
	base_aco::base_aco(int n, int nc, int nic):base(n,n,1,nc,nic,0.){
	}

//the default behaviour is to set to 1 all the values corresponding to values inside the bounds and 0 elsewhere
void base_aco::set_heuristic_information_matrix() {
	//allocates the memory for eta.
	create_heuristic_information_matrix();  //that MUST be called at the begining of each set_heuristic_information_matrix implementation!

	for(std::vector<std::vector<std::vector<fitness_vector> > >::size_type k = 0; k < m_eta.size(); ++k) {
		for(std::vector<std::vector<fitness_vector> >::size_type i=0; i < m_eta[0].size(); ++i) {
			for(std::vector<fitness_vector>::size_type  j = 0; j < m_eta[0][0].size(); ++j) {
				if(i <= get_ub()[k] - get_lb()[k] && j <= get_ub()[k+1] - get_lb()[k+1]) {
					m_eta[k][i][j][0] = 1;
				}
				else {
					//CR - This is just not a number as there is no edge in the graph!!!!
					m_eta[k][i][j][0] = std::numeric_limits<double>::quiet_NaN();
				}
			}
		}
	}
}

const std::vector<std::vector<std::vector<fitness_vector> > > &base_aco::get_heuristic_information_matrix() const {
	return m_eta;
}

void base_aco::create_heuristic_information_matrix() {
	double max_size = 0;
	for(problem::base::size_type i = 0; i < get_i_dimension(); ++i) {
		if(max_size < (get_ub()[i] - get_lb()[i])) {
			max_size = get_ub()[i] - get_lb()[i];
		}
	}
	const std::vector<decision_vector>::size_type nComponents = boost::numeric_cast<std::vector<decision_vector>::size_type>(max_size) + 1; 

	fitness_vector tempA(get_f_dimension(),0);	//used for initialisation purpouses
	std::vector<fitness_vector> tempB(nComponents,tempA); //used for initialisation purpouses
	std::vector<std::vector<fitness_vector> > tempC(nComponents,tempB); //used for initialisation purpouses
	std::vector<std::vector<std::vector<fitness_vector> > > eta(get_i_dimension(), tempC); //heuristic information matrix 
	m_eta = eta;
}



bool base_aco::check_partial_feasibility(const decision_vector &x) const{
	(void)x; //to avoid the  unused parameter ‘x’ warning by compiler
	return true;
}

}} //namespaces
