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

#include "ChooseBestMigrationSelectionPolicy.h"
#include <vector>
#include <algorithm>

// 09/03/2009: Initial version by Marek Rucinski.

namespace pagmo
{

std::vector<individual> ChooseBestMigrationSelectionPolicy::selectForMigration(const population& population)
{
	int migrationRate = getNumberOfIndividualsToMigrate(population);

	//Create a temporary array of individuals
	std::vector<individual> result(population.begin(), population.end());

	/*std::cout << "Before sorting:" << std::endl;
	for(std::vector<individual>::const_iterator it = result.begin(); it != result.end(); ++it) {
		std::cout << it->get_fitness() << " ";
	}
	std::cout << std::endl;*/

	//Sort the individuals (best go first)
	std::sort(result.begin(), result.end(), individual::compare_by_fitness);

	/*std::cout << "After sorting:" << std::endl;
	for(std::vector<individual>::const_iterator it = result.begin(); it != result.end(); ++it) {
		std::cout << it->get_fitness() << " ";
	}
	std::cout << std::endl;*/

	//Leave only desired number of elements in the result
	result.erase(result.begin() + migrationRate, result.end());

	/*std::cout << "After erease:" << std::endl;
	for(std::vector<individual>::const_iterator it = result.begin(); it != result.end(); ++it) {
		std::cout << it->getFitness() << " ";
	}
	std::cout << std::endl;*/

	return result;
}

}
