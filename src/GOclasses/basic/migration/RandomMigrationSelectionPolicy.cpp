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

#include "RandomMigrationSelectionPolicy.h"
#include <vector>
#include <algorithm>

// 09/03/2009: Initial version by Marek Rucinski.

namespace pagmo
{

std::vector<individual> RandomMigrationSelectionPolicy::selectForMigration(const population& population)
{
	int migrationRate = getNumberOfIndividualsToMigrate(population);

	//Create a temporary array of individuals
	std::vector<individual> result(population.begin(), population.end());

	//Permute the temporary array randomly (only first [migrationRate] elements)
	for (int i = 0; i < migrationRate; i++) {
		//Draw the next individual index
		int nextIndividualIndex = i + (rng() % (population.size() - i));

		//Swap the individual at the current position with the drawn one
		if (nextIndividualIndex != i) {
			std::swap(result[i], result[nextIndividualIndex]);
		}
	}

	//Remove remaining elements from the array
	result.erase(result.begin() + migrationRate, result.end());

	return result;
}

}
