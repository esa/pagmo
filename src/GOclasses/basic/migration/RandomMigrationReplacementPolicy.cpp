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

#include <algorithm>

#include "RandomMigrationReplacementPolicy.h"
#include "../../../exceptions.h"

// 09/03/2009: Initial version by Marek Rucinski.

namespace pagmo
{

std::list<std::pair<int, int> > RandomMigrationReplacementPolicy::selectForReplacement(const std::vector<individual>& incomingpopulation, const population& destinationpopulation)
{
	int migrationRateLimit = std::min<int>(getMaxMigrationRate(destinationpopulation), (int)incomingpopulation.size());

	std::vector<int> incomingpopulationIndices(migrationRateLimit);
	std::vector<int> destinationpopulationIndices(destinationpopulation.size());

	//Fill in the arrays of indices
	std::generate(incomingpopulationIndices.begin(), incomingpopulationIndices.end(), IndexGenerator());
	std::generate(destinationpopulationIndices.begin(), destinationpopulationIndices.end(), IndexGenerator());

	//Permute the indices (incoming population)
	for (int i = 0; i < migrationRateLimit; i++) {
		int nextIncomingpopulationIndex = i + (rng() % (migrationRateLimit - i));
		if (nextIncomingpopulationIndex != i) {
			std::swap(incomingpopulationIndices[i], incomingpopulationIndices[nextIncomingpopulationIndex]);
		}
	}

	//Permute the indices (destination population)
	for (int i = 0; i < migrationRateLimit; i++) {
		int nextDestinationpopulationIndex = i + (rng() % (destinationpopulation.size() - i));
		if (nextDestinationpopulationIndex != i) {
			std::swap(destinationpopulationIndices[i], destinationpopulationIndices[nextDestinationpopulationIndex]);
		}
	}

	// Create the result
	std::list<std::pair<int, int> > result;

	for (int i = 0; i < migrationRateLimit; i++) {
		result.push_back(std::make_pair(destinationpopulationIndices[i], incomingpopulationIndices[i]));
	}

	//std::cout << "Result size: " << result.size() << std::endl;

	return result;
}

}
