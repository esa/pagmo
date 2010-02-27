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

#include "BestReplaceWorstIfBetterMigrationReplacementPolicy.h"
#include "../../../exceptions.h"
#include <algorithm>
#include <iostream>

namespace pagmo
{

std::ostream& print_vector_of_individuals_indirect(std::ostream& str, const std::vector<individual>& vect_indiv, const std::vector<int>& vect_index)
{
	for (std::vector<int>::const_iterator it = vect_index.begin(); it != vect_index.end(); ++it) {
		str << vect_indiv[*it].get_fitness() << " ";
	}
	return str;
}

// 09/03/2009: Initial version by Marek Rucinski.

std::list<std::pair<int, int> > BestReplaceWorstIfBetterMigrationReplacementPolicy::selectForReplacement(const std::vector<individual>& incomingpopulation, const population& destinationpopulation)
{
	int migrationRateLimit = std::min<int>(getMaxMigrationRate(destinationpopulation), (int)incomingpopulation.size());

	// Temporary vectors to store sorted poopulations
	std::vector<int> incomingpopulationIndices(incomingpopulation.size());
	std::vector<int> destinationpopulationIndices(destinationpopulation.size());

	//Fill in the arrays of indices
	std::generate(incomingpopulationIndices.begin(), incomingpopulationIndices.end(), IndexGenerator());
	std::generate(destinationpopulationIndices.begin(), destinationpopulationIndices.end(), IndexGenerator());

	/*
	std::cout << "Initial individuals orders:" << std::endl;
	print_vector_of_individuals_indirect(std::cout, incomingpopulation, incomingpopulationIndices) << std::endl;
	print_vector_of_individuals_indirect(std::cout, destinationpopulation.toVector(), destinationpopulationIndices) << std::endl;
	*/

	//Sort the arrays of indices
	//From best to worst
	std::sort(incomingpopulationIndices.begin(), incomingpopulationIndices.end(), IndirectIndividualSorter(incomingpopulation));
	//From worst to best
	std::sort(destinationpopulationIndices.begin(), destinationpopulationIndices.end(), IndirectIndividualSorter(destinationpopulation.toVector()));
	std::reverse(destinationpopulationIndices.begin(), destinationpopulationIndices.end());

	/*
	std::cout << "Sorted orders:" << std::endl;
	print_vector_of_individuals_indirect(std::cout, incomingpopulation, incomingpopulationIndices) << std::endl;
	print_vector_of_individuals_indirect(std::cout, destinationpopulation.toVector(), destinationpopulationIndices) << std::endl;
	*/

	// Create the result
	std::list<std::pair<int, int> > result;

	//std::cout << "Assigments:" << std::endl;

	for (int i = 0; i < migrationRateLimit; i++) {
		//Add the assigment only if the new one is better
		if (individual::compare_by_fitness(incomingpopulation[incomingpopulationIndices[i]], destinationpopulation[destinationpopulationIndices[i]])) {
			result.push_back(std::make_pair(destinationpopulationIndices[i], incomingpopulationIndices[i]));

			//std::cout << destinationpopulationIndices[i] << " <- " << incomingpopulationIndices[i] << std::endl;
		}
	}

	return result;
}

}
