/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

#include "RandomMigrationReplacementPolicy.h"
#include "../../exceptions.h"
#include <algorithm>

// 09/03/2009: Initial version by Marek Rucinski.

std::list<std::pair<int, int> > RandomMigrationReplacementPolicy::selectForReplacement(const std::vector<Individual>& incomingPopulation, const Population& destinationPopulation)
{
	int migrationRateLimit = std::min<int>(getMaxMigrationRate(destinationPopulation), (int)incomingPopulation.size());
	
	std::vector<int> incomingPopulationIndices(migrationRateLimit);
	std::vector<int> destinationPopulationIndices(destinationPopulation.size());
	
	//Fill in the arrays of indices
	std::generate(incomingPopulationIndices.begin(), incomingPopulationIndices.end(), IndexGenerator());
	std::generate(destinationPopulationIndices.begin(), destinationPopulationIndices.end(), IndexGenerator());
	
	//Permute the indices (incoming population)
	for(int i = 0; i < migrationRateLimit; i++) {
		int nextIncomingPopulationIndex = i + (rng() % (migrationRateLimit - i));
		if(nextIncomingPopulationIndex != i) {
			std::swap(incomingPopulationIndices[i], incomingPopulationIndices[nextIncomingPopulationIndex]);
		}				
	}
	
	//Permute the indices (destination population)
	for(int i = 0; i < migrationRateLimit; i++) {
		int nextDestinationPopulationIndex = i + (rng() % (destinationPopulation.size() - i));
		if(nextDestinationPopulationIndex != i) {
			std::swap(destinationPopulationIndices[i], destinationPopulationIndices[nextDestinationPopulationIndex]);
		}				
	}
	
	// Create the result
	std::list<std::pair<int, int> > result;
	
	for(int i = 0; i < migrationRateLimit; i++) {
		result.push_back(std::make_pair(destinationPopulationIndices[i], incomingPopulationIndices[i]));
	}
	
	//std::cout << "Result size: " << result.size() << std::endl;
	
	return result;
}
