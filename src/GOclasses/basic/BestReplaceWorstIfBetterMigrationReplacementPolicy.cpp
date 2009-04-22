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

#include "BestReplaceWorstIfBetterMigrationReplacementPolicy.h"
#include "../../exceptions.h"
#include <algorithm>
#include <iostream>

std::ostream& print_vector_of_individuals_indirect(std::ostream& str, const std::vector<Individual>& vect_indiv, const std::vector<int>& vect_index) {
	for(std::vector<int>::const_iterator it = vect_index.begin(); it != vect_index.end(); ++it) {
		str << vect_indiv[*it].getFitness() << " ";
	}
	return str;
}

// 09/03/2009: Initial version by Marek Rucinski.

std::list<std::pair<int, int> > BestReplaceWorstIfBetterMigrationReplacementPolicy::selectForReplacement(const std::vector<Individual>& incomingPopulation, const Population& destinationPopulation)
{
	int migrationRateLimit = std::min(getMaxMigrationRate(destinationPopulation), (int)incomingPopulation.size());
	
	// Temporary vectors to store sorted poopulations
	std::vector<int> incomingPopulationIndices(incomingPopulation.size());
	std::vector<int> destinationPopulationIndices(destinationPopulation.size());
	
	//Fill in the arrays of indices
	std::generate(incomingPopulationIndices.begin(), incomingPopulationIndices.end(), IndexGenerator());
	std::generate(destinationPopulationIndices.begin(), destinationPopulationIndices.end(), IndexGenerator());
	
	/*
	std::cout << "Initial individuals orders:" << std::endl;
	print_vector_of_individuals_indirect(std::cout, incomingPopulation, incomingPopulationIndices) << std::endl;
	print_vector_of_individuals_indirect(std::cout, destinationPopulation.toVector(), destinationPopulationIndices) << std::endl;
	*/
	
	//Sort the arrays of indices
	//From best to worst
	std::sort(incomingPopulationIndices.begin(), incomingPopulationIndices.end(), IndirectIndividualSorter(incomingPopulation));
	//From worst to best
	std::sort(destinationPopulationIndices.begin(), destinationPopulationIndices.end(), IndirectIndividualSorter(destinationPopulation.toVector()));
	std::reverse(destinationPopulationIndices.begin(), destinationPopulationIndices.end());
	
	/*
	std::cout << "Sorted orders:" << std::endl;
	print_vector_of_individuals_indirect(std::cout, incomingPopulation, incomingPopulationIndices) << std::endl;
	print_vector_of_individuals_indirect(std::cout, destinationPopulation.toVector(), destinationPopulationIndices) << std::endl;
	*/
	
	// Create the result
	std::list<std::pair<int, int> > result;
	
	//std::cout << "Assigments:" << std::endl;
	
	for(int i = 0; i < migrationRateLimit; i++) {
		//Add the assigment only if the new one is better
		if(Individual::compare_by_fitness(incomingPopulation[incomingPopulationIndices[i]], destinationPopulation[destinationPopulationIndices[i]])) {
			result.push_back(std::make_pair(destinationPopulationIndices[i], incomingPopulationIndices[i]));
			
			//std::cout << destinationPopulationIndices[i] << " <- " << incomingPopulationIndices[i] << std::endl;
		}
	}
	
	return result;
}
