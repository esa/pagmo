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

#include "ConstInRateMigrationScheme.h"

// 09/03/2009: Initial version by Marek Rucinski.

void ConstInRateMigrationScheme::preEvolutionCallback(island& _island)
{
	// Choose a random incoming neighbour and perform the migration.
	lock_type lock(topology_mutex);
	
	const std::vector<size_t>& neighbours = topology->get_neighbours_in(_island.id());
	
	if(neighbours.size() > 0) { //the island must have neighbours
			
		//Draw a neighbour
		size_t chosenNeighbourIndex = rng() % neighbours.size();
		size_t chosenNeighbour = neighbours[chosenNeighbourIndex];
		
		//Get neighbour's immigrants
		std::vector<Individual> immigrants = immigrantsDB[chosenNeighbour];
		
		if(immigrants.size() > 0) { //if there's anything to migrate, do it
			_island.acceptMigratingIndividuals(immigrants);
		}	
	}
}

void ConstInRateMigrationScheme::postEvolutionCallback(island& _island)
{
	//Update _island's database with it's migrating individuals.
	lock_type lock(topology_mutex);
	
	//Get individuals from the islands
	std::vector<Individual> immigrants = _island.getMigratingIndividuals();
	
	//Replace the contents of the database for the island
	immigrantsDB[_island.id()].swap(immigrants);
}
