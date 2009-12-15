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

#include "MigrationScheme.h"

// 09/03/2009: Initial version by Marek Rucinski.

namespace pagmo
{

void MigrationScheme::preEvolutionCallback(island& _island)
{
	lock_type lock(topology_mutex);

	//std::cout << "preEvolutionCallback( " << _island.id() << " )" << std::endl;

	//1. Obtain immigrants
	std::vector<individual> immigrants;

	if (!migrationDirection) {
		// For migrationDirection == 0, immigrantsByIsland are just "inboxes" of particular islands.
		immigrants.swap(immigrantsByIsland[_island.id()]);

		//std::cout << immigrants.size() << " individuals were available in the \"inbox\"!" << std::endl;
	} else {
		// For migrationDirection == 1, immigrantsByIsland behave like "outboxes", i.e. each is a "database of best individuals" for corresponding island

		const std::vector<size_t> neighbours = topology->get_neighbours_in(_island.id());

		//std::cout << "The island has " << neighbours.size() << " neigbours." << std::endl;

		if (neighbours.size() > 0) {
			if (!distributionType) {
				// For one-to-one migration choose a random neighbour island and fetch immigrants
				size_t chosenNeighbourIndex = (size_t)(rng() * neighbours.size());
				size_t chosenNeighbour = neighbours[chosenNeighbourIndex];

				//std::cout << "Randomly chosen neighbour id: " << chosenNeighbour << std::endl;

				immigrants = immigrantsByIsland[chosenNeighbour];

				//std::cout << "The neighbour provides " << immigrants.size() << " individuals" << std::endl;
			} else {
				// For broadcast migration fetch immigrants from all neighbour islands' databases
				for (size_t i = 0; i < neighbours.size(); i++) {
					//std::cout << "Neighbour island with id " << neighbours[i] << " provides " << immigrantsByIsland[neighbours[i]].size() << " individuals" << std::endl;
					immigrants.insert(immigrants.end(), immigrantsByIsland[neighbours[i]].begin(), immigrantsByIsland[neighbours[i]].end());
				}

				//std::cout << "The neighbours provide in total " << immigrants.size() << " individuals" << std::endl;
			}
		}
	}

	//2. Insert immigrants into population
	if (immigrants.size() > 0) { //if there's anything to migrate, do it
		//but take into account the island's migration probability

		if (rng() < _island.getMigrationProbability()) {
			_island.acceptMigratingIndividuals(immigrants);
		}
	}
}

void MigrationScheme::postEvolutionCallback(island& _island)
{
	lock_type lock(topology_mutex);

	//std::cout << "postEvolutionCallback( " << _island.id() << " )" << std::endl;

	//1. Obtain individuals migrating from the island
	std::vector<individual> immigrants;

	//2. Store them in an appropriate place
	if (!migrationDirection) {
		// For migrationDirection == 0, immigrantsByIsland are "inboxes" of particular islands.

		const std::vector<size_t> neighbours = topology->get_neighbours_out(_island.id());

		//std::cout << "The island has " << neighbours.size() << " neigbours." << std::endl;

		if (neighbours.size() > 0) {
			immigrants = _island.getMigratingIndividuals();

			//std::cout << "The island migrates " << immigrants.size() << " individuals" << std::endl;

			if (immigrants.size() > 0) {
				if (!distributionType) {
					// For one-to-one migration choose a random neighbour island and put immigrants to its inbox
					size_t chosenNeighbourIndex = (size_t)(rng() * neighbours.size());
					size_t chosenNeighbour = neighbours[chosenNeighbourIndex];

					//std::cout << "Randomly chosen neighbour id: " << chosenNeighbour << std::endl;

					immigrantsByIsland[chosenNeighbour].insert(immigrantsByIsland[chosenNeighbour].end(), immigrants.begin(), immigrants.end());

					//std::cout << "The size of the neighbour's inbox is now "  << immigrantsByIsland[chosenNeighbour].size() << std::endl;
				} else {
					// For broadcast migration put immigrants to all neighbour islands' inboxes
					for (size_t i = 0; i < neighbours.size(); i++) {
						immigrantsByIsland[neighbours[i]].insert(immigrantsByIsland[neighbours[i]].end(), immigrants.begin(), immigrants.end());

						//std::cout << "The size of the inbox of a neighbour with id " << neighbours[i] << " is now " << immigrantsByIsland[neighbours[i]].size() << std::endl;
					}
				}
			}
		}
	} else {
		// For migrationDirection == 1, immigrantsByIsland behave like "outboxes", i.e. each is a "database of best individuals" for corresponding island
		immigrants = _island.getMigratingIndividuals();

		//std::cout << "Storing " << immigrants.size() << " individuals int the DB of the island " << _island.id() << std::endl;

		immigrantsByIsland[_island.id()].swap(immigrants);
	}
}

std::ostream &operator<<(std::ostream &s, const MigrationScheme& ms)
{
	s << "Migration algorithm: " << (ms.distributionType ? "broadcast, " : "one-to-one, ") << (ms.migrationDirection ? "initiated by destination" : "initiated by source") << std::endl;

	s << "Topology:            ";
	if (ms.topology) {
		s << std::endl << *(ms.topology) << std::endl;
	} else {
		s << "none" << std::endl;
	}

	return s;
}

}
