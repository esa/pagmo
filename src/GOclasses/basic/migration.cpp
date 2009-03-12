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

#include "migration.h"

// 09/03/2009: Initial version by Marek Rucinski.

// ==== DummySelectionPolicy ====

boost::shared_ptr<Population> DummySelectionPolicy::selectForMigration(const Population& population)
{
	boost::shared_ptr<Population> result(new Population(population.problem())); //See shared pointers best practices
	return result;
}




// ==== DummyReplacementPolicy ====

boost::shared_ptr<std::list<std::pair<int, int> > > DummyReplacementPolicy::selectForReplacement(const Population& incomingPopulation, const Population& destinationPopulation)
{
	boost::shared_ptr<std::list<std::pair<int, int> > > result(new std::list<std::pair<int, int> >()); //See shared pointers best practices
	return result;
}
