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

#include "MigrationPolicy.h"
#include "../../../exceptions.h"

// 09/03/2009: Initial version by Marek Rucinski.

namespace pagmo
{

void MigrationPolicy::setMigrationProbability(const double _migrationProbability)
{
	if ((_migrationProbability < 0.0) || (_migrationProbability > 1.0)) {
		pagmo_throw(value_error, "Migration probability must be within the [0.0, 1.0] range.");
	}
	migrationProbability = _migrationProbability;
}

std::ostream &operator<<(std::ostream &s, const MigrationPolicy& msp)
{
	s << "Migration probability: " << msp.migrationProbability << std::endl;
	s << "Selection policy:      " << std::endl << *(msp.selectionPolicy) << std::endl;
	s << "Replacement policy:    " << std::endl << *(msp.replacementPolicy) << std::endl;
	return s;
}

}
