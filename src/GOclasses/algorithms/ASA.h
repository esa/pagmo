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

// 16/05/2008: Initial version by Dario Izzo.

#ifndef PAGMO_ASA_H
#define PAGMO_ASA_H

#include <iostream>

#include "../../config.h"
#include "../basic/population.h"
#include "go_algorithm.h"

class __PAGMO_VISIBLE ASAalgorithm: public go_algorithm {
	public:
		//This method initialise the SA-AN algorithm starting and final temperature setting deafult values for
		//the StartStep, the niterTemp and the niterRange. Tcoeff is evaluated accordingly.
		ASAalgorithm(int, const double &, const double &);
		virtual Population evolve(const Population &) const;
		virtual ASAalgorithm *clone() const {return new ASAalgorithm(*this);}
	private:
		virtual void log(std::ostream &) const;
		size_t niterTot;
		size_t niterTemp;
		size_t niterRange;
		double Ts;
		double Tf;
		double StartStep;
};

#endif
