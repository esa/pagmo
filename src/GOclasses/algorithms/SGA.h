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

// Created by Dario Izzo on 10/05/08.

#ifndef PAGMO_SGA_H
#define PAGMO_SGA_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../../config.h"
#include "../basic/population.h"
#include "go_algorithm.h"

class __PAGMO_VISIBLE SGAalgorithm: public go_algorithm {
	public:
		SGAalgorithm(int, const double &, const double &, int insert_bestInit);
		virtual Population evolve(const Population &) const;
		virtual SGAalgorithm *clone() const {return new SGAalgorithm(*this);}
		virtual std::string id_object() const { return id_name(); }
	private:
		virtual void log(std::ostream &) const;
		const size_t	generations;
		const double 	CR;		//crossover
		const double	M;		//mutation
		const int		insert_best;
		mutable			rng_uint32 rng;
};

#endif
