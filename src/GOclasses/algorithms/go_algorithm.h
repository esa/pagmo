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

#ifndef PAGMO_GO_ALGORITHM_H
#define PAGMO_GO_ALGORITHM_H

#include <string>
#include <typeinfo>

#include "../../Functions/rng/rng.h"
#include "../../config.h"
#include "../basic/population.h"

class __PAGMO_VISIBLE go_algorithm
{
	public:
		go_algorithm();
		go_algorithm(const go_algorithm &);
		go_algorithm &operator=(const go_algorithm &);
		virtual Population evolve(const Population &) const = 0;
		virtual go_algorithm *clone() const = 0;
        virtual ~go_algorithm() {}
		std::string id_name() const {return typeid(*this).name();}
	protected:
		mutable rng_double drng;
};

#endif
