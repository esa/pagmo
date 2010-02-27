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

// 07/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_PY_CONTAINER_UTILS_H
#define PAGMO_PY_CONTAINER_UTILS_H

#include "../../exceptions.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace pagmo
{

template <class Derived>
class py_container_utils
{
	protected:
		template <class ConstIterator>
		ConstIterator it_from_index(int n) const {
			const size_t index = get_ra_index(n);
			ConstIterator it = derived_const_cast->begin();
			for (size_t i = 0; i < index; ++i, ++it) {}
			return it;
		}
		template <class Iterator>
		Iterator it_from_index(int n) {
			const size_t index = get_ra_index(n);
			Iterator it = derived_cast->begin();
			for (size_t i = 0; i < index; ++i, ++it) {}
			return it;
		}
		size_t get_ra_index(int n) const {
			const size_t size = derived_const_cast->size();
			if (n >= 0) {
				if ((size_t)n >= size) {
					pagmo_throw(index_error,"container index out of range");
				}
				return n;
			} else {
				if (size < (size_t)(-n)) {
					pagmo_throw(index_error,"container index out of range");
				}
				return (size - (size_t)(-n));
			}
		}
};

}

#undef derived_cast
#undef derived_const_cast

#endif
