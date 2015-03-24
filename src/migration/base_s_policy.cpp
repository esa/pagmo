/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include "base_s_policy.h"
#include "base.h"

namespace pagmo { namespace migration {

/// Constructor from migration rate and type.
/**
 * By default, one single individual will be selected for migration.
 *
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 *
 * @see migration::base::base().
 */
base_s_policy::base_s_policy(const double &rate, rate_type type):base(rate,type)
{}

/// Trivial destructor.
/**
 * No side effects.
 */
base_s_policy::~base_s_policy() {}

}}
