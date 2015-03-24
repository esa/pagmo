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

#include <gsl/gsl_errno.h>
#include <iostream>

#include "gsl_init.h"

namespace pagmo
{

gsl_init::gsl_init():m_init(false)
{
	gsl_set_error_handler(&custom_gsl_error_handler);
	m_init = true;
}

void custom_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno)
{
	std::cout << "WARNING: GSL reported an error, details follow.\n\n";
	std::cout << "Reason: " << reason << '\n';
	std::cout << "File: " << file << '\n';
	std::cout << "Line: " << line << '\n';
	std::cout << "GSL error code: " << gsl_errno << '\n';
}

}
