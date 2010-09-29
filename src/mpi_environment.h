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

#ifndef PAGMO_MPI_ENVIRONMENT_H
#define PAGMO_MPI_ENVIRONMENT_H

#include <boost/mpi/environment.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/utility.hpp>

#include "config.h"

/*!
\page mpi_support MPI support in PaGMO
PaGMO can be configured to use MPI to dispatch parallel optimisation instances to computers participating in a cluster.
\section mpi_requirements Requirements
The requirements to enable and use MPI support in PaGMO are the following:

- foo
- bar


*/

namespace pagmo
{

/// MPI environment class.
class __PAGMO_VISIBLE mpi_environment: private boost::noncopyable
{
	public:
		mpi_environment();
		~mpi_environment();
		int size() const;
	private:
		void listen();
	private:
		boost::scoped_ptr<boost::mpi::environment> m_environment;
};

}

#endif
