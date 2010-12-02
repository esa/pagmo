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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/utility.hpp>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

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
		static bool is_multithread();
		static int get_size();
		static int get_rank();
		template <class T>
		static void recv(T &retval, int source)
		{
			MPI_Status status;
			// First receive the size.
			int size;
			MPI_Recv(static_cast<void *>(&size),1,MPI_INT,source,0,MPI_COMM_WORLD,&status);
			// Prepare the vector of chars.
			std::vector<char> buffer_char(boost::numeric_cast<std::vector<char>::size_type>(size),0);
			// Receive the payload.
			MPI_Recv(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,source,1,MPI_COMM_WORLD,&status);
			// Build the string from the vector.
			const std::string buffer_str(buffer_char.begin(),buffer_char.end());
			// Unpickle the payload.
			std::stringstream ss(buffer_str);
			boost::archive::text_iarchive ia(ss);
			ia >> retval;
		}
		template <class T>
		static void send(const T &payload, int destination)
		{
			std::stringstream ss;
			boost::archive::text_oarchive oa(ss);
			oa << payload;
			const std::string buffer_str(ss.str());
			std::vector<char> buffer_char(buffer_str.begin(),buffer_str.end());
			// Send the size.
			int size = boost::numeric_cast<int>(buffer_char.size());
			MPI_Send(static_cast<void *>(&size),1,MPI_INT,destination,0,MPI_COMM_WORLD);
			// Send the string.
			MPI_Send(static_cast<void *>(&buffer_char[0]),size,MPI_CHAR,destination,1,MPI_COMM_WORLD);
		}
		static bool iprobe(int);
	private:
		void listen();
};

}

#endif
