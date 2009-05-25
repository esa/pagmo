#ifndef CARTESIAN_CS_H
#define CARTESIAN_CS_H

#include <boost/shared_ptr.hpp>

#include "coordinate_system.h"

struct cartesian_cs: public coordinate_system {
	boost::shared_ptr<coordinate_system> clone() const
	{
		return boost::shared_ptr<coordinate_system>(new cartesian_cs());
	}
};

#endif
