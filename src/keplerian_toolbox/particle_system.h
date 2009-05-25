#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <boost/shared_ptr.hpp>
#include <vector>

#include "orbit_propagator.h"
#include "particle.h"

class particle_system {

	private:
		std::vector<particle>			m_container;
		boost::shared_ptr<orbit_propagator>	m_op;
};

#endif
