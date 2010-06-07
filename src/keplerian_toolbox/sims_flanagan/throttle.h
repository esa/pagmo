#ifndef THROTTLE_H
#define THROTTLE_H

#include "../epoch.h"
#include "../astro_constants.h"
#include <numeric>

namespace kep_toolbox { namespace sims_flanagan{

/// A single throttle
/**
 * This class models a single throttle in the Sims-Flanagan model. It essentialy contains the cartesian
 * components of one throttle (non dimensional impulse)
 *impulse
 * @author David di Lorenzo
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class throttle {
public:
	throttle() {}

	throttle(epoch _start, epoch _end, const array3D& _value)
		: m_start(_start), m_end(_end), m_value(_value) {}

	epoch get_start() const {
		return m_start;
	}

	epoch get_end() const {
		return m_end;
	}

	const array3D& get_value() const {
		return m_value;
	}

	double get_norm() const {
		return std::sqrt(std::inner_product(m_value.begin(), m_value.end(), m_value.begin(), 0.));
	}


private:
	epoch m_start;
	epoch m_end;
	array3D m_value;
};

}} //namespaces

#endif
