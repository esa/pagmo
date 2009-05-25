#ifndef COORDINATE_SYSTEM_H
#define COORDINATE_SYSTEM_H

#include <boost/shared_ptr.hpp>
#include <typeinfo>
#include <vector>

struct coordinate_system {
	virtual void to_cartesian(std::vector<double> &) const {}
	virtual void from_cartesian(std::vector<double> &) const {}
	virtual boost::shared_ptr<coordinate_system> clone() const = 0;
	virtual bool operator==(const coordinate_system &other) const
	{
		return (typeid(*this) == typeid(other));
	}
	virtual ~coordinate_system() {}
};

#endif
