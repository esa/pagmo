#ifndef PAGMO_EXCEPTIONS_H
#define PAGMO_EXCEPTIONS_H

#include <iostream>
#include <string>

class base_exception
{
	public:
		base_exception(const std::string &s): m_what(s) {}
		const std::string &what() const {
			return m_what;
		}
	protected:
		std::string m_what;
};

struct index_error: public base_exception {
	index_error(const std::string &s): base_exception(s) {}
};

#endif
