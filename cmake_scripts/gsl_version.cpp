#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_version.h>
#include <iostream>
#include <string>
#include <vector>

int main()
{
	std::vector<std::string> v;
	boost::split(v,GSL_VERSION,boost::algorithm::is_any_of("."));
	std::cout << "\n\tDetected GSL version: " << GSL_VERSION;
	if (v.size() < 2) {
		return -1;
	}
	if (boost::lexical_cast<int>(v[0]) < 1) {
		return -1;
	} else if (boost::lexical_cast<int>(v[0]) == 1 && boost::lexical_cast<int>(v[1]) < 15) {
		return -1;
	}
	return 0;
}
