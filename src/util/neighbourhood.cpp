# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

# include "neighbourhood.h"

using namespace std;
namespace pagmo{ namespace util {namespace neighbourhood {

/**
 * Compute the neighbourhood graph. At the end of the call retval[i][j] will contain the j-th closest vector
 * (according to the euclidian distance) to the i-th vector.
 * @param[out] retval a matrix representing the neigborhood graph
 * @param[in]  weights the vector of real vectors
 */
void euclidian::compute_neighbours(std::vector<std::vector<pagmo::population::size_type> > &retval, const std::vector<std::vector<double> > &weights) {
	for(unsigned int i = 0; i < weights.size(); ++i) {
		std::vector<double> distances;
		for(unsigned int j = 0; j < weights.size(); ++j) {
			distances.push_back(distance(weights[i],weights[j]));
		}
		retval.push_back(order(distances));
	}
}

/**
 * Compute the euclidian distance between two real vectors
 * @param a first vector
 * @param b second vector
 * @return euclidian distance between the vectors
 */
double euclidian::distance(const std::vector<double> &a, const std::vector<double> &b) {
	double rtr = 0.0;
	for(std::vector<double>::size_type i = 0; i < a.size(); ++i) {
		rtr += pow(a[i]-b[i], 2);
	}
	return sqrt(rtr);
}

}}} //namespaces
