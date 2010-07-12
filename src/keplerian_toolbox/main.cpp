#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "keplerian_toolbox.h"



using namespace std;
using namespace kep_toolbox;
int main() {
	double s,c,tof,a,p,vr1,vt1,vr2,vt2;
	int lw = 0, iter;
	s = 3.0/2.0;
	c = 1;
	tof = 0.01;
	iter = lambert_2d(s,c,tof, lw, a, p, vr1, vt1, vr2, vt2);
	std::cout << "a: " << a << std::endl;
	std::cout << "p: " << p << std::endl;
	std::cout << "vr1: " << vr1 << std::endl;
	std::cout << "vt1: " << vt1 << std::endl;
	std::cout << "vr2: " << vr2 << std::endl;
	std::cout << "vt2: " << vt2 << std::endl;
	std::cout << "iter: " << iter << std::endl;

}

