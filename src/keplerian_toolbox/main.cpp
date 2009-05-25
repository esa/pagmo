#include "cartesian_cs.h"
#include "particle.h"
#include "spherical_cs.h"

int main()
{
	double state[] = {1,2,3,4,5,6};
	particle p(1,state,cartesian_cs());
	std::cout << p << '\n';
	p.set_cs(spherical_cs());
	std::cout << p << '\n';
	p.set_cs(cartesian_cs());
	std::cout << p << '\n';

//	pv_state<double> s(std::vector<double>(3,1),std::vector<double>(3,1));
//	array_d3 v = {{1,2,3}}, p = {{1,2,3}};
//	s.set_velocity(v);
//	s.set_position(v);
//	std::cout << s << '\n';
//	std::cout << s.set_coordinate_system(spherical_coordinate_system<double>()) << '\n';
//	std::cout << s.set_coordinate_system(spherical_coordinate_system<double>()) << '\n';
//	std::cout << s.set_coordinate_system(cartesian_coordinate_system<double>()) << '\n';
//	std::cout << pv_state<double>(s).set_coordinate_system(cartesian_coordinate_system<double>()) << '\n';
//	std::cout << sizeof(pv_state<double>) << '\n';
//
//	dynamical_system<double> d;

	//for (size_t i = 0; i < 1000000; ++i) {
	//	s.propagate(i);
	//}
}
