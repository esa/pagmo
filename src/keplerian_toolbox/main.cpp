#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random.hpp>
#include "keplerian_toolbox.h"
#include "../AstroToolbox/lambert.h"



using namespace std;
using namespace kep_toolbox;
int main() {
	array3D r1,r2,v1,v2,v_1,v_2;
	double a,p,th;

	int Ntrials = 100,lw;
	double tof,err1=0,err2=0;
	boost::mt19937 rng;
	boost::uniform_int<> dist(0, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_bit(rng, dist);
	boost::uniform_real<> dist1(-2,2);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > drng(rng, dist1);

	for (int i = 0; i<Ntrials; ++i){
		r1[0] = drng() * 2; r1[1] = drng() * 2; r1[2] = drng() * 2;
		r2[0] = drng() * 2; r2[1] = drng() * 2; r2[2] = drng() * 2;
		tof = (drng () + 2.1 ) * 10;

		lw = rand_bit();
				std::cout << std::setprecision(15) << r1 << r2 << tof << " " << lw << std::endl;
		int retval = lambert_3d(v1,v2,a,p,r1,r2,tof,1,lw);
		LambertI (&r1[0],&r2[0],tof,1,lw,&v_1[0],&v_2[0],a,p,th,retval);
		err1 += sqrt((v1[0]-v_1[0])*(v1[0]-v_1[0]) + (v1[1]-v_1[1])*(v1[1]-v_1[1]) + (v1[2]-v_1[2])*(v1[2]-v_1[2]));
		err2 += sqrt((v2[0]-v_2[0])*(v2[0]-v_2[0]) + (v2[1]-v_2[1])*(v2[1]-v_2[1]) + (v2[2]-v_2[2])*(v2[2]-v_2[2]));
	}
	std::cout << err1/Ntrials <<std::endl;
	std::cout << err2/Ntrials <<std::endl;

}
