#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "epoch.h"
#include "core_functions/propagate_lagrangian.h"
#include "planet.h"
#include "astro_constants.h"
#include "sims_flanagan/leg.h"
#include "spacecraft.h"
#include "sims_flanagan/sc_state.h"
#include "sims_flanagan/fb_traj.h"
#include "codings.h"



using namespace std;
using namespace kep_toolbox;
int main() {
	long hrs,min,sec;
	double mjd2000 = 10 / 24 + 12.0 / 1440 + 30.123456789/ 86400;
	hrs = (long)(mjd2000 * 24);
	min = (long) ( (mjd2000*24-hrs) * 60);
	sec = (long) ( ( ( (mjd2000*24-hrs) * 60) - min ) * 60 );
	double dblfsec = ( ( ( (mjd2000*24-hrs) * 60) - min ) * 60 ) - sec;
	cout << endl << endl << "DEBUG OF CLASS EPOCH" << endl;
	cout << hrs << " " << min << " " << sec << " " <<  setiosflags(ios::fixed) << dblfsec << endl;

	epoch start(mjd2000);
	boost::posix_time::ptime time(boost::gregorian::date(2002,1,10), 
				      boost::posix_time::time_duration(1,2,3));
	

	planet earth(planet::EARTH);
	
	epoch date1(0.);
	epoch date2(200);
	array3D r,v;
	earth.get_eph(date1,r,v);
	cout << endl << endl <<  "DEBUG OF CLASS PLANET" << endl;
	cout << "Earth in " << date1 <<endl;
	cout << "Pos: " << r[0] << " " << r[1] << " " << r[2] << endl;
	cout << "Vel: " << v[0] << " " << v[1] << " " << v[2] << endl;
	earth.get_eph(date2,r,v);
	cout << "Earth in " << date2 <<endl;
	cout << "Pos: " << r[0] << " " << r[1] << " " << r[2] << endl;
	cout << "Vel: " << v[0] << " " << v[1] << " " << v[2] << endl;

	earth.get_eph(date1,r,v);
	propagate_lagrangian(r,v, 100*ASTRO_DAY2SEC,ASTRO_MU_SUN);
	propagate_lagrangian(r,v, 100*ASTRO_DAY2SEC,ASTRO_MU_SUN);
	cout << "Earth in " << date2 <<endl;
	cout << "Pos: " << r[0] << " " << r[1] << " " << r[2] << endl;
	cout << "Vel: " << v[0] << " " << v[1] << " " << v[2] << endl;
	
	
	cout << endl << endl <<  "DEBUG OF CLASS LEG" << endl;
	int n_seg = 10;
	spacecraft sc(2000,0.05,2000);

	cout << endl << endl <<  "DEBUG OF CLASS FB_TRAJ" << endl;
	std::vector<planet> sequence; sequence.push_back(earth); sequence.push_back(earth); sequence.push_back(earth);
	sims_flanagan::fb_traj ee_lt(sequence,n_seg,sc);
	std::vector<double> x;
	//Departure date
	x.push_back(0);

	//First Leg
	//Starting mass
	//x.push_back(sc.get_mass());
	//Starting Vinf
	x.push_back(0); x.push_back(0); x.push_back(0);
	//Throttles
	for(int i =0;i < n_seg;i++){
		x.push_back(0); x.push_back(0); x.push_back(0);
	}
	//Final Vinf
	x.push_back(0); x.push_back(0); x.push_back(0);
	//Final  mass
	//x.push_back(sc.get_mass());
	//Transfer time
	x.push_back(200);// * ASTRO_DAY2SEC);

	//Second leg
	//Starting mass
	//x.push_back(sc.get_mass());
	//Starting Vinf
	x.push_back(0); x.push_back(0); x.push_back(0);
	//Throttles
	for(int i =0;i < n_seg;i++){
		x.push_back(0); x.push_back(0); x.push_back(0);
	}
	//Final Vinf
	x.push_back(0); x.push_back(0); x.push_back(0);
	//Final  mass
	//x.push_back(sc.get_mass());
	//Transfer time
	x.push_back(400);

	//Init trajectory
	base_format coding(2, n_seg, sc.get_mass());
	ee_lt.init_from_full_vector(x.begin(), x.end(), coding);
	//Print trajectory
	std::cout << ee_lt << std::endl;

}

