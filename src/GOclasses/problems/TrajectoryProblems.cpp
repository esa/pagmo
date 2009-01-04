/*
 *  GOProblem.cpp
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <vector>

#include "TrajectoryProblems.h"
#include "GOproblem.h"
#include "trajobjfuns.h"

//***********************************************************************************
//Trajectory problems
//***********************************************************************************

const double messengerfullProb::lb[26] = {1900, 3,    0, 0, 100, 100, 100, 100, 100, 100, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
const double messengerfullProb::ub[26] = {2200, 4.05, 1, 1, 500, 500, 500, 500, 500, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};
const int messengerfullProb::sequence[7] = {3, 2, 2, 1, 1, 1, 1};

messengerfullProb::messengerfullProb():GOProblem(26,lb,ub),mgadsm(orbit_insertion,sequence,7,0,0,0,0.704,2440 + 200) {
	//Bounds shrinked further with the c3=25
	//double lb[26] = {2000, 4, 0, 0.6, 428, 210, 210, 250, 345, 520, 0.01, 0.01, 0.5, 0.5, 0.6, 0.8, 1.1, 1.1, 1.05, 1.05, 1.05, 0, 1, -M_PI, 0, 0};
	//double ub[26] = {2070, 5, 0.65, 1  , 468,    230,    230,    270,    365,    550,    0.99,    0.99,    0.99,    0.7,    0.9,    0.95,   6,   1.5,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  3};
	
	//Bounds shrinked with the c3=25
	//double lb[26] = {2000, 4, 0, 0, 428, 210, 210, 250, 340, 520, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
	//double ub[26] = {2100, 5, 1, 1, 468, 240, 240, 280, 370, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};

    //Bounds shrinked to those of the actual mission (main difference is that thi forces the used resonances at mercury)
	//double lb[26] = {1900, 3,    0, 0, 400, 200, 200, 200, 300, 400, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
	//double ub[26] = {2200, 4.05, 1, 1, 500, 300, 300, 300, 400, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};
};

double messengerfullProb::objfun(const std::vector<double>& x) const {
   	double obj = 0.0;
	MGA_DSM(x, mgadsm,
			obj);
	return obj;
}


const double messengerProb::lb[18] = {1000, 1, 0, 0, 200, 30,  30,  30,  0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.1, -M_PI, -M_PI, -M_PI};
const double messengerProb::ub[18] = {1200, 5, 1, 1, 400, 400, 400, 400, 0.99, 0.99, 0.99, 0.99, 6,   6,   6,    M_PI,  M_PI,  M_PI};
const int messengerProb::sequence[5] = {3, 3, 2, 2, 1};

messengerProb::messengerProb():GOProblem(18,lb,ub),mgadsm(total_DV_rndv,sequence,5,0,0,0,0,0) {};

double messengerProb::objfun(const std::vector<double>& x) const {
   	double obj = 0.0;
	MGA_DSM(x, mgadsm,
			obj);
	return obj;
}

const double tandemProb::lb[18] = {5475, 2.5, 0, 0, 20   , 20  ,  20 , 20  , 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI};
const double tandemProb::ub[18] = {9132, 4.9, 1, 1, 2500 , 2500, 2500, 2500, 0.99, 0.99, 0.99, 0.99,    6,    6,    6,  M_PI,  M_PI,  M_PI};
const int tandemProb::sequence[5] = {3, 2, 3, 3, 6};

tandemProb::tandemProb():GOProblem(18,lb,ub),mgadsm(orbit_insertion,sequence,5,0,0,0,0.98531407996358,80330.0) {
  /*double lb[18] = { 5475, 2.938860000000, 0.405043000000, 0.315992000000, 94.289800000000, 342.079000000000, 1008.410000000000, 1200.380000000000, 0.028803000000, 0.340659000000, 0.060535100000, 0.112446000000, 1.051750000000, 1.078730000000, 1.053800000000, -1.870500000000, -2.225840000000, -1.886840000000};
    double ub[18] = {9132, 4.093870000000, 0.586676000000, 0.707644000000, 365, 600.040000000000, 1200.320000000000, 2000, 0.897733000000, 0.804798000000, 0.805629000000, 0.419066000000, 2.370910000000, 1.422330000000, 2.156690000000, -1.169580000000, -1.127050000000, -0.917483000000};*/

  /*double lb[18] = { 6844.451830612043, 3.067699723709, 0.477284845440, 0.417781690881, 108.507200707471, 511.976881365260, 1004.050120849356, 1784.974706859021, 0.081005224956, 0.364221914356, 0.102296747234, 0.117776456111, 1.052250007616, 1.217251823179, 1.050000000000, -1.927995228102, -1.742401844681, -1.991359970502};
    double ub[18] = { 9063.429408756601, 3.640972901933, 0.541267874511, 0.666894655742, 1284.558367334185, 1393.888844450796, 1193.773380479754, 2314.454268212400, 0.837754757113, 0.803612856456, 0.559199754928, 0.340380612138, 1.578888906483, 1.414129557493, 1.473108613200, -1.207415304290, -1.162377264313, -1.075256620654};*/

   /*double lb[18] = { 7006.238587842692, 3.097717073232, 0.487979915911, 0.420687094540, 109.771414690208, 940.854505183490, 1084.392758236149, 1884.930756684131, 0.203810643630, 0.524034320427, 0.061286469418, 0.094306544631, 1.083463972542, 1.236225186517, 1.050000000000, -1.937841604629, -1.743932449066, -1.832730790569};
     double ub[18] = { 9060.215094956218, 3.611338874179, 0.527214806521, 0.657358693786, 1197.610002758771, 1394.575620926851, 1196.061120842642, 2282.478189832147, 0.802342883580, 0.888447799151, 0.686859262454, 0.297661986078, 1.543881588983, 1.456965103045, 1.469911762813, -1.262189135419, -1.184228350040, -1.053263771491};*/
};

double tandemProb::objfun(const std::vector<double>& x) const {
    double obj = 0;

	MGA_DSM(x, mgadsm, obj);

	//evaluating the mass from the dvs
	double rE[3];
	double vE[3];
	Planet_Ephemerides_Analytical (x[0],3,rE,vE);
	double VINFE = x[1];
	double udir = x[2];
	double vdir = x[3];
	double vtemp[3];
	vtemp[0]= rE[1]*vE[2]-rE[2]*vE[1];
	vtemp[1]= rE[2]*vE[0]-rE[0]*vE[2];
	vtemp[2]= rE[0]*vE[1]-rE[1]*vE[0];
	double iP1[3];
	double normvE=sqrt(vE[0]*vE[0]+vE[1]*vE[1]+vE[2]*vE[2]);
	iP1[0]=	vE[0]/normvE;
	iP1[1]=	vE[1]/normvE;
	iP1[2]=	vE[2]/normvE;
	double zP1[3];
	double normvtemp=sqrt(vtemp[0]*vtemp[0]+vtemp[1]*vtemp[1]+vtemp[2]*vtemp[2]);
	zP1[0]= vtemp[0]/normvtemp;
	zP1[1]= vtemp[1]/normvtemp;
	zP1[2]= vtemp[2]/normvtemp;
	double jP1[3];
	jP1[0]= zP1[1]*iP1[2]-zP1[2]*iP1[1];
	jP1[1]= zP1[2]*iP1[0]-zP1[0]*iP1[2];
	jP1[2]= zP1[0]*iP1[1]-zP1[1]*iP1[0];
	double theta=2*M_PI*udir; 		//See Picking a Point on a Sphere
	double phi=acos(2*vdir-1)-M_PI/2; //In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation
	double vinf[3];
	vinf[0]=VINFE*(cos(theta)*cos(phi)*iP1[0]+sin(theta)*cos(phi)*jP1[0]+sin(phi)*zP1[0]);
	vinf[1]=VINFE*(cos(theta)*cos(phi)*iP1[1]+sin(theta)*cos(phi)*jP1[1]+sin(phi)*zP1[1]);
	vinf[2]=VINFE*(cos(theta)*cos(phi)*iP1[2]+sin(theta)*cos(phi)*jP1[2]+sin(phi)*zP1[2]);
	//We rotate it to the equatorial plane
	ecl2equ(vinf,vinf);
	//And we find the declination in degrees
	double normvinf=sqrt(vinf[0]*vinf[0]+vinf[1]*vinf[1]+vinf[2]*vinf[2]);
	double sindelta = vinf[2] / normvinf;
	double declination = asin(sindelta)/M_PI*180;

	//double m_initial = SoyuzFregat(VINFE,declination);
	double m_initial = Atlas501(VINFE,declination);

	//We evaluate the final mass
	double Isp = 312;
	double g0 = 9.80665;
	double sumDVvec=0;
	//double totaltime=x[4]+x[5]+x[6]+x[7];
	for(unsigned int i=1;i<=5;i++) {
		sumDVvec=sumDVvec+mgadsm.DV[i];
	}
	double m_final;
	sumDVvec=sumDVvec+0.165; //losses for 3 swgbys + insertion

	m_final = m_initial * exp(-sumDVvec/Isp/g0*1000);

	return -log(m_final);
}


const double cassini1Prob::lb[6] = {-1000, 30,100,30,400,1000};
const double cassini1Prob::ub[6] = {0,400,470,400,2000,6000};

cassini1Prob::cassini1Prob():GOProblem(6,lb,ub) {};

double cassini1Prob::objfun(const std::vector<double>& x) const {
    return cassini1(x);
}

const double gtoc1Prob::lb[8] = {3000,14,14,14,14,100,366,300};
const double gtoc1Prob::ub[8] = {10000,2000,2000,2000,2000,9000,9000,9000};

gtoc1Prob::gtoc1Prob():GOProblem(8,lb,ub) {};

double gtoc1Prob::objfun(const std::vector<double>& x) const {
    return gtoc1(x);
}

const double cassini2Prob::lb[22] = {-750, 3, 0, 0, 100, 100, 30, 400, 800, 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.15, 1.7, -M_PI, -M_PI, -M_PI, -M_PI};
const double cassini2Prob::ub[22] = {780,  5, 1, 1, 400, 500, 300, 1600, 2200, 0.9, 0.9, 0.9, 0.9, 0.9, 6, 6, 6.5, 291, M_PI, M_PI, M_PI, M_PI};
const int cassini2Prob::sequence[6] = {3, 2, 2, 3, 5, 6};

cassini2Prob::cassini2Prob():GOProblem(22,lb,ub),mgadsm(total_DV_rndv,sequence,6,0,0,0,0,0) {};

double cassini2Prob::objfun(const std::vector<double>& x) const {
   	double obj = 0;
	MGA_DSM(x, mgadsm,
			obj);
	return obj;
}

const double rosettaProb::lb[22] = {1460, 3, 0, 0, 300, 150, 150, 300, 700, 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI};
const double rosettaProb::ub[22] = {780,  5, 1, 1, 500, 800, 800, 800, 1850, 0.9, 0.9, 0.9, 0.9, 0.9, 9, 9, 9, 9, M_PI, M_PI, M_PI, M_PI};
const int rosettaProb::sequence[6] = {3, 3, 4, 3, 3, 10};

rosettaProb::rosettaProb():GOProblem(22,lb,ub),mgadsm(rndv,sequence,6,0,0,0,0,0) {
	//MGA_DSM stuff
	mgadsm.asteroid.keplerian[0] = 3.50294972836275;
	mgadsm.asteroid.keplerian[1] = 0.6319356;
	mgadsm.asteroid.keplerian[2] = 7.12723;
	mgadsm.asteroid.keplerian[3] = 50.92302;
	mgadsm.asteroid.keplerian[4] = 11.36788;
	mgadsm.asteroid.keplerian[5] = 0.0;
	mgadsm.asteroid.epoch = 52504.23754000012;
	mgadsm.asteroid.mu = 0.0;
};

double rosettaProb::objfun(const std::vector<double>& x) const {
	double obj = 0;
	MGA_DSM(x, mgadsm,
			obj);
	return obj;
}

const double sagasProb::lb[12] = {7000, 0, 0, 0, 50, 300, 0.01, 0.01, 1.05, 8, -M_PI, -M_PI};
const double sagasProb::ub[12] = {9100, 7, 1, 1, 2000, 2000, 0.9, 0.9 ,7 ,500 ,M_PI,M_PI};
const int sagasProb::sequence[3] = {3,3,5};

sagasProb::sagasProb():GOProblem(12,lb,ub),mgadsm(time2AUs,sequence,3,50.0,6.782,1.782,0,0) {};

double sagasProb::objfun(const std::vector<double>& x)  const{
   	double obj = 0;
	MGA_DSM(x, mgadsm,
			obj);
	return obj;

}
