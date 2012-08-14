// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //


#include "misc4Tandem.h"
#include "math.h"
#include <vector>

//The following data refer to the Launcher Atlas501 performances as given by NASA
 static const double x_atl[9] = {2.5,3,3.5,4,4.5,5,5.5,5.75,6};
 static const double y_atl[15] = {-90,-40, -30,-29,-28.5, -20, -10, 0, 10, 20,28.5,29,30, 40, 90};
 static const double data_atl[15][9] = {
	{1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1},
	{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
	{10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0},
	{1160,1100,1010,930,830,740,630,590,550},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{2335.0	,2195.0,2035.0,1865.0,1675.0,1480.0,1275.0,1175.0,1075.0},
	{1160,1100,1010,930,830,740,630,590,550},
	{10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0},
	{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
	{1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1}
};

//The following data refer to the Launcher Soyuz-Fregat performances as given by ESOC. The data here and consider
// an elaborated five impulse strategy to exploit the launcher performances as much as possible.
static const double x_sf[5] = {1,2,3,4,5};
static const double y_sf[15] = {-90, -65, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 65, 90};
static const double data_sf[15][5] = {
	{1e-3,1e-3,1e-3,1e-3,1e-3},
	{100.00000,100.00000,100.00000,100.00000,100.00000},
	{1830.50000,1815.90000,1737.70000,1588.00000,1344.30000},
	{1910.80000,1901.90000,1819.00000,1636.40000,1369.30000},
	{2001.80000,1995.30000,1891.30000,1673.90000,1391.90000},
	{2108.80000,2088.60000,1947.90000,1708.00000,1409.50000},
	{2204.00000,2167.30000,1995.50000,1734.50000,1419.60000},
	{2270.80000,2205.80000,2013.60000,1745.10000,1435.20000},
	{2204.70000,2133.60000,1965.40000,1712.80000,1413.60000},
	{2087.90000,2060.60000,1917.70000,1681.10000,1392.50000},
	{1979.17000,1975.40000,1866.50000,1649.00000,1371.70000},
	{1886.90000,1882.20000,1801.00000,1614.60000,1350.50000},
	{1805.90000,1796.00000,1722.70000,1571.60000,1327.60000},
	{100.00000,100.00000,100.00000,100.00000,100.00000},
	{1e-3,1e-3,1e-3,1e-3,1e-3}
	};

int xant(const double &x) {
	int i;
	for(i=1; i<4; i++) {
		if (x_sf[i]>x) break;
	}
	return i-1;
}

int yant(const double &y) {
	int i;
	for(i=1; i<14; i++) {
		if (y_sf[i]>y) break;
	}
	return i-1;
}

int xantA5(const double &x) {
	int i;
	for(i=1; i<8; i++) {
		if (x_atl[i]>x) break;
	}
	return i-1;
}

int yantA5(const double &y) {
	int i;
	for(i=1; i<14; i++) {
		if (y_atl[i]>y) break;
	}
	return i-1;
}



double interp2SF(const double &VINF, const double &declination) {

	double retval;
	int v_indx=xant(VINF);
	int dec_indx=yant(declination);
	if (fabs(declination)>=90) return 1e-3;
	if ((VINF<1)||(VINF>5)) return 1e-3;

	double dx = x_sf[v_indx+1]-x_sf[v_indx];
    double dydx = dx * (y_sf[dec_indx+1]-y_sf[dec_indx]);
	retval =  data_sf[dec_indx][v_indx]        /  dydx * (x_sf[v_indx+1]-VINF) * (y_sf[dec_indx+1]-declination);
    retval += data_sf[dec_indx][v_indx+1]      /  dydx * (VINF - x_sf[v_indx]) * (y_sf[dec_indx+1]-declination);
    retval += data_sf[dec_indx+1][v_indx]      /  dydx * (x_sf[v_indx+1]-VINF) * (declination - y_sf[dec_indx]);
	retval += data_sf[dec_indx+1][v_indx+1]    /  dydx * (VINF - x_sf[v_indx]) * (declination - y_sf[dec_indx]);
	return retval;
}


double interp2A5(const double &VINF, const double &declination) {

    double retval;
    int v_indx=xantA5(VINF);
    int dec_indx=yantA5(declination);
    if ((VINF<2.5)||(VINF>6)) return 1e-1;
    if (fabs(declination)>=90) return 1e-1;

    double dx = x_atl[v_indx+1]-x_atl[v_indx];
    double dydx = dx * (y_atl[dec_indx+1]-y_atl[dec_indx]);
	retval =  data_atl[dec_indx][v_indx]        /  dydx * (x_atl[v_indx+1]-VINF) * (y_atl[dec_indx+1]-declination);
    retval += data_atl[dec_indx][v_indx+1]      /  dydx * (VINF - x_atl[v_indx]) * (y_atl[dec_indx+1]-declination);
    retval += data_atl[dec_indx+1][v_indx]      /  dydx * (x_atl[v_indx+1]-VINF) * (declination - y_atl[dec_indx]);
	retval += data_atl[dec_indx+1][v_indx+1]    /  dydx * (VINF - x_atl[v_indx]) * (declination - y_atl[dec_indx]);
	return retval;
}


double SoyuzFregat (const double &VINF, const double &declination) {
//This function returns the mass that a Soyuz-Fregat launcher can inject
//into a given escape velocity and asymptote declination. The data here
//refer to ESOC WP-521 and consider an elaborated five impulse strategy to
//exploit the launcher performances as much as possible.
	return interp2SF(VINF,declination);
}


double Atlas501 (const double &VINF, const double &declination) {
//This function returns the mass that a Atlas501 Launcher
	return interp2A5(VINF,declination);
}
void ecl2equ (double (&ecl)[3],double (&equ)[3]){
	static const double incl=0.409072975;
	double temp[3];
	temp[0]=ecl[0];
	temp[1]=ecl[1];
	temp[2]=ecl[2];
    equ[0]=temp[0];
	equ[1]=temp[1]*cos(incl)-temp[2]*sin(incl);
	equ[2]=temp[1]*sin(incl)+temp[2]*cos(incl);
}


