/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<cmath>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base_stochastic.h"
#include "spheres.h"

static const int nr_input = 8;
static const int nr_output = 3;
static const int nr_spheres = 3;
static const int nr_eq = 9;
static const double target_distance2 = 0.25;
static const double target_distance = 0.5;


double norm(double v[3]) {
	return(std::sqrt(v[0]*v[0] + v[1]*v[1] +v[2]*v[2]));
}

double norm2(double v[3]) {
	return(v[0]*v[0] + v[1]*v[1] +v[2]*v[2]);
}

void versify(double v[3]) {
	double V = norm(v);
	if  (V<=0) printf("Warning zero norm vector passed to versify");
	v[0] = v[0] / V;
	v[1] = v[1] / V;
	v[2] = v[2] / V;
}

void cross_vector(double vout[3], const double vin1[3], const double vin2[3]) {
	vout[0] = vin1[1]*vin2[2] - vin1[2]*vin2[1];
	vout[1] = vin1[2]*vin2[0] - vin1[0]*vin2[2];
	vout[2] = vin1[0]*vin2[1] - vin1[1]*vin2[0];
}

void q2C(double Cout[3][3], const double q[4]) {
	Cout[0][0] = 1. - 2. * (q[1]*q[1]+q[2]*q[2]);
	Cout[1][1] = 1. - 2. * (q[0]*q[0]+q[2]*q[2]);
	Cout[2][2] = 1. - 2. * (q[1]*q[1]+q[0]*q[0]);
	Cout[0][1] = 2. * (q[0]*q[1] + q[2]*q[3]);
	Cout[0][2] = 2. * (q[0]*q[2] - q[1]*q[3]);
	Cout[1][0] = 2. * (q[0]*q[1] - q[2]*q[3]);
	Cout[1][2] = 2. * (q[1]*q[2] + q[0]*q[3]);
	Cout[2][0] = 2. * (q[2]*q[0] + q[1]*q[3]);
	Cout[2][1] = 2. * (q[2]*q[1] - q[0]*q[3]);
}

void matrix_transformation(double vout[3], const double R[3][3]) {
	int i;
	double vout_cp[3];
	vout_cp[0] = vout[0]; vout_cp[1] = vout[1]; vout_cp[2] = vout[2];
	for (i=0;i<3;++i) {
		vout[i] = R[i][0] * vout_cp[0] + R[i][1] * vout_cp[1] + R[i][2] * vout_cp[2];
	}
}

void matrix_inv_transformation(double vout[3], const double R[3][3]) {
	int i;
	double vout_cp[3];
	vout_cp[0] = vout[0]; vout_cp[1] = vout[1]; vout_cp[2] = vout[2];
	for (i=0;i<3;++i) {
		vout[i] = R[0][i] * vout_cp[0] + R[1][i] * vout_cp[1] + R[2][i] * vout_cp[2];
	}
}

void quat_dyn(double dq[4], const double q[4], const double w[3]) {
	dq[0] = 0.5 * (q[1]*w[2] - q[2]*w[1] + q[3]*w[0]);
	dq[1] = 0.5 * (q[2]*w[0] - q[0]*w[2] + q[3]*w[1]);
	dq[2] = 0.5 * (q[0]*w[1] - q[1]*w[0] + q[3]*w[2]);
	dq[3] = -0.5 * (q[0]*w[0] + q[1]*w[1] + q[2]*w[2]);
}



namespace pagmo { namespace problem {

spheres::spheres(int n_evaluations, int n_hidden_neurons,
		 double numerical_precision, unsigned int seed) :
	base_stochastic((nr_input + 1) * n_hidden_neurons + (n_hidden_neurons + 1) * nr_output, seed),
	m_ffnn(nr_input,n_hidden_neurons,nr_output), m_n_evaluations(n_evaluations),
	m_n_hidden_neurons(n_hidden_neurons), m_numerical_precision(numerical_precision),
	m_ic(nr_eq) {
	// Here we set the bounds for the problem decision vector, i.e. the nn weights
	set_lb(-1);
	set_ub(1);
	// We then instantiate the ode integrator system using gsl
	m_sys = {ode_func,NULL,nr_eq,&m_ffnn};
	m_gsl_drv_pntr = gsl_odeiv2_driver_alloc_y_new(&m_sys, gsl_odeiv2_step_rk8pd, 1e-6,m_numerical_precision,0.0);
}

spheres::spheres(const spheres &other):
	base_stochastic(other),
	m_ffnn(other.m_ffnn),
	m_n_evaluations(other.m_n_evaluations),m_n_hidden_neurons(other.m_n_hidden_neurons),
	m_numerical_precision(other.m_numerical_precision),m_ic(other.m_ic)
{
	// Here we set the bounds for the problem decision vector, i.e. the nn weights
	m_sys = {ode_func,NULL,nr_eq,&m_ffnn};
	m_gsl_drv_pntr = gsl_odeiv2_driver_alloc_y_new(&m_sys, gsl_odeiv2_step_rk8pd, 1e-6,m_numerical_precision,0.0);
}

spheres::~spheres(){
	gsl_odeiv2_driver_free(m_gsl_drv_pntr);
}

/// Clone method.
base_ptr spheres::clone() const
{
	return base_ptr(new spheres(*this));
}

//This function evaluates the fitness of a given spheres configuration ....
double spheres::single_fitness( const std::vector<double> &y, const ffnn& neural_net) const {
	double	fit = 0.0;
	double	context[8], vel_f[3];
	int	k;

	// for each sphere



	for( int i = 0; i < nr_spheres; i++ ){	// i - is the sphere counter 0 .. 1 .. 2 ..
		k = 0;

		// we now load in context the perceived data (as decoded from the world state y)
		for( int n = 1; n <= nr_spheres - 1; n++ ){				// consider the vector from each other sphere
			for( int j = 0; j < 3; j++ ){				// consider each component from the vectors
				context[k++] = y[i*3 + j] - y[ (i*3 + j + n*3) % 9 ];
			}
		}

		// context now contains the relative position vectors (6 components) in the absolute frame
		// we write, on the last two components of context, the norms of these relative positions
		context[6] = context[0]*context[0] + context[1]*context[1] + context[2]*context[2];
		context[7] = context[3]*context[3] + context[4]*context[4] + context[5]*context[5];

		//We evaluate the output from the neural net
		neural_net.eval(vel_f, context);
		vel_f[0] = vel_f[0] * 2 * 0.3 - 0.3;
		vel_f[1] = vel_f[1] * 2 * 0.3 - 0.3;
		vel_f[2] = vel_f[2] * 2 * 0.3 - 0.3;

		//first add |L² - R²| to the fitness
		double temp = std::abs(target_distance2 - context[6]);
		fit += temp;
		temp = std::abs(target_distance2 - context[7]);
		fit += temp;

		// and then the final velocity
		fit += norm2(vel_f);
	}
	return (fit / 2);
}

int spheres::ode_func( double t, const double y[], double f[], void *params ) {

	// Here we recover the neural network
	ffnn	*ptr_ffnn = (ffnn*)params;

	// The fixed-size vector context represent the sensory data perceived from each sphere. These are
	// the body axis components of the relative positions of the other spheres, and their modules
	double  context[nr_input];
	double  out[nr_output];

	// Here are some counters
	int  k;


	for( int i = 0; i < nr_spheres; i++ ){	// i - is the sphere counter 0 .. 1 .. 2 ..
		k = 0;

		// we now load in context the perceived data (as decoded from the world state y)
		for( int n = 1; n <= nr_spheres - 1; n++ ){				// consider the vector from each other sphere
			for( int j = 0; j < 3; j++ ){				// consider each component from the vectors
				context[k++] = y[i*3 + j] - y[ (i*3 + j + n*3) % 9 ];
			}
		}

		// context now contains the relative position vectors (6 components) in the absolute frame
		// we write, on the last two components of context, the norms of these relative positions
		context[6] = context[0]*context[0] + context[1]*context[1] + context[2]*context[2];
		context[7] = context[3]*context[3] + context[4]*context[4] + context[5]*context[5];

		//We evaluate the output from the neural net
		ptr_ffnn->eval(out, context);

		//Here we set the dynamics transforming the nn output [0,1] in desired velocities [-0/3,0.3]
		f[i*3] = out[0] * 0.3*2 - 0.3;
		f[i*3+1] = out[1] * 0.3*2 - 0.3;
		f[i*3+2] = out[2] * 0.3*2 - 0.3;

	}
	return GSL_SUCCESS;
}

spheres::ffnn::ffnn(const unsigned int n_inputs, const unsigned int n_hidden,const unsigned int n_outputs) :
	m_n_inputs(n_inputs), m_n_hidden(n_hidden), m_n_outputs(n_outputs),
	m_weights((n_inputs + 1) * n_hidden + (n_hidden + 1) * n_outputs), m_hidden(n_hidden)
{}

void spheres::ffnn::set_weights(const std::vector<double> &weights) {
	m_weights = weights;
}

void spheres::ffnn::eval(double out[], const double in[]) const {
	// Offset for the weights to the output nodes
	unsigned int offset = m_n_hidden * (m_n_inputs + 1);

	// -- PROCESS CONTEXT USING THE NEURAL NETWORK --
	for( int i = 0; i < m_n_hidden; i++ ){
		// Set the bias (the first weight to the i'th hidden node)
		m_hidden[i] = m_weights[i * (m_n_inputs + 1)];

		for( int j = 0; j < m_n_inputs; j++ ){
			// Compute the weight number
			int ji = i * (m_n_inputs + 1) + (j + 1);
			// Add the weighted input
			m_hidden[i] += m_weights[ji] * in[j];
		}

		// Apply the transfer function (a sigmoid with output in [0,1])
		m_hidden[i] = 1.0 / ( 1 + std::exp( -m_hidden[i] ));
	}

	// generate values for the output nodes
	for( int i = 0; i < m_n_outputs; i++ ){
		// add the bias (weighted by the first weight to the i^th output node
		out[i] = m_weights[offset + i * (m_n_hidden + 1)];

		for( int j = 0; j < m_n_hidden; j++ ){
			// compute the weight number
			int ji = offset + i * (m_n_hidden + 1) + (j + 1);
			// add the weighted input
			out[i] += m_weights[ji] * m_hidden[j];
		}

		out[i] = 1.0 / ( 1 + std::exp( -out[i] ));
	}
}

void spheres::objfun_impl(fitness_vector &f, const decision_vector &x) const {
	f[0]=0;
	// Make sure the pseudorandom sequence will always be the same
	m_drng.seed(m_seed);
	// Set the ffnn weights
	m_ffnn.set_weights(x);
	// Loop over the number of repetitions
	for (int count=0;count<m_n_evaluations;++count) {
		// Creates the initial conditions at random
		// Position starts in a [-1,1] box
		for (int i=0; i<9; ++i) {
			m_ic[i] = (m_drng()*2 - 1);
		}
		// Integrate the system
		double t0 = 0.0;
		double tf = 50.0;
		//gsl_odeiv2_driver_set_hmin (m_gsl_drv_pntr, 1e-6);
		int status = gsl_odeiv2_driver_apply( m_gsl_drv_pntr, &t0, tf, &m_ic[0] );
		// Not sure if this help or what it does ....
		gsl_odeiv2_driver_reset (m_gsl_drv_pntr);
		if( status != GSL_SUCCESS ){
			printf ("ERROR: gsl_odeiv2_driver_apply returned value = %d\n", status);
			break;
		}
		f[0] += single_fitness(m_ic,m_ffnn);

	}
	f[0] /= m_n_evaluations;
}

bool my_sort_function (std::vector<double> i,std::vector<double> j) { return (i[9] < j[9]); }

std::vector<std::vector<double> > spheres::post_evaluate(const decision_vector & x, int N, unsigned int seed) const {
	std::vector<double> one_row(10,0.0);
	std::vector<std::vector<double> > ret(N,one_row);
	// Make sure the pseudorandom sequence will always be the same
	m_drng.seed(seed);
	// Set the ffnn weights
	m_ffnn.set_weights(x);
	// Loop over the number of repetitions
	for (int count=0;count<N;++count) {
		// Creates the initial conditions at random
		// Position starts in a [-1,1] box
		for (int i=0; i<9; ++i) {
			m_ic[i] = (m_drng()*2 - 1);
			one_row[i] = m_ic[i];
		}

		// Integrate the system
		double t0 = 0.0;
		double tf = 50.0;
		//gsl_odeiv2_driver_set_hmin (m_gsl_drv_pntr, 1e-6);
		int status = gsl_odeiv2_driver_apply( m_gsl_drv_pntr, &t0, tf, &m_ic[0] );
		// Not sure if this help or what it does ....
		//gsl_odeiv2_driver_reset (m_gsl_drv_pntr);
		if( status != GSL_SUCCESS ){
			printf ("ERROR: gsl_odeiv2_driver_apply returned value = %d\n", status);
			break;
		}
		one_row[9] = single_fitness(m_ic,m_ffnn);
		ret[count] = one_row;
	}
	// sorting by fitness
	std::sort (ret.begin(), ret.end(), my_sort_function);
	return ( ret );
}

std::vector<std::vector<double> > spheres::simulate(const decision_vector &x, const std::vector<double> &ic, int N) const {
	std::vector<double> y0(ic);
	std::vector<double> one_row(10,0.0);
	std::vector<std::vector<double> > ret;
	// Set the ffnn weights
	m_ffnn.set_weights(x);
	// Integrate the system
	double ti, t0=0;
	double tf = 70.0;

	// pushing back the initial conditions
	one_row[0] = 0.0;
	std::copy(y0.begin(),y0.end(),one_row.begin()+1);
	ret.push_back(one_row);

	for( int i = 1; i <= N; i++ ){
		ti = i * tf / N;
		int status = gsl_odeiv2_driver_apply( m_gsl_drv_pntr, &t0, ti, &y0[0] );
		//pushing_back the result
		one_row[0] = ti;
		std::copy(y0.begin(),y0.end(),one_row.begin()+1);
		ret.push_back(one_row);
		if( status != GSL_SUCCESS ){
			printf ("ERROR: gsl_odeiv2_driver_apply returned value = %d\n", status);
			break;
		}
	}
	//Not sure if this help or what it does ....
	//gsl_odeiv2_driver_reset (m_gsl_drv_pntr);
	return ( ret );
}
} //namespace problem
} //namespace pagmo

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::spheres);
