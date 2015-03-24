/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
#include<algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base_stochastic.h"
#include "spheres.h"

static const int nr_input = 8;
static const int nr_output = 3;
static const int nr_spheres = 3;
static const int nr_eq = 9;

static double norm2(double v[3]) {
	return(v[0]*v[0] + v[1]*v[1] +v[2]*v[2]);
}

namespace pagmo { namespace problem {

spheres::spheres(int n_evaluations, int n_hidden_neurons,
		 double numerical_precision, unsigned int seed, bool symmetric, double sim_time, const std::vector<double>& sides) :
	base_stochastic((nr_input/(int(symmetric)+1) + 1) * n_hidden_neurons + (n_hidden_neurons + 1) * nr_output, seed),
	m_ffnn(nr_input,n_hidden_neurons,nr_output), m_n_evaluations(n_evaluations),
	m_n_hidden_neurons(n_hidden_neurons), m_numerical_precision(numerical_precision),
	m_ic(nr_eq), m_symm(symmetric), m_sim_time(sim_time), m_sides(sides) {
	// Here we set the bounds for the problem decision vector, i.e. the nn weights
	set_lb(-1);
	set_ub(1);
	// We then instantiate the ode integrator system using gsl
	gsl_odeiv2_system sys = {ode_func,NULL,nr_eq,&m_ffnn};
	m_sys = sys;
	m_gsl_drv_pntr = gsl_odeiv2_driver_alloc_y_new(&m_sys, gsl_odeiv2_step_rk8pd, 1e-6,m_numerical_precision,0.0);
	// And make sure the three sides are ordered and squared here
	std::sort(m_sides.begin(),m_sides.end());
	m_sides[0]*=m_sides[0];	m_sides[1]*=m_sides[1];	m_sides[2]*=m_sides[2];
}

spheres::spheres(const spheres &other):
	base_stochastic(other),
	m_ffnn(other.m_ffnn),
	m_n_evaluations(other.m_n_evaluations),m_n_hidden_neurons(other.m_n_hidden_neurons),
	m_numerical_precision(other.m_numerical_precision),m_ic(other.m_ic), m_symm(other.m_symm), m_sim_time(other.m_sim_time),m_sides(other.m_sides)
{
	// Here we set the bounds for the problem decision vector, i.e. the nn weights
	gsl_odeiv2_system sys = {ode_func,NULL,nr_eq,&m_ffnn};
	m_sys = sys;
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
		for( int n = 1; n <= nr_spheres - 1; n++ ){		// consider the vector from each other sphere
			for( int j = 0; j < 3; j++ ){			// consider each component from the vectors
				context[k++] = y[i*3 + j] - y[ (i*3 + j + n*3) % 9 ];
			}
		}

		// context now contains the relative position vectors (6 components) in the absolute frame
		// we write, on the last two components of context, the norms2 of these relative positions
		context[6] = context[0]*context[0] + context[1]*context[1] + context[2]*context[2];
		context[7] = context[3]*context[3] + context[4]*context[4] + context[5]*context[5];

		//We evaluate the output from the neural net
		neural_net.eval(vel_f, context);
		vel_f[0] = vel_f[0] * 2 * 0.3 - 0.3;
		vel_f[1] = vel_f[1] * 2 * 0.3 - 0.3;
		vel_f[2] = vel_f[2] * 2 * 0.3 - 0.3;

		//and we add the final velocity violation (two times since we will divide by two later)
		fit += 2 *(norm2(vel_f));

		//and we keep them within the box [-1 1]
		//if (std::abs(y[i*3]) > 1)   fit += std::abs(y[i*3]);
		//if (std::abs(y[i*3+1]) > 1) fit += std::abs(y[i*3+1]);
		//if (std::abs(y[i*3+2]) > 1) fit += std::abs(y[i*3+2]);
	}
	fit = fit/2;
	
	// Here we compute the violation from the desired velocity
#define SRT_PAIR(a,b) if ((b)<(a)) {(tmp=a);(a=b);(b=tmp);}
	double tmp = 0;
	double a2 = (y[0]-y[3])*(y[0]-y[3]) + (y[1]-y[4])*(y[1]-y[4]) + (y[2]-y[5])*(y[2]-y[5]);
	double b2 = (y[0]-y[6])*(y[0]-y[6]) + (y[1]-y[7])*(y[1]-y[7]) + (y[2]-y[8])*(y[2]-y[8]);
	double c2 = (y[6]-y[3])*(y[6]-y[3]) + (y[7]-y[4])*(y[7]-y[4]) + (y[8]-y[5])*(y[8]-y[5]);

	SRT_PAIR(a2,b2)
	SRT_PAIR(b2,c2)
	SRT_PAIR(a2,b2)
#undef SRT_PAIR

	fit += (a2-m_sides[0])*(a2-m_sides[0])+(b2-m_sides[1])*(b2-m_sides[1])+(c2-m_sides[2])*(c2-m_sides[2]);
	//fit += std::abs(a2-m_sides[0])+std::abs(b2-m_sides[1])+std::abs(c2-m_sides[2]);
	return fit;
}

int spheres::ode_func( double t, const double y[], double f[], void *params ) {
	(void)t;
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
		f[i*3] = out[0] * 0.3 * 2 - 0.3;
		f[i*3+1] = out[1] * 0.3 * 2 - 0.3;
		f[i*3+2] = out[2] * 0.3 * 2 - 0.3;

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
	for( unsigned int i = 0; i < m_n_hidden; i++ ){
		// Set the bias (the first weight to the i'th hidden node)
		m_hidden[i] = m_weights[i * (m_n_inputs + 1)];

		for( unsigned int  j = 0; j < m_n_inputs; j++ ){
			// Compute the weight number
			int ji = i * (m_n_inputs + 1) + (j + 1);
			// Add the weighted input
			m_hidden[i] += m_weights[ji] * in[j];
		}

		// Apply the transfer function (a sigmoid with output in [0,1])
		m_hidden[i] = 1.0 / ( 1 + std::exp( -m_hidden[i] ));
	}

	// generate values for the output nodes
	for( unsigned int  i = 0; i < m_n_outputs; i++ ){
		// add the bias (weighted by the first weight to the i^th output node
		out[i] = m_weights[offset + i * (m_n_hidden + 1)];

		for( unsigned int  j = 0; j < m_n_hidden; j++ ){
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
	// Set the ffnn weights from x, by accounting for symmetries in neurons weights
	set_nn_weights(x);
	// Loop over the number of repetitions
	for (int count=0;count<m_n_evaluations;++count) {
		// Creates the initial conditions at random
		// Positions starts in a [-1,1] box
		for (int i=0; i<6; ++i) {
			m_ic[i] = (m_drng()*2 - 1);
		}

		// Centered around the origin
		m_ic[6] = - (m_ic[0] + m_ic[3]);
		m_ic[7] = - (m_ic[1] + m_ic[4]);
		m_ic[8] = - (m_ic[2] + m_ic[5]);

		// Integrate the system
		double t0 = 0.0;
		double tf = m_sim_time;
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

static bool my_sort_function (std::vector<double> i,std::vector<double> j) { return (i[9] < j[9]); }

std::vector<std::vector<double> > spheres::post_evaluate(const decision_vector & x, int N, unsigned int seed) const {
	std::vector<double> one_row(10,0.0);
	std::vector<std::vector<double> > ret(N,one_row);
	// Make sure the pseudorandom sequence will always be the same
	m_drng.seed(seed);
	// Set the ffnn weights
	set_nn_weights(x);
	// Loop over the number of repetitions
	for (int count=0;count<N;++count) {
		// Creates the initial conditions at random

		// Position starts in a [-1,1] box (evolution is in [-2,2])
		for (int i=0; i<6; ++i) {
		m_ic[i] = (m_drng()*2 - 1);
		}
		// Centered around the origin
		m_ic[6] = - (m_ic[0] + m_ic[3]);
		m_ic[7] = - (m_ic[1] + m_ic[4]);
		m_ic[8] = - (m_ic[2] + m_ic[5]);

		for (int i=0; i<9; ++i) {
			one_row[i] = m_ic[i];
		}

		// Integrate the system
		double t0 = 0.0;
		double tf = m_sim_time;
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

void spheres::set_nn_weights(const decision_vector &x) const {
	if (m_symm) { //symmetric weigths activated
		int w = 0;
		for(unsigned int h = 0; h < m_ffnn.m_n_hidden; h++)
		{
			int start_index = h * 5; // (nr_input/2+1)
			// bias, dx1, dy1, dz1
			for(int j = 0; j < 4; j++)
			{
				m_ffnn.m_weights[w] = x[start_index+j]; w++;
			}
			// dx2, dy2, dz2
			for(int j = 1; j <= 3; j++)
			{
				m_ffnn.m_weights[w] = x[start_index+j]; w++;
			}
			// distance 1
			m_ffnn.m_weights[w] = x[start_index+4]; w++;
			// distance 2
			m_ffnn.m_weights[w] = x[start_index+4]; w++;
		}
		int ind = 0;
		for(unsigned int ww = w; ww < m_ffnn.m_weights.size(); ww++)
		{
			m_ffnn.m_weights[ww] = x[(nr_input/2+1)*m_ffnn.m_n_hidden+ind];
			ind++;
		}
	} else {//no symmetric weights
		m_ffnn.m_weights = x;
	}
}

std::vector<std::vector<double> > spheres::simulate(const decision_vector &x, const std::vector<double> &ic, int N) const {
	std::vector<double> y0(ic);
	std::vector<double> one_row(10,0.0);
	std::vector<std::vector<double> > ret;
	// Set the ffnn weights
	set_nn_weights(x);
	// Integrate the system
	double ti, t0=0;
	double tf = m_sim_time;

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

std::vector<double> spheres::get_nn_weights(decision_vector x) const {
	set_nn_weights(x);
	return m_ffnn.m_weights;
}

std::string spheres::get_name() const
{
	return "MIT SPHERES - Neurocontroller Evolution";
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight.
 */
std::string spheres::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\n\tSample Size: " << m_n_evaluations << '\n';
	oss << "\tHidden Neurons: " << m_n_hidden_neurons << '\n';
	oss << "\tODE precision: " << m_numerical_precision << '\n';
	oss << "\tSeed: " << m_seed << '\n';
	oss << "\tSymmetric Weights: " << m_symm << '\n';
	oss << "\tSimulation time: " << m_sim_time << '\n';
	oss << "\tTriangle sides (squared): " << m_sides << '\n';
	return oss.str();
}

} //namespace problem
} //namespace pagmo

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::spheres)
