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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/round.hpp>
#include <string>
#include <vector>
#include <algorithm>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base.h"
#include "sga_gray.h"
#include "../problem/base_stochastic.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability (of each individual)
 * @param[in] m Mutation probability (of each encoded bit)
 * @param[in] elitism: the best individual is reinserted in the population each elitism generations.
 * If a elitism of 0 is set, then the elitism is avoided.
 * @param[in] mut Mutation type. sga_gray::mutation::UNIFORM
 * @param[in] sel Selection type. One of sga_gray::selection::BEST20, sga_gray::selection::ROULETTE
 * @param[in] cro Crossover type. One of sga_gray::crossover::SINGLE_POINT
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1]\f$, mutation probability is not \f$ \in [0,1]\f$,
 * elitism is <=0
 *
 */
sga_gray::sga_gray(int gen, const double &cr, const double &m, int elitism, mutation::type mut, selection::type sel, crossover::type cro)
	:base(),m_gen(gen),m_cr(cr),m_m(m),m_elitism(elitism),m_mut(mut),m_sel(sel),m_cro(cro)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (cr > 1 || cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1] range");
	}
	if (m < 0 || m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (elitism < 0) {
		pagmo_throw(value_error,"elitisim must be greater or equal than zero");
	}

	m_bit_encoding = 25;
	m_max_encoding_integer = int(std::pow(2., m_bit_encoding));
}

/// Clone method.
base_ptr sga_gray::clone() const
{
	return base_ptr(new sga_gray(*this));
}

/// Evolve implementation.
/**
 * Run the simple genetic algorithm for the number of generations specified in the constructors.
 * At each improvment velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void sga_gray::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();

	//We perform some checks to determine wether the problem/population are suitable for sga_gray
	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and sga_gray is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and sga_gray is not suitable to solve it");
	}

	if (NP < 5) {
		pagmo_throw(value_error,"for sga_gray at least 5 individuals in the population are needed");
	}

	if (prob.get_i_dimension() > 1) {
		pagmo_throw(value_error,"This version of SGA gray is not compatble with discrete problems.");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy), Xnew(NP,dummy);

	std::vector< std::vector<int> > chrom_vect(NP);

	std::vector<fitness_vector > fit(NP);		//fitness

	fitness_vector bestfit;
	decision_vector bestX(D,0);

	std::vector<int> selection(NP);

	// Initialise the phenotypes and their fitness to that of the initial deme
	for (pagmo::population::size_type i = 0; i<NP; i++ ) {
		X[i] = pop.get_individual(i).cur_x;
		fit[i] = pop.get_individual(i).cur_f;
	}

	// Find the best member and store in bestX and bestfit
	double bestidx = pop.get_best_idx();
	bestX = pop.get_individual(bestidx).cur_x;
	bestfit = pop.get_individual(bestidx).cur_f;

	// Main sga_gray loop
	for(int j=0; j<m_gen; j++) {

		selection = this->selection(fit,prob);

		// Xnew stores the new selected generation of genotypes
		for(population::size_type i = 0; i < NP; i++) {
			chrom_vect[i] = this->encode(X.at(selection.at(i)), lb, ub);
		}

		this->crossover(chrom_vect);
		this->mutate(chrom_vect);

		//Xnew stores the new selected generation of genotypes
		for(population::size_type i=0; i<NP; i++) {
			Xnew[i] = this->decode(chrom_vect.at(i), lb, ub);
		}

		// If the problem is a stochastic optimization chage the seed and re-evaluate taking care to update also best and local bests
		try
		{
			//4 - Evaluate the new population (stochastic problem)
			dynamic_cast<const pagmo::problem::base_stochastic &>(prob).set_seed(m_urng());
			pop.clear(); // Removes memory based on different seeds (champion and best_x, best_f, best_c)
			
			// We re-evaluate the best individual (for elitism)
			prob.objfun(bestfit,bestX);
			// Re-evaluate wrt new seed the particle position and memory
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				// We evaluate here the new individual fitness
				prob.objfun(fit[i],Xnew[i]);
				// We update the velocity (in case coupling with PSO via archipelago)
				//dummy = Xnew[i];
				//std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				///We now set the cleared pop. cur_x is the best_x, re-evaluated with new seed.
				pop.push_back(Xnew[i]);
				//pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		catch (const std::bad_cast& e)
		{
			//4 - Evaluate the new population (deterministic problem)
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				prob.objfun(fit[i],Xnew[i]);
				dummy = Xnew[i];
				std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				//updates x and v (cache avoids to recompute the objective function and constraints)
				pop.set_x(i,Xnew[i]);
				pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		
		//5 - Reinsert best individual every m_elitism generations
		if ((m_elitism != 0) && (j % m_elitism == 0) ) {
			int worst=0;
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if ( prob.compare_fitness(fit[worst],fit[i]) ) worst=i;
			}
			Xnew[worst] = bestX;
			fit[worst] = bestfit;
			dummy = Xnew[worst];
			std::transform(dummy.begin(), dummy.end(), pop.get_individual(worst).cur_x.begin(), dummy.begin(),std::minus<double>());
			//updates x and v (cache avoids to recompute the objective function)
			pop.set_x(worst,Xnew[worst]);
			pop.set_v(worst,dummy);
		}
		X = Xnew;
	} // end of main sga_gray loop
}

/// Algorithm name
std::string sga_gray::get_name() const
{
	return "A Simple Genetic Algorithm with Gray encoding";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sga_gray::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "CR:" << m_cr << ' ';
	s << "M:" << m_m << ' ';
	s << "elitism:" << m_elitism << ' ';
	s << "mutation:";
	switch (m_mut) {
	case mutation::UNIFORM: {
		s << "UNIFORM ";
		break;
	}
	}
	s << "selection:";
	switch (m_sel) {
	case selection::BEST20: {
		s << "BEST20 ";
		break;
	}
	case selection::ROULETTE: {
		s << "ROULETTE ";
		break;
	}
	}
	s << "crossover:";
	switch (m_cro) {
	case crossover::SINGLE_POINT: {
		s << "SINGLE_POINT ";
		break;
	}
	}

	return s.str();
}

/// Select the individuals.
/**
 * Returns a vector of selected individuals positions.
 *
 * @param[in] pop_f: population fitness.
 * @param[in] prob: problem used to compare the individuals fitness.
 * @param[out] vector output of positions of selected individuals.
 */
std::vector<int> sga_gray::selection(const std::vector<fitness_vector> &pop_f, const problem::base &prob) const
{
	population::size_type NP = pop_f.size();

	std::vector<int> selection(NP);

	switch (m_sel) {
	case selection::BEST20: {
		//selects the best 20% and puts multiple copies in Xnew
		// in practice, for performance reasons, the selection, fitnessID should not be created each time
		// thus we might prefer to use a class of operators
		int tempID;
		std::vector<int> fitnessID(NP);

		for (population::size_type i=0; i<NP; i++) {
			fitnessID[i]=i;
		}
		for (population::size_type i=0; i < (NP-1); ++i) {
			for (population::size_type j=i+1; j<NP; ++j) {
				if( prob.compare_fitness(pop_f[fitnessID[j]],pop_f[fitnessID[i]]) ) {
					//if ( pop_f[fitnessID[j]] < pop_f[fitnessID[i]]) {
					//swap fitness values
					// fit[i].swap(fit[j]);
					//swap id's
					tempID = fitnessID[i];
					fitnessID[i] = fitnessID[j];
					fitnessID[j] = tempID;
				}
			}
		}

		// fitnessID now contains the position of individuals ranked from best to worst
		int best20 = NP/5;
		for (pagmo::population::size_type i=0; i<NP; ++i) {
			selection[i] = fitnessID[i % best20]; // multiple copies
		}
		break;
	}
	case selection::ROULETTE: { // ROULETTE
		std::vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);

		//We scale all fitness values from 0 (worst) to absolute value of the best fitness
		fitness_vector worstfit=pop_f[0];
		for (population::size_type i = 1; i < NP;i++) {
			if (prob.compare_fitness(worstfit,pop_f[i])) {
				worstfit=pop_f[i];
			}
		}

		for (population::size_type i = 0; i < NP; i++) {
			selectionfitness[i] = fabs(worstfit[0] - pop_f[i][0]);
		}

		// We build and normalise the cumulative sum
		cumsumTemp[0] = selectionfitness[0];
		for (population::size_type i = 1; i< NP; i++) {
			cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
		}
		for (population::size_type i = 0; i < NP; i++) {
			cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
		}

		//we throw a dice and pick up the corresponding index
		double r2;
		for (population::size_type i = 0; i < NP; i++) {
			r2 = m_drng();
			for (population::size_type j = 0; j < NP; j++) {
				if (cumsum[j] > r2) {
					selection[i]=j;
					break;
				}
			}
		}
		break;
	}
	}
	return selection;
}

/// Crossover the individuals.
/**
 * Crossover the individuals chromosomes.
 *
 * @param[in/out] pop_x: vector of chromosomes to crossover.
 */
void sga_gray::crossover(std::vector< std::vector<int> > &pop_x) const
{
	population::size_type NP = pop_x.size();

	// chromosome dimension
	fitness_vector::size_type D = pop_x.at(0).size();

	std::vector<population::size_type> mating_pool(0);
	// creates the mating pool
	for(population::size_type i=0; i<NP; i++) {
		if (m_drng() < m_m) {
			mating_pool.push_back(i);
		}
	}
	population::size_type mating_pool_size = mating_pool.size();

	if (mating_pool_size % 2 != 0) {
		// we remove or add one individual
		if(m_drng() < 0.5) {
			// add one randomly selected
			mating_pool.push_back(boost::uniform_int<int>(0,mating_pool_size - 1)(m_urng));
		} else {
			mating_pool.erase(mating_pool.begin()+boost::uniform_int<int>(0,mating_pool_size - 1)(m_urng));
		}
	}
	// update the mating pool size
	mating_pool_size = mating_pool.size();

	if(mating_pool_size == 0) {
		return;
	}

	switch (m_cro) {
	//0 - single point crossover
	case crossover::SINGLE_POINT:
	{
		// random mating of the individuals
		for (population::size_type i=0; i<mating_pool_size/2; i++) {
			// we randomly select the individuals
			std::vector<int> &member1 = pop_x[mating_pool[i*2]];
			std::vector<int> &member2 = pop_x[mating_pool[i*2 + 1]];

			// we mate them at a random position
			int position = boost::uniform_int<int>(1, D-1)(m_urng);

			std::swap_ranges(member1.begin(), member1.begin()+position, member2.begin());
		}
		break;
	}
	}
}

/// Mutate the individuals.
/**
 * Mutate the individuals chromosomes.
 *
 * @param[in/out] pop_x: vector of chromosomes to mutate.
 */
void sga_gray::mutate(std::vector< std::vector<int> > &pop_x) const
{
	const problem::base::size_type D = pop_x.at(0).size();
	const population::size_type NP = pop_x.size();

	switch (m_mut) {
	case mutation::UNIFORM: {
		for (population::size_type i=0; i<NP; i++) {
			for (problem::base::size_type j=0; j<D; j++) {
				if (m_drng() < m_m) {
					pop_x[i][j] = pop_x[i][j] == 0 ? 1:0;//boost::uniform_int<int>(0,1)(m_urng);
				}
			}
		}
		break;
	}
	}
}

/// Convert a number to its binary representation.
/**
 *
 * @param[in] number: number to convert in binary representation.
 * @param[in] lb: lower bound for the encoding.
 * @param[in] ub: upper bound for the encoding.
 * @param[out] chromosome containing the binary representation of the number.
 */
std::vector<int> sga_gray::double_to_binary(const double &number, const double &lb, const double &ub) const
{
	// initialize the vector of size m_bit_encoding
	std::vector<int> binary(m_bit_encoding, 0);

	// convert the current number into its binary representation considering the domain
	// available
	int temp_number = (number - lb) * (m_max_encoding_integer - 1) / (ub - lb) ;

	// store the binary number
	int position = 0;
	while (temp_number!=0)
	{
		if( (temp_number % 2) == 0 ) {
			binary[position] = 0;
		} else {
			binary[position] = 1;
		}
		temp_number = (int)std::floor(temp_number/2);
		position++;
	}
	// reverse the order as this algorithm provides the reverse binary reprentation
	std::reverse(binary.begin(),binary.end());

	return binary;
}

/// Convert a binary to its double representation.
/**
 * @param[in] binary: binary chromosome to convert in double.
 * @param[in] lb: lower bound for the encoding.
 * @param[in] ub: upper bound for the encoding.
 * @param[out] chromosome containing the binary representation of the number.
 */
double sga_gray::binary_to_double(const std::vector<int> &binary, const double &lb, const double &ub) const
{
	// find the representation of the binary number in the integer domain
	int temp_number = 0;
	for(int i=0; i<m_bit_encoding; i++) {
		temp_number += binary.at(m_bit_encoding - 1 - i) * std::pow(2.,i);
	}

	// rescaling back into the domain double domain

	return temp_number * (ub - lb) / (m_max_encoding_integer - 1) + lb;
}

/// Convert a gray chromosome to its binary representation.
/**
 * @param[in] gray: gray chromosome to convert in its binar representation.
 * @param[out] chromosome containing the binary representation of the gray chromosome.
 */
std::vector<int> sga_gray::gray_to_binary(const std::vector<int> &gray) const
{
	// uses the XOR table
	int length = gray.size();

	std::vector<int> binary = gray;

	// the MSB is the same
	for(int i=1; i<length; i++) {
		binary[i] = (binary[i-1]+gray[i]) % 2;
	}

	return binary;
}

/// Convert a binary chromosome to its gray representation.
/**
 * @param[in] binary: binary chromosome to convert in its gray representation.
 * @param[out] chromosome containing the gray representation of the binary chromosome.
 */
std::vector<int> sga_gray::binary_to_gray(const std::vector<int> &binary) const
{
	int length = binary.size();
	std::vector<int> gray = binary;

	// the MSB is the same
	// uses the XOR table
	for(int i = 1; i < length; i++) {
		gray[i] = binary[i]^binary[i-1];
	}

	return gray;
}

/// Encode a decision vector in its a gray representation chromosome.
/**
 * @param[in] decision_vector x to convert in its gray representation.
 * @param[in] lb: lower bounds for the encoding.
 * @param[in] ub: upper bounds for the encoding.
 * @param[out] chromosome containing the gray representation of the decision vector.
 */
std::vector<int> sga_gray::encode(const decision_vector &x, const decision_vector &lb, const decision_vector &ub) const
{
	std::vector<int> encoded_x(x.size() * m_bit_encoding, 0);

	for(decision_vector::size_type i=0; i<x.size(); i++) {
		std::vector<int> encoded_gene = binary_to_gray(double_to_binary(x.at(i),lb.at(i),ub.at(i)));
		// copy the gene at the right location
		for(int j=0; j<m_bit_encoding; j++) {
			encoded_x[i*m_bit_encoding + j] = encoded_gene.at(j);
		}
	}

	return encoded_x;
}

/// Decode a gray encoded chromosome in its a decision vector representation.
/**
 * @param[in] chrom: chromosome to convert in its vector of double (decision vector) representation.
 * @param[in] lb: lower bounds for the encoding.
 * @param[in] ub: upper bounds for the encoding.
 * @param[out] chromosome containing the decision vector representation of the chromosome.
 */
decision_vector sga_gray::decode(const std::vector<int> &chrom, const decision_vector &lb, const decision_vector &ub) const
{
	decision_vector decoded_x(chrom.size() / m_bit_encoding, 0.);

	for(decision_vector::size_type i=0; i<decoded_x.size(); i++) {
		// extract the part that need to be decoded
		std::vector<int> encoded_gene(m_bit_encoding);

		for(int j=0; j<m_bit_encoding; j++) {
			encoded_gene[j] = chrom.at(i*m_bit_encoding + j);
		}

		decoded_x[i] = binary_to_double(gray_to_binary(encoded_gene),lb.at(i),ub.at(i));
	}

	return decoded_x;
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::sga_gray)
