#include "de.h"


namespace pagmo { namespace algorithm {

/// Advanced constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] m_f rate of choosing from memory.
 * @param[in] cr minimum pitch adjustment rate.
 * @param[in] strategy maximum pitch adjustment rate.
 * @throws value_error if m_f,cr are not in the [0,1] interval, strategy is not one of 1 .. 10, gen is negative
 */
de::de(int gen, double m_f, double cr, int strategy):m_gen(gen),m_f(m_f),m_cr(cr),m_strategy(strategy) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (strategy < 1 || strategy > 10) {
		pagmo_throw(value_error,"strategy index must be one of 1 ... 10");
	}
	if (cr < 0 || m_f < 0 || cr > 1 || m_f > 1) {
		pagmo_throw(value_error,"the m_f and m_cr parameters must be in the [0,1] range");
	}
}

/// Clone method.
base_ptr de::clone() const
{
	return base_ptr(new de(*this));
}

void de::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;


	//We perform some checks to determine wether the problem/population are suitable for DE
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for DE to optimise");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and DE is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and DE is not suitable to solve it");
	}

	if (NP < 6) {
		pagmo_throw(value_error,"for DE at least 6 individuals in the population are needed");
	}


	// Some vectors used during evolution are allocated here. As the algorithm does not know about NP and D
	// it does not help too much to have these as mutable members of the algorithm as they would anyway need to be resized
	// here. Resizing introdues the risk that the same algorithm, called on different problems would
	// mess up everything as resizing vector of vectors is tricky.

	decision_vector dummy(D), tmp(D); //dummy is used for initialisation purposes, tmp to contain the mutated candidate
	std::vector<decision_vector> popold(NP,dummy), popnew(NP,dummy), popswap(NP,dummy);
	std::vector<double> fit(NP);
	decision_vector gbX(D),gbIter(D);
	double newfitness;		//new fitness of the mutaded candidate
	double gbfit;			//global best fitness

	//We extract from pop the chromosomes and fitness associated
	for (std::vector<double>::size_type i = 0; i < NP; ++i) {
		popold[i] = pop.get_individual(i).cur_x;
		popnew[i] = popold[i];
		fit[i] = pop.get_individual(i).cur_f[0];
	}

	// Initialise the global bests
	gbX=popold[0];
	gbfit=fit[0];

	for (size_t i = 1; i < NP; ++i) {		//the int i = 1 jumps the first member as it is already set as the best
		if (fit[i] < gbfit) {
			gbfit = fit[i];			// save best member ever
			gbX = popold[i];
		}
	}
	gbIter = gbX;				// save best member of generation

	// Main DE iterations
	size_t r1,r2,r3,r4,r5;	//indexes to the selected population members
	for (int gen = 0; gen < m_gen; ++gen) {
		//Start of the loop through the deme
		for (size_t i = 0; i < NP; ++i) {
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 2 !!!     */
				r1 = (size_t)(m_drng()*NP);
			} while (r1==i);

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 3 !!!     */
				r2 = (size_t)(m_drng()*NP);
			} while ((r2==i) || (r2==r1));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 4 !!!     */
				r3 = (size_t)(m_drng()*NP);
			} while ((r3==i) || (r3==r1) || (r3==r2));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 5 !!!     */
				r4 = (size_t)(m_drng()*NP);
			} while ((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 6 !!!     */
				r5 = (size_t)(m_drng()*NP);
			} while ((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (m_strategy == 1) { /* strategy DE0 (not in our paper) */
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc), L = 0;
				do {
					tmp[n] = gbIter[n] + m_f*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"gbIter[]"-schemes experience misconvergence. Try e.g. m_f=0.7 and m_cr=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (m_strategy == 2) { /* strategy DE1 in the techreport */
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc), L = 0;
				do {
					tmp[n] = popold[r1][n] + m_f*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try m_f=0.85 and m_cr=1.------*/
			/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (m_strategy == 3) { /* similiar to DE2 but generally better */
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc), L = 0;
				do {
					tmp[n] = tmp[n] + m_f*(gbIter[n] - tmp[n]) + m_f*(popold[r1][n]-popold[r2][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (m_strategy == 4) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc), L = 0;
				do {
					tmp[n] = gbIter[n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (m_strategy == 5) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc), L = 0;
				do {
					tmp[n] = popold[r5][n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (m_strategy == 6) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc);
				for (size_t L = 0; L < Dc; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] + m_f*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (m_strategy == 7) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc);
				for (size_t L = 0; L < Dc; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r1][n] + m_f*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (m_strategy == 8) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc);
				for (size_t L = 0; L < Dc; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = tmp[n] + m_f*(gbIter[n] - tmp[n]) + m_f*(popold[r1][n]-popold[r2][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (m_strategy == 9) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc);
				for (size_t L = 0; L < Dc; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (m_strategy == 10) {
				tmp = popold[i];
				size_t n = (size_t)(m_drng()*Dc);
				for (size_t L = 0; L < Dc; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r5][n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					}
					n = (n+1)%Dc;
				}
			}


			/*=======Trial mutation now in tmp[]. force feasibility and how good this choice really was.==================*/
			// a) feasibility
			size_t i2 = 0;
			while (i2<Dc) {
				if ((tmp[i2] < lb[i2]) || (tmp[i2] > ub[i2]))
					tmp[i2] = m_drng()*(ub[i2]-lb[i2]) + lb[i2];
				++i2;
			}
			newfitness = prob.objfun(tmp)[0];    /* Evaluate new vector in tmp[] */
			//b) how good?
			if (newfitness <= fit[i]) {  /* improved objective function value ? */
				fit[i]=newfitness;
				popnew[i] = tmp;
				if (newfitness<gbfit) {        /* Was this a new minimum for the deme? */
					/* if so...*/
					gbfit=newfitness;          /* reset gbfit to new low...*/
					gbX=tmp;
				}
			} else {
				popnew[i] = popold[i];
			}
			/* swap population arrays. New generation becomes old one */

		}//End of the loop through the deme

		/* Save best population member of current iteration */
		gbIter = gbX;

		/* swap population arrays. New generation becomes old one */
		for (size_t i = 0; i < NP; ++i) {
			popswap[i] = popold[i];
			popold[i] = popnew[i];
			popnew[i] = popswap[i];
		}

	}//end main DE iterations

	//we end by constructing the object population containing the final results
	//WASTE OF OBJECTIVE FUNCTIONS EVALUATIONS!!!! CACHE NEEDS TO REMEMBER AT LEAST NP
	for (population::size_type i=0;i<NP;++i){
		pop.set_x(i,popnew[i]);
	}

}



/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string de::human_readable_extra() const
{
	std::ostringstream s;
	s << "\tGenerations:\t" << m_gen << '\n';
	s << "\tWeight parameter (m_f):\t\t" << m_f << '\n';
	s << "\tCrossover parameter (m_cr):\t" << m_cr << '\n';
	s << "\tStrategy selected:\t" << m_strategy << '\n';
	return s.str();
}

}} //namespaces
