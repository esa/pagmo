#ifndef PAGMO_ALGORITHM_DE_H
#define PAGMO_ALGORITHM_DE_H

#include "../config.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// Differential Evolution Algorithm
/**
 *
 * \image html de.jpg "Differential Evolution block diagram"
 * \image latex de.jpg "Differential Evolution block diagram" width=5cm
 *
 * Differential Evolution is an heuristic optimizer developed by Rainer Storn and Kenneth Price.
 *
 * ''A breakthrough happened, when Ken came up with the idea of using vector differences for perturbing
 * the vector population. Since this seminal idea a lively discussion between Ken and Rainer and endless
 * ruminations and computer simulations on both parts yielded many substantial improvements which
 * make DE the versatile and robust tool it is today'' (from the official web pages....)
 *
 * The implementation provided for PaGMO derives from the code provided in the official
 * DE web site and is suitable for box-constrained continuous optimization.
 *
 * At each call of the evolve method a number of function evaluations equal to m_gen * pop.size()
 * is performed.
 *
 * @see http://www.icsi.berkeley.edu/~storn/code.html for the official DE web site
 * @see http://www.springerlink.com/content/x555692233083677/ for the paper that introduces Differential Evolution
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 */
		
class __PAGMO_VISIBLE de: public base
{
public:
	de(int gen, double f, double cr, int strategy);
	base_ptr clone() const;
	void evolve(population &) const;
protected:
	std::string human_readable_extra() const;
private:
	// Number of generations.
	const int m_gen;
	// Weighting factor
	const double m_f;
	// Crossover probability
	const double m_cr;
	// Startegy
	const int m_strategy;
};

}} //namespaces

#endif // DE_H
