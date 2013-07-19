#include "racing.h"
#include "../rng.h"
#include "../problem/base_stochastic.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>

#include <utility>

namespace pagmo{ namespace util{
	
namespace racing{
//! @cond Doxygen skips the following

/// Special population tailored to the needs of racing
/**
 * This is a special type of population which allows direct manipulation
 * of m_container with fitness vectors and constraint vectors.
 *
 * param[in] pop Population to be copied over
 **/
racing_population::racing_population(const population &pop): population(pop)
{
}

/// Update decision_vector without invoking objective function
/**
 * One of the most bizarre things you could in the world of PaGMO -- setting a
 * decision vector of an individual without evaluating its objective function.
 * This causes the fitness and constraint vectors of that individual to be
 * completely invalid, so as the best_x, best_f, etc. Only use this if you know
 * what you are doing.
 *
 * (Essentially, this is the first halve of the canonical set_x)
 **/
void racing_population::set_x_noeval(const size_type idx, const decision_vector &x)
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}
	if (!problem().verify_x(x)) {
		pagmo_throw(value_error,"decision vector is not compatible with problem");

	}
	// Set decision vector.
	m_container[idx].cur_x = x;
}

/// Update directly fitness and constraint
/**
 * One of the most bizarre things you could do in the world of PaGMO --
 * directly setting fitness and constraint vectors. It only does what it says
 * -- set_fc -- meaning cur_x and best_x might become invalid as set_fc simply
 *  ignores and does not check their validity. Only use this if you know what
 *  you are doing.
 *
 *  (Essentially, this is the last halve of the canonical set_x)
 **/
void racing_population::set_fc(const size_type idx, const fitness_vector &f, const constraint_vector &c)
{
	if (idx >= size()) {
		pagmo_throw(index_error, "Invalid individual position in set_fc");
	}
	if (f.size() != problem().get_f_dimension()) {
		pagmo_throw(value_error, "Incompatible fitness dimension in set_fc");
	}
	if (c.size() != problem().get_c_dimension()) {
		pagmo_throw(value_error, "Incompatible constraint dimension in set_fc");
	}
	m_container[idx].cur_f = f;
	m_container[idx].cur_c = c;
	if (!m_container[idx].best_x.size() ||
		problem().compare_fc(m_container[idx].cur_f,m_container[idx].cur_c,m_container[idx].best_f,m_container[idx].best_c))
	{
		m_container[idx].best_x = m_container[idx].cur_x;
		m_container[idx].best_f = m_container[idx].cur_f;
		m_container[idx].best_c = m_container[idx].cur_c;
	}
	update_dom(idx);
}

/// Friedman rank assignment (before every racing iteration)
/**
 * Updates racers with the friedman ranks, assuming that the required individuals
 * in racing_pop have been re-evaluated with the newest seed. Ranking is the same
 * as the one returned by get_best_idx(N), but in case of ties an average rank
 * will be assigned among those who tied.
 *
 * @param[out] racers Data strcture storing the racing data which will be updated
 * @param[in] racing_pop Population on which racing will run
 *
**/
void f_race_assign_ranks(std::vector<racer_type>& racers, const racing_population& racing_pop)
{
	// Update ranking to be used for stat test
	// Note that the ranking returned by get_best_idx() requires some post-processing,
	// as individuals not currently in race should be ignored.

	typedef population::size_type size_type;

	std::vector<size_type> raw_order = racing_pop.get_best_idx(racing_pop.size());
	int cur_rank = 1;
	std::vector<size_type> ordered_idx_active;
	for(size_type i = 0; i < raw_order.size(); i++){
		int ind_idx = raw_order[i];
		// Note that some individuals would have been dropped along the way and
		// become inactive. They should do not affect the latest rankings in any way.
		if(racers[ind_idx].active){
			racers[ind_idx].m_hist.push_back(cur_rank);
			cur_rank++;
			ordered_idx_active.push_back(ind_idx);
		}
	}
	
	// --Adjust ranking to cater for ties--
	// 1. Check consecutively ranked individuals whether they are tied.
	std::vector<bool> tied(ordered_idx_active.size() - 1, false);
	if(racing_pop.problem().get_f_dimension()==1){
		// Single-objective case
		population::trivial_comparison_operator comparator(racing_pop);
		for(size_type i = 0; i < ordered_idx_active.size() - 1; i++){	
			size_type idx1 = ordered_idx_active[i];
			size_type idx2 = ordered_idx_active[i+1];
			if (!comparator(idx1, idx2) && !comparator(idx2, idx1)){
				tied[i] = true;
			}
		}
	}	
	else{
		// Multi-objective case
		// TODO: Possibe big jumps in observation data due to Pareto ranks
		population::crowded_comparison_operator comparator(racing_pop);
		for(size_type i = 0; i < ordered_idx_active.size() - 1; i++){	
			size_type idx1 = ordered_idx_active[i];
			size_type idx2 = ordered_idx_active[i+1];
			if (!comparator(idx1, idx2) && !comparator(idx2, idx1)){
				tied[i] = true;
			}
		}
	}
	// std::cout << "Tied stats: "; for(size_type i = 0; i < tied.size(); i++) std::cout << tied[i] <<" "; std::cout<<std::endl;

	// 2. For all the individuals who are tied, modify their rankings to
	// be the average rankings in case of no tie.
	size_type cur_pos = 0;
	size_type begin_avg_pos;
	while(cur_pos < tied.size()){
		begin_avg_pos = cur_pos;
		while(tied[cur_pos]==1){
			cur_pos++;
			if(cur_pos >= tied.size()){
				break;
			}
		}

		double avg_rank = 0;
		// std::cout << "Ties between: ";
		for(size_type i = begin_avg_pos; i <= cur_pos; i++){
			avg_rank += racers[ordered_idx_active[i]].m_hist.back() / ((double)cur_pos - begin_avg_pos + 1);
			// std::cout << ordered_idx_active[i] << " (r=" << racers[ordered_idx_active[i]].m_hist.back() << ") ";
		}
		// std::cout << std::endl;
		// if(cur_pos - begin_avg_pos + 1 > 1)
			// std::cout << "Setting tied ranks between " << cur_pos - begin_avg_pos + 1 << " individuals to be " << avg_rank   << std::endl;
		for(size_type i = begin_avg_pos; i <= cur_pos; i++){
			racers[ordered_idx_active[i]].m_hist.back() = avg_rank;
		}

		// If no tie at all for this begin pos
		if(begin_avg_pos == cur_pos){
			cur_pos++;
		}
	}

	// Update mean of rank (which also reflects the sum of rank, useful for later
	// pair-wise test)
	for(size_type i = 0; i < ordered_idx_active.size(); i++){
		racer_type& cur_racer = racers[ordered_idx_active[i]];
		cur_racer.m_mean = 0;
		for(unsigned int i = 0; i < cur_racer.length(); i++){
			cur_racer.m_mean += (cur_racer.m_hist[i]) / (double)cur_racer.length();
		}
	}

	//std::cout << "Adjusted ranking: "; for(size_type i = 0; i < ordered_idx_active.size(); i++) std::cout << "(" << ordered_idx_active[i] << ")-" << racers[ordered_idx_active[i]].m_hist.back() << " "; std::cout << std::endl;
}

/// Rank adjustment (after every racing iteration)
/**
 * Once some individuals are removed (de-activated) from the racing pool,
 * adjust the previous ranks of the remaining individuals. The resulted ranks
 * are free of influence from those just-deleted individuals.
 *
 * @param[out] racers Data structure for storing racing data which will be updated
 * @param[in] deleted_racers Indices of the individuals who have just be de-activated
 *
**/
void f_race_adjust_ranks(std::vector<racer_type>& racers, const std::vector<population::size_type>& deleted_racers)
{
	for(unsigned int i = 0; i < racers.size(); i++){
		if(!racers[i].active) continue;
		for(unsigned int j = 0; j < racers[i].length(); j++){
			int adjustment = 0;
			for(unsigned int k = 0; k < deleted_racers.size(); k++){
				if(racers[i].m_hist[j] > racers[deleted_racers[k]].m_hist[j]){
					adjustment++;
				}
			}
			racers[i].m_hist[j] -= adjustment;
		}
		//std::cout << "After deleting, rank of [" << i << "] adjusted as: " << racers[i].m_hist << std::endl;
	}
	//TODO: Tie cases not handled yet, worth it?
}

/// Perform a Friedman test
/**
 * Friedman test has the following procedures:
 * (1) Check if the null hypothesis that all rakings of the treatments are
 *     equivalent can be rejected. This involves the use of a chi-squared
 *     distribution.
 * (2) If so, perform pair-wise comparison betwen the treatments based on
 *     rank sum. This can be achieved by using a t-distribution.
 *
 * @param[in] X Observation data, each element (vector) represents the measurements for each "treatment"
 * @param[in] delta Confidence level for the statistical test
 *
 * @return Result of the statistical test 
 */
stat_test_result friedman_test(const std::vector<std::vector<double> >& X, double delta)
{	
	// TODO: throw when X is empty
	
	unsigned int N = X.size(); // # of different configurations
	unsigned int B = X[0].size(); // # of different instances

	// std::cout << "N = " << N << " B = " << B << std::endl;

	// (Now obtained the rankings, done with problems and pop.)

	// ----------- Stat. tests  starts-------------

	// Compute mean rank
	std::vector<double> X_mean(N, 0);
	for(unsigned int i = 0; i < N; i++){
		for(unsigned int j = 0; j < X[i].size(); j++){
			X_mean[i] += X[i][j];
		}
		X_mean[i] /= (double)X[i].size();
	}

	// Fill in R and T
	std::vector<double> R(N, 0);
	double A1 = 0;
	double C1 = B * N * (N+1) * (N+1) / 4.0;

	for(unsigned int i = 0; i < N; i++){
		for(unsigned int j = 0; j < B; j++){
			R[i] += X[i][j];
			A1 += (X[i][j])*(X[i][j]);
		}
	}

	double T1 = 0;
	for(unsigned int i = 0; i < N; i++){
		T1 += ((R[i] - B*(N+1)/2.0) * (R[i] - B*(N+1)/2.0));
	}
	T1 *= (N - 1) / (A1 - C1);
	
	using boost::math::chi_squared;
	using boost::math::students_t;
	using boost::math::quantile;

	chi_squared chi_squared_dist(N - 1);

	double delta_quantile = quantile(chi_squared_dist, 1 - delta);

	//std::cout << "T1: "<< T1 << "; delta_quantile = " << delta_quantile << std::endl;

	// Null hypothesis 0: All the ranks observed are equally likely
	bool null_hypothesis_0 = (boost::math::isnan(T1) || T1 < delta_quantile);

	stat_test_result res;

	// True if null hypothesis is rejected -- some differences between
	// the observations will be significant
	res.trivial = null_hypothesis_0;

	res.is_better = std::vector<std::vector<bool> >(N, std::vector<bool>(N, false));

	if(!null_hypothesis_0){	
		std::vector<std::vector<bool> >& is_better = res.is_better;

		students_t students_t_dist((N-1)*(B-1));
		double t_delta2_quantile = quantile(students_t_dist, 1 - delta/2.0);
		double Q = sqrt(((A1-C1) * 2.0 * B / ((B-1) * (N-1))) * (1 - T1 / (B*(N-1))));
		for(unsigned int i = 0; i < N; i++){
			for(unsigned int j = i + 1; j < N; j++){
				double diff_r = fabs(R[i] - R[j]);
				// std::cout<< "diff_r = " << diff_r << std::endl;
				// Check if a pair is statistically significantly different
				if(diff_r > t_delta2_quantile * Q){
					if(X_mean[i] < X_mean[j]){
						is_better[i][j] = true;
					}
					if(X_mean[j] < X_mean[i]){
						is_better[j][i] = true;
					}
				}
			}
		}
	}

	return res;

}
//! @endcond Doxygen comments the following

}}}
