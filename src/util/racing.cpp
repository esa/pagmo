#include "racing.h"
#include "../rng.h"
#include "../problem/base_stochastic.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>

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

/// Constructor based on a problem
/**
 * Creates an empty population which can be further manipulated using special
 * routines in racing_population
 *
 * param[in] prob Problem of interest
 **/
racing_population::racing_population(const problem::base &prob): population(prob)
{
}

/// Update decision_vector without invoking objective function
/**
 * One of the most bizarre things you could in the world of PaGMO -- setting a
 * decision vector of an individual without evaluating its objective function.
 * This causes the fitness and constraint vectors of that individual to be
 * completely invalid, so as the best_x, best_f, etc. Use with caution.
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
 *  ignores and does not check their validity. Another note is that best_f and
 *  best_c will always mirror cur_f and cur_c. Use with caution.
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
	// NOTE: As update_dom() uses best_f and best_c when computing Pareto ranks
	// and hence, racing_population can be used a way to by pass this in order
	// to respect more the concept of racing
	m_container[idx].best_f = f;
	m_container[idx].best_c = c;
	update_dom(idx);
}


/// Append individual with given decision vector without invoking objective function
/**
 * Only allocates spaces for the incoming decision vector in m_container,
 * m_dom_list, and m_dom_count. If the fitnesses and contraints will be used,
 * make sure to call set_fc() prior to using them, or bear the consequences.
 *
 * NOTE: Differences compared with the real push_back:
 * - No evaluation (no set_x)
 * - No initialization of velocities (v's are all zeros) 
 * - Champion will not be allocated correctly subsequently (requires real set_x)
 **/
void racing_population::push_back_noeval(const decision_vector &x)
{
	// No checking on the validity of x:
	// Accept a dummy x, as when using racing_population, the main focus is on
	// the fitness and constraint vectors. The main purpose of this function is
	// to allocate the spaces but skip the evaluation.

	// Store sizes temporarily.
	const fitness_vector::size_type f_size = problem().get_f_dimension();
	const constraint_vector::size_type c_size = problem().get_c_dimension();
	const decision_vector::size_type p_size = problem().get_dimension();
	// Push back an empty individual.
	m_container.push_back(individual_type());
	m_dom_list.push_back(std::vector<size_type>());
	m_dom_count.push_back(0);
	// Resize individual's elements.
	m_container.back().cur_x.resize(p_size);
	m_container.back().cur_v.resize(p_size);
	m_container.back().cur_c.resize(c_size);
	m_container.back().cur_f.resize(f_size);

	// As we bypass set_x which will allocate spaces for bests, they must be
	// explicitly allocated here.
	m_container.back().best_f.resize(f_size);
	m_container.back().best_c.resize(c_size);
	
	// Set the individual.
	set_x_noeval(m_container.size() - 1, x);
}


/// Returns the ``rankings" of the individuals.
/**
 * The rankings are obtained using get_best_idx() in the base population class.
 * The rankings are further processed to cater for possible ties, which is a
 * step typically required by ranking-based statistical testing.
 **/
std::vector<double> racing_population::get_rankings() const
{
	std::vector<double> rankings(size());
	std::vector<size_type> raw_order = get_best_idx(size());

	int cur_rank = 1;
	for(size_type i = 0; i < raw_order.size(); i++){
		int ind_idx = raw_order[i];
		rankings[ind_idx] = cur_rank;
		cur_rank++;
	}
	
	// --Adjust ranking to cater for ties--
	// 1. Check consecutively ranked individuals whether they are tied.
	std::vector<bool> tied(size() - 1, false);
	if(problem().get_f_dimension() == 1){
		// Single-objective case
		population::trivial_comparison_operator comparator(*this);
		for(size_type i = 0; i < size() - 1; i++){	
			if (!comparator(i, i+1) && !comparator(i+1, i)){
				tied[i] = true;
			}
		}
	}	
	else{
		// Multi-objective case
		population::crowded_comparison_operator comparator(*this);
		for(size_type i = 0; i < size() - 1; i++){	
			if (!comparator(i, i+1) && !comparator(i+1, i)){
				tied[i] = true;
			}
		}
	}

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
		for(size_type i = begin_avg_pos; i <= cur_pos; i++){
			avg_rank += rankings[i] / ((double)cur_pos - begin_avg_pos + 1);
		}
		for(size_type i = begin_avg_pos; i <= cur_pos; i++){
			rankings[i] = avg_rank;
		}

		// If no tie at all for this begin pos
		if(begin_avg_pos == cur_pos){
			cur_pos++;
		}
	}
	return rankings;
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
void f_race_assign_ranks(std::vector<racer_type>& racers, const racing_population& racing_pop_full)
{
	// Update ranking to be used for stat test
	// Note that the ranking returned by get_best_idx() requires some post-processing,
	// as individuals not currently in race should be ignored.

	typedef population::size_type size_type;

	// Construct a more condensed population with only active individuals
	racing_population racing_pop(racing_pop_full.problem());
	decision_vector dummy_x(racing_pop_full.problem().get_dimension(), 0);
	std::vector<size_type> idx_mapping;
	for(size_type i = 0; i < racers.size(); i++){
		if(racers[i].active){
			racing_pop.push_back_noeval(dummy_x);
			racing_pop.set_fc(racing_pop.size()-1, racing_pop_full.get_individual(i).cur_f, racing_pop_full.get_individual(i).cur_c);
			idx_mapping.push_back(i);
		}
	}
	
	// Get the rankings in the sense of satistical testing
	std::vector<double> rankings = racing_pop.get_rankings();
	for(unsigned int i = 0; i < racing_pop.size(); i++){
		racers[idx_mapping[i]].m_hist.push_back(rankings[i]);
	}

	// Update mean of rank (which also reflects the sum of rank, useful for
	// later pair-wise test)
	for(size_type i = 0; i < rankings.size(); i++){
		racer_type& cur_racer = racers[idx_mapping[i]];
		cur_racer.m_mean = 0;
		for(unsigned int j = 0; j < cur_racer.length(); j++){
			cur_racer.m_mean += (cur_racer.m_hist[j]) / (double)cur_racer.length();
		}
	}
	
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
	}
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
 *
 */
stat_test_result core_friedman_test(const std::vector<std::vector<double> >& X, double delta)
{	
	pagmo_assert(X.size() > 0);
	
	unsigned int N = X.size(); // # of different configurations
	unsigned int B = X[0].size(); // # of different instances

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

/// Returns the pair-wise statistical testing results based on Friedman Test
/**
 * @param[in] racers List of racers which will be filled up with rank data
 * @param[in] pop Population storing the updated fitness and constraint of each
 * active individual in race.
 * @return A structure containing the statistical testing results 
 **/
stat_test_result friedman_test(std::vector<racer_type> &racers, const std::vector<population::size_type> &in_race, const racing_population &pop, double delta)
{
	f_race_assign_ranks(racers, pop);

	// Observation data (TODO: is this necessary ? ... a lot of memory
	// allocation gets done here and we already have in memory all we need.
	// could we not pass by reference directly racers and in_race to the
	// friedman test?)
	std::vector<std::vector<double> > X;
	for(unsigned int i = 0; i < in_race.size(); i++){
		X.push_back(racers[in_race[i]].m_hist);
	}

	// Friedman Test
	stat_test_result ss_result = core_friedman_test(X, delta);
	return ss_result;
}

/// Returns the pair-wise statistical testing results based on Wilcoxon test
/**
 * @param[in] racers Racers structure whose m_hist fields will be updated
 * @param[in] in_race Indices of the two active individuals
 * @param[in] pop Population storing full evaluation history of the two active
 * individuals
 * @return A structure containing the statistical testing results 
 **/
stat_test_result wilcoxon_ranksum_test(std::vector<racer_type> &racers, const std::vector<population::size_type> &in_race, const racing_population& wilcoxon_pop, double delta)
{
	pagmo_assert(in_race.size() == 2);

	// Get the rankings of all the data points
	std::vector<double> rankings = wilcoxon_pop.get_rankings();
	
	// Need to update racers in case no statistical significant results can be
	// found, which will then default to selecting the one with best mean. Two
	// specific individuals in the wilcoxon_pop correspond to the newest two
	// evaluated points.
	racers[in_race[0]].m_hist.push_back(rankings[wilcoxon_pop.size()/2 - 1]);
	racers[in_race[1]].m_hist.push_back(rankings[wilcoxon_pop.size() - 1]);

	// Compute the mean ranks
	for(std::vector<population::size_type>::const_iterator it = in_race.begin(); it != in_race.end(); it++){
		racer_type& cur_racer = racers[*it];
		cur_racer.m_mean = 0;
		for(unsigned int j = 0; j < cur_racer.length(); j++){
			cur_racer.m_mean += (cur_racer.m_hist[j]) / (double)cur_racer.length();
		}
	}

	std::vector<std::vector<double> > X(2);
	unsigned int n_samples = wilcoxon_pop.size() / 2;
	for(unsigned int i = 0; i < 2; i++){
		X[i].resize(n_samples);
	}
	for(unsigned int i = 0; i < wilcoxon_pop.size(); i++){
		X[i%2][i/2] = rankings[i];
	}
	return core_wilcoxon_ranksum_test(X, delta);
}

// Helper routine for Wilcoxon test
std::size_t wilcoxon_faculty( std::size_t n )
{
	if( n == 1 )
		return( 1 );

	return( n * wilcoxon_faculty( n-1 ) );
}

// Helper routine for Wilcoxon test
double wilcoxon_frequency( double u, int sampleSizeA, int sampleSizeB )
{
	if( u < 0. || sampleSizeA < 0 || sampleSizeB < 0 )
		return( 0. );
	
	if( u == 0 && sampleSizeA >= 0 && sampleSizeB >= 0 )
		return( 1 );

	return( wilcoxon_frequency( u - sampleSizeB, sampleSizeA - 1, sampleSizeB ) + wilcoxon_frequency( u, sampleSizeA, sampleSizeB - 1 ) );
}

// Performs Wilcoxon Rank Sum Test (a.k.a.Wilcoxon–Mann–Whitney test)
// Intended to be used when the number of active racers is only two, as it has
// been reported that under such circumstances this test is more data efficient
// than Friedman test.
stat_test_result core_wilcoxon_ranksum_test(const std::vector<std::vector<double> > &X, double delta)
{
	if(X.size() != 2){
		pagmo_throw(value_error, "Wilcoxon rank-sum test is applicable only when the group size is 2");
	}
	unsigned int N = 2;

	// Compute sum of ranks
	std::vector<double> rank_sum(N, 0);
	for(unsigned int i = 0; i < N; i++){
		for(unsigned int j = 0; j < X[i].size(); j++){
			rank_sum[i] += X[i][j];
		}
	}

	// Fill in pair-wise comparison results
	stat_test_result res(N);

	int sizeA = X[0].size();
	int sizeB = X[1].size();
	if(sizeA < 12 && sizeB < 12){
		// Cannot use normal approximation for the rank-sum statistic when
		// sample size is small. Boost does not have the required Wilcoxon
		// rank-sum distribution (currently in their todo list)... This piece
		// of code is adapted from the Shark machine learning library.
		int wA = rank_sum[0];
		//int wB = rank_sum[1];
		double uA = wA - sizeA * ( sizeA + 1 ) / 2.;
		//double uB = wB - sizeB * ( sizeB + 1 ) / 2.;
		double pA = (double)( wilcoxon_faculty( sizeA ) * wilcoxon_faculty( sizeB ) ) / (double)( wilcoxon_faculty( sizeA + sizeB ) ) * wilcoxon_frequency( uA, sizeA, sizeB );
		//double pB = (double)( wilcoxon_faculty( sizeA ) * wilcoxon_faculty( sizeB ) ) / (double)( wilcoxon_faculty( sizeA + sizeB ) ) * wilcoxon_frequency( uB, sizeA, sizeB );
		if(pA < delta){
			res.trivial = false;
			if(rank_sum[0] > rank_sum[1]){
				res.is_better[1][0] = true;
			}
			else{
				res.is_better[0][1] = true;
			}
		}
	}	
	else{
		// Use the normal approximation
		using boost::math::normal;
		int n1 = X[0].size();
		int n2 = X[1].size();
		normal normal_dist(n1*n2/2, sqrt(n1*n2*(n1+n2+1.0)/12.0));
		double delta_quantile_upp = quantile(normal_dist, 1 - delta);
		if(rank_sum[0] > rank_sum[1] && rank_sum[0] > delta_quantile_upp){
			res.trivial = false;
			res.is_better[1][0] = true;
		}
		else if(rank_sum[1] > rank_sum[0] && rank_sum[1] > delta_quantile_upp){
			res.trivial = false;
			res.is_better[0][1] = true;
		}
	}

	return res;
}

//! @endcond Doxygen comments the following

}}}
