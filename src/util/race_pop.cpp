#include "race_pop.h"
#include "../problem/ackley.h"
#include "../problem/base_stochastic.h"

#include <map>
#include <utility>

namespace pagmo { namespace util { namespace racing {

/// Constructor
/**
 * Construct a race_pop object from an external population and a seed. The seed
 * will determine all racing conditions.
 *
 * @param[in] pop population containing the individuals to race
 * @param[in] seed seed of the race
 */
race_pop::race_pop(const population& pop, unsigned int seed): m_race_seed(seed), m_pop(pop), m_pop_wilcoxon(pop), m_seeds(), m_seeder(seed), m_use_caching(true), m_cache_data(pop.size()), m_cache_averaged_data(pop.size())
{
	register_population(pop);
}

/// Constructor
/**
 * Construct a race_pop object from a seed. The seed will determine all racing
 * conditions. The exact population on which race will run can be (and must be)
 * supplied later via register_population().
 *
 * @param[in] seed seed of the race
 */
race_pop::race_pop(unsigned int seed): m_race_seed(seed), m_pop(population(problem::ackley())), m_pop_wilcoxon(population(problem::ackley())), m_pop_registered(false), m_seeds(), m_seeder(seed), m_use_caching(true), m_cache_data(0), m_cache_averaged_data(0)
{
}

/// Update the population on which the race will run
/**
 * This also re-allocate the spaces required to the cache entries.
 *
 * @param[in] pop The new population
 **/
void race_pop::register_population(const population &pop)
{
	m_pop = pop;
	// This is merely to set up the problem in wilcoxon pop
	m_pop_wilcoxon = pop;
	reset_cache();
	if(m_cache_data.size() != pop.size()){
		m_cache_data.resize(pop.size());
		m_cache_averaged_data.resize(pop.size());
	}
	cache_register_signatures(pop);
	m_pop_registered = true;
}

/// Get the number of individuals in the registered population
population::size_type race_pop::size() const
{
	return m_pop.size();
}

// Check if the provided active_set is valid.
void race_pop::_validate_active_set(const std::vector<population::size_type>& active_set, unsigned int pop_size) const
{
	if(active_set.size() == 0)
		return;
	std::vector<bool> hit(pop_size, 0);
	for(unsigned int i = 0; i < active_set.size(); i++){
		if(active_set[i] >= pop_size){
			pagmo_throw(index_error, "Racing: Active set contains invalid (out of bound) index.");
		}
		if(hit[active_set[i]] == true){
			pagmo_throw(index_error, "Racing: Active set contains repeated indices.");
		}
		hit[active_set[i]] = true;
	}
}

// Check if the problem is stochastic
void race_pop::_validate_problem_stochastic(const problem::base& prob) const
{
	try
	{
		dynamic_cast<const pagmo::problem::base_stochastic &>(prob).get_seed();
	}
	catch (const std::bad_cast& e)
	{
		pagmo_throw(type_error, "Attempt to call racing routines on a non-stochastic problem, use get_best_idx() instead");
	}
}

// Check if the other parameters for racing are sensible
void race_pop::_validate_racing_params(const population& pop, const population::size_type n_final, double delta) const
{
	if(n_final > pop.size()){
		pagmo_throw(value_error, "Number of intended winner is too large");
	}
	if(delta < 0 || delta > 1){
		pagmo_throw(value_error, "Confidence level should be a small value greater than zero");
	}
}

// Check if the provided budget is sensible
void race_pop::_validate_budget(const unsigned int min_trials, const unsigned int max_f_evals, const std::vector<population::size_type>& in_race) const
{
	// The fevals required to complete the first iteration
	unsigned int min_required_1 = compute_required_fevals(in_race, 1);
	if (max_f_evals < min_required_1) {
		pagmo_throw(value_error, "Maximum number of function evaluations is smaller than the number of racers");
	}
	// The fevals required to complete first min_trials-th iterations
	unsigned int min_required_2 = 0;
	for(unsigned int i = 1; i <= min_trials; i++){
		min_required_2 += compute_required_fevals(in_race, i);
	}
	if (max_f_evals < min_required_2) {
		pagmo_throw(value_error, "You are asking for a mimimum amount of trials which cannot be made with the allowed function evaluation budget");
	}
}

// Construct the output (winner / loser) list based on the book kept data
// Two possible scenarios:
// (1) Goal is BEST: Output will be decided + best ones from in_race (if required).
// (2) Goal is WORST: Output will be discarded + worst ones from in_race (if required).
std::vector<population::size_type> race_pop::construct_output_list(
		const std::vector<racer_type>& racers,
		const std::vector<population::size_type>& decided,
		const std::vector<population::size_type>& in_race,
		const std::vector<population::size_type>& discarded,
		const population::size_type n_final,
		const bool race_best)
{
	typedef population::size_type size_type;
	// To be sorted either in ascending (BEST) or descending (WORST)
	std::vector<std::pair<double, size_type> > argsort;
	for(std::vector<size_type>::const_iterator it = in_race.begin(); it != in_race.end(); ++it){
		argsort.push_back(std::make_pair(racers[*it].m_mean, *it));
	}
	std::vector<size_type> output;
	if(race_best){
		output = decided;
		std::sort(argsort.begin(), argsort.end(), std::less<std::pair<double, size_type> >());
	}
	else{
		output = discarded;
		std::sort(argsort.begin(), argsort.end(), std::greater<std::pair<double, size_type> >());
	}
	int sorted_idx = 0;
	while(output.size() < n_final){
		output.push_back(argsort[sorted_idx++].second);
	}	

	return output;
}

// Update m_pop with the evaluation data w.r.t current seed for Friedman test
//
// The resulting population is aligned with the racers, i.e. m_pop[0]
// corresponds to racers[0], storing the newest fitness and constraint vector
// evaluated under the new seed. Evaluation data can come from cache or fresh
// computation.
// 
// @return The number of objective function calls made
unsigned int race_pop::prepare_population_friedman(const std::vector<population::size_type>& in_race, unsigned int count_iter)
{
	unsigned int count_nfes = 0;
	// Perform re-evaluation on necessary individuals under current seed
	for(std::vector<population::size_type>::const_iterator it = in_race.begin(); it != in_race.end(); ++it) {
		// Case 1: Current racer has previous data that can be reused, no
		// need to be evaluated with this seed
		if(m_use_caching && cache_data_exist(*it, count_iter-1)){
			const eval_data& cached_data = cache_get_entry(*it, count_iter-1);
			m_pop.set_fc(*it, cached_data.f, cached_data.c);
		}
		// Case 2: No previous data can be reused, perform actual
		// re-evaluation and update the cache
		else{
			count_nfes++;
			const population::individual_type &ind = m_pop.get_individual(*it);
			fitness_vector f_vec = m_pop.problem().objfun(ind.cur_x);
			constraint_vector c_vec = m_pop.problem().compute_constraints(ind.cur_x);
			m_pop.set_fc(*it, f_vec, c_vec);
			if(m_use_caching)
				cache_insert_data(*it, f_vec, c_vec);
		}
	}
	return count_nfes;
}

/// Update m_pop_wilcoxon to contain evaluation data required for Wilcoxon test
/**
 * In the end, m_pop_wilcoxon should hold all the fc vectors evaluated so far
 * for the two active individuals. This is to facilitate the rankings required
 * by Wilcoxon rank-sum test.
 **/
unsigned int race_pop::prepare_population_wilcoxon(const std::vector<population::size_type>& in_race, unsigned int count_iter)
{
	unsigned int count_nfes = 0;
	if(in_race.size() != 2){
		pagmo_throw(value_error, "Wilcoxon rank sum test is only applicable when there are two active individuals");
	}	

	// Need to bootstrap by pulling all the previous evaluation data into
	// wilcoxon_pop for the first time when race is left with the two
	// individuals. Subsequent racing iterations can just re-use this
	// wilcoxon_pop memory structure and pad in necessary fc data for each
	// particular iteration.
	unsigned int start_count_iter;
	if(m_pop_wilcoxon.size() == 0){
		start_count_iter = 1;
	}
	else{
		start_count_iter = count_iter;
	}
	for(std::vector<population::size_type>::const_iterator it = in_race.begin(); it != in_race.end(); ++it) {
		decision_vector dummy_x;
		for(unsigned int i = start_count_iter; i <= count_iter; i++){
			// Case 1: Current racer has previous data that can be reused, no
			// need to be evaluated with this seed
			if(m_use_caching && cache_data_exist(*it, i-1)){
				m_pop_wilcoxon.push_back_noeval(dummy_x);
				const eval_data& cached_data = cache_get_entry(*it, i-1);
				m_pop_wilcoxon.set_fc(m_pop_wilcoxon.size()-1, cached_data.f, cached_data.c);
			}
			// Case 2: No previous data can be reused, perform actual
			// re-evaluation and update the cache
			else{
				count_nfes++;
				const population::individual_type &ind = m_pop.get_individual(*it);
				fitness_vector f_vec = m_pop.problem().objfun(ind.cur_x);
				constraint_vector c_vec = m_pop.problem().compute_constraints(ind.cur_x);
				m_pop_wilcoxon.push_back_noeval(dummy_x);
				m_pop_wilcoxon.set_fc(m_pop_wilcoxon.size()-1, f_vec, c_vec);
				if(m_use_caching)
					cache_insert_data(*it, f_vec, c_vec);
			}
		}
	}
	return count_nfes;
}

/// Computes the required number of actual fevals to complete the current iteration
/*
 * This function takes into account the existence of cache. For example, if the
 * cache has already been filled with many data, a call to run race with zero
 * budget is actually possible.
 *
 * @param[in] racers List of racers
 * @param[in] num_iter Index of current iteration (starts from 1)
 */
unsigned int race_pop::compute_required_fevals(const std::vector<population::size_type>& in_race, unsigned int num_iter) const
{
	unsigned int required_fevals = 0;
	for(unsigned int i = 0; i < in_race.size(); i++){
		if(!cache_data_exist(in_race[i], num_iter-1)){
			required_fevals++;	
		}
	}
	return required_fevals;
}

/// Races some individuals in a population
/**
 * Performs an F-Race among certain individuals in a population.
 * F-Races are based on the Friedman test, thus
 * a routine to determine the ranking of individuals w.r.t.
 * any given rng seed is needed. Here, population::get_best_idx() is used
 * extending the whole racing concept to multi-objective and constrained
 * optimization.
 *
 * Specifically, racing contains the following steps:
 *
 * (1) Re-evaluate the active individuals with the newest random seed
 * (2) Assign ranks to the active individuals and append to the observation data
 * (3) Perform statistical test based on Friedman test (thus obtain pair-wise
 *     comparison result with statistical significance)
 * (4) Update three individual index lists: decided, discarded, in_race
 * 	   - decided: No need to be further evaluated, clearly superior
 * 	   - discarded: No need to be further evaluated, clearly inferior
 * 	   - in_race: Undecided, more evaluations needed on them (a.k.a active)
 * (5) Repeat (1) until termination condidtion met
 * (6) Return decided. If too few individuals are in decided, append it
 *     with individuals from in_race based on their rank sum (smaller the better).
 *
 * @param[in] n_final Desired number of winners.
 * @param[in] min_trials Minimum number of trials to be executed before dropping individuals.
 * @param[in] max_f_evals Maximum number of objective function evaluation before the race ends.
 * @param[in] delta Confidence level for statistical testing.
 * @param[in] active_set Indices of individuals that should participate in the race. If empty, race on the whole population.
 * @param[in] term_cond Termination condition (MAX_BUDGET, MAX_DATA_COUNT)
 * @param[in] race_best When true races for the best individuals
 * @param[in] screen_output Whether to log racing status on the console output.
 *
 * @return std::pair first: the indices of the individuals that remain in the
 * race in the end, a.k.a the winners, second: the function evaluations actually needed
 *
 * @throws type_error if the underlying problem is not stochastic
 * @throws index_error if active_set is invalid (out of bound / repeated indices)
 * @throws value_error if other specified racing parameters are not sensible
 *
 * @see Birattari, M., Stützle, T., Paquete, L., & Varrentrapp, K. (2002). A Racing Algorithm for Configuring Metaheuristics. GECCO ’02 Proceedings of the Genetic and Evolutionary Computation Conference (pp. 11–18). Morgan Kaufmann Publishers Inc.
 * @see Heidrich-Meisner, Verena, & Christian Igel (2009). Hoeffding and Bernstein Races for Selecting Policies in Evolutionary Direct Policy Search. Proceedings of the 26th Annual International Conference on Machine Learning, pp. 401-408. ACM Press.
 */
std::pair<std::vector<population::size_type>, unsigned int> race_pop::run(
			 const population::size_type n_final, 
			 const unsigned int min_trials, 
			 const unsigned int max_f_evals, 
			 const double delta, 
			 const std::vector<population::size_type>& active_set, 
			 termination_condition term_cond, 
			 const bool race_best, 
			 const bool screen_output
			 )
{
	// First check whether the a population has been properly registered
	if(!m_pop_registered){
		pagmo_throw(value_error, "Attempt to run race but no population is registered");
	}
	// We start validating the inputs:
	// a - Problem has to be stochastic
	_validate_problem_stochastic(m_pop.problem());
	// b - active_set has to contain valid indexes
	_validate_active_set(active_set, m_pop.size());
	// c - Other parameters have to be sane
	_validate_racing_params(m_pop, n_final, delta);

	typedef population::size_type size_type;
	
	// Temporary: Consider a fresh start every time race() is called
	std::vector<racer_type> racers(m_pop.size(), racer_type());

	// If active_set is empty, default to race all individuals
	if(active_set.size() == 0){
		for(size_type i = 0; i < m_pop.size(); i++){
			racers[i].active = true;
		}
	}
	else{
		for(size_type i = 0; i < active_set.size(); i++){
			racers[active_set[i]].active = true;
		}
	}
	
	// Indices of racers who are currently active in the pop's sense
	// Examples:
	// (1) in_race.size() == 5 ---> Only 5 individuals are still being raced.
	// (2) pop.get_individual(in_race[0]) gives a reference to an actual
	//     individual which is active in the race.
	std::vector<size_type> in_race;
	std::vector<size_type> decided;
	std::vector<size_type> discarded;

	for(size_type i = 0; i < racers.size(); i++){
		if(racers[i].active){
			in_race.push_back(i);
		}
	}

	if(term_cond == MAX_BUDGET){
		// d - Check if the given budget is too small
		_validate_budget(min_trials, max_f_evals, in_race);		
	}

	size_type N_begin = in_race.size();
	size_type n_final_best;

	// Complimentary relationship between race-for-best and race-for-worst
	if(race_best){
		n_final_best = n_final;
	}
	else{
		n_final_best = N_begin - n_final;
	}

	// Reset data holder for wilcoxon test
	m_pop_wilcoxon.clear();
	bool use_wilcoxon = false;

	unsigned int count_iter = 0;
	unsigned int count_nfes = 0;

	// The stochastic problem's seed will be changed using a pre-determined sequence
	unsigned int seed_idx = 0;

	// Start of the main loop. It will stop as soon as we have decided enough winners or
	// discarded enough losers
	while(decided.size() < n_final_best && decided.size() + in_race.size() > n_final_best){
		
		count_iter++;

		if(screen_output){
			std::cout << "\n-----Iteration: " << count_iter << ", evaluation count = " << count_nfes << std::endl;
			std::cout << "Decided: " << decided << std::endl;
			std::cout << "In-race: " << in_race << std::endl;
			std::cout << "Discarded: " << discarded << std::endl;
			std::cout << "Mean ranks: ";
			for(unsigned int i = 0; i < in_race.size(); i++){
				if(racers[in_race[i]].active){
					std::cout <<  "(" << in_race[i] << "): " << racers[in_race[i]].m_mean << " ";
				}
			}
			std::cout << std::endl;
			print_cache_stats(in_race);
		}
		
		if(term_cond == MAX_BUDGET){
			// Check if there is enough budget for evaluating the individuals in the race 
			unsigned int required_fevals = compute_required_fevals(in_race, count_iter);
			if(count_nfes + required_fevals > max_f_evals){
				break;
			}
		}
		else{
			// Here max_f_evals has the meaning of "maximum number of data
			// points" to be considered for each individual
			if(count_iter > max_f_evals){
				break;
			}
		}

		unsigned int cur_seed = get_current_seed(seed_idx++);
		dynamic_cast<const pagmo::problem::base_stochastic &>(m_pop.problem()).set_seed(cur_seed);

		// NOTE: Here after resetting to a new seed, we do not perform
		// re-evaluation of the whole population, as this defeats the purpose
		// of doing race! Only the required individuals (i.e. those still
		// active in racing) shall be re-evaluated. A direct consequence is
		// that the champion of the population is not valid anymore nor the
		// individuals best_x and best_f -- they do not correspond to the
		// latest seed.  This is OK, as this only affects the local copy the
		// population, which will not be accessed elsewhere, and the champion
		// information is not used during racing.

		// Update m_pop with re-evaluation results or possibly data from cache,
		// and invoke statistical testing routines.
		stat_test_result ss_result;
		if(use_wilcoxon && in_race.size() == 2){
			// Perform Wilcoxon rank-sum test
			count_nfes += prepare_population_wilcoxon(in_race, count_iter);
			ss_result = wilcoxon_ranksum_test(racers, in_race, m_pop_wilcoxon, delta);
		}
		else{
			// Perform Friedman test
			count_nfes += prepare_population_friedman(in_race, count_iter);
			ss_result = friedman_test(racers, in_race, m_pop, delta);

		}

		if(count_iter < min_trials)
			continue;

		if(!ss_result.trivial){
			// Inside here some pairs must be statistically different, let's find them out
		
			const std::vector<std::vector<bool> >& is_better = ss_result.is_better;

			std::vector<size_type> out_of_race;

			std::vector<bool> to_decide(in_race.size(), false), to_discard(in_race.size(), false);

			for(unsigned int i = 0; i < in_race.size(); i++){
				unsigned int vote_decide = 0;
				unsigned int vote_discard = 0;
				for(unsigned int j = 0; j < in_race.size(); j++){
					if(i == j || to_decide[j] || to_discard[j])
						continue;
					// Check if a pair is statistically significantly different
					if(is_better[i][j]){
						vote_decide++;
					}
					else if(is_better[j][i]){
						vote_discard++;
					}
				}

				if(vote_decide >= N_begin - n_final_best - discarded.size()){
					to_decide[i] = true;
				}
				else if(vote_discard >= n_final_best - decided.size()){
					to_discard[i] = true;
				}
			}

			std::vector<size_type> new_in_race;
			for(unsigned int i = 0; i < in_race.size(); i++){
				if(to_decide[i]){
					decided.push_back(in_race[i]);
					out_of_race.push_back(in_race[i]);
					racers[in_race[i]].active = false;
				}
				else if(to_discard[i]){
					discarded.push_back(in_race[i]);
					out_of_race.push_back(in_race[i]);
					racers[in_race[i]].active = false;
				}
				else{
					new_in_race.push_back(in_race[i]);
				}
			}

			in_race = new_in_race;

			// Check if this is that important
			if(!use_wilcoxon || in_race.size() > 2){
				f_race_adjust_ranks(racers, out_of_race);
			}
		}

	};

	// Note that n_final (instead of n_final_best) is required here
	std::vector<size_type> winners =
		construct_output_list(racers, decided, in_race, discarded, n_final, race_best);

	if(screen_output){
		std::cout << "\nRace ends after " << count_iter << " iterations, incurred nfes = " << count_nfes << std::endl;
		std::cout << "Returning winners: " << std::vector<size_type>(winners.begin(), winners.end()) << std::endl;
	}

	return std::make_pair(winners, count_nfes);
}

/// Returns mean fitness of the individuals based on past evaluation data
/**
 * @param[in] ind_list The indices of the individuals whose mean fitness
 * vectors are to be extracted. If this is empty, mean data of all the
 * individuals will be returned.
 * 
 * @throws value_error if any of the requested individuals has not been raced before.
 *
 * @return Mean fitness vectors of the individuals in ind_list
 **/
std::vector<fitness_vector> race_pop::get_mean_fitness(const std::vector<population::size_type> &ind_list) const
{
	_validate_active_set(ind_list, m_pop.size());

	std::vector<population::size_type> active_set;
	// If empty list is given then assume all individuals are of interest
	if(ind_list.size() == 0){
		active_set.resize(size());
		for(population::size_type i = 0; i < size(); i++){
			active_set[i] = i;
		}
	}
	else{
		active_set = ind_list;
	}

	std::vector<fitness_vector> mean_fitness(active_set.size());
	for(unsigned int i = 0; i < active_set.size(); i++){
		if(m_cache_data[active_set[i]].size() == 0){
			pagmo_throw(value_error, "Request the mean fitness of an individual which has not been raced before");
		}
		mean_fitness[i] = m_cache_averaged_data[active_set[i]].f;
	}
	return mean_fitness;
}

/// Clear all the cache
void race_pop::reset_cache()
{
	for(unsigned int i = 0; i < m_cache_data.size(); i++){
		m_cache_data[i].clear();
		m_cache_averaged_data[i].f.clear();
		m_cache_averaged_data[i].c.clear();
	}
}

/// Insert a data_point
/**
 * @param[in] key_idx The key is just the position (index) of the individual
 * @param[in] f Fitness vector to be inserted
 * @param[in] c Constraint vector to be inserted
 **/
void race_pop::cache_insert_data(unsigned int key_idx, const fitness_vector &f, const constraint_vector &c)
{
	if(key_idx >= m_cache_data.size()){
		pagmo_throw(index_error, "cache_insert_data: Invalid key index");
	}
	m_cache_data[key_idx].push_back(eval_data(f,c));
	// Update the averaged data to be returned upon each race call
	if(m_cache_data[key_idx].size() == 1){
		m_cache_averaged_data[key_idx] = m_cache_data[key_idx].back();
	}
	else{
		unsigned int len = m_cache_data[key_idx].size();
		// Average for each fitness dimension
		for(unsigned int i = 0; i < m_cache_averaged_data[key_idx].f.size(); i++){
			m_cache_averaged_data[key_idx].f[i] =
				(m_cache_averaged_data[key_idx].f[i]*(len-1) +
				 m_cache_data[key_idx].back().f[i]) / (double)len;
		}
		// Average for each constraint dimension
		for(unsigned int i = 0; i < m_cache_averaged_data[key_idx].c.size(); i++){
			m_cache_averaged_data[key_idx].c[i] =
				(m_cache_averaged_data[key_idx].c[i]*(len-1) +
				 m_cache_data[key_idx].back().c[i]) / (double)len;
		}
	}
}

/// Delete the data associated with the key index corresponding to the
/// concerned decision vector
void race_pop::cache_delete_entry(unsigned int key_idx)
{
	m_cache_data.erase(m_cache_data.begin() + key_idx);
}

/// Check if the data point exist in the current cache for a particular key index
/**
 * This would be used to check if there is enough data points in the cache for
 * a decision vector. It is assumed that the key index supplied will be valid.
 **/
bool race_pop::cache_data_exist(unsigned int key_idx, unsigned int data_location) const
{
	if(key_idx >= m_cache_data.size()){
		pagmo_throw(index_error, "cache_data_exist: Invalid key index");
	}
	if(data_location < m_cache_data[key_idx].size()){
		return true;
	}
	return false;
}

/// Get a const reference to a data point
const race_pop::eval_data &race_pop::cache_get_entry(unsigned int key_idx, unsigned int data_location) const
{
	if(key_idx >= m_cache_data.size()){
		pagmo_throw(index_error, "cache_get_entry: Invalid key index");
	}
	if(data_location >= m_cache_data[key_idx].size()){
		pagmo_throw(index_error, "cache_get_netry: Invalid data location");
	}
	return m_cache_data[key_idx][data_location];
}

// Each cache entry is dedicated to an individual in the population, and it can
// be associated with a signature based on the decision vector of the
// individual. This can be useful to identify a match when trying to inherit
// memory from another race_pop structure to maximize information reuse.
void race_pop::cache_register_signatures(const population& pop)
{
	m_cache_signatures.clear();	
	for(population::size_type i = 0; i < pop.size(); i++){
		m_cache_signatures.push_back(pop.get_individual(i).cur_x);
	}
}

/// Inherits the memory of another race_pop object
/** If compatible, inherits past evaluation data from another race_pop object.
 * Useful in scenarios when racing individuals in a cross
 * generation setting.
*/
void race_pop::inherit_memory(const race_pop& src)
{
	// If seeds are different, no memory transfer is possible
	if(src.m_race_seed != m_race_seed){
		pagmo_throw(value_error, "Incompatible seed in inherit_memory");
	}
	std::map<decision_vector, unsigned int> src_cache_locations;
	for(unsigned int i = 0; i < src.m_cache_data.size(); i++){
		src_cache_locations.insert(std::make_pair(src.m_cache_signatures[i], i));
	}
	int cnt_transferred = 0;
	for(unsigned int i = 0; i < m_cache_data.size(); i++){
		std::map<decision_vector, unsigned int>::iterator it
			= src_cache_locations.find(m_cache_signatures[i]);
		if(it != src_cache_locations.end()){
			if(src.m_cache_data[it->second].size() > m_cache_data[i].size()){
				m_cache_data[i] = src.m_cache_data[it->second];
				m_cache_averaged_data[i] = src.m_cache_averaged_data[it->second];
				cnt_transferred++;
			}
		}
	}
}


/// Print some stats about the cache, for debugging purposes
void race_pop::print_cache_stats(const std::vector<population::size_type> &in_race) const
{
	for(std::vector<population::size_type>::const_iterator it = in_race.begin(); it != in_race.end(); it++){
		std::cout << "Cache of ind#" << *it << ": length = " << m_cache_data[*it].size() << std::endl;
	}
}


/// Set a new racing seed.
/**
 * Internally, the list of seed is cleared, all cache entry get invalidated,
 * and a new list of seeds will of be generated as required when racing
 * starts subsequently.
 *
 * @param[in] seed New seed to be set
 */
void race_pop::set_seed(unsigned int seed)
{
	m_race_seed = seed;
	m_seeder.seed(seed);
	m_seeds.clear();
	reset_cache();
}

// Produce new seeds and append to the list of seeds
void race_pop::generate_seeds(unsigned int num_seeds)
{
	for(unsigned int i = 0; i < num_seeds; i++){
		m_seeds.push_back(m_seeder());
	}
}

// Get the n-th seed to be used for evaluation of the n-th data point. With
// this we can ensure that all the aligned data points are generated using the
// same rng seed.
unsigned int race_pop::get_current_seed(unsigned int seed_idx)
{
	if(seed_idx >= m_seeds.size()){
		const unsigned int expanding_length = 500;
		generate_seeds(expanding_length);
	}
	return m_seeds[seed_idx]; 
}

}}}
