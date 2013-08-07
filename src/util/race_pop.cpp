#include "race_pop.h"
#include "../problem/base_stochastic.h"

namespace pagmo { namespace util { namespace racing {

/// Constructor
/**
 * Construct a race_pop object from an external population and a seed. The seed
 * will determine all racing conditions.
 *
 * @param[in] pop population containing the individuals to race
 * @param[in] seed seed of the race
 */
race_pop::race_pop(const population& pop, unsigned int seed): m_pop(pop), m_seeds(), m_seeder(seed), m_cache_data(pop.size())
{
}

// Check if the provided active_set is valid.
void race_pop::_validate_active_set(const std::vector<population::size_type>& active_set, unsigned int pop_size)
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
void race_pop::_validate_problem_stochastic(const problem::base& prob)
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
void race_pop::_validate_racing_params(const population& pop, const population::size_type n_final, const unsigned int min_trials, const unsigned int max_f_evals, double delta, unsigned int active_set_size)
{
	if(n_final > pop.size()){
		pagmo_throw(value_error, "Number of intended winner is too large");
	}
	if(delta < 0 || delta > 1){
		pagmo_throw(value_error, "Confidence level should be a small value greater than zero");
	}
	if (max_f_evals < active_set_size) {
		pagmo_throw(value_error, "Maximum number of function evaluations is smaller than the number of racers");
	}
	if (max_f_evals < min_trials*active_set_size) {
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
 * @param[in] race_goal Whether to extract the best or the worst individuals.
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
std::pair<std::vector<population::size_type>, unsigned int> race_pop::run(const population::size_type n_final, const unsigned int min_trials, const unsigned int max_f_evals, const double delta, const std::vector<population::size_type>& active_set, const bool race_best, const bool screen_output)
{
	// We start validating the inputs:
	// a - Problem has to be stochastic
	_validate_problem_stochastic(m_pop.problem());
	// b - active_set has to contain valid indexes
	_validate_active_set(active_set, m_pop.size());
	// c - Other parameters have to be sane
	_validate_racing_params(m_pop, n_final, min_trials, max_f_evals, delta, active_set.size());

	typedef population::size_type size_type;

	// The stochastic problem's seed will be changed using a pre-determined sequence
	unsigned int seed_idx = 0;

	// The race will terminate when a given number of function evaluations will be exceeded
	//race_termination_condition::type term_cond = race_termination_condition::EVAL_COUNT;
	
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

	size_type N_begin = in_race.size();
	size_type n_final_best;

	// Complimentary relationship between race-for-best and race-for-worst
	if(race_best){
		n_final_best = n_final;
	}
	else{
		n_final_best = N_begin - n_final;
	}

	unsigned int count_iter = 0;
	unsigned int count_nfes = 0;

	// Start of the main loop. It will stop as soon as we have decided enough winners or
	// discarded enough losers
	while(decided.size() < n_final_best && decided.size() + in_race.size() > n_final_best){
		
		count_iter++;

		if(screen_output){
			std::cout << "\n-----Iteration: " << count_iter << ", evaluation count = " << count_nfes << std::endl;
			std::cout << "Decided: " << decided << std::endl;
			std::cout << "In-race: " << in_race << std::endl;
			std::cout << "Discarded: " << discarded << std::endl;
		}

		// Check if there is enough budget for evaluating the individuals in the race 
		// TODO: Need to discount for caching mechanism
		if(count_nfes + in_race.size() > max_f_evals){
			break;
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

		// Perform re-evaluation on necessary individuals under current seed
		for(std::vector<size_type>::iterator it = in_race.begin(); it != in_race.end(); ++it) {
			// Case 1: Current racer has previous data that can be reused, no
			// need to be evaluated with this seed
			if(cache_data_exist(*it, count_iter-1)){
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
				cache_insert_data(*it, f_vec, c_vec);
			}
		}

		f_race_assign_ranks(racers, m_pop);

		// Enforce a minimum required number of trials
		if(count_iter < min_trials)
			continue;

		// Observation data (TODO: is this necessary ? ... a lot of memory allocation gets done here and we
		// already have in memory all we need. could we not pass by reference directly racers and in_race 
		// to the friedman test?)
		std::vector<std::vector<double> > X;
		for(unsigned int i = 0; i < in_race.size(); i++){
			X.push_back(racers[in_race[i]].m_hist);
		}

		// Friedman Test
		stat_test_result ss_result = friedman_test(X, delta);

		if(!ss_result.trivial){
			// Inside here some pairs must be statistically different, let's find them out
		
			const std::vector<std::vector<bool> >& is_better = ss_result.is_better;

			// std::cout << "Null hypothesis 0 rejected!" << std::endl;

			// std::vector<size_type> out_of_race;

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
				// std::cout << "[" << *it_i << "]: vote_decide = " << vote_decide << ", vote_discard = " << vote_discard << std::endl;
				if(vote_decide >= N_begin - n_final_best - discarded.size()){
					to_decide[i] = true;
				}
				else if(vote_discard >= n_final_best - decided.size()){
				//else if(vote_discard >= 1){ // Equivalent to the previous more aggressive approach
					to_discard[i] = true;
				}
			}

			std::vector<size_type> new_in_race;
			for(unsigned int i = 0; i < in_race.size(); i++){
				if(to_decide[i]){
					decided.push_back(in_race[i]);
					racers[in_race[i]].active = false;
				}
				else if(to_discard[i]){
					discarded.push_back(in_race[i]);
					racers[in_race[i]].active = false;
				}
				else{
					new_in_race.push_back(in_race[i]);
				}
			}

			in_race = new_in_race;

			// Check if this is that important
			// f_race_adjust_ranks(racers, out_of_race);
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

/// Clear all the cache
void race_pop::reset_cache()
{
	m_cache_data.clear();
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

// Produce new seeds and append to the list of seeds
void race_pop::generate_seeds(unsigned int num_seeds){
	for(unsigned int i = 0; i < num_seeds; i++){
		m_seeds.push_back(m_seeder());
	}
}

// Get the n-th seed to be used for evaluation of the n-th data point. With
// this we can ensure that all the aligned data points are generated using the
// same rng seed.
unsigned int race_pop::get_current_seed(unsigned int seed_idx){
	if(seed_idx >= m_seeds.size()){
		const unsigned int expanding_length = 500;
		generate_seeds(expanding_length);
	}
	return m_seeds[seed_idx]; 
}

}}}
