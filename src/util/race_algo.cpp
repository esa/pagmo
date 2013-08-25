#include "race_algo.h" 
#include "race_pop.h"
#include "../problem/base_stochastic.h"

namespace pagmo { namespace util { namespace racing {

namespace metrics_algos {

/**
 * This class implements the mechanism to assign a quality measure to each
 * member of a set of algorithms, while being treated as "a population of
 * algorithms".  Required by internally by race_algo to transform the racing of
 * algorithm into racing of individuals.
 *
 * The decision vector in this case is a single integer, representing which
 * algorithm is to be invoked to evolve the population, e.g. if x = [k], this
 * problem returns some measure of how well the k-th algorithm could evolve a
 * population.
 *
 * The standard way of assigning fitness to an algorithm contains two cases:
 * (1) Single problem: Algo's fitness assigned to be the fitness of the evolved
 * champion on the problem, evolved by the specified algorithm.
 * (2) Multiple problem: Randomly sample a single problem (based on current
 * seed) from the pool of problems, then proceed as (1).
 */
class standard : public problem::base_stochastic
{
	public:
		standard(const std::vector<problem::base_ptr> &probs = std::vector<problem::base_ptr>(), const std::vector<algorithm::base_ptr> &algos = std::vector<algorithm::base_ptr>(), unsigned int seed = 0, unsigned int pop_size = 100);

	protected:
		//copy constructor
		standard(const standard &);
		problem::base_ptr clone() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		//TODO: Now only handles unconstrained single objective
		//void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
	
	private:

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_stochastic>(*this);
			ar & m_algos;
			ar & m_probs;
			ar & m_pop_size;
		}

		void setup(const std::vector<problem::base_ptr> &probs, const std::vector<algorithm::base_ptr> &algos);

		std::vector<algorithm::base_ptr> m_algos;
		std::vector<problem::base_ptr> m_probs;
		unsigned int m_pop_size;
};

/// Constructor of the performance metrics of algorithms on a single problem
/**
 * The performance of each algorithm is encoded in the fitness implementation
 * of this meta-problem. For example, to obtain the performance of the first
 * algorithm, one can do (in fact only race_algo should know about this):
 *
 * std::vector<unsigned int> idx = {0};
 * fitness_values f = this->objfun(idx);
 *
 * @param[in] problem Problem against which the algorithms will be evaluated
 * @param[in] algos The set of algorithms
 * @param[in] seed Seed to be used internally as a stochastic problem
 * @param[in] pop_size Size of the population to be evolved
 *
 * @throws value_error if there are incompatible problems in the supplied set. All the problems need to have the same fitness and constraint dimension.
 *
 */
standard::standard(const std::vector<problem::base_ptr> &probs, const std::vector<algorithm::base_ptr> &algos, unsigned int seed, unsigned int pop_size): base_stochastic(1, 1, probs.front()->get_f_dimension(), probs.front()->get_c_dimension(), probs.front()->get_ic_dimension(), 0, seed), m_pop_size(pop_size)
{
	setup(probs, algos);
}

/// Set up the internal structures of target problem and algorithm sets
/**
 * Check the sanity of the supplied problems. All fitness and constraint
 * dimensions have to be the same across all the problems. This is one of the
 * assumption made in the current implementation of race_algo.
 */
void standard::setup(const std::vector<problem::base_ptr> &probs, const std::vector<algorithm::base_ptr> &algos)
{
	// Sanity check on the algos and probs
	if(algos.size() == 0){
		pagmo_throw(value_error, "Empty algorithm set in race_algo");
	}
	problem::base::f_size_type fdim = probs[0]->get_f_dimension();
	problem::base::c_size_type cdim = probs[0]->get_c_dimension();
	for(unsigned int i = 1; i < probs.size(); i++){
		if(probs[i]->get_f_dimension() != fdim){
			pagmo_throw(value_error, "Incompatible fitness dimension among the problems");
		}
		if(probs[i]->get_c_dimension() != cdim){
			pagmo_throw(value_error, "Incompatible constraint dimension among the problems");
		}
	}
	
	// Take snapshots of the supplied algos and problems
	for(unsigned int i = 0; i < algos.size(); i++){
		m_algos.push_back(algos[i]->clone());
	}
	for(unsigned int i = 0; i < probs.size(); i++){
		m_probs.push_back(probs[i]->clone());
	}
	
	// Set bounds. Decision variable is simply the index of the selected algorithm
	set_bounds(0, m_algos.size()-1);
}

/// Copy Constructor. Performs a deep copy
standard::standard(const standard &standard_copy):
	// TODO: Can this be simplified? Required often when writing meta-problems
	base_stochastic(1, 1, standard_copy.get_f_dimension(),
			standard_copy.get_c_dimension(),
			standard_copy.get_ic_dimension(), 0, standard_copy.m_seed),
	m_algos(standard_copy.m_algos),
	m_probs(standard_copy.m_probs),
	m_pop_size(standard_copy.m_pop_size)
{
	set_bounds(standard_copy.get_lb(), standard_copy.get_ub());
}

/// Clone method.
problem::base_ptr standard::clone() const
{
	return problem::base_ptr(new standard(*this));
}

/// The performance of an algorithm encoded in the fitness function
/**
 * The performance of an algorithm is defined as the champion a population,
 * after being evolved by the algorithm.
 */
void standard::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// Bounds check (TODO: No need? already done by problem if bounds set correctly)
	if(x[0] < 0 || x[0] >= m_algos.size()){
		pagmo_throw(value_error, "Out of bound algorithm index");
	}

	// Seeding control
	for(unsigned int i = 0; i < m_algos.size(); i++){
		m_algos[i]->reset_rngs(m_seed);
	}
	m_drng.seed(m_seed);

	// Randomly sample a problem if required
	unsigned int prob_idx;
	if(m_probs.size() == 1){
		prob_idx = 0;
	}
	else{
		prob_idx = (unsigned int)(m_drng() * 100000) % m_probs.size();
	}

	// Fitness defined as the quality of the champion in the evolved
	// population, evolved by the selected algorithm
	population pop(*m_probs[prob_idx], m_pop_size, m_seed);
	m_algos[x[0]]->evolve(pop);
	f = pop.champion().f;
}

//TODO: Computation of fitness and constraints are decoupled -- a single pop
//has to be evolved twice? How to cache?
//void compute_constraints_impl(constraint_vector &, const decision_vector &) const
//{
//}

}


/// Constructor of the racing mechanism for algorithms
/**
 * Construct for a single problem
 *
 * @param[in] algos The set of algorithms to be raced
 * @param[in] prob The problem to be considered
 * @param[in] The size of the population that the algorithms will be evolving
 * @param[in] seed Seed to be used in racing mechanisms
 */
race_algo::race_algo(const std::vector<algorithm::base_ptr> &algos, const problem::base &prob, unsigned int pop_size, unsigned int seed): m_pop_size(pop_size), m_seed(seed)
{
	for(unsigned int i = 0; i < algos.size(); i++){
		m_algos.push_back(algos[i]->clone());
	}
	m_probs.push_back(prob.clone());
}

/// Constructor of the racing mechanism for algorithms
/**
 * Construct for a set of problems
 *
 * @param[in] algos The set of algorithms to be raced
 * @param[in] probs The set of problems to be considered
 * @param[in] The size of the population that the algorithms will be evolving
 * @param[in] seed Seed to be used in racing mechanisms
 */
race_algo::race_algo(const std::vector<algorithm::base_ptr> &algos, const std::vector<problem::base_ptr> &probs, unsigned int pop_size, unsigned int seed): m_pop_size(pop_size), m_seed(seed)
{
	for(unsigned int i = 0; i < algos.size(); i++){
		m_algos.push_back(algos[i]->clone());
	}
	for(unsigned int i = 0; i < probs.size(); i++){
		m_probs.push_back(probs[i]->clone());
	}
}

/// Juice of racing mechanisms for algorithms
/**
 * The interface of race_algo mirrors race_pop.
 *
 * @param[in] n_final Desired number of winners.
 * @param[in] min_trials Minimum number of trials to be executed before dropping algorithms.
 * @param[in] max_f_evals Maximum number of objective function evaluation before the race ends.
 * @param[in] delta Confidence level for statistical testing.
 * @param[in] active_set Indices of individuals that should participate in the race. If empty, race on the whole algorithm set.
 * @param[in] race_goal Whether to extract the best or the worst algorithm.
 * @param[in] screen_output Whether to log racing status on the console output.
 *
 * @see Refer to util::racing::race_pop for the details of the racing mechanisms.
 */
std::pair<std::vector<unsigned int>, unsigned int> race_algo::run(
	const unsigned int n_final,
	const unsigned int min_trials,
	const unsigned int max_count,
	double delta,
	const std::vector<unsigned int> &active_set,
	const bool race_best,
	const bool screen_output)
{
	/**	
	 * NOTE: The race of algorithms is in essence very similar to the race of
	 * individuals. A direct re-implementation of the racing routines for
	 * algorithms will cause unnecessary code duplication. To avoid that,
	 * internally, here the task of race_algo is reformulated into its ``dual
	 * form" in pop_race. The results of race_algo is the reconstructed from
	 * the results of pop_race. If Hoeffding race or Bernstein race is desired,
	 * they can be implemented in racing.cpp and be invoked in race_pop.cpp.
	 * This way, all the existing / to-be-implemented racing mechanisms can be
	 * shared by both race_pop and race_algo.
	 */	

	// Construct an internal population, such that the winners of the race in
	// this population corresponds to the winning algorithm
	metrics_algos::standard metrics(m_probs, m_algos, m_seed, m_pop_size);
	population algos_pop(metrics);
	for(unsigned int i = 0; i < m_algos.size(); i++){
		decision_vector algo_idx(1);
		algo_idx[0] = i;
		algos_pop.push_back(algo_idx);
	}

	// Conversion to types that pop_race is familiar with
	std::vector<population::size_type> pop_race_active_set(active_set.size());
	for(unsigned int i = 0; i < active_set.size(); i++){
		pop_race_active_set[i] = active_set[i];
	}

	// Run the actual race
	std::pair<std::vector<population::size_type>, unsigned int> res =
		algos_pop.race(n_final, min_trials, max_count, delta, pop_race_active_set,
				race_best, screen_output);

	// Convert the result to the algo's context
	std::pair<std::vector<unsigned int>, unsigned int> res_algo_race;
	res_algo_race.first.clear();
	for(unsigned int i = 0; i < res.first.size(); i++){
		res_algo_race.first.push_back(res.first[i]);
	}
	res_algo_race.second = res.second;

	return res_algo_race;
}

}}}
