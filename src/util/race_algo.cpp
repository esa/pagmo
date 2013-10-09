#include "race_algo.h" 
#include "race_pop.h"
#include "../problem/base_stochastic.h"

#include <algorithm>

namespace pagmo { namespace util { namespace racing {

namespace metrics_algos {

/// Returns the maximum constraint dimension of from set of problems
static problem::base::c_size_type get_max_c_dimension(const std::vector<problem::base_ptr> &probs)
{
	problem::base::c_size_type max_c_dim = 0;
	for(std::vector<problem::base_ptr>::size_type i = 0; i < probs.size(); i++){
		max_c_dim = std::max(max_c_dim, probs[i]->get_c_dimension());
	}
	return max_c_dim;
}

/// Returns the maximum inequality constraint dimension from a set of problems
static problem::base::c_size_type get_max_ic_dimension(const std::vector<problem::base_ptr> &probs)
{
	problem::base::c_size_type max_ic_dim = 0;
	for(std::vector<problem::base_ptr>::size_type i = 0; i < probs.size(); i++){
		max_ic_dim = std::max(max_ic_dim, probs[i]->get_ic_dimension());
	}
	return max_ic_dim;
}

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
 * (2) Multiple problems: Randomly sample a single problem (based on current
 * seed) from the pool of problems, then proceed as (1).
 *
 * Currently supports box constrained and equality / inequality constrained
 * single-objective problems
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
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
	
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
		constraint_vector zero_pad_constraint(const constraint_vector&, problem::base::c_size_type) const;	
		
		void evaluate_algorithm(unsigned int) const;

		std::vector<algorithm::base_ptr> m_algos;
		std::vector<problem::base_ptr> m_probs;
		unsigned int m_pop_size;

		// To avoid inefficiency resulted from the decoupled fitness and
		// constraint computation
		mutable std::vector<bool> m_is_first_evaluation;
		mutable std::vector<unsigned int> m_database_seed;
		mutable std::vector<fitness_vector> m_database_f;
		mutable std::vector<constraint_vector> m_database_c;
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
 * @param[in] probs std::vector of pagmo::problem::base_ptr against which the algorithms will be evaluated
 * @param[in] algos std::vector of pagmo::algorithm::base_ptr
 * @param[in] seed Seed to be used internally as a stochastic problem
 * @param[in] pop_size Size of the population to be evolved
 *
 * @throws value_error if there are incompatible algorithms or problems in the supplied sets (multi-objective algorithms not supported yet).
 *
 */
standard::standard(const std::vector<problem::base_ptr> &probs, const std::vector<algorithm::base_ptr> &algos, unsigned int seed, unsigned int pop_size): base_stochastic(1, 1, probs.front()->get_f_dimension(), get_max_c_dimension(probs), get_max_ic_dimension(probs), 0, seed), m_pop_size(pop_size), m_is_first_evaluation(algos.size(), true), m_database_seed(algos.size()), m_database_f(algos.size()), m_database_c(algos.size())
{
	setup(probs, algos);
}

/// Set up the internal structures of target problem and algorithm sets
/**
 * Check the sanity of the supplied problems. Currently, racing of algorithms
 * over multi-objective problems is not yet supported.
 *
 */
void standard::setup(const std::vector<problem::base_ptr> &probs, const std::vector<algorithm::base_ptr> &algos)
{
	// Sanity check on the algos and probs
	if(algos.size() == 0){
		pagmo_throw(value_error, "Empty algorithm set in race_algo");
	}
	for(unsigned int i = 0; i < probs.size(); i++){
		if(probs[i]->get_f_dimension() > 1){
			pagmo_throw(value_error, "Racing of multi-objective algorithms is not supported yet");
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
	base_stochastic(1, 1, standard_copy.get_f_dimension(),
			standard_copy.get_c_dimension(),
			standard_copy.get_ic_dimension(), 0, standard_copy.m_seed),
	m_algos(standard_copy.m_algos),
	m_probs(standard_copy.m_probs),
	m_pop_size(standard_copy.m_pop_size),
	m_is_first_evaluation(standard_copy.m_is_first_evaluation),
	m_database_seed(standard_copy.m_database_seed),
	m_database_f(standard_copy.m_database_f),
	m_database_c(standard_copy.m_database_c)
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
	evaluate_algorithm(x[0]);
	f = m_database_f[x[0]];
}

/// The performance of an algorithm in terms of a constraint vector
/**
 * Similar to objfun_impl(), the constraint violation is determined by the
 * champion of the evolved population.
 */
void standard::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	evaluate_algorithm(x[0]);
	c = m_database_c[x[0]];
}

/// Pad constraint vector with non-violating values (i.e. 0)
/**
 * This is necessary when the underlying problems have different dimension. The
 * contraint dimensions of this meta-problem is set as the maximum constraint
 * dimension of all the underlying problems. Constraint vector having less than
 * that dimension will be padded here, up to the necessary dimension.
 */
constraint_vector standard::zero_pad_constraint(const constraint_vector& c,  problem::base::c_size_type ic_dim) const
{
	constraint_vector c_padded(get_c_dimension(), 0);
	c_size_type c_dim = c.size();
	for(problem::base::c_size_type i = 0; i < c_dim - ic_dim; i++){
		c_padded[i] = c[i];
	}
	// Note that this->get_ic_dimension() must be larger than ic_dim by definition
	for(problem::base::c_size_type i = c_dim - ic_dim; i < c_dim; i++){
		c_padded[i + get_ic_dimension() - ic_dim] = c[i];
	}
	return c_padded;
}

/// Evaluate the target algorithm and stores the evaluation results
/**
 * If objfun_impl() and compute_constraints_impl() are decoupled, this
 * meta-problem needs to evolve an identical population twice to get
 * respectively the fitness and constraint vectors. To avoid this, this
 * function serves as a proxy to the actual evaluation job, and it will take
 * care of skipping unnecessary evolution of population.
 *
 */
void standard::evaluate_algorithm(unsigned int algo_idx) const
{
	if(algo_idx >= m_algos.size()){
		pagmo_throw(value_error, "Out of bound algorithm index");
	}

	// The requested data is ready, nothing to do
	if(!m_is_first_evaluation[algo_idx] && m_database_seed[algo_idx] == m_seed){
		return;
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
	m_algos[algo_idx]->evolve(pop);

	// Store the data, to be retrieved by objfun_impl() or compute_constraints_impl()
	m_is_first_evaluation[algo_idx] = false;
	m_database_seed[algo_idx] = m_seed;
	m_database_f[algo_idx] = pop.champion().f;
	m_database_c[algo_idx] = zero_pad_constraint(pop.champion().c, m_probs[prob_idx]->get_ic_dimension());
}


}


/// Constructor of the racing mechanism for algorithms
/**
 * Construct for a single problem
 *
 * @param[in] algos The set of algorithms to be raced
 * @param[in] prob The problem to be considered
 * @param[in] pop_size The size of the population that the algorithms will be evolving
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
 * @param[in] pop_size The size of the population that the algorithms will be evolving
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
 * @param[in] max_count Maximum number of objective function evaluation before the race ends.
 * @param[in] delta Confidence level for statistical testing.
 * @param[in] active_set Indices of individuals that should participate in the race. If empty, race on the whole algorithm set.
 * @param[in] race_best Whetn true extracts the best, otherwise the worst algorithm.
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
	 * form" in pop_race. The results of race_algo is then reconstructed from
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
	    algos_pop.race(n_final, min_trials, max_count, delta,
	                   pop_race_active_set, race_best, screen_output);

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
