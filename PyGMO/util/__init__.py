from _util import *

def _race_pop_ctor(self, population=None, seed=0):
	"""
	Constructs a racing object responsible for racing individuals in a population
	
	USAGE: race_pop(population, seed=0)
	
	* population: The population containint the individuals to be raced
	* seed: Seed of the racing object
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	if(population != None):
		arg_list.append(population)
	arg_list.append(seed)
	self._orig_init(*arg_list)

race_pop._orig_init = race_pop.__init__
race_pop.__init__ = _race_pop_ctor

# enum
_util.race_pop.termination_condition = _util._termination_condition

def _race_pop_run(self, n_final, min_trials=0, max_count=500, delta=0.05, active_set=[], term_cond=race_pop.termination_condition.MAX_BUDGET, race_best=True, screen_output=False):
	"""
	Start a race among the individuals
	
	USAGE: race_pop.run(n_final, min_trials=0, max_count=500, delta=0.05, active_set=[], term_cond=termination_condition.MAX_BUDGET, race_best=True, screen_output=False):
	
	* n_final: The desired number of winners when race ends
	* min_trials: Each individuals be evaluated at least this number of times before being dropped
	* max_count: The allow number of function evaluations (MAX_BUDGET), or the maximum number of data points to be considered for each individual (MAX_DATA_COUNT)
	* delta: Confidence level of the statistical test
	* active_set: Indices of the individuals to be raced, empty means race all individuals
	* term_cond: Can be race_pop.termination_condition.MAX_BUDGET or race_pop.termination_condition.MAX_DATA_COUNT
	* race_best: When True winners are the best, otherwise winners are the worst
	* screen_output: Log racing stats at each iteration onto the screen
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(n_final)
	arg_list.append(min_trials)
	arg_list.append(max_count)
	arg_list.append(delta)
	arg_list.append(active_set)
	arg_list.append(term_cond)
	arg_list.append(race_best)
	arg_list.append(screen_output)
	return self._orig_run(*arg_list)

race_pop._orig_run = race_pop.run
race_pop.run = _race_pop_run


def _race_algo_ctor(self, algo_list, probs, pop_size=100, seed=0):
	"""
	Construct the racing objective responsible for racing algorithms
	
	USAGE: race_algo(algo_list, probs, pop_size=100, seed=0)

	* algo_list: The algorithms to be raced
	* probs: Can be a single PyGMO problem or a list of them
	* pop_size: All the algorithms will be evolving internally some random population of this size
	* seed: Seed of the race
	"""
	# We set the defaults or the kwargs
	arg_list=[]

	algo_vec = vector_of_algorithm_base_ptr()
	algo_vec.extend(algo_list)
	arg_list.append(algo_vec)

	try:
		l = len(probs);
		prob_vec = vector_of_problem_base_ptr()
		prob_vec.extend(probs)
		arg_list.append(prob_vec)
	except TypeError:
		arg_list.append(probs)

	arg_list.append(pop_size)
	arg_list.append(seed)

	self._orig_init(*arg_list)

race_algo._orig_init = race_algo.__init__
race_algo.__init__ = _race_algo_ctor

def _race_algo_run(self, n_final, min_trials=0, max_count=500, delta=0.05, active_set=[], race_best=True, screen_output=False):

	"""
	Start a race among several algorithms
	
	USAGE: race_algo.run(n_final, min_trials=0, max_count=500, delta=0.05, active_set=[], race_best=True, screen_output=False):
	
	* n_final: The desired number of winners when race ends
	* min_trials: Each algorithms be evaluated at least this number of times before being dropped
	* max_count: The allow number of algorithm performance evaluation (i.e. number of calls to evolve)
	* delta: Confidence level of the statistical test
	* active_set: Indices of the individuals to be raced, empty means race all algorithms
	* race_best: When True winners are the best, otherwise winners are the worst
	* screen_output: Log racing stats at each iteration onto the screen
	"""
	# We set the defaults or the kwargs
	arg_list=[]
	arg_list.append(n_final)
	arg_list.append(min_trials)
	arg_list.append(max_count)
	arg_list.append(delta)
	arg_list.append(active_set)
	arg_list.append(race_best)
	arg_list.append(screen_output)
	return self._orig_run(*arg_list)

race_algo._orig_run = race_algo.run
race_algo.run = _race_algo_run
