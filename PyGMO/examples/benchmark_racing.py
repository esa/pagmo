from PyGMO import *
import numpy as np
import matplotlib.pyplot as plt

def evolve_and_evaluate(algo, prob, pop_size, prob_orig, seed):

    algo.reset_rngs(seed)

    # Seed used to ensure both algorithm evolves an identical population
    pop = population(prob, pop_size, seed)

    pop = algo.evolve(pop)

    #winner_idx = pop.race(1)[0][0];
    #print "race winner", winner_idx, "vs champion idx", pop.get_best_idx()
    #champion_true_fitness = prob_orig.objfun(pop[winner_idx].cur_x)

    champion_true_fitness = prob_orig.objfun(pop.champion.x)
    #print 'Final champion =', champion_true_fitness

    return champion_true_fitness

def start_experiment(num_trials=20, pop_size=40, fevals_max=100000, nr_eval_per_x=40):

    # 1. Set up the problem
    prob_orig = problem.ackley(10)
    #prob_orig = problem.griewank(10)
    noise_level = 0.05
    prob_single_eval = problem.noisy(prob_orig, trials = 1, param_second = noise_level)
    prob_regular = problem.noisy(prob_orig, trials = nr_eval_per_x, param_second = noise_level)


    """
    Notes for SAME_PROB: Both algorithm will operate on the same version of the
    problem (same n_eval). pso_gen will evolve for fevals_max/2*pop_size times;
    pso_gen_racing will evolve until fevals_max is hit. In this case a single
    feval count in is referred to a single call to objfun() of the problem
    (with n_eval as 10).
    """
    SAME_PROB = True

    # 2A. Set up pso_gen algorithm without racing:
    # Each generation of pso_gen requires 2*pop_size*nr_eval_per_x
    # evaluations. Ignoring the cost of initialization here.
    # NOTE: No need to scale down if both algo has the same version of problem
    if SAME_PROB:
        gen_budget = fevals_max/(2*pop_size)
    else:
        gen_budget = fevals_max/(2*pop_size*nr_eval_per_x)
    print 'Non-racing pso gen will evolve for %d generations' % gen_budget
    algo = algorithm.pso_gen(gen_budget,0.7298,2.05,2.05,0.05,5,2,4)
    
    # 2B. Set up pso_gen algorithm with racing:
    # Setting gen number to be an arbitrarily large number, let fevals
    # decide when to terminate.
    nr_eval_per_x_racing = nr_eval_per_x
    algo_with_racing = algorithm.pso_gen_racing(1000000,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x_racing,fevals_max)
    # TODO: Use below to check the sanity of racing in factoring out the effect of exceeded fevals
    # algo_with_racing = algorithm.pso_gen_racing(gen_budget,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x_racing,999999999)

    # 3. Run both algorithms and record their performance
    fitnesses = []
    fitnesses_racing = []

    seed = 1234
    for i in range(num_trials):
        print '::: Trial #%d :::' % i
        seed += 177
        fitnesses.append(evolve_and_evaluate(algo, prob_regular, pop_size, prob_orig, seed)[0])
        if SAME_PROB:
            fitnesses_racing.append(evolve_and_evaluate(algo_with_racing, prob_regular, pop_size, prob_orig, seed)[0])
        else:
            fitnesses_racing.append(evolve_and_evaluate(algo_with_racing, prob_single_eval, pop_size, prob_orig, seed)[0])
        print (fitnesses[-1], fitnesses_racing[-1])

    avg_fitness = np.mean(fitnesses)
    avg_fitness_racing = np.mean(fitnesses_racing)

    print '----------------------------------------------'
    print 'Final averaged actual fitness over %d trials:' % num_trials
    print 'pso_gen without racing: %f' % avg_fitness
    print 'pso_gen with racing: %f' % avg_fitness_racing
    print '----------------------------------------------'

    return (avg_fitness, avg_fitness_racing)

if __name__ == '__main__':
    start_experiment(num_trials=20, pop_size=20, nr_eval_per_x=5, fevals_max=500000)
