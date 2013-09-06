from PyGMO import *
import numpy as np
import matplotlib.pyplot as plt


prob_orig = problem.ackley(10)

nr_eval_per_x = 3
noise_level = 0.5
prob_racing = problem.noisy(prob_orig, trials = 1, param_second = noise_level)
prob_regular = problem.noisy(prob_orig, trials = nr_eval_per_x, param_second = noise_level)

fevals_max = 100000
#fevals_max = 50000
pop_size = 20

num_trials = 20

def get_pso_gen():
    print 'pso_gen without racing'
    # Each generation of pso_gen requires 2*pop_size*nr_eval_per_x
    # evaluations. Ignoring the cost of initialization here.
    gen_budget = fevals_max/(2*pop_size*nr_eval_per_x)
    print 'Gen_budget is', gen_budget
    algo = algorithm.pso_gen(gen_budget,0.7298,2.05,2.05,0.05,5,2,4)
    return algo

def get_pso_gen_racing():
    print 'pso_gen with racing'
    # Setting gen number to be an arbitrarily large number, let fevals
    # decide when to terminate.
    gen_budget = 1000000
    algo = algorithm.pso_gen_racing(gen_budget,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x,fevals_max)
    return algo

def run_algo(algo, prob, seed):

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

if __name__ == '__main__':
    fitnesses = []
    fitnesses_racing = []

    algo = get_pso_gen()
    algo_with_racing = get_pso_gen_racing()

    seed = 1234
    for i in range(num_trials):
        print '::: Trial #%d :::' % i
        seed += 177
        fitnesses.append(run_algo(algo, prob_regular, seed)[0])
        fitnesses_racing.append(run_algo(algo_with_racing, prob_racing, seed)[0])

    avg_fitness = np.mean(fitnesses)
    avg_fitness_racing = np.mean(fitnesses_racing)

    print '----------------------------------------------'
    print 'Final averaged actual fitness over %d trials:' % num_trials
    print 'pso_gen without racing: %f' % avg_fitness
    print 'pso_gen with racing: %f' % avg_fitness_racing
    print '----------------------------------------------'
