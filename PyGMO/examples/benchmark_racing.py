from PyGMO import *
import numpy as np
import matplotlib.pyplot as plt

prob_orig = problem.ackley(10)
prob = problem.noisy(prob_orig, param_second = 0.5)
fevals_max = 100000
#fevals_max = 50000
pop_size = 40
nr_eval_per_x = 5

num_trials = 20

def get_pso_gen():
    print 'pso_gen without racing'
    algo = algorithm.pso_gen(100000,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x,False,fevals_max)
    return algo

def get_pso_gen_racing():
    print 'pso_gen with racing'
    algo = algorithm.pso_gen(100000,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x,True,fevals_max)
    return algo

def run_algo(algo, seed):

    algo.reset_rngs(seed)

    pop = population(prob, pop_size, seed)

    pop = algo.evolve(pop)

    winner_idx = pop.race(1)[0][0];
    #print "race winner", winner_idx, "vs champion idx", pop.get_best_idx()
    champion_true_fitness = prob_orig.objfun(pop[winner_idx].cur_x)

    #champion_true_fitness = prob_orig.objfun(pop.champion.x)
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
        fitnesses.append(run_algo(algo, seed)[0])
        fitnesses_racing.append(run_algo(algo_with_racing, seed)[0])

    avg_fitness = np.mean(fitnesses)
    avg_fitness_racing = np.mean(fitnesses_racing)

    print '----------------------------------------------'
    print 'Final averaged actual fitness over %d trials:' % num_trials
    print 'pso_gen without racing: %f' % avg_fitness
    print 'pso_gen with racing: %f' % avg_fitness_racing
    print '----------------------------------------------'
