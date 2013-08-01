from PyGMO import *
import numpy as np
import matplotlib.pyplot as plt

prob_orig = problem.ackley(10)
prob = problem.noisy(prob_orig, param_second = 0.5)
fevals_max = 200000
pop_size = 20
pop_seed = 1234

num_trials = 30

def get_pso_gen():
    print 'pso_gen without racing'
    algo = algorithm.pso_gen(100000,0.7298,2.05,2.05,0.05,5,2,4,False,fevals_max)
    return algo

def get_pso_gen_racing():
    print 'pso_gen with racing'
    algo = algorithm.pso_gen(100000,0.7298,2.05,2.05,0.05,5,2,4,True,fevals_max)
    return algo

def run_algo():

    pop = population(prob, pop_size, pop_seed)

    pop = algo.evolve(pop)

    champion_fitness = prob_orig.objfun(pop.champion.x)

    print 'Final champion =', champion_fitness

    return champion_fitness

if __name__ == '__main__':
    fitnesses = []

    algo = get_pso_gen()
    #algo = get_pso_gen_racing()

    for i in range(num_trials):
        print '::: Trial #%d :::' % i
        fitnesses.append(run_algo()[0])
        pop_seed += 177

    avg_fitness = np.mean(fitnesses)
    print 'Final averaged over %d trials =' % num_trials, avg_fitness
