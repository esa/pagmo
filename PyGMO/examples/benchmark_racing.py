from PyGMO import *
import numpy as np
import matplotlib.pyplot as plt

import copy

# stochastic_type = 'NOISY'
stochastic_type = 'ROBUST'
base_problem = problem.ackley(10)


class post_eval:

    """
    Obtain the post-evaluated fitness via repeated averaing over different seeds.
    """

    def __init__(self, post_eval_prob, post_eval_n=500, seed=5):
        self.post_eval_prob = post_eval_prob
        self.post_eval_n = post_eval_n
        self.seed = seed

    def objfun(self, x):
        post_f = 0
        np.random.seed(self.seed)
        for i in range(self.post_eval_n):
            self.post_eval_prob.seed = np.random.randint(1000000)
            post_f += self.post_eval_prob.objfun(x)[0] / \
                float(self.post_eval_n)
        return (post_f,)


def start_experiment(
        num_trials=20,
        pop_size=40,
        fevals_max=100000,
        nr_eval_per_x=40,
        noise_level=0.05,
        seed=123):

    # 1. Set up the problem
    if(stochastic_type == 'NOISY'):
        prob_single_eval = problem.noisy(
            base_problem,
            trials=1,
            param_second=noise_level,
            noise_type=problem.noisy.noise_distribution.UNIFORM)
        prob_regular = problem.noisy(
            base_problem,
            trials=nr_eval_per_x,
            param_second=noise_level,
            noise_type=problem.noisy.noise_distribution.UNIFORM)
    else:
        prob_single_eval = problem.robust(
            base_problem, trials=1, rho=noise_level)
        prob_regular = problem.robust(
            base_problem, trials=nr_eval_per_x, rho=noise_level)

    prob_post_eval = post_eval(prob_single_eval, post_eval_n=500)

    """
    Notes for SAME_PROB: Both algorithm will operate on the same version of the
    problem (same n_eval). pso_gen will evolve for fevals_max/2*pop_size times;
    pso_gen_racing will evolve until fevals_max is hit. In this case a single
    feval count in is referred to a single call to objfun() of the problem
    (with n_eval as 10).
    """
    SAME_PROB = False

    # 2A. Set up pso_gen algorithm without racing:
    # Each generation of pso_gen requires 2*pop_size*nr_eval_per_x
    # evaluations. Ignoring the cost of initialization here.
    # NOTE: No need to scale down if both algo has the same version of problem
    if SAME_PROB:
        gen_budget = fevals_max / (2 * pop_size)
    else:
        gen_budget = fevals_max / (2 * pop_size * nr_eval_per_x)
    print('Non-racing pso gen will evolve for %d generations' % gen_budget)
    algo_psogen = algorithm.pso_gen(
        gen_budget, 0.7298, 2.05, 2.05, 0.05, 5, 2, 4)

    # 2B. Set up pso_gen algorithm with racing:
    # Setting gen number to be an arbitrarily large number, let fevals
    # decide when to terminate.
    nr_eval_per_x_racing = nr_eval_per_x
    algo_psogen_racing = algorithm.pso_gen_racing(
        1000000,
        0.7298,
        2.05,
        2.05,
        0.05,
        5,
        2,
        4,
        nr_eval_per_x_racing,
        fevals_max)
    # TODO: Use below to check the sanity of racing in factoring out the effect of exceeded fevals
    # algo_with_racing = algorithm.pso_gen_racing(gen_budget,0.7298,2.05,2.05,0.05,5,2,4,nr_eval_per_x_racing,999999999)

    # 3. Run both algorithms and record their performance
    if SAME_PROB:
        algo_prob_pairs = [
            (algo_psogen, prob_regular), (algo_psogen_racing, prob_regular)]
    else:
        algo_prob_pairs = [
            (algo_psogen,
             prob_regular),
            (algo_psogen_racing,
             prob_single_eval)]

    post_evaluated_fitnesses = []

    np.random.seed(seed)
    for i in range(num_trials):
        print('::: Trial #%d :::' % i)
        results = []
        seed += np.random.randint(100000)
        for algo, prob in algo_prob_pairs:
            algo.reset_rngs(seed)
            # Seed used to ensure both algorithm evolves an identical
            # population
            pop = population(prob, pop_size, seed)
            pop = algo.evolve(pop)
            # winner_idx = pop.race(1)[0][0];
            # print("race winner", winner_idx, "vs champion idx", pop.get_best_idx())
            # champion_true_fitness = prob_orig.objfun(pop[winner_idx].cur_x)
            champion_true_fitness = prob_post_eval.objfun(pop.champion.x)[0]
            # print('Final champion =', champion_true_fitness)
            results.append(champion_true_fitness)
        print(results)
        post_evaluated_fitnesses.append(results)

    post_evaluated_fitnesses = list(zip(*post_evaluated_fitnesses))

    averaged_no_racing = np.mean(post_evaluated_fitnesses[0])
    averaged_racing = np.mean(post_evaluated_fitnesses[1])

    print('----------------------------------------------')
    print('Final averaged actual fitness over %d trials:' % num_trials)
    print('pso_gen without racing: %f' % averaged_no_racing)
    print('pso_gen with racing: %f' % averaged_racing)
    print('----------------------------------------------')

    return (averaged_no_racing, averaged_racing)


def vary_nr_eval_per_x(default_params):

    pars = copy.deepcopy(default_params)

    param_list = list(range(3, 20, 2))
    f_no_racing_list = []
    f_racing_list = []
    for n in param_list:
        pars['nr_eval_per_x'] = n
        f_no_racing, f_racing = start_experiment(**pars)
        f_no_racing_list.append(f_no_racing)
        f_racing_list.append(f_racing)
    plt.ion()
    plt.figure()
    plt.plot(param_list, f_racing_list, '-o')
    plt.plot(param_list, f_no_racing_list, '-s')
    plt.legend(['PSO racing', 'PSO without racing'])
    plt.xlabel('nr_eval_per_x')
    plt.ylabel('Post-evaluated fitness')
    prob_stat = '%s-%s' % (stochastic_type, base_problem.get_name())
    plt.title('%s\nPSO: With/without racing (fevals=%d) (%d trials)' %
              (prob_stat, pars['fevals_max'], pars['num_trials']))
    # plt.savefig('%s-psogenracing-nr_eval_per_x.png' % prob_stat)


def vary_neighbourhood_size(default_params):

    pars = copy.deepcopy(default_params)

    param_list = np.linspace(0.01, 0.2, num=20)
    f_no_racing_list = []
    f_racing_list = []
    for p in param_list:
        pars['noise_level'] = p
        f_no_racing, f_racing = start_experiment(**pars)
        f_no_racing_list.append(f_no_racing)
        f_racing_list.append(f_racing)
    plt.ion()
    plt.figure()
    plt.plot(param_list, f_racing_list, '-o')
    plt.plot(param_list, f_no_racing_list, '-s')
    plt.legend(['PSO racing', 'PSO without racing'], loc='best')
    plt.xlabel('Robust\'s neighbourhood size')
    plt.ylabel('Post-evaluated fitness')
    prob_stat = '%s-%s' % (stochastic_type, base_problem.get_name())
    plt.title('%s\nPSO: With/without racing (fevals=%d) (%d trials)' %
              (prob_stat, pars['fevals_max'], pars['num_trials']))
    # plt.savefig('%s-psogenracing-robust_neighbourhood_small.png' % prob_stat)


def vary_fevals_budget(num_trials=20, nr_eval_per_x=10, nb_size=0.5):

    pars = copy.deepcopy(default_params)

    param_list = list(range(10000, 200000, 20000))
    f_no_racing_list = []
    f_racing_list = []
    for fevals_max in param_list:
        pars['fevals_max'] = fevals_max
        f_no_racing, f_racing = start_experiment(**pars)
        f_no_racing_list.append(f_no_racing)
        f_racing_list.append(f_racing)
    plt.ion()
    plt.figure()
    plt.plot(param_list, f_racing_list, '-o')
    plt.plot(param_list, f_no_racing_list, '-s')
    plt.legend(['PSO racing', 'PSO without racing'], loc='best')
    plt.xlabel('Evaluation budget (# of fevals)')
    plt.ylabel('Post-evaluated fitness')
    prob_stat = '%s-%s' % (stochastic_type, base_problem.get_name())
    plt.title(
        '%s\nPSO: With/without racing (neighbourhood size = %.2f) (%d trials)' %
        (prob_stat, pars['noise_level'], pars['num_trials']))
    # plt.savefig('%s-psogenracing-robust_fevals.png' % prob_stat)

if __name__ == '__main__':
    # start_experiment(num_trials=20, pop_size=20, nr_eval_per_x=20, fevals_max=200000)

    default_params = dict(
        num_trials=10,
        pop_size=20,
        nr_eval_per_x=10,
        fevals_max=100000,
        noise_level=0.3)

    vary_nr_eval_per_x(default_params)
    vary_neighbourhood_size(default_params)
    # vary_fevals_budget(default_params)
