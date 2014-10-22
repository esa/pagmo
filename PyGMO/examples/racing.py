from PyGMO import *
import random
import numpy as np
import matplotlib.pyplot as plt
import copy


def brute_force_average_f(pop_noisy, num_winner, eval_budget):
    """
    Allocate evenly the evaluation budget to all the individuals.
    The winners will be determined by the averaged objective values.
    (Note: Only applicable to single objective, non-constrained problems)
    """
    pop_n = len(pop_noisy)
    f_mean = [0.0] * pop_n
    cnt = [0] * pop_n
    eval_total = 0

    f_prev = [-999.0] * pop_n
    pop_noisy_local = copy.deepcopy(pop_noisy)

    while eval_total < eval_budget:
        cur_seed = random.randint(0, 100000000)
        pop_noisy_local._problem_reference.seed = cur_seed
        for (i, p) in enumerate(pop_noisy_local):
            newest_f = pop_noisy_local.problem.objfun(p.cur_x)[0]
            f_mean[i] = f_mean[i] + newest_f
            cnt[i] = cnt[i] + 1
            eval_total = eval_total + 1
            if eval_total >= eval_budget:
                break

    f_mean = [f_mean[i] / cnt[i] for i in range(0, pop_n)]

    winners = np.argsort(f_mean)[:num_winner]

    return winners


def brute_force_average_rank(pop_noisy, num_winner, eval_budget):
    """
    Allocate evenly the evaluation budget to all the individuals.
    The winners will be determined by the averaged ranking.
    """
    pop_n = len(pop_noisy)
    rank_sum = [0.0] * pop_n
    eval_total = 0

    pop_noisy_local = copy.deepcopy(pop_noisy)

    while eval_total < eval_budget:
        cur_seed = random.randint(0, 100000000)
        pop_noisy_local._problem_reference.seed = cur_seed
        for (i, p) in enumerate(pop_noisy_local):
            pop_noisy_local.set_x(i, p.cur_x)
        cur_winners = pop_noisy_local.get_best_idx(pop_n)
        for (rank, ind_idx) in enumerate(cur_winners):
            rank_sum[ind_idx] += rank
        eval_total += pop_n

    winners = np.argsort(rank_sum)[:num_winner]

    return winners


def get_capture_rates(
        prob_orig=problem.ackley(10),
        noise=0.3,
        pop_n=20,
        num_winner=4,
        eval_budget=200):
    """
    Returns a list containing the capture rates of different methods.

    Capture rate: Percentage of matching with the ground truth ordering.

    Example: ground_truth = [1,3,5,7,9], winners = [1,2,3,4,5], capture
    rate of winners = 3 / 5 = 0.6
    """

    prob_noisy = problem.noisy(prob_orig, 1, 0, noise)

    pop_orig = population(prob_orig, pop_n)
    pop_noisy = population(prob_noisy)
    for p in pop_orig:
        pop_noisy.push_back(p.cur_x)

    rates = []

    # Ground truth
    winners_orig = pop_orig.get_best_idx(num_winner)

    # Results from different methods:
    winners_racing, fevals = pop_noisy.race(
        num_winner, 0, eval_budget, 0.05, [])
    capture_race = 100.0 * sum([int(p in winners_orig)
                               for p in winners_racing]) / num_winner
    rates.append(capture_race)

    winners_brute_force_rank = brute_force_average_rank(
        pop_noisy, num_winner, eval_budget)
    capture_bf_rank = 100.0 * \
        sum([int(p in winners_orig)
            for p in winners_brute_force_rank]) / num_winner
    rates.append(capture_bf_rank)

    if(prob_noisy.f_dimension == 1 and prob_noisy.c_dimension == 0):
        winners_brute_force_f = brute_force_average_f(
            pop_noisy, num_winner, eval_budget)
        capture_bf_f = 100.0 * \
            sum([int(p in winners_orig)
                for p in winners_brute_force_f]) / num_winner
        rates.append(capture_bf_f)

    return rates


def get_rank_sum_errors(
        prob_orig=problem.ackley(10),
        noise=0.3,
        pop_n=20,
        num_winner=4,
        eval_budget=200):
    """
    Returns a list containing the rank-sum error rates of different methods.

    Rank-sum error rates: The same metric used in the C++ tests.

    The error is defined as the difference between the true ranks sum and the
    returned ranks sum. For example if race returns [1,3,5], then the error is
    (1+3+5) - (0+1+2) = 6. The allowed error is the size of the returned list,
    corresponding, to allow just losing the winner (i.e. [1,2,3] is still valid, but
    [1,2,4] not)
    """

    prob_noisy = problem.noisy(prob_orig, 1, 0, noise)

    pop_orig = population(prob_orig, pop_n)

    # Increase the level of difficulty
    # algo = algorithm.pso(gen=50)
    # pop_orig = algo.evolve(pop_orig)

    # True ordering of the individuals
    winners_orig = pop_orig.get_best_idx(pop_n)

# Individuals are already sorted by their true quality in pop_noisy
    pop_noisy = population(prob_noisy)
    for ind_idx in winners_orig:
        pop_noisy.push_back(pop_orig[ind_idx].cur_x)

# In perfect case, rank sum should be sum([0,1,2,....,num_winner-1])
    ground_truth = ((num_winner - 1) + 1) * (num_winner - 1) / 2

    # Worse possible error
    # normalize_factor = (((pop_n-1)+1)*(pop_n-1) / 2  - ground_truth) - ground_truth

    errors = []

    # Results from different methods:
    winners_racing, fevals = pop_noisy.race(
        num_winner, 0, eval_budget, 0.05, [])
    error_race = sum(winners_racing) - ground_truth
    errors.append(error_race)

    winners_brute_force_rank = brute_force_average_rank(
        pop_noisy, num_winner, eval_budget)
    error_bf_rank = sum(winners_brute_force_rank) - ground_truth
    errors.append(error_bf_rank)

    if(prob_noisy.f_dimension == 1 and prob_noisy.c_dimension == 0):
        winners_brute_force_f = brute_force_average_f(
            pop_noisy, num_winner, eval_budget)
        error_bf_f = sum(winners_brute_force_f) - ground_truth
        errors.append(error_bf_f)

    return errors

# Setting some common parameters for the experimentations
num_trials = 200
final_n = 5
default_noise = 0.5

"""
metric_fn = get_capture_rates;
metric_name = 'Capture rate (%)'
"""

metric_fn = get_rank_sum_errors
metric_name = 'Rank-sum error'


def repeat_and_average(fn, *args):
    s = []
    for i in range(num_trials):
        s.append(fn(*args))
    s = np.array(s)
    return np.mean(s, 0)


def run_varying_noise(prob_orig):
    # --- Set-up A: Test with different noise levels ---
    plt.close()
    noise_levels = [
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]

    start_n = 100
    eval_budget = start_n * 5
    success_rates = []

    for noise in noise_levels:
        print('Setup A: noise = %lf' % noise)
        averaged_rates = repeat_and_average(
            metric_fn,
            prob_orig,
            noise,
            start_n,
            final_n,
            eval_budget)
        success_rates.append(averaged_rates)

    success_rates = np.array(success_rates).T

    plt.ion()
    plt.hold(True)
    for (i, succ) in enumerate(success_rates):
        plt.plot(noise_levels, succ)
    plt.xlabel('Noise level (sigma)')
    plt.ylabel(metric_name)
    plt.title(
        '%s: Varying noise levels\n# of inds: %d -> %d, eval. budget = %d' %
        (prob_orig.get_name(), start_n, final_n, eval_budget))
    plt.legend(
        ('Racing',
         'Brute-force (averaged rank)',
         'Brute-force (averaged f)'),
        loc='best')
    plt.savefig('%s-racing-varying-noise' %
                prob_orig.get_name().replace(' ', ''), formant='png')
    # --- End of set-up A ----


def run_varying_initial_size(prob_orig):
    # --- Set-up B: Test with different initial pop sizes ---
    plt.close()
    start_n_list = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
    success_rates = []

    for start in start_n_list:
        print('Setup B: initial popsize = %d' % start)
        averaged_rates = repeat_and_average(metric_fn, prob_orig,
                                            default_noise, start, final_n,
                                            start * 5)
        success_rates.append(averaged_rates)

    success_rates = np.array(success_rates).T

    plt.ion()
    plt.hold(True)
    for (i, succ) in enumerate(success_rates):
        plt.plot(start_n_list, succ)
    plt.xlabel('Initial population size')
    plt.ylabel(metric_name)
    plt.title(
        '%s: Varying initial pop size\n# of winners: %d, eval. budget = popsize * 5' %
        (prob_orig.get_name(), final_n))
    plt.legend(
        ('Racing',
         'Brute-force (averaged rank)',
         'Brute-force (averaged f)'),
        loc='best')
    plt.savefig('%s-racing-varying-initialpopsize' %
                prob_orig.get_name().replace(' ', ''), formant='png')
    # --- End of set-up B ---


def run_varying_eval_budget(prob_orig):
    # --- Set-up C: Test with different evaluation budget ---
    plt.close()
    eval_budget_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    start_n = 100
    success_rates = []

    for eval_budget in eval_budget_list:
        print('Setup C: evaluation budget = %d' % eval_budget)
        averaged_rates = repeat_and_average(metric_fn, prob_orig,
                                            default_noise, start_n, final_n,
                                            eval_budget)
        success_rates.append(averaged_rates)

    success_rates = np.array(success_rates).T

    plt.ion()
    plt.hold(True)
    for (i, succ) in enumerate(success_rates):
        plt.plot(eval_budget_list, succ)
    plt.xlabel('Allowed evaluation budget')
    plt.ylabel(metric_name)
    plt.title('%s: Varying eval. budget\n# of inds: %d -> %d' %
              (prob_orig.get_name(), start_n, final_n))
    plt.legend(
        ('Racing',
         'Brute-force (averaged rank)',
         'Brute-force (averaged f)'),
        loc='best')
    plt.savefig('%s-racing-varying-budget' %
                prob_orig.get_name().replace(' ', ''), formant='png')
    # --- End of set-up C ---

if __name__ == '__main__':
    prob_orig = problem.ackley(10)
    # prob_orig = problem.cec2006(5)
    # prob_orig = problem.zdt(1,10)

    run_varying_noise(prob_orig)
    run_varying_initial_size(prob_orig)
    run_varying_eval_budget(prob_orig)
