from PyGMO.algorithm import base


class py_cross_entropy(base):

    """
    Cross-Entropy algorithm (Python)
    """

    def __init__(
            self,
            gen=500,
            elite=0.5,
            scale=0.3,
            variant=1,
            screen_output=False):
        """
        Constructs a Cross-Entropy Algorithm (Python)

        USAGE: algorithm.py_cross_entropy(gen = 1, elite = 0.5, scale = 0.2, variant=1, screen_output = False))

        NOTE: A multivariate normal distribution is used.
              The first sample is centered around the population champion.
              Covariance matrix and mean is evaluated using ind.best_x

        * gen: number of generations
        * elite: fraction of the population considered as elite (in (0,1])
        * scale: scaling factor for the estimated covariance matrix
        * variant: algoritmic variant to use (one of [1,2])
                 1. 'Canonical' - Covariance Matrix is evaluated as sum (x_(i+1)-mu_i)^T (x_(i+1)-mu_i)
                 2. 'Dario's' - Covariance Matrix is evaluated as   sum (x_(i+1)-mu_i^T)^T (x_(i+1)-mu_i^T)
        * screen_output: activates screen_output (output at each generation)
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "This algorithm needs numpy to run. Is numpy installed?")

        base.__init__(self)
        self.__gen = gen
        self.__elite = elite
        self.__scale = scale
        self.__screen_output = screen_output
        self.__weights = []
        self.__variant = variant
        np.random.seed()

    def evolve(self, pop):
        from numpy import matrix, array, log, diag
        from numpy.random import multivariate_normal, random, normal
        from numpy.linalg import norm, cholesky, LinAlgError, eig
        import matplotlib.pyplot as pl

        # Let's rename some variables
        prob = pop.problem
        lb = prob.lb
        ub = prob.ub
        dim, cont_dim, int_dim, c_dim = prob.dimension, prob.dimension - \
            prob.i_dimension, prob.i_dimension, prob.c_dimension

        # And perform checks on the problem type
        if cont_dim == 0:
            raise ValueError(
                "There is no continuous dimension for cross_entropy to optimise!!")

        if c_dim > 0:
            raise ValueError(
                "This version of cross_entropy is not suitable for constrained optimisation")

        if int_dim > 0:
            raise ValueError(
                "The chromosome has an integer part .... this version of cross_entropy is not able to deal with it")

        # We then check that the elite is not empty
        n_elite = int(len(pop) * self.__elite)
        if n_elite == 0:
            raise ValueError(
                "Elite contains no individuals ..... maybe increase the elite parameter?")

        # If the incoming population is empty ... do nothing
        np = len(pop)
        if np == 0:
            return population

        # Let's start the algorithm
        mu = matrix(pop.champion.x)
        C = matrix([[0] * n_elite] * n_elite)
        variation = array([[0.0] * dim] * np)
        newpop = array([[0.0] * dim] * np)

        self.__weights = [log(n_elite + 0.5) - log(i + 1)
                          for i in range(n_elite)]  # recombination weights
        # normalize recombination weights array
        self.__weights = [w / sum(self.__weights) for w in self.__weights]

        for gen in range(self.__gen):

            # 1 - We extract the elite from this generation (NOTE: we use
            # best_f to rank)
            elite = [matrix(pop[idx].best_x)
                     for idx in pop.get_best_idx(n_elite)]

            pl.plot(0, 0.1, 'og')
            for ind in elite:
                pl.plot(ind[0, 0], ind[0, 1], 'or')
            pl.show()
            input()

            # 2 - We evaluate the Covariance Matrix
            if self.__variant == 1:
                # as least square estimator of the elite (with mean mu)
                C = (elite[0] - mu).T * (elite[0] - mu) * self.__weights[0]
                for i in range(1, n_elite):
                    C = C + (elite[i] - mu).T * \
                        (elite[i] - mu) * self.__weights[i]
                    1 / 0
            if self.__variant == 2:
                # using Dario's method
                mu = mu.T
                C = (elite[0] - mu).T * (elite[0] - mu) * self.__weights[0]
                for i in range(1, n_elite):
                    C = C + (elite[i] - mu).T * \
                        (elite[i] - mu) * self.__weights[i]
            # C = C / n_elite

            # 3 - We compute the new elite mean
            mu = elite[0] * self.__weights[0]
            for i in range(1, n_elite):
                mu = mu + elite[i] * self.__weights[i]
            pl.plot(mu[0, 0], mu[0, 1], 'ob')
            input()

            # 4 - We generate the new sample
            variation = multivariate_normal([0] * dim, C, [np])
            # eigen decomposition, B==normalized eigenvectors, O(N**3)
            D, B = eig(C)
            D = [d ** 0.5 for d in D]  # D contains standard deviations now
            variation = [B * diag(D) * normal(0, 1, [dim, 1])
                         for i in range(np)]
            variation = [[j[0, 0] for j in matr] for matr in variation]

            for i, d_mu in enumerate(variation):
                newpop[i] = mu + d_mu * self.__scale
                pl.plot(newpop[i][0], newpop[i][1], 'ok')
            pl.show()
            input()

            # 5 - We fix it within the bounds
            for row in range(newpop.shape[0]):
                for col in range(newpop.shape[1]):
                    if newpop[row, col] > ub[col]:
                        newpop[row, col] = lb[col] + \
                            random() * (ub[col] - lb[col])
                    elif newpop[row, col] < lb[col]:
                        newpop[row, col] = lb[col] + \
                            random() * (ub[col] - lb[col])

            # 6 - And perform reinsertion
            for i in range(np):
                pop.set_x(i, newpop[i])

            # 7 - We print to screen if necessary
            if self.__screen_output:
                if not(gen % 20):
                    print("\nGen.\tChampion\tHighest\t\tLowest\t\tVariation")
                print(
                    "%d\t%e\t%e\t%e\t%e" %
                    (gen, pop.champion.f[0], max(
                        [ind.cur_f[0] for ind in pop]), min(
                        [ind.cur_f[0] for ind in pop]), norm(d_mu)))
        return pop

    def get_name(self):
        return "Cross Entropy (Python)"

    def human_readable_extra(self):
        return "gen=" + str(self.__gen) + " elite fraction=" + str(self.__elite) + \
            " covariance scaling=" + str(self.__scale) + " variant=" + str(self.__variant)
