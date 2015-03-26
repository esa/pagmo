from PyGMO.algorithm import base


class py_cmaes(base):

    """
    Covariance Matrix Adaptation Evolutionary Strategy (Python)
    """

    def __init__(
            self,
            gen=500,
            cc=-1,
            cs=-1,
            c1=-1,
            cmu=-1,
            sigma0=0.5,
            ftol=1e-6,
            xtol=1e-6,
            memory=False,
            screen_output=False):
        """
        Constructs a Covariance Matrix Adaptation Evolutionary Strategy (Python)

        USAGE: algorithm.py_cmaes(gen = 500, cc = -1, cs = -1, c1 = -1, cmu = -1, sigma0=0.5, ftol = 1e-6, xtol = 1e-6, memory = False, screen_output = False)

        NOTE: In our variant of the algorithm, particle memory is used to extract the elite and reinsertion
        is made aggressively ..... getting rid of the worst guy). Also, the bounds of the problem
        are enforced, as to allow PaGMO machinery to work. Fine control on each iteration can be achieved
        by calling the algo with gen=1 (algo state is stored, cmaes will continue at next call ... without
        initializing again all its state!!)

        * gen: number of generations
        * cc: time constant for C cumulation (in [0,1]) if -1 automatic values are set
        * cs: time constant for sigma cumulation (in [0,1]) if -1 automatic values are set
        * c1: learning rate for rank-1 update (in [0,1]) if -1 automatic values are set
        * cmu: learning rate for rank-mu update (in [0,1]) if -1 automatic values are set
        * sigma0: starting step (std)
        * xtol: stopping criteria on the x tolerance
        * ftol: stopping criteria on the f tolerance
        * memory: when True the algorithm preserves memory of covariance, step and more between successive runs
        * screen_output: activates screen_output (output at each generation)
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "This algorithm needs numpy to run. Is numpy installed?")

        if (gen <= 0):
            raise ValueError("gen needs to be > 0")

        if ((cc < 0 or cc > 1) and not cc == -1):
            raise ValueError("cc needs to be in [0,1] or -1 for auto value")

        if ((cs < 0 or cs > 1) and not cc == -1):
            raise ValueError("cs needs to be in [0,1] or -1 for auto value")

        if ((c1 < 0 or c1 > 1) and not cc == -1):
            raise ValueError("c1 needs to be in [0,1] or -1 for auto value")

        if ((cmu < 0 or cmu > 1) and not cc == -1):
            raise ValueError("cmu needs to be in [0,1] or -1 for auto value")

        base.__init__(self)

        # Data members
        self.__cc = cc
        self.__cs = cs
        self.__c1 = c1
        self.__cmu = cmu
        self.__gen = gen
        self.__xtol = xtol
        self.__ftol = ftol
        self.__sigma0 = sigma0
        self.__memory = memory
        self.screen_output = screen_output

        # Algorithm memory
        self.__mean = 0
        self.__variation = 0
        self.__newpop = np.matrix([[1]])
        self.__B = 0
        self.__D = 0
        self.__C = 0
        self.__invsqrtC = 0
        self.__pc = 0
        self.__ps = 0
        self.__counteval = 0
        self.__eigeneval = 0

        np.random.seed()

    def de_operator(self, genotype):
        from numpy import matrix
        from random import randint, random, shuffle
        from copy import deepcopy
        from numpy.random import normal
        F = 0.7
        CR = 0.9
        dim, popsize = genotype.shape
        retval = matrix([[0.0] * dim] * popsize).T
        for i in range(popsize):
            idxs = range(popsize)
            shuffle(idxs)
            r1, r2, r3 = idxs[:3]
            retval[:, i] = deepcopy(genotype[:, i])
            n = randint(0, dim - 1)
            L = 0
            while ((random() < CR) and (L < dim)):
                retval[n, i] = genotype[n, r3] + F * (genotype[n, r1] - genotype[n, r2])
                n = (n + 1) % dim
                L = L + 1
            #retval[:, i] = normal(0.0, 1.0, [dim, 1])
        return retval

    def evolve(self, pop):
        from numpy import matrix, array, log, diag, eye, sqrt, exp, ones
        from numpy.random import normal, random
        from numpy.linalg import norm, eig, inv

        # Let's rename some variables
        prob = pop.problem
        lb = prob.lb
        ub = prob.ub
        dim, cont_dim, int_dim, c_dim, f_dim = prob.dimension, prob.dimension - \
            prob.i_dimension, prob.i_dimension, prob.c_dimension, prob.f_dimension

        # And perform checks on the problem type
        if cont_dim == 0:
            raise ValueError(
                "There is no continuous dimension for CMAES to optimise!!")

        if c_dim > 0:
            raise ValueError(
                "This version of CMAES is not suitable for constrained optimisation")

        if int_dim > 0:
            raise ValueError(
                "The chromosome has an integer part .... this version of CMAES is not able to deal with it")
        if f_dim > 1:
            raise ValueError(
                "The problem is not single objective and CMAES is not suitable to solve it")

        if len(pop) < 5:
            raise ValueError(
                "for CMAES at least 5 individuals in the population are required")

        # Setting sizes .....
        N = dim
        lam = len(pop)
        mu = lam / 2

        # Setting coefficients for Selection
        weights = [log(mu + 0.5) - log(i + 1) for i in range(mu)]
        sumW = sum(weights)
        weights = [w / sumW for w in weights]
        # weights for weighted recombination
        # variance-effectiveness of sum w_i x_i
        mueff = 1.0 / sum(w ** 2 for w in weights)

        # Setting coefficients for Adaptation automatically or to user defined
        # data
        cc = self.__cc
        cs = self.__cs
        c1 = self.__c1
        cmu = self.__cmu
        if self.__cc == -1:
            cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N)
            # t-const for cumulation for C
        if self.__cs == -1:
            cs = (mueff + 2) / (N + mueff + 5)
            # t-const for cumulation for sigma control
        if self.__c1 == -1:
            c1 = 2 / ((N + 1.3) ** 2 + mueff)
            # learning rate for rank-one update of C
        if self.__cmu == -1:
            cmu = 2 * (mueff - 2 + 1 / mueff) / ((N + 2) ** 2 + mueff)
            # and for rank-mu update

        damps = 1 + 2 * max(0, sqrt((mueff - 1) / (N + 1)) - 1) + cs
        # damping for sigma
        # expectation of ||N(0,I)|| == norm(randn(N,1))
        chiN = N ** 0.5 * (1 - 1.0 / (4 * N) + 1.0 / (21 * N ** 2))

        # Initializing and allocating
        if (self.__newpop.shape == (N, lam)) and (self.__memory):
            mean = self.__mean
            variation = self.__variation
            newpop = self.__newpop
            B = self.__B
            D = self.__D
            C = self.__C
            invsqrtC = self.__invsqrtC
            pc = self.__pc
            ps = self.__ps
            counteval = self.__counteval
            eigeneval = self.__eigeneval
        else:
            mean = matrix(pop.champion.x).T
            variation = array([[0.0] * N] * lam)
            newpop = matrix([[0.0] * lam] * N)
            B = matrix(eye(N, N))
            # B defines the coordinate system
            D = ones(N)
            # diagonal D defines the scaling
            C = matrix(eye(N, N))
            # covariance matrix C
            invsqrtC = matrix(eye(N, N))
            # inverse of sqrt(C)
            pc = matrix([[0]] * N)
            ps = matrix([[0]] * N)
            counteval = 0
            eigeneval = 0

        sigma = self.__sigma0

        if self.screen_output:
            print("CMAES 4 PaGMO (Python)\n")
            print("mu: " + str(mu) + " - lambda: " + str(lam) +
                  " - N: " + str(N) + " - muef: " + str(mueff) + "\n")
            print(
                "cc: " +
                str(cc) +
                " - cs: " +
                str(cs) +
                " - c1: " +
                str(c1) +
                " - cmu: " +
                str(cmu) +
                " - sigma: " +
                str(sigma) +
                " - damps: " +
                str(damps) +
                " - chiN: " +
                str(chiN) +
                "\n")

        genotypes = matrix([[0.0] * dim] * lam).T
        # Let's start the algorithm
        for gen in range(self.__gen):

            inversemap = inv(B * diag(D))

            # We transform the chromosomes in the genotype space
            phenotypes = [ind.best_x for ind in pop]
            for i, x in enumerate(phenotypes):
                genotypes[:, i] = inversemap * (matrix(x).T - mean) / sigma

            variation = self.de_operator(genotypes)
            variation = [B * diag(D) * variation[:, i] for i in range(lam)]
            variation = [[j[0, 0] for j in matr] for matr in variation]

            # 1 - We generate and evaluate lam new individuals
            #variation = [B * diag(D) * normal(0, 1, [dim, 1])
            #             for i in range(lam)]
            #variation = [[j[0, 0] for j in matr] for matr in variation]
            for i, d_mu in enumerate(variation):
                newpop[:, i] = mean + sigma * matrix(d_mu).T

            # fixing the bounds
            for row in range(newpop.shape[0]):
                for col in range(newpop.shape[1]):
                    if newpop[row, col] > ub[row]:
                        newpop[row, col] = lb[row] + \
                            random() * (ub[row] - lb[row])
                    elif newpop[row, col] < lb[row]:
                        newpop[row, col] = lb[row] + \
                            random() * (ub[row] - lb[row])

            # insert in population
            for i in range(lam):
                #idx = pop.get_worst_idx()
                pop.set_x(i, [newpop[j, i] for j in range(N)])
            counteval += lam

            # 2 - We extract the elite from this generation
            # a = sorted(pop,lambda x,y: cmp(x.cur_f,y.cur_f))
            elite = [matrix(pop[idx].best_x).T for idx in pop.get_best_idx(mu)]
            # elite = [matrix(ind.cur_x).T for ind in a]
            # elite = elite[:mu]

            # 3 - Compute the new elite mean storing the old one
            meanold = mean
            mean = elite[0] * weights[0]
            for i in range(1, mu):
                mean += elite[i] * weights[i]

            # 4 - Update evolution paths
            ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * \
                invsqrtC * (mean - meanold) / sigma
            hsig = ((ps.T * ps)[0, 0] /
                    (1 - (1 - cs) ** (2.0 * counteval / lam)) / N) < (2.0 + 4.0 /
                                                                      (N + 1))
            hsig = int(hsig)
            pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * \
                (mean - meanold) / sigma

            # 5 - Adapt Covariance Matrix
            Cold = C
            C = (elite[0] - meanold) * (elite[0] - meanold).T * weights[0]
            for i in range(1, mu):
                C += (elite[i] - meanold) * (elite[i] - meanold).T * weights[i]
            C /= sigma ** 2
            C = (1 - c1 - cmu) * Cold + cmu * C + c1 * \
                ((pc * pc.T) + (1 - hsig) * cc * (2 - cc) * Cold)

            # 6 - Adapt sigma
            sigma *= exp((cs / damps) * (norm(ps) / chiN - 1))

            # 7 - Perform eigen-decomposition of C
            # achieve O(N^2)
            if ((counteval - eigeneval) > (lam / (c1 + cmu) / N / 10)):
                eigeneval = counteval
                C = (C + C.T) / 2  # enforce symmetry
                D, B = eig(C)
                # eigen decomposition, B==normalized eigenvectors
                D = [s ** 0.5 for s in D]  # D contains standard deviations now
            # if not (0 in D):                                #Avoids numerical
            # nans skipping evaluation of invsqrtC
                invsqrtC = B * diag([1 / d for d in D]) * B.T

            # 8 - Print to screen if necessary
            if self.screen_output:
                if not(gen % 20):
                    print(
                        "\nGen.\tChampion\tHighest\t\tLowest\t\tVariation\t\tStep")
                print(
                    "%d\t%e\t%e\t%e\t%e\t%e" %
                    (gen, pop.champion.f[0], max(
                        [ind.cur_f[0] for ind in pop]), min(
                        [ind.cur_f[0] for ind in pop]), norm(d_mu), sigma))

            # 9 - Check the exit conditions (every 40 generations)
            if not(gen % 40):
                if (norm(d_mu) < self.__xtol):
                    if self.screen_output:
                        print("Exit condition -- xtol < " + str(self.__xtol))
                    return pop

                tmp = abs(pop[pop.get_worst_idx()].best_f[0] -
                          pop[pop.get_best_idx()].best_f[0])

                if (tmp < self.__ftol):
                    if self.screen_output:
                        print("Exit condition -- ftol < " + str(self.__ftol))
                    return pop

        # Update algorithm memory
        if self.__memory:
            self.__mean = mean
            self.__variation = variation
            self.__newpop = newpop
            self.__B = B
            self.__D = D
            self.__C = C
            self.__invsqrtC = invsqrtC
            self.__pc = pc
            self.__ps = ps
            self.__counteval = counteval
            self.__eigeneval = eigeneval
            self.__sigma0 = sigma

        if self.screen_output:
            print("Exit condition -- iteration > " + str(self.__gen))
        return pop

    def get_name(self):
        return "CMAES (Python)"

    def human_readable_extra(self):
        return "gen=" + str(self.__gen) + " cc=" + str(self.__cc) + " cs=" + str(self.__cs) + " c1=" + str(self.__c1) + \
            " cmu=" + str(self.__cmu) + " sigma0=" + str(self.__sigma0) + " xtol=" + str(self.__xtol) + " ftol=" + str(self.__ftol)
