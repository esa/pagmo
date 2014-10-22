from PyGMO import *

# We instantiate the optimization problem as usual. We will choose a 4-dimensional 2-objective
# problem.
prob = problem.kur(4)

# Now we instantiate the analysis class. In this case we will sample 500 points via sobol method
# and choose output to file.
inspector = util.analysis(prob, 500, method='sobol', output_to_file=True)

# **********************************************************************************************
pop = population(prob, 500)
isl = island(algorithm.nsga_II(), pop)
isl.evolve(20)
pop_after_20_gen = isl.population

inspector_kur_20 = util.analysis(pop_after_20_gen, output_to_file=True)

# **********************************************************************************************
inspector.f_distribution()

inspector.f_linearity_convexity()

# We will compute all regression properties for linear, quadratic and cubic regression, with
# and without interaction.
inspector.f_regression(degree=[1, 1, 2, 2, 3, 3], interaction=[False, True, False, True, False, True])

inspector.f_sensitivity()

# The levelset test is very slow, we do not recommend launching it out of development purposes.
# Anyway this is an example so here it is. We will launch it using percentile 45 as threshold
# instead of the default 50, and only with 2-fold crossvalidation for tuning of svm.
# inspector.levelset(threshold=45, k_tune=2)

# For the local searches, we will choose Nelder-Mead simplex algorithm. As it is a multi-objective
# problem, we will also choose a decomposition method, for instance weighted with uniform weights.
# To keep it short, we will only show properties of the best 3 clusters instead of the default 10.
# We will also generate a scatter plot for dimensions x1 and x2.
inspector.local_search(clusters_to_show=3, scatter_plot_dimensions=[1, 2],
                       algo=algorithm.gsl_nm(), decomposition_method='weighted', weights='uniform')

# Being only a 2-objective problem, we probably will not obtain dimensionality reduction, but we
# will launch this test all the same to look into the f-correlation matrix.
inspector.f_correlation()

# **********************************************************************************************
prob = problem.cec2006(5)
inspector = util.analysis(prob, 1000, 'lhs', output_to_file=True)

inspector.c_feasibility()
inspector.c_linearity()
inspector.c_regression(degree=[1, 1, 2, 2], interaction=[False, True, False, True])
inspector.c_sensitivity()
