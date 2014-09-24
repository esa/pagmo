Dr. PyGMO - Analysis Module
===========================

A Quick Look
------------

Optimization problems are often provided as black box functions. A number of exploratory techniques based on sampling of the search space are 
available to help to gain knowledge about the problem. In particular we are interested about problem characteristics that can be critical for 
the performance of the algorithms or those that can help the user to reformulate the problem in a more efficient way.
You may follow the :ref:`landscape_analysis_with_DrPyGMO` tutorial for a better insight on its use. Hereafter the list of test available in the 
current version of the module *PyGMO.util.analysis*.

**NOTE:** The packages, pandas and numpy (for plots), scikit-learn (for levelset analysis) need to be installed for having the complete set of functionalities.

Tests
^^^^^
======================================= =============================================================== ==========================================================
Common Name                          	Name in PyGMO							Comments
======================================= =============================================================== ==========================================================
Sample                              	:py:func:`PyGMO.util.analysis.sample`			        Sampling of the search space
F-distribution	                     	:py:func:`PyGMO.util.analysis.f_distribution`                   Distribution of sampled fitness values
F-linearity and convexity  	    	:py:func:`PyGMO.util.analysis.f_linearity_convexity`	        Probability of linearity and convexity of fitness function
F-regression			     	:py:func:`PyGMO.util.analysis.f_regression`			Polynomial regression on fitness function
F-correlation			     	:py:func:`PyGMO.util.analysis.f_correlation`			Fitness dimensionality reduction via PCA
F-sensitivity			     	:py:func:`PyGMO.util.analysis.f_sensitivity`			Jacobian and Hessian of fitness function
Levelset			     	:py:func:`PyGMO.util.analysis.levelset`				SVM binary classification
Local search			     	:py:func:`PyGMO.util.analysis.local_search`			Clustering of local minima
C-feasibility				:py:func:`PyGMO.util.analysis.c_feasibility`			Constraint effectiveness and redundancy
C-linearity				:py:func:`PyGMO.util.analysis.c_linearity`			Probability of linearity of constraints function
C-regression				:py:func:`PyGMO.util.analysis.c_regression`			Polynomial regression on constraints function
C-sensitivity				:py:func:`PyGMO.util.analysis.c_sensitivity`			Jacobian of constraints function
======================================= =============================================================== ==========================================================

Detailed Documentation
----------------------

.. autoclass:: PyGMO.util.analysis

----------------------

   .. automethod:: PyGMO.util.analysis.__init__()

----------------------

   .. automethod:: PyGMO.util.analysis.sample()

----------------------

   .. automethod:: PyGMO.util.analysis.f_distribution()

----------------------

   .. automethod:: PyGMO.util.analysis.f_linearity_convexity()

----------------------

   .. automethod:: PyGMO.util.analysis.f_regression()

----------------------

   .. automethod:: PyGMO.util.analysis.f_correlation()

----------------------

   .. automethod:: PyGMO.util.analysis.f_sensitivity()

----------------------

   .. automethod:: PyGMO.util.analysis.levelset()

----------------------

   .. automethod:: PyGMO.util.analysis.local_search()

----------------------

   .. automethod:: PyGMO.util.analysis.c_feasibility()

----------------------

   .. automethod:: PyGMO.util.analysis.c_linearity()

----------------------

   .. automethod:: PyGMO.util.analysis.c_regression()

----------------------

   .. automethod:: PyGMO.util.analysis.c_sensitivity()

----------------------

   .. automethod:: PyGMO.util.analysis._scale_sample()

----------------------

   .. automethod:: PyGMO.util.analysis._skew()

----------------------

   .. automethod:: PyGMO.util.analysis._kurtosis()

----------------------

   .. automethod:: PyGMO.util.analysis._mean()

----------------------

   .. automethod:: PyGMO.util.analysis._var()

----------------------

   .. automethod:: PyGMO.util.analysis._std()

----------------------

   .. automethod:: PyGMO.util.analysis._ptp()

----------------------

   .. automethod:: PyGMO.util.analysis._percentile()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_f_distr()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_x_pcp()

----------------------

   .. automethod:: PyGMO.util.analysis._n_peaks_f()

----------------------

   .. automethod:: PyGMO.util.analysis._p_lin_conv()

----------------------

   .. automethod:: PyGMO.util.analysis._regression_coefficients()

----------------------

   .. automethod:: PyGMO.util.analysis._regression_properties()

----------------------

   .. automethod:: PyGMO.util.analysis._regression_press()

----------------------

   .. automethod:: PyGMO.util.analysis._build_polynomial()

----------------------

   .. automethod:: PyGMO.util.analysis._regression_predict()

----------------------

   .. automethod:: PyGMO.util.analysis._f_correlation()

----------------------

   .. automethod:: PyGMO.util.analysis._perform_f_pca()

----------------------

   .. automethod:: PyGMO.util.analysis._get_gradient()

----------------------

   .. automethod:: PyGMO.util.analysis._richardson_gradient()

----------------------

   .. automethod:: PyGMO.util.analysis._get_hessian()

----------------------

   .. automethod:: PyGMO.util.analysis._richardson_hessian()

----------------------

   .. automethod:: PyGMO.util.analysis._grad_properties()

----------------------

   .. automethod:: PyGMO.util.analysis._hess_properties()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_gradient_sparsity()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_gradient_pcp()

----------------------

   .. automethod:: PyGMO.util.analysis._get_local_extrema()

----------------------

   .. automethod:: PyGMO.util.analysis._cluster_local_extrema()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_local_cluster_pcp()

----------------------

   .. automethod:: PyGMO.util.analysis.plot_local_cluster_scatter()

----------------------

   .. automethod:: PyGMO.util.analysis._svm()

----------------------

   .. automethod:: PyGMO.util.analysis._svm_p_values()

----------------------

   .. automethod:: PyGMO.util.analysis._c_lin()

----------------------

   .. automethod:: PyGMO.util.analysis._compute_constraints()

----------------------

   .. automethod:: PyGMO.util.analysis._c_effectiveness()

----------------------

   .. automethod:: PyGMO.util.analysis._ic_redundancy()
 
