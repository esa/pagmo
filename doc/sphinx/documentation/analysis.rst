Analysis
========

Tests
^^^^^
.. _analysis_tests:
======================================= =============================================================== =========================================
Common Name                          	Name in PyGMO							Comments
======================================= =============================================================== =========================================
Sample                              	:py:func:`PyGMO.util.analysis.sample`			
F-distribution	                     	:py:func:`PyGMO.util.analysis.f_distribution`            
F-linearity and convexity  	    	:py:func:`PyGMO.util.analysis.f_linearity_convexity`	
F-regression			     	:py:func:`PyGMO.util.analysis.f_regression`			Polynomial regression on F
F-sensitivity			     	:py:func:`PyGMO.util.analysis.f_sensitivity`			Jacobian and Hessian of F
Levelset			     	:py:func:`PyGMO.util.analysis.levelset`				SVM binary classification
Local search			     	:py:func:`PyGMO.util.analysis.local_search`			Clustering of local extrema
F-correlation			     	:py:func:`PyGMO.util.analysis.f_correlation`			Fitness dimensionality reduction via PCA
C-feasibility				:py:func:`PyGMO.util.analysis.c_feasibility`			Constraint effectiveness and redundancy
C-linearity				:py:func:`PyGMO.util.analysis.c_linearity`		
C-regression				:py:func:`PyGMO.util.analysis.c_regression`			Polynomial regression on C
C-sensitivity				:py:func:`PyGMO.util.analysis.c_sensitivity`			Jacobian of C
======================================= =============================================================== =========================================

Detailed Documentation
----------------------

.. autoclass:: PyGMO.util.analysis

----------------------

   .. automethod:: PyGMO.util.analysis.__init__

----------------------

   .. automethod:: PyGMO.util.analysis.sample

----------------------

   .. automethod:: PyGMO.util.analysis.f_distribution

----------------------

   .. automethod:: PyGMO.util.analysis.f_linearity_convexity

----------------------

   .. automethod:: PyGMO.util.analysis.f_regression

----------------------

   .. automethod:: PyGMO.util.analysis.f_sensitivity

----------------------

   .. automethod:: PyGMO.util.analysis.levelset

----------------------

   .. automethod:: PyGMO.util.analysis.local_search

----------------------

   .. automethod:: PyGMO.util.analysis.f_correlation

----------------------

   .. automethod:: PyGMO.util.analysis.c_feasibility

----------------------

   .. automethod:: PyGMO.util.analysis.c_linearity

----------------------

   .. automethod:: PyGMO.util.analysis.c_regression

----------------------

   .. automethod:: PyGMO.util.analysis.c_sensitivity

----------------------

   .. automethod:: PyGMO.util.analysis._scale_sample

----------------------

   .. automethod:: PyGMO.util.analysis._skew

----------------------

   .. automethod:: PyGMO.util.analysis._kurtosis

----------------------

   .. automethod:: PyGMO.util.analysis._mean

----------------------

   .. automethod:: PyGMO.util.analysis._var

----------------------

   .. automethod:: PyGMO.util.analysis._std

----------------------

   .. automethod:: PyGMO.util.analysis._ptp

----------------------

   .. automethod:: PyGMO.util.analysis._percentile

----------------------

   .. automethod:: PyGMO.util.analysis.plot_f_distr

----------------------

   .. automethod:: PyGMO.util.analysis.plot_x_pcp

----------------------

   .. automethod:: PyGMO.util.analysis._n_peaks_f

----------------------

   .. automethod:: PyGMO.util.analysis._p_lin_conv

----------------------

   .. automethod:: PyGMO.util.analysis._regression_coefficients

----------------------

   .. automethod:: PyGMO.util.analysis._regression_properties

----------------------

   .. automethod:: PyGMO.util.analysis._regression_press

----------------------

   .. automethod:: PyGMO.util.analysis._build_polynomial

----------------------

   .. automethod:: PyGMO.util.analysis._regression_predict

----------------------

   .. automethod:: PyGMO.util.analysis._f_correlation

----------------------

   .. automethod:: PyGMO.util.analysis._perform_f_pca

----------------------

   .. automethod:: PyGMO.util.analysis._get_gradient

----------------------

   .. automethod:: PyGMO.util.analysis._richardson_gradient

----------------------

   .. automethod:: PyGMO.util.analysis._get_hessian

----------------------

   .. automethod:: PyGMO.util.analysis._richardson_hessian

----------------------

   .. automethod:: PyGMO.util.analysis._grad_properties

----------------------

   .. automethod:: PyGMO.util.analysis._hess_properties

----------------------

   .. automethod:: PyGMO.util.analysis.plot_gradient_sparsity

----------------------

   .. automethod:: PyGMO.util.analysis.plot_gradient_pcp

----------------------

   .. automethod:: PyGMO.util.analysis._get_local_extrema

----------------------

   .. automethod:: PyGMO.util.analysis._cluster_local_extrema

----------------------

   .. automethod:: PyGMO.util.analysis.plot_local_cluster_pcp

----------------------

   .. automethod:: PyGMO.util.analysis.plot_local_cluster_scatter

----------------------

   .. automethod:: PyGMO.util.analysis._svm

----------------------

   .. automethod:: PyGMO.util.analysis._svm_p_values

----------------------

   .. automethod:: PyGMO.util.analysis._c_lin

----------------------

   .. automethod:: PyGMO.util.analysis._compute_constraints

----------------------

   .. automethod:: PyGMO.util.analysis._c_effectiveness

----------------------

   .. automethod:: PyGMO.util.analysis._ic_redundancy

----------------------
 
