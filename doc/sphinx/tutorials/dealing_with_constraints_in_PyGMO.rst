.. _dealing_with_constraints_in_PyGMO:

================================================================
Dealing with constraints in PyGMO
================================================================

.. toctree::
   :maxdepth: 2

   death_penalty
   multi_objective_transformation
   self_adaptive_penalty
   co_evolution_penalty_method
   immune_system
   repair_methods

Evolutionary Optimization techniques have been successfully 
applied in a wide range of real-life optimization problems and 
a variety of them have been already integrated into the 
PaGMO/PyGMO collection of heuristic optimization methods. 
Several constraint handling approaches for evolutionary 
algorithms have been developed in the recent years and divided 
in five main classes (survey of techniques in [1]): Penalty 
Methods, Special Representations and Operators, Repair 
Algorithms, Separation of Objective and Constraints and Hybrid 
Methods. Among them few approaches have been identified as 
promising techniques and have been integrated into the existing 
PaGMO/PyGMO infrastructure. 

In particular:

* Four Penalty Methods: Death penalty (infeasible individual are rejected in the selection procedure regardless their level of infeasibility, [2]); Self Adaptive Fitness Formulation (the fitness of the infeasible individual is increased as to favor those solutions which are nearly feasible and also have a good objective function value, [3]); and Co-evolutionary Penalty (collaborative evolution of two populations, the second one encoding the penalty factors to be used in the fitness evaluation of the individuals of the first population, [4]).
* Separation of Objective and Constraints: Multi-objective Optimization (the constrained optimization problem is transformed into a multiobjective optimization problem where each constraint violation represents an additional objective [5]).
* Hybrid Methods: Immune System (part of the feasible population is transformed into antigen population, the infeasible solutions generate antibodies those evolving, with particular immune operators, until they become sufficiently similar to the antigen population and then recombined to continue with the evolution [6]); Local Repair Techniques (the infeasible solutions are locally repaired with gradient based unconstrained optimization techniques minimizing the sum of the constraints violations [7]).

All of these techniques have been implmented in PaGMO/PyGMO to 
effectively solve constrained optimization problems and to have a 
consistent framework to assess the performances of the different 
constraints handling techniques. In fact, from [1], it is important 
to be able to establish under what conditions certain constraints 
handling techniques are more efficient and more convenient than other 
ones (that might depends the number of linear/non-linear constraints, 
number of local optima, continuity of the objective function...). 
From the previous studies, only inconclusive evidences of techniques 
from benchmarks could be drawn. This is due to the fact that most of 
the studies were incomplete (i.e run only on one test case) or 
inconsistent, making it difficult to retrieve information for best 
practices in constrained stochastic optimization.

Each approach has a different implementation strategy:

* Death Penalty and Multi-objective Optimization redefine statically the optimization problem (i.e. the new optimization problem is the same during all the evolutions), two meta-problems have been created to cope with these techniques. The meta-problem takes the original constrained problem as input and generate a new unconstrained problem (respectively single and multi-objective).
* Self Adaptive redefines dynamically the optimization problem (information are extrapolated from the current population to assign the new fitness function), a meta-algorithm that internally invokes the evolution of the population on a modified problem has been implemented. The meta-algorithm is constructed given an heuristic algorithm and the constrained problem.
* Co-Evolutionary Penalty and Immune System are meta-algorithms as well. They are constructed as before with an heuristic algorithm and a constrained problem. Internally the method performs a preprocessing on the population and defines a modified problem suitable for the evolution.
* Local Repair Technique is a meta-algorithm that takes as argument two algorithms (an heuristic and a local one), it instantiates a new population (the infeasible solutions) with an unconstrained problem defined by the sum of the violation of the constraints and invokes the evolve method of the local technique to solve it.

References:

[1] Coello Coello C.A. Theoretical and Numerical Constraints-Handling Techniques 
Used with Evolutionary Algorithms: A Survey of the State of the Art, Comput. 
Methods Appl. Mech Engrg. 191, pp. 1245-1287, 2000.

[2] Bäck, T., Hoffmeister, F. and Schwell, H.P. A Survey of evolution strategies, 
Proceedings of the Fourth International Conference on Genetic Algorithms, Morgan 
Kaufmann, 2-9, 1991.

[3] Farmani, R. and Wright, J. Self-adaptive fitness formulation for constrained 
optimization, IEEE Transactions on Evolutionary Computation, 7 (5), pp. 445- 455, 2003.

[4] Coello Coello C.A. Use of Self-Adaptive Penalty Approach for Engineering 
Optimization Problems, Comput. Ind. 41(2), pp. 113-127, 2000.

[5] Coello Coello C.A. Treating Constraints as Objectives for Single-Objective 
Evolutionary Optimization, Engrg. Optim., 32(3), pp. 275-308, 2000.

[6] Hajela P. and Lee J. Constrained Genetic Search via Schema Adaptation. An 
Immune Network Solution, Proceedings of the First World Congress of Structural 
and Multidisciplinary Optimization, pp. 915-920, 1995.

[7] Belur S.V. CORE: Constrained Optimization by Random Evolution, Late Breaking 
Papers at the Genetic Programming Conference, pp. 280-286, 1997.

[8] J. J. Liang, T. P. Runarsson, E. Mezura-Montes, M. Clerc, P. N. Suganthan, C. A. 
Coello Coello, K. Deb, Problem Definitions and Evaluation Criteria for the 
CEC 2006, Special Session on Constrained Real-Parameter Optimization, Technical 
Report, Nanyang Technological University, Singapore, 2006.

