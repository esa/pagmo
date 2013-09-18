.. _migration_based_on_hypervolumes:

================================================================
Migration based on hypervolumes
================================================================

In this tutorial we will cover some migration strategies that are based on the hypervolume computation.
There are in total 4 migration policies which are based on the hypervolume feature.

#. hv_greedy_s_policy
#. hv_greedy_r_policy
#. hv_best_s_policy
#. hv_fair_r_policy

We can divide these into two groups, by the strategy with which we apply the ordering to given set of individuals:

#. Successive contribution
#. Current contribution

The distinction between these originates from the following observation:

Given a set of points for which we compute the exclusive contribution to the hypervolume, we can order them descending.
This suggests that an individual with small contribution are less "interesting", and in general less favored by the evolutionary strategy.
It turns out that this reasoning does not necessarily have to be true, e.g. two individuals located not far apart will appear to be weak, only due to the fact that each of these contributes (exclusively) very small amount of volume. On the other hand, any individual that may be relatively weak (according to some golden-standard of "best" individual), will turn out to be ranked higher.
This may lead to worsening the state of the population.

One solution is computing the least (or greatest, depending on the type of policy) contributor successively, e.g. compute the least contributor (store it for later) and remove it from the set of points.
In such case, any sibling that was nearby the previous individual (point) is now valued properly.


Successive contribution
===============================================================

Selection policy (hv_greedy_s_policy)
-------------------------------------

Replacement policy (hv_greedy_r_policy)
---------------------------------------

Single contribution
=======================================================

Selection policy (hv_best_s_policy)
-------------------------------------

Replacement policy (hv_fair_r_policy)
-------------------------------------
