from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
import unittest


class ArchipelagoTests(unittest.TestCase):
    """Test archipelago properties"""

    def test_simple_probs(self):
        """ Testing whether migration history matches the expected output for each migration_direction and distribution_type """

        for migr_dir in [migration_direction.source, migration_direction.destination]:
            for dist_type in [distribution_type.point_to_point, distribution_type.broadcast]:
                prob = problem.rosenbrock(10)
                alg = algorithm.jde(20)
                archi = archipelago(alg, prob, 3, 20, migration_direction=migr_dir, distribution_type=dist_type)
                top = topology.ring(3)
                top.set_weight(0, 1, 0.0)
                top.set_weight(0, 2, 0.0)
                top.set_weight(1, 2, 0.0)
                top.set_weight(1, 0, 1.0)
                top.set_weight(2, 0, 1.0)
                top.set_weight(2, 1, 1.0)
                archi.topology = top
                archi.evolve(200)
                migr_hist = archi.dump_migr_history()
                # Below: After 200 evaluations, there should be some migrants from 1->0, 2->0 and 2->1
                # There should be no migrants from 1->0, 2->0 and 2->1
                self.assertTrue("(1,0,1)" not in migr_hist)
                self.assertTrue("(1,0,2)" not in migr_hist)
                self.assertTrue("(1,1,2)" not in migr_hist)

    def do_test_migr_setup(self, pop_xs, out_pop_xs, top, n_evolves):
        """ Generic procedure for testing whether the state of populations in 'pop_xs',
        after performing 'n_evolve' migration steps is equal to the expected 'out_pop_xs', given topology 'top'. """
        prob = problem.identity()
        alg = algorithm.null()
        pops = []
        for xs in pop_xs:
            pop = population(prob)
            for x in xs:
                pop.push_back(x)
            pops.append(pop)

        archi = archipelago(distribution_type=distribution_type.broadcast, migration_direction=migration_direction.destination)
        for pop in pops:
            archi.push_back(island(alg, pop, s_policy=migration.best_s_policy(), r_policy=migration.fair_r_policy()))
        archi.topology = top
        archi.evolve_batch(n_evolves, 1, False)
        out_xs = []
        for ii, isl in enumerate(archi, 1):
            out_xs.append(tuple(sorted([i.cur_f for i in isl.population])))
        out_xs = tuple(out_xs)
        self.assertEqual(out_xs, out_pop_xs)

    def test_one_way_ring_null_alg(self):
        pop_xs = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((2.0, ), (2.0, ), (2.0, )),
            ((3.0, ), (3.0, ), (3.0, )),
            ((4.0, ), (4.0, ), (4.0, )),
        )
        # Since we evolve the populations in index-order, the champion of pop1 ([1, ])
        # will travel along the ring until pop4
        out_pop_xs = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((1.0, ), (2.0, ), (2.0, )),
            ((1.0, ), (3.0, ), (3.0, )),
            ((1.0, ), (4.0, ), (4.0, )),
        )
        top = topology.one_way_ring(4)
        top.set_weight(1.0)
        self.do_test_migr_setup(pop_xs, out_pop_xs, top, 1)

        # We set the probability to 0.0 between islands 0 and 1. During migration,
        # it is the the champion of pop2 - (2.0, ) which travels along the ring
        top.set_weight(0, 1, 0.0)
        out_pop_xs_2 = (
            ((1.0, ), (2.0, ), (3.0, )),
            ((2.0, ), (2.0, ), (2.0, )),
            ((2.0, ), (3.0, ), (3.0, )),
            ((2.0, ), (4.0, ), (4.0, )),
        )
        self.do_test_migr_setup(pop_xs, out_pop_xs_2, top, 1)

        # We set the probability to 0.0, not migration should happen
        top.set_weight(0.0)
        self.do_test_migr_setup(pop_xs, pop_xs, top, 1)

    def test_distribution_type(self):
        """Testing whether the distribution_type property works"""

        # Ensure the default distribution type is point_to_point
        a = archipelago()
        self.assertEqual(a.distribution_type, distribution_type.point_to_point)

        for dist_type in [distribution_type.point_to_point, distribution_type.broadcast]:
            a = archipelago()
            a.distribution_type = dist_type

            self.assertEqual(a.distribution_type, dist_type)


def get_archipelago_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(ArchipelagoTests))
    return suite
