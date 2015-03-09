from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
import unittest


class ArchipelagoTests(unittest.TestCase):
    """ Test archipelago migration with different migration probabilities """

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
                # Although this test does not seem to fail, it is not deterministic
                self.assertTrue("(1,1,0)" in migr_hist)
                self.assertTrue("(1,2,0)" in migr_hist)
                self.assertTrue("(1,2,1)" in migr_hist)
                # Below: There should be no migrants from 1->0, 2->0 and 2->1
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


class TopologyTests(unittest.TestCase):
    """ Tests for topology object and migration probabability """

    def test_overloaded_add_edge(self):
        """ Tests the add_edge(n, m, weight) setter """

        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1, 0.1)
        t.add_edge(2, 0)  # 1.0 by default
        t.add_edge(0, 2, 0.3)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 0.3)

    def test_incorrect_probability(self):
        """ Tests the setter for probabilities outside [0, 1] """

        top = topology.custom()
        top.push_back()
        top.push_back()

        # Add edges with valid probability
        top.add_edge(0, 1, 0.1)
        top.add_edge(1, 0, 0.1)
        top.set_weight(0.1)
        top.set_weight(0, 0.1)
        top.set_weight(0, 1, 0.1)

        self.assertRaises(ValueError, top.set_weight, -0.1)
        self.assertRaises(ValueError, top.set_weight, 0, -0.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1, -0.1)

        self.assertRaises(ValueError, top.set_weight, 1.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1.1)
        self.assertRaises(ValueError, top.set_weight, 0, 1, 1.1)

        # Add new vertex and try to add_edge with incorrect probability
        top.push_back()
        self.assertRaises(ValueError, top.add_edge, 1, 2, -0.1)
        self.assertRaises(ValueError, top.add_edge, 2, 1, 1.1)

    def test_basic_probs(self):
        """ Tests the basic getter/setter for custom migration probability """
        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1)
        t.add_edge(2, 0)
        t.add_edge(0, 2)
        t.set_weight(2, 0, 0.1)
        t.set_weight(0, 2, 0.2)
        self.assertEqual(t.get_weight(0, 1), 1.0)  # This one should be set by default
        self.assertEqual(t.get_weight(2, 0), 0.1)
        self.assertEqual(t.get_weight(0, 2), 0.2)

    def test_ring_topo(self):
        """ Tests the getter/setter for the ring topology """

        t = topology.ring()
        t.push_back()
        t.push_back()

        # By default the weights are 1.0
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 1.0)

        # Once you alter one weight, the remaining are still 1.0
        t.set_weight(0, 1, 0.1)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 1.0)

        # Pushing new node should not change anything to the previous ones
        t.push_back()
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 1), 1.0)

        # Test the all-custom weights for ring topology
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.set_weight(0, 2, 0.3)
        t.set_weight(2, 0, 0.4)
        t.set_weight(1, 2, 0.5)
        t.set_weight(2, 1, 0.6)
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 0.2)
        self.assertEqual(t.get_weight(0, 2), 0.3)
        self.assertEqual(t.get_weight(2, 0), 0.4)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        self.assertEqual(t.get_weight(2, 1), 0.6)

    def test_mixed_migr_representation(self):
        """ Testing the mixing of different edge weight setters """

        t = topology.custom()
        t.push_back()
        t.push_back()
        t.push_back()
        t.push_back()
        t.add_edge(0, 1)
        t.add_edge(1, 0)
        t.add_edge(1, 2)
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 1.0)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        t.set_weight(1, 0.5)  # Set weight for each out-going edge of vertex 1
        self.assertEqual(t.get_weight(0, 1), 1.0)
        self.assertEqual(t.get_weight(1, 0), 0.5)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        t.add_edge(1, 3)  # Weight of a new out-going edge of vertex 1 should still be 1.0 by default
        self.assertEqual(t.get_weight(1, 0), 0.5)
        self.assertEqual(t.get_weight(1, 2), 0.5)
        self.assertEqual(t.get_weight(1, 3), 1.0)

        t.set_weight(0.3)  # Set the weight to all edges
        self.assertEqual(t.get_weight(0, 1), 0.3)
        self.assertEqual(t.get_weight(1, 0), 0.3)
        self.assertEqual(t.get_weight(1, 2), 0.3)
        self.assertEqual(t.get_weight(1, 3), 0.3)
        t.add_edge(2, 3)  # Weight of a new edge should still be 1.0 by default

        self.assertEqual(t.get_weight(0, 1), 0.3)
        self.assertEqual(t.get_weight(1, 0), 0.3)
        self.assertEqual(t.get_weight(1, 2), 0.3)
        self.assertEqual(t.get_weight(1, 3), 0.3)
        self.assertEqual(t.get_weight(2, 3), 1.0)

    def test_chain_break(self):
        """ Testing 'chain breaking' after push_back in one-way-ring """

        t = topology.ring(2)
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.push_back()
        # In two-way ring, the weights between the first and second island should be the same after push_back
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 0), 0.2)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 1), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)
        self.assertEqual(t.get_weight(0, 2), 1.0)

        # In one-way-ring, the edge 1->0 will be broken
        # and replaced by 1->2 and 2->0 with weights 1.0 each
        t = topology.one_way_ring(2)
        t.set_weight(0, 1, 0.1)
        t.set_weight(1, 0, 0.2)
        t.push_back()
        self.assertEqual(t.get_weight(0, 1), 0.1)
        self.assertEqual(t.get_weight(1, 2), 1.0)
        self.assertEqual(t.get_weight(2, 0), 1.0)

    def test_topology_serialize(self):
        """ Testing whether the weights are retained after serialization of an archipelago """

        prob = problem.rosenbrock(10)
        alg = algorithm.jde(20)
        archi = archipelago(alg, prob, 4, 20)
        top = topology.ring(4)
        top.set_weight(0, 1, 0.01)
        top.set_weight(0, 3, 0.03)
        top.set_weight(1, 2, 0.12)
        top.set_weight(2, 1, 0.21)
        archi.topology = top
        archi.evolve(5)
        import pickle
        pickle.loads(pickle.dumps(archi))
        self.assertEqual(archi.topology.get_weight(0, 1), 0.01)
        self.assertEqual(archi.topology.get_weight(0, 3), 0.03)
        self.assertEqual(archi.topology.get_weight(1, 2), 0.12)
        self.assertEqual(archi.topology.get_weight(2, 1), 0.21)

    def test_incorrect_weights(self):
        """ Testing ValueError for getter/setter """

        top = topology.ring(4)
        self.assertRaises(ValueError, top.get_weight, 0, 2)
        self.assertRaises(ValueError, top.set_weight, 0, 2, 0.5)

    def test_average_path_length(self):
        """ Testing the average path length algorithm """

        top = topology.ring(2)
        self.assertEqual(top.get_average_shortest_path_length(), 1.0)
        top = topology.ring(3)
        self.assertEqual(top.get_average_shortest_path_length(), 1.0)
        top = topology.ring(5)
        self.assertEqual(top.get_average_shortest_path_length(), 1.5)
        top = topology.ring(9)
        self.assertEqual(top.get_average_shortest_path_length(), 2.5)
        top = topology.custom()
        top.push_back()
        top.push_back()
        top.push_back()
        top.push_back()
        top.add_edge(0, 1)
        top.add_edge(1, 2)
        top.add_edge(2, 3)
        top.add_edge(3, 0)
        top.add_edge(0, 2)
        self.assertEqual(top.get_average_shortest_path_length(), 1.75)


def get_topology_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(TopologyTests))
    suite.addTests(unittest.makeSuite(ArchipelagoTests))
    return suite
