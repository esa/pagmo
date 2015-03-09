from PyGMO import topology, algorithm, problem, archipelago, distribution_type, migration_direction, island, population, migration
import unittest



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
    return suite
