from PyGMO import archipelago, distribution_type
import unittest


class ArchipelagoTests(unittest.TestCase):
    """Test archipelago properties"""

    def test_distribution_type(self):
        """Testing whether the distribution_type property works"""

        # Ensure the default distribution type is point_to_point
        a = archipelago()
        self.assertEqual(a.distribution_type, distribution_type.point_to_point)

        for dist_type in [distribution_type.point_to_point, distribution_type.broadcast]:
            a = archipelago()
            a.distribution_type = dist_type

            self.assertEqual(a.distribution_type, dist_type)


def get_property_test_suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(ArchipelagoTests))
    return suite
