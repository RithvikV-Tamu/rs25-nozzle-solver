import unittest

from rs25_nozzle_solver import area_mach_ratio, solve_supersonic_mach, solve_nozzle


class AreaMachTests(unittest.TestCase):
    def test_area_ratio_is_one_at_sonic_point(self):
        for gamma in (1.2, 1.4):
            self.assertAlmostEqual(area_mach_ratio(1.0, gamma), 1.0, places=12)

    def test_known_value_for_air(self):
        # For gamma = 1.4 and M = 2: 0.5 * (1.5)**3 = 1.6875 exactly.
        self.assertAlmostEqual(area_mach_ratio(2.0, 1.4), 1.6875, places=12)

    def test_root_finder_recovers_mach_two_for_air(self):
        self.assertAlmostEqual(solve_supersonic_mach(1.6875, 1.4), 2.0, places=8)

    def test_root_finder_round_trip(self):
        for gamma in (1.2, 1.3, 1.4):
            for mach in (1.05, 1.5, 2.5, 4.0, 4.6, 6.0):
                area_ratio = area_mach_ratio(mach, gamma)
                self.assertAlmostEqual(solve_supersonic_mach(area_ratio, gamma), mach, places=8)

    def test_throat_returns_mach_one(self):
        self.assertEqual(solve_supersonic_mach(1.0, 1.2), 1.0)

    def test_area_ratio_below_one_is_rejected(self):
        with self.assertRaises(ValueError):
            solve_supersonic_mach(0.9, 1.2)

    def test_nozzle_mach_increases_from_throat_to_exit(self):
        result = solve_nozzle()
        self.assertEqual(result["mach"][0], 1.0)
        self.assertTrue((result["mach"][1:] > 1.0).all())
        self.assertTrue((result["mach"][1:] >= result["mach"][:-1]).all())


if __name__ == "__main__":
    unittest.main()
