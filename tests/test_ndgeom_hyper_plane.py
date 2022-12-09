import unittest

from math import isclose
from pysfc.linalg import unit
from pysfc.ndgeom.hyper_plane import HyperPlane

class HyperPlaneTest(unittest.TestCase):
    def test_one(self):
        hp = HyperPlane.from_points([[0,0], [10,0]])
        self.assertEqual(hp.signed_distance(hp.through), 0)

    def test_two(self):
        pts = [[0, 10], [10, 0]]
        hp = HyperPlane.from_points(pts)
        for pt in pts:
            self.assertAlmostEqual(hp.signed_distance(pt), 0)
        self.assertAlmostEqual(hp.signed_distance(hp.through), 0)

    def test_three(self):
        # if we construct the plane, then we get a point through it
        pt = [-5, -5]
        hp = HyperPlane.from_normal_and_point(unit([1, -1]), pt)
        self.assertAlmostEqual(hp.signed_distance(pt), 0)
        self.assertAlmostEqual(hp.signed_distance(hp.through), 0)

if __name__ == "__main__":
    unittest.main()