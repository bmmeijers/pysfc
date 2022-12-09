import unittest
import math
from pysfc.linalg import determinant
from pysfc.linalg import add, sub, mul, div, angle

class TestDeterminant(unittest.TestCase):
    # def test_2d(self):
    #     self.assertEqual(determinant([[1000., -1.], [1000., +1.]]), 0.0)
    # def test_2d_other(self):
    #     self.assertEqual(determinant([[1000., +1.], [1000., -1.]]), 0.0)

    def test_orient2d(self):
        """consistent with orient2d of predicates"""
        # assert orient2d( (0, 0), (0, 10), (-10, 10)) == 100.
        # assert orient2d( (0, 0), (0, 10), (10, 10)) == -100.

        # orientation test in 2d , raised to 3d:
        self.assertEqual(determinant([[0, 0, 1], [0, 10, 1], [-10, 10, 1]]), 100.0)
        self.assertEqual(determinant([[0, 0, 1], [0, 10, 1], [10, 10, 1]]), -100.0)

    def test_3d(self):
        # 3D test case
        self.assertEqual(determinant([[1, 3, 5], [2, 0, 4], [4, 2, 7]]), 18.0)
        self.assertEqual(
            determinant([[5, 3, 58], [-4, 23, 11], [34, 2, -67]]), -53317.0
        )

    def test_4d(self):
        # 4D test case
        self.assertEqual(
            determinant([[2, 1, 3, 0], [1, 0, 2, 3], [3, 2, 0, 1], [2, 0, 1, 3]]),
            -24,
        )


class TestVectorOps(unittest.TestCase):
    def test(self):
        v1, v2 = list(map(float, list(range(0, 3)))), list(map(float, list(range(4, 7))))
        # v1 = [0,1,2]
        # v2 = [4,5,6]

        assert sub(v2, v1) == (4, 4, 4), sub(v2, v1)
        assert sub(v1, v2) == (-4, -4, -4), sub(v1, v2)
        assert sub(v1, 5) == (-5, -4, -3), sub(v1, 5)

        assert add(v1, v2) == (4, 6, 8)
        assert add(v1, 10) == (10, 11, 12)

        assert div((4, 4, 4), 4) == (1, 1, 1)
        assert div(v1, v2) == (0, 1. / 5., 2. / 6.), div(v1, v2)

        assert mul(v1, v2) == (0, 1 * 5, 2 * 6), mul(v1, v2)
        assert mul(v2, 10) == (40, 50, 60)

        assert angle((1, 1), (-1, 1)) == math.pi * 0.5


if __name__ == "__main__":
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=5))