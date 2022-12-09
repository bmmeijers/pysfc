import unittest
from pysfc.ndgeom.ndbox import ndbox

class nDBoxTest(unittest.TestCase):
    """from nd box create hyperplanes
    """
    def test_1d_case(self):
        # 1d box
        box = ndbox([-1], [1])
        # creates two hyperplanes
        hyperplanes = box.as_hyperplanes()
        # check their relation to the origin, and point on the left and right
        origin = [0]
        negative = [-2]
        positive = [+2]
        #
        left = hyperplanes[0]
        self.assertEqual(left.signed_distance([-1]), 0.0)
        self.assertEqual(left.signed_distance(origin), -1)
        self.assertEqual(left.signed_distance(negative), 1)
        self.assertEqual(left.signed_distance(positive), -3)
        #
        right = hyperplanes[1]
        self.assertEqual(right.signed_distance([1]), 0.0)
        self.assertEqual(right.signed_distance(origin), -1)
        self.assertEqual(right.signed_distance(negative), -3)
        self.assertEqual(right.signed_distance(positive), 1)

    def test_nd_case(self):
        for d in range(1, 15):
            box = ndbox([-1] * d, [1] * d)
            hyperplanes = box.as_hyperplanes()
            origin = [0] * d
            for hp in hyperplanes:
                self.assertEqual(hp.signed_distance(origin), -1)


if __name__ == "__main__":
    unittest.main()