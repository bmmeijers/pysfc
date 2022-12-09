import unittest
from pysfc.encode_decode import nenc, ndec, henc, hdec


# 
# def do_at_depth(depth, n_dims):
# the morton curve makes jumps, where the amount of bits needed
# to represent the curve so far is not sufficient enough any more
#     x = (2 ** n_dims) ** depth
#     around = 20  # number of values before and after the 'jump'
#     for init_val in range(max(0, x - around), x + around):
#         # sfc.append((val, ndec(val, 2)))
#         coord = ndec(init_val, n_dims)
#         _ = nenc(coord=coord)
#         print()
#
# do_at_depth(25, 2)
#

class TestBug(unittest.TestCase):
    def test_1_nsmall(self):
        "First 9 steps of the N-order curve in 2D"
        # test first 9 steps of the n-order curve in 2D
        tests = [
            ((0, 0), 0),
            ((0, 1), 1),
            ((1, 0), 2),
            ((1, 1), 3), 
            ((0, 2), 4),
            ((0, 3), 5),
            ((1, 2), 6),
            ((1, 3), 7),
            ((2, 0), 8)]
        for c, outcome in tests:
            self.assertEqual(nenc(c), outcome)
            self.assertEqual(c, ndec(outcome, 2))

    def test_1_hsmall(self):
        "First 9 steps of the H-order curve in 2D"
        # test first 9 steps of the h-order curve in 2D
        tests = [
            ((0, 0), 0),
            ((1, 0), 1),
            ((1, 1), 2),
            ((0, 1), 3), 
            ((0, 2), 4),
            ((0, 3), 5),
            ((1, 3), 6),
            ((1, 2), 7),
            ((2, 2), 8)]
        for c, outcome in tests:
            self.assertEqual(henc(c), outcome)
            self.assertEqual(c, hdec(outcome, 2))

    def test_2_bug_dims(self):
        correct_sfc_keys = [
            1125899906842623,
            1125899906842624,
            1125899906842625,
            1125899906842626,
            1125899906842627,
            1125899906842628,
        ]
        coords = [(33554431, 33554431), (0, 33554432), (0, 33554433), (1, 33554432), (1, 33554433), (0, 33554434)]
        #for enc, dec, in [(nenc, ndec), (henc, hdec)]:
        for coord, correct_sfc_key in zip(coords, correct_sfc_keys):
            # assert nenc(coord) == correct_sfc_key
            self.assertEqual(nenc(coord), correct_sfc_key)
            # assert ndec(correct_sfc_key, 2) == coord
            self.assertEqual(ndec(correct_sfc_key, 2), coord)

    def _do_at_depth(self, encode, decode, depth, n_dims):
        # the morton sfc curve makes jumps, where the amount of bits needed
        # to represent the curve so far is not sufficient enough any more
        x = (2 ** n_dims) ** depth
        around = 5  # number of keys before and after the 'jump'
        for init_val in range(max(0, x - around), x + around):
            coord = decode(init_val, n_dims)
            key = encode(coord=coord)
            self.assertEqual(key, init_val)

    def test_jumps(self):
        """The Morton SFC makes larger jumps, where the amount of bits needed
        to represent the curve so far is not sufficient enough any more.

        This test test both the Morton and Hilbert code to see if this `jump' in
        the needed number of bits is handled correctly (in earlier code there was
        a bug).
        """
        for n_dims in range(1, 10):
            for level in range(64):
                self._do_at_depth(nenc, ndec, level, n_dims)
                self._do_at_depth(henc, hdec, level, n_dims)


if __name__ == "__main__":
    unittest.main()
