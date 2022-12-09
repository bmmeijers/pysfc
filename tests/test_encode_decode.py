import unittest

from pysfc.encode_decode import henc, hdec
from pysfc.encode_decode import nenc, ndec

# maxside = 2**2
maxside = 2**3
# maxside = 2**4
# maxside = 2**5


class TestEncodingDecoding(unittest.TestCase):

    def test_hilbert_lengthy(self):
        "H-order enc and dec functions have to agree"
        self._do_lengthy_upto_including_5d(maxside, henc, hdec)

    def test_norder_lengthy(self):
        "N-order enc and dec functions have to agree"
        self._do_lengthy_upto_including_5d(maxside, nenc, ndec)

    def _do_lengthy_upto_including_5d(self, maxside, encode, decode):
        """Test that encode and decode functions agree with each other,
        also in higher dims"""
        # 2D
        for x in range(maxside):
            for y in range(maxside):
                c = (x, y)
                e = encode(c)
                d = decode(e, 2)
                self.assertEqual(c, d)
        # 3D
        for x in range(maxside):
            for y in range(maxside):
                for z in range(maxside):
                    c = (x, y, z)
                    e = encode(c)
                    d = decode(e, 3)
                    self.assertEqual(c, d)
        # 4D
        for x in range(maxside):
            for y in range(maxside):
                for z in range(maxside):
                    for t in range(maxside):
                        c = (x, y, z, t)
                        e = encode(c)
                        d = decode(e, 4)
                        self.assertEqual(c, d)
        # 5D
        for x in range(maxside):
            for y in range(maxside):
                for z in range(maxside):
                    for t in range(maxside):
                        for s in range(maxside):
                            c = (x, y, z, t, s)
                            e = encode(c)
                            d = decode(e, 5)
                            self.assertEqual(c, d)


    def test_nsmall(self):
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
            assert nenc(c) == outcome
            assert c == ndec(outcome, 2)

    def test_hsmall(self):
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
            assert henc(c) == outcome
            assert c == hdec(outcome, 2)


    def test_3d(self):
        "One 3D coordinate forth and back (H- and N-order)"
        assert ndec(1903, 3) == (7, 9, 15)
        assert nenc((7, 9, 15)) == 1903
        assert henc((7,9,15)) == 2049
        assert hdec(2049, 3) == (7, 9, 15)


    def test_4d(self):
        key = 72057594037927936
        n_dims = 4
        assert ndec(key, n_dims) != (0,0,0,0)

    def test_manual(self):
        # key = 72057594037927936-1
        key = 72057594037927936
        key = 72057594037927936+2
        n_dims = 4
        assert ndec(key, n_dims) != (0,0,0,0)

#    _test_nsmall()
#    _test_hsmall()
##    _test_lengthy(2**4, nenc, ndec)
#    _test_lengthy(2**4, henc, hdec)

if __name__ == "__main__":
    #test_manual()
    unittest.main()

