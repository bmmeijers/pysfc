import unittest
from relate import ndbox
from newapi import hquery
from newapi import nquery
from query_brute import brute
from range_pack import filter_ranges_by_gapsize

def is_sorted(seq):
    "Returns whether a sequence is sorted"
    if len(seq) > 1:
        return all(a < b for a, b in zip(seq, seq[1:]))
    # we consider a sequence of 0 and 1 element to be always sorted
    return True


class TestSorted(unittest.TestCase):
    "Test for is_sorted (used inside tests below)"

    def test_sorted(self):
        "Increasing sequence is sorted"
        for i in range(10):
            seq = range(i)
            assert is_sorted(seq) == True

    def test_sorted_rev(self):
        "Decreasing sequence is not sorted"
        for i in range(2, 10):
            seq = range(i)
            seq.reverse()
            assert is_sorted(seq) == False


boxes2d = [
    # cube2x2; the 4 main boxes, 1x1
    ndbox([0, 0], [1, 1]),
    ndbox([0, 1], [1, 2]),
    ndbox([1, 0], [2, 1]),
    ndbox([1, 1], [2, 2]),

    # cube4x4; the 4 main boxes, 2x2
    ndbox([0, 0], [2, 2]),
    ndbox([2, 2], [4, 4]),
    ndbox([0, 2], [2, 4]),
    ndbox([2, 0], [4, 2]),

    # cube4x4; a center box, 2x2
    ndbox([1, 1], [3, 3]),

    # cube4x4; halves of the cube, 4x2 | 2x4
    ndbox([0, 0], [4, 2]),
    ndbox([0, 2], [4, 4]),
    ndbox([0, 0], [2, 4]),
    ndbox([2, 0], [4, 4]),

    # cube4x4; 2x3 | 3x2
    ndbox([0, 0], [2, 3]),
    ndbox([2, 0], [4, 3]),
    ndbox([0, 1], [2, 4]),
    ndbox([2, 1], [4, 4]),
    ndbox([0, 0], [3, 2]),
    ndbox([0, 2], [3, 4]),
    ndbox([1, 0], [4, 2]),
    ndbox([1, 2], [4, 4]),
    ndbox([0, 1], [3, 3]),
    ndbox([1, 0], [4, 3]),

    # cube4x4; arbitrary boxes
    ndbox([1, 1], [4, 4]),
    ndbox([0, 1], [2, 4]),
    ndbox([2, 0], [4, 3]),

    # cube8x8; center
    ndbox([3, 3], [5, 5]),
    # cube8x8; left
    ndbox([1, 1], [3, 7]),
    # cube8x8; right
    ndbox([5, 2], [7, 4]),
]


class TestAgainstBrute2D(unittest.TestCase):

    def setUp(self):
        self.boxes = boxes2d

    def test_brute(self):
        "Compare against brute-force computation"
        for tp, func in (('h', hquery), ('n', nquery)):
            for box in self.boxes:
                expected = brute(box, tp)
                obtained = func(box)
                # make the ranges list sorted
                expected.sort()
                obtained.sort()
                # glue adjacent + consecutive ranges
                expected = filter_ranges_by_gapsize(expected)
                obtained = filter_ranges_by_gapsize(obtained)
                assert expected == obtained, "{} {} {} {}".format(box, expected, obtained, func)


class TestQueryFramework2D(unittest.TestCase):

    def setUp(self):
        self.boxes = boxes2d

    def test_sorted(self):
        "Ranges sorted by their start"
        for query in hquery, nquery:
            for box in self.boxes:
                result = query(box)
                assert is_sorted([start for start, end in result]), "{} {} {}".format(result, query, box)

    def test_have_length(self):
        "Ranges have length (their end > start)"
        for query in hquery, nquery:
            for box in self.boxes:
                result = query(box)
                assert all([end > start for start, end in result]), "{} {} {}".format(result, query, box)


#    def test_given_queries(self):
#        "Some Hilbert queries with answers"
#        hilbert_queries = [
#            (ndbox([1, 1], [4, 4]), [(2, 3), (6, 7), (7, 8), (8, 12), (12, 13), (13, 14)]),
#            (ndbox([0, 1], [2, 4]), [(2, 3), (3, 4), (4, 8)]),
#            (ndbox([2, 0], [4, 3]), [(8, 9), (11, 12), (12, 16)]),
#        ]

#        for box, expected in hilbert_queries:
#            obtained = query(box)
#            assert expected == obtained, "exp {} != obt {}".format(expected, obtained)


if __name__ == "__main__":
    unittest.main()
