import unittest
from pysfc.ndgeom.hyper_plane import _visualize_hyperplane_2d
from pysfc.ndgeom.ndbox import ndbox
from pysfc.query.nary_tree import hquery, nquery, range_generation
from pysfc.query.range_glue import glue_ranges_by_gapsize
from pysfc.encode_decode import henc as hencode, nenc as nencode
from pysfc.ndgeom.ndbox import ndbox

def make_boxes(mbits=4, ndims=3):
    # 2**M = max side of the cube
    assert mbits > 2
    assert ndims > 2
    boxes = []

    # generate a set of nd-boxes for querying
    # the boxes are positioned around the center of the cube
    # taking a full slice of 1 high in 1 dimension,
    # while using all other dimensions to their full extents

    for n in range(2, ndims):
        for m in range(2, mbits):
            val = 2**m
            half = (val - 1) // 2
            full_range = (0, val)
            below_half_range = (half, half + 1)
            over_half_range = (half + 1, half + 2)
            for where in below_half_range, over_half_range:
                for m in range(0, n):
                    lo, hi = [], []
                    for _ in range(n):
                        if _ != m:
                            r = full_range
                        else:
                            r = where
                        lo.append(r[0])
                        hi.append(r[1])
                    box = ndbox(lo, hi)
                    boxes.append(box)
    return boxes


from itertools import product  # for cartesian product


def compute_deltas(query):
    """Returns how large every dimension is
    """
    dims = query.dims
    diff = [None] * dims
    for d in range(dims):
        diff[d] = int(query.hi[d] - query.lo[d])
    return diff


def brute(query, tp = 'h'):
    """A way of specifying ranges to search with in a brute force way

    Note, this function is only intended to test whether a more elegant
    way of querying performs in a correct way.
    """
    # the type of the SFC (h: hilbert, n: n-order, z: z-order)
    if tp == 'h':
        enc = hencode
    elif tp == 'n':
        enc = nencode
    else:
        raise ValueError("Unknown encoding given")

    dims = query.dims
    diff = compute_deltas(query)
    # materialize the values that are there along each dimension
    along_dims = []
    for d in range(dims):
        start_of_dim = query.lo[d]
        along_d = []
        for val in range(diff[d]):
            along_d.append(start_of_dim + val)
        along_dims.append(along_d)

    # we obtain the Cartesian product of all the dimensions
    # and for each value we map the coordinate to the sfc key value
    ranges =[]
    for c in product(*along_dims):
        start = enc(c)
        ranges.append((start, start+1))

    # now sort the ranges we obtained (necessary for hilbert)
    ranges.sort()

    return ranges


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
            seq = list(range(i))
            assert is_sorted(seq) == True

    def test_sorted_rev(self):
        "Decreasing sequence is not sorted"
        for i in range(2, 10):
            seq = list(range(i))
            seq.reverse()
            assert is_sorted(seq) == False


boxes2d = [
    # cube2d 2x2; the 4 main boxes, 1x1
    ndbox([0, 0], [1, 1]),
    ndbox([0, 1], [1, 2]),
    ndbox([1, 0], [2, 1]),
    ndbox([1, 1], [2, 2]),

    # cube2d 4x4; the 4 main boxes, 2x2
    ndbox([0, 0], [2, 2]),
    ndbox([2, 2], [4, 4]),
    ndbox([0, 2], [2, 4]),
    ndbox([2, 0], [4, 2]),

    # cube2d 4x4; a center box, 2x2
    ndbox([1, 1], [3, 3]),

    # cube2d 4x4; halves of the cube, 4x2 | 2x4
    ndbox([0, 0], [4, 2]),
    ndbox([0, 2], [4, 4]),
    ndbox([0, 0], [2, 4]),
    ndbox([2, 0], [4, 4]),

    # cube2d 4x4; 2x3 | 3x2
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

    # cube2d 4x4; arbitrary boxes
    ndbox([1, 1], [4, 4]),
    ndbox([0, 1], [2, 4]),
    ndbox([2, 0], [4, 3]),

    # cube2d 8x8; center
    ndbox([3, 3], [5, 5]),
    # cube2d 8x8; left
    ndbox([1, 1], [3, 7]),
    # cube2d 8x8; right
    ndbox([5, 2], [7, 4]),
]

boxes3d = [
    # cube3d; the 8 main boxes, 1x1x1
    ndbox([0, 0, 0], [1, 1, 1]),
    ndbox([0, 0, 1], [1, 1, 2]),
    ndbox([0, 1, 0], [1, 2, 1]),
    ndbox([0, 1, 1], [1, 2, 2]),
    ndbox([1, 0, 0], [2, 1, 1]),
    ndbox([1, 0, 1], [2, 1, 2]),
    ndbox([1, 1, 0], [2, 2, 1]),
    ndbox([1, 1, 1], [2, 2, 2]),

    # cube3d; boxes of 2x2x2
    ndbox([0, 0, 0], [2, 2, 2]),
    ndbox([0, 0, 2], [2, 2, 4]),
    ndbox([2, 2, 0], [4, 4, 2]),
    ndbox([2, 2, 2], [4, 4, 4]),
    ndbox([0, 2, 0], [2, 4, 2]),
    ndbox([0, 2, 2], [2, 4, 4]),
    ndbox([2, 0, 0], [4, 2, 2]),
    ndbox([2, 0, 2], [4, 2, 4]),

    # cube3d; center box of 2x2x2
    ndbox([1, 1, 1], [3, 3, 3]),

    # cube3d; horizontal 
    ndbox([0, 0, 1], [4, 4, 2]),
    ndbox([0, 0, 2], [4, 4, 3]),
]

boxes4d = [
    # cube4d; the 16 small boxes, 1x1x1x1
    ndbox([0, 0, 0, 0], [1, 1, 1, 1]),
    ndbox([0, 0, 1, 1], [1, 1, 2, 2]),
    ndbox([0, 0, 0, 0], [1, 1, 1, 1]),
    ndbox([0, 0, 1, 1], [1, 1, 2, 2]),
    ndbox([0, 1, 0, 0], [1, 2, 1, 1]),
    ndbox([0, 1, 1, 1], [1, 2, 2, 2]),
    ndbox([0, 1, 0, 0], [1, 2, 1, 1]),
    ndbox([0, 1, 1, 1], [1, 2, 2, 2]),
    ndbox([1, 0, 0, 0], [2, 1, 1, 1]),
    ndbox([1, 0, 1, 1], [2, 1, 2, 2]),
    ndbox([1, 0, 0, 0], [2, 1, 1, 1]),
    ndbox([1, 0, 1, 1], [2, 1, 2, 2]),
    ndbox([1, 1, 0, 0], [2, 2, 1, 1]),
    ndbox([1, 1, 1, 1], [2, 2, 2, 2]),
    ndbox([1, 1, 0, 0], [2, 2, 1, 1]),
    ndbox([1, 1, 1, 1], [2, 2, 2, 2]),

    # cube4d; boxes of 2x2x2
    ndbox([0, 0, 0, 0], [2, 2, 2, 2]),
    ndbox([0, 0, 2, 2], [2, 2, 4, 4]),
    ndbox([0, 0, 0, 0], [2, 2, 2, 2]),
    ndbox([0, 0, 2, 2], [2, 2, 4, 4]),
    ndbox([2, 2, 0, 0], [4, 4, 2, 2]),
    ndbox([2, 2, 2, 2], [4, 4, 4, 4]),
    ndbox([2, 2, 0, 0], [4, 4, 2, 2]),
    ndbox([2, 2, 2, 2], [4, 4, 4, 4]),
    ndbox([0, 2, 0, 0], [2, 4, 2, 2]),
    ndbox([0, 2, 2, 2], [2, 4, 4, 4]),
    ndbox([0, 2, 0, 0], [2, 4, 2, 2]),
    ndbox([0, 2, 2, 2], [2, 4, 4, 4]),
    ndbox([2, 0, 0, 0], [4, 2, 2, 2]),
    ndbox([2, 0, 0, 2], [4, 2, 2, 4]),
    ndbox([2, 0, 2, 0], [4, 2, 4, 2]),
    ndbox([2, 0, 2, 2], [4, 2, 4, 4]),

    # cube4d; center box of 2x2x2x2
    ndbox([1, 1, 1, 1], [3, 3, 3, 3]),
]

boxes5d = [
    # cube5d; center box of 2x2x2x2x2
    ndbox([1, 1, 1, 1, 1], [3, 3, 3, 3, 3]),
]

boxes6d = [
    # cube6d; center box of 2x2x2x2x2x2
    ndbox([1, 1, 1, 1, 1, 1], [3, 3, 3, 3, 3, 3]),
]

boxes7d = [
    # cube7d; center box of 2x2x2x2x2x2
    ndbox([1, 1, 1, 1, 1, 1, 1], [3, 3, 3, 3, 3, 3, 3]),
]

boxesnd = make_boxes(mbits=4, ndims=6)


class TestQueryFrameworkND(unittest.TestCase):
    "Testing the query functionality"

    def setUp(self):
        self.boxes = boxes2d + boxes3d + boxes4d + boxes5d + boxes6d + boxes7d + boxesnd

    def test_against_brute_as_hyperplanes(self):
        "Ranges same as obtained by brute-force computation"
        for tp, func in (('h', hquery), ('n', nquery), ):
            for box in self.boxes:
                expected = brute(box, tp)
                # expected = []
                obtained = func(box.as_hyperplanes())
                # make the ranges list sorted
                expected.sort()
                # the list from our query procedure should be sorted
                # obtained.sort()
                # glue adjacent + consecutive ranges
                expected = glue_ranges_by_gapsize(expected)
                obtained = glue_ranges_by_gapsize(obtained)
                # assert expected == obtained, "\n{}\n expected ranges: {}\n obtained ranges: {}\n function: {}".format(box, expected, obtained, func)
                # after glue'ing the ranges should be equivalent
                self.assertEqual(expected, obtained)

    def test_against_brute_as_box(self):
        "Ranges same as obtained by brute-force computation"
        for tp, func in (('h', hquery), ('n', nquery), ):
            for box in self.boxes:
                expected = brute(box, tp)
                # expected = []
                obtained = func(box)
                # make the ranges list sorted
                expected.sort()
                # the list from our query procedure should be sorted
                # obtained.sort()
                # glue adjacent + consecutive ranges
                expected = glue_ranges_by_gapsize(expected)
                obtained = glue_ranges_by_gapsize(obtained)
                # assert expected == obtained, "\n{}\n expected ranges: {}\n obtained ranges: {}\n function: {}".format(box, expected, obtained, func)
                # after glue'ing the ranges should be equivalent
                self.assertEqual(expected, obtained)

    def test_sorted(self):
        "Ranges are sorted by their start"
        for query in hquery, nquery:
            for box in self.boxes:
                result = query(box)
                assert is_sorted([start for start, end in result]), "{} {} {}".format(result, query, box)

    def test_have_length(self):
        "Ranges do have length (their end > start)"
        for query in hquery, nquery:
            for box in self.boxes:
                result = query(box)
                assert all([end > start for start, end in result]), "{} {} {}".format(result, query, box)

# def test_manual():
#     boxes = boxes2d
#     tp = 'h'
#     for box in boxes:
#         expected = brute(box, tp)
#         obtained = hquery(box)# .as_hyperplanes())
#         expected.sort()
#         expected = glue_ranges_by_gapsize(expected)
#         obtained = glue_ranges_by_gapsize(obtained)
#         print(expected)
#         print(obtained)
#         assert expected == obtained

# def test_manual2():
#     boxes = boxes2d
#     is_hilbert = True
#     for box in boxes:
#         query_planes = box.as_hyperplanes()
#         ndims = box.dims
#         mbits = 25
#         max_level = mbits
#         result = range_generation(query_planes, ndims, mbits, max_level, is_hilbert)

#         with open("/tmp/query_planes_2d.csv", "w") as fh:
#             print("geometry;what", file=fh)
#             for hs in query_planes:
#                 _visualize_hyperplane_2d(hs, fh, size=2 ** mbits)

#         with open("/tmp/query_result.csv", "w") as fh:
#             print("interaction;geometry;sfc_start;sfc_end;mask", file=fh)
#             for node, interaction, mask in reversed(result):
#                 if is_hilbert:
#                     range_gen = node.as_hilbert_range
#                 else:
#                     range_gen = node.as_morton_range
#                 line = "{0};{1};{2[0]};{2[1]};{3}".format(
#                     interaction, node.as_ndbox().as_wkt_2d(), range_gen(), mask
#                 )
#                 print(line, file=fh)
#         input('paused')

if __name__ == "__main__":
    # test_manual2()
    unittest.main()
    
