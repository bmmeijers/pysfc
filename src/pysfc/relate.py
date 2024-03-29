"""
Determine the spatial relation between two nD-boxes

Based on http://sfclib.github.io
"""
# from itertools import product
from pysfc.ndgeom import ndbox, bitreverse


def relate(rect, qrt):
    """
    Spatial relationship between two nd-boxes.

    Outcome can be:
        0: equal
        1: contains
        2: intersects
        -1: no overlap
    """
    # how many dimensions to check?
    dims = rect.dims

    # equal, all coordinates are equal
    ncmp = 1
    for d in range(dims):
        ncmp &= rect.lo[d] == qrt.lo[d] and rect.hi[d] == qrt.hi[d]
    if ncmp:
        return 0

    # fully contains, rect fully contains qrt
    ncmp = 1
    for d in range(dims):
        ncmp &= rect.lo[d] <= qrt.lo[d] and rect.hi[d] >= qrt.hi[d]
    if ncmp:
        return 1

    # intersects, the two nd-boxes interact 
    # (either on the boundary or internally)
    ncmp = 1
    for d in range(dims):
        ncmp &= rect.lo[d] < qrt.hi[d] and rect.hi[d] > qrt.lo[d]
    if ncmp:
        return 2

    # no overlap
    return -1


def _test():
    # test that equal boxes do result in equal (0)
    one = ndbox((0, 0, 0), (10, 10, 10))
    other = ndbox((0, 0, 0), (10, 10, 10))
    assert relate(one, one) == 0
    assert relate(one, other) == 0

    # test that non-equal boxes do not result in equal (not 0)
    one = ndbox((0, 0, 0), (10, 10, 10))
    other = ndbox((0, 0, 0), (1, 10, 10))
    assert relate(one, other) != 0
    other = ndbox((0, 0, 0), (10, 1, 10))
    assert relate(one, other) != 0
    other = ndbox((0, 0, 0), (10, 10, 1))
    assert relate(one, other) != 0
    other = ndbox((1, 0, 0), (10, 10, 10))
    assert relate(one, other) != 0
    other = ndbox((0, 1, 0), (10, 10, 10))
    assert relate(one, other) != 0
    other = ndbox((0, 0, 1), (10, 10, 10))
    assert relate(one, other) != 0

    # test that this box has 8 corners
    corners = list(pt for pt in one.vertices)
    assert len(corners) == 8
    # the following vertices (not present in the box definition)
    # should be there in the generated vertices
    assert (0, 0, 10) in corners
    assert (0, 10, 0) in corners
    assert (10, 0, 0) in corners
    # these vertices are there in the constructor
    assert (0, 0, 0) in corners
    assert (10, 10, 10) in corners


if __name__ == "__main__":
    _test()
