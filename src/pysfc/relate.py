"""
Determine the spatial relation between two nD-boxes

Based on http://sfclib.github.io
"""

from bitreverse import bitreverse


class ndbox(object):
    """A nD box object"""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi):
        """lo and hi attributes are lists with equal number of coordinates"""
        assert len(lo) == len(hi)
        self.lo = tuple(lo)
        self.hi = tuple(hi)

    @property
    def dims(self):
        """How many dimensions does this box have (i.e. how big is n in nD)"""
        return len(self.lo)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "ndbox({}, {})".format(repr(self.lo), repr(self.hi))

    def nth_child_norder(self, n):
        """
        calculate coordinates of the nth-child given coordinates of the parent nd-box
        where n is given as norder integer 
        (to be able to derive which of the nary childs you get)

        An example for 2D

        Note that 0 = at low side, and 1 = at high side of dimension
        for n=0 -> z=0, y=0, x=0 means the child is located with respect to center:
        (at low side, at low side, at low side)

        n-order:

            n     xy
            -  -----
            0     00
            1     01
            2     10
            3     11

        z-order (note: flipped xy on top):

            n     yx
            -  -----
            0     00
            1     01
            2     10
            3     11

        """
        c = [l + ((h-l) / 2) for h, l in zip(self.hi, self.lo)]
        l = list(self.lo[:])
        h = list(self.hi[:])
        dim = self.dims
        for d in range(dim): # for 3D: 2 = x, y = 1, z = 0
    #        print n, ">>", d
            rshifted = (bitreverse(n, dim) >> d)  # bit shift the child index number to the right by d bits AFTER THE BITS OF N ARE REVERSED
            at_high_side_for_dim_d = bool((rshifted) & 1)
            #        print "", names[d],  ":", binstr(rshifted, dim), "&", binstr(1, dim), "=", at_high_side_for_dim_d
            if at_high_side_for_dim_d:
                # keep the high, shift the low ordinate to the center
                l[d] = c[d]
            else:
                # keep the low, shift the high ordinate to the center
                h[d] = c[d]
        #    print " >>>", l, h
        # return the coordinates of this child box
        return ndbox(l, h)


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


if __name__ == "__main__":
    _test()
