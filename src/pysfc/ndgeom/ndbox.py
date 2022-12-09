from itertools import product
from pysfc.ndgeom.hyper_plane import HyperPlane

def bitreverse(n, bits):
    N = 1 << bits           # find N: shift left 1 by the number of bits
    nrev = n                # nrev will store the bit-reversed pattern
    for i in range(1, bits):
        n >>= 1
        nrev <<= 1
        nrev |= n & 1       # give LSB of n to nrev
    nrev &= N - 1           # clear all bits more significant than N-1
    return nrev


class ndbox(object):
    """A nD box object"""
    __slots__ = ('lo', 'hi', 'dims')

    def __init__(self, lo, hi):
        """lo and hi attributes are lists with equal number of coordinates, giving extents in lo and hi direction of specific dimension"""
        assert len(lo) == len(hi)
        assert all([lo[d] <= hi[d] for d in range(len(lo))])
        self.lo = tuple(lo)
        self.hi = tuple(hi)
        # How many dimensions does this box have (i.e. how big is n in nD)
        self.dims = len(self.lo)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "ndbox({}, {})".format(repr(self.lo), repr(self.hi))

    @property
    def center(self):
        half_size = self.half_size
        center = [self.lo[d] + half_size[d] for d in range(self.dims)]
        return center

    @property
    def half_size(self):
        half_size = [0.5 * (self.hi[d] - self.lo[d]) for d in range(self.dims)]
        return half_size

    @property
    def vertices(self):
        """Returns an iterator of the corner vertices of this hypercube """
        dims = ((self.lo[n], self.hi[n]) for n in range(self.dims))
        return product(*dims)

    def as_wkt_2d(self):
        # x = 0
        # y = 1
        pts = [
            (self.lo[0], self.lo[1]), # lo [x] - lo [y]
            (self.hi[0], self.lo[1]), # hi-lo 
            (self.hi[0], self.hi[1]), # hi-hi
            (self.lo[0], self.hi[1]), # lo-hi
            (self.lo[0], self.lo[1]), # lo-lo
        ]
        c = ", ".join(["{0[0]} {0[1]}".format(pt) for pt in pts])
        return f"POLYGON (({c}))"

    def __getitem__(self, i):
        if i == 0:
            return self.lo
        elif i == 1:
            return self.hi
        else:
            raise ValueError(f'no such axis {i}')

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
            # print ("", d, ":", binstr(rshifted, dim), "&", binstr(1, dim), "=", at_high_side_for_dim_d)
            if at_high_side_for_dim_d:
                # keep the high, shift the low ordinate to the center
                l[d] = c[d]
            else:
                # keep the low, shift the high ordinate to the center
                h[d] = c[d]
        #    print " >>>", l, h
        # return the coordinates of this child box
        return ndbox(l, h)

    def as_hyperplanes(self):
        """
        Given a ndbox, return a list of HyperPlanes that bound the same space
        as the box
        """
        planes = []
        dims = [0.0] * len(self.lo)
        # lo dimensions, normal negative, towards origin (assuming lo > 0)

        # FIXME: does this also work with the box on both sides of the origin
        # (i.e. having negative coordinates)
        for dim, _ in enumerate(self.lo):
            w = dims[:]
            w[dim] = -1.0
            planes.append(HyperPlane.from_normal_and_point(w, self.lo))
        # hi dimensions, negative normal
        for dim, _ in enumerate(self.hi):
            w = dims[:]
            w[dim] = +1.0
            planes.append(HyperPlane.from_normal_and_point(w, self.hi))
        return planes