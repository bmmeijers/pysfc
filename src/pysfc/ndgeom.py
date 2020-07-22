import sys
from pysfc.vectorops import mul, dot, norm, cross, unit, sub, add
from itertools import product

# [*] Figure out 'handedness' of planes <--> what should be positive side/negative side?
##
# [*] make 3d work
##
# not arbitrary planes, but planes between min and max value (orthogonal to axis) -> i.e. nd-box selection 
# [*] with all dims bounded 
# [*] with 1 dimension unbounded in one direction
## 
# Test better:
# [ ] 4d
# [ ] 5d
# [ ] 6d
##
# [ ] convenience methods for conversion of regular convex geometry to hyperplanes 
#     - e.g. square, hexagon, combined with other dims, like time

# [*] hyperplane set -> test a box against more hyperplanes --> relate(hyperplanes, box)
#       [*] nd-box -> obtain corner points <--> SFC code

# [ ] Decide on whether to always normalize HyperPlane representation in constructor

# [ ] Move tests to own folder / under test framework

# [ ] Move 'playground' functions out of this module

# [ ] Integrate in DBMS using PL/Python

# we can generate all corner points on an ndbox by using itertools product using [min,max] values per dimension:
# >>> for pt in product([1,2],[3,4],[5,6],[7,8]):
#...     print(pt)
#


class HyperPlane(object):

    def __init__(self, w, b):
        """

        Arguments:

            w   list of float (ω)
            b   float (β)

        "Any hyperplane can be written as the set of points x satisfying
            ω ⋅ x - β = 0,
        where ω is the (not necessarily normalized) normal vector to the hyperplane.
        This is much like Hesse normal form, except that ω  is not necessarily a unit vector." -- https://en.wikipedia.org/wiki/Support_vector_machine

        "Any hyperplane of a Euclidean space has exactly two unit normal vectors." -- https://en.wikipedia.org/wiki/Hyperplane
        """

        self.w = tuple(map(float, w))   # normal vector, not necessarily normalized
        self._b = float(b)              # do use 'offset' property, instead of _b directly
        self.is_normalized = False
        # there should be at least 1 component that is non-zero in w
        assert len(tuple(filter(lambda x: x!= 0.0, self.w))) > 0
        # TODO: should we always normalize in the constructor?
        # this way we can remove the is_normalized property

    def __repr__(self):
        return "HyperPlane({}, {})".format(repr(self.w), repr(self.offset))

    def normalize(self):
        """
        Normalize the representation of the halfplane, 
        i.e. normalize its normal vector and the offset
        """
        self._b = self.offset        # offset property will do the correct thing, also when the plane is already normalized!
        self.w = unit(self.w)        # make w a unit vector by dividing by its own length
        self.is_normalized = True

    @property
    def offset(self):
        """
        Offset (distance) from the origin, along the direction of the normal vector ω

        "The parameter: 
            β / ‖ ω ‖ 
        determines the offset of the hyperplane from the 
        origin *along* the normal vector ω" -- https://en.wikipedia.org/wiki/Support_vector_machine

        """
        if self.is_normalized:
            return self._b
        else:
            return self._b / norm(self.w)

    @property
    def through(self):
        """
        Returns nD-point through which this hyperplane passes
        """
        if self.is_normalized:
            return mul(self.w, self.offset)
        else:
            # normalize the normal vector
            normalized = unit(self.w)
            # multiply with the offset
            return mul(normalized, self.offset)


def signed_distance(point, hyperplane):
    """
    Returns signed distance from point to hyperplane
    """
    return dot(hyperplane.w, point) - hyperplane.offset


def lies_on_closed_upper_halfspace(point, hyperplane): 
    """

    This hyperplane splits the affine space in 2, so we can check whether the given point is on the positive half of these 2 spaces (the upper part)
    As the check is >= (note '=') we also consider the halfplane itself to be included.

    "In the case of a real affine space, in other words when the coordinates are real numbers, this affine space separates the space into two half-spaces, which are the connected components of the complement of the hyperplane, and are given by the inequalities" -- https://en.wikipedia.org/wiki/Hyperplane

    "In geometry, a half-space is either of the two parts into which a plane divides the three-dimensional Euclidean space. More generally, a half-space is either of the two parts into which a hyperplane divides an affine space. That is, the points that are not incident to the hyperplane are partitioned into two convex sets (i.e., half-spaces), such that any subspace connecting a point in one set to a point in the other must intersect the hyperplane. A half-space can be either open or closed. An open half-space is either of the two open sets produced by the subtraction of a hyperplane from the affine space. A closed half-space is the union of an open half-space and the hyperplane that defines it." -- https://en.wikipedia.org/wiki/Half-space_%28geometry%29

    This method returns true when point is above[*] the plane
    (i.e. has a positive distance to the hyperplane)

    [*] MM: FIXME: not sure what I would call this method: above / below | positive/negative distance in nD -- should we test with < 0.0 or with > 0.0
    """
    # TODO: == rename method??? ==
    # this could be a method on a nD point class: point.is_covered_by(plane)

    # |> points with a positive distance to the plane 
    # |> in 2D True if point on the side that the normal ω points towards
    # |> in 3D True if point above the plane (on the side that the normal ω points towards)
    d = signed_distance(point, hyperplane)
    return d >= 0.0 ## closed upper halfspace (including halfplane space itself)


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

    @property
    def vertices(self):
        """Returns an iterator of the corner vertices of this hypercube """
        dims = ((self.lo[n], self.hi[n]) for n in range(self.dims))
        return product(*dims)

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


class ndsphere(object):
    """
    nd sphere defined by center point (nD) and radius (scalar)
    """
    __slots__ = ('center', 'radius')

    def __init__(self, center, radius):
        assert radius >= 0
        self.center = tuple(center)
        self.radius = float(radius)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "ndsphere({}, {})".format(repr(self.center), repr(self.radius))


def as_ndsphere(ndbox):
    """
    Given a box, return a sphere that goes through the box its vertices
    """
    delta = sub(ndbox.hi, ndbox.lo)
    center = add(ndbox.lo, mul(delta, 0.5))
    radius = norm(sub(center, ndbox.lo))
    return ndsphere(center, radius)


def as_hyperplanes(ndbox):
    """
    Given a box, return a list of HyperPlanes that bound the same space
    as the box
    """
    planes = []
    dims = [0.0] * len(ndbox.lo)
    # lo dimensions, positive normal
    for dim, lo in enumerate(ndbox.lo):
        w = dims[:]
        w[dim] = +1.0
        b = lo
        planes += HyperPlane(w, b),
    # hi dimensions, negative normal
    for dim, hi in enumerate(ndbox.hi):
        w = dims[:]
        w[dim] = -1.0
        b = -hi
        planes += HyperPlane(w, b),
    return planes


def covers(hp, sphere):
    """
    The closed, upper halfspace (i.e. including the 'boundary' of the halfspace)
    totally covers the sphere (including the sphere its boundary)
    """
    return signed_distance(sphere.center, hp) - sphere.radius >= 0.0


# TODO: remove following code
#def anyinteract_with_closed_upper_halfspace(sphere, hyperplane): 
#    """
#    This hyperplane splits the affine space in 2 (positive/negative),
#    so we can check whether the given sphere has any interaction with positive part of the space
#    """
#    distance_center_to_plane = signed_distance(sphere.center, hyperplane)
#    return (distance_center_to_plane + sphere.radius) >= 0.0


def intersects(hp, sphere):
    """
    The closed, upper halfspace intersects the sphere
    (i.e. there exists a spatial relation between the two)
    """
    return signed_distance(sphere.center, hp) + sphere.radius >= 0.0

# -----------------------------------------------------------------------------
# -- spatial relate
# -----------------------------------------------------------------------------


#    Outcome of original relate test is:
#
#        0: equal
#        1: contains
#        2: intersects
#        -1: no overlap
#


NO_INTERACT = 0     # no interaction
INTERACT = 1        # possibly some interaction
CONTAINED = 2       # interaction, for sure


def relate_sphere(planes, sphere):
    """
    spatial relation between hyperplanes and nd-sphere

    note, hyperplanes should form a minimal description, 
    i.e. no redundant planes allowed
    """
    for plane in planes:
        if not intersects(plane, sphere):
            return NO_INTERACT
    return INTERACT


### 'PSEUDOCODE':
###def relate(planes, ndbox):
###    outcomes = []
###    for plane in planes:
###        # simple test (*all* corner points of ndbox on wrong side of half plane)
###        # *report* box as outside, when *all* points are at wrong side (i.e. all tests need to be False) -- -1:no overlap
###    # *report* as inside (all corner points at right side of all half planes)
###    # ...
###    # *report* as unknown / to be refined
###    # ...
def relate_box__box_based(planes, box):
    """
    Spatial relation between hyperplanes and nd-box

    note, hyperplanes should form a minimal description, 
    i.e. no redundant planes allowed
    """

#       [ ] small optimization: 3^d - 2^d (may not be needed when using sphere around box?)
#            -- test only new corners for 'kids' boxes in traversal algo of n-ary tree
#           2D =  9 -  4 =  5 refined
#           3D = 27 -  8 = 19 refined
#           4D = 81 - 16 = 65 refined

    # -- test if the box is outside of one of the planes
    # if so, for sure no interaction
    for plane in planes:
        are_covered = [lies_on_closed_upper_halfspace(pt, plane) for pt in box.vertices]
        if True not in are_covered:
            return NO_INTERACT
    # -- if we end up here, figure out interaction
    # check full containment versus partial overlap
    gen = (lies_on_closed_upper_halfspace(pt, plane) for pt in box.vertices for plane in planes)
    if all(gen):
        # all points on correct side of all planes 
        # -> box is fullly contained by planes
        return CONTAINED
    else:
        # some interaction
        return INTERACT


def relate_box__sphere_based(planes, box):
    """
    Spatial relation between hyperplanes and nd-box

    note, hyperplanes should form a minimal description, 
    i.e. no redundant planes allowed

    In 2D, this method is ~2x as fast as the box based method,
    but less accurate when box is close with the hyperplane

    """
    # derive a sphere from the box
    sphere = as_ndsphere(box)
    # -- test if the spere is outside of one of the planes
    # if so, for sure no interaction
    for plane in planes:
        if not intersects(plane, sphere):
            return NO_INTERACT
    # -- if we end up here, figure out interaction
    # check full containment versus partial overlap
    gen = (covers(plane, sphere) for plane in planes)
    if all(gen):
        # sphere fully covered by all planes -> box is fullly contained by planes
        return CONTAINED
    else:
        # some interaction
        # note, box may still be fully contained (but not sure due to approximation by sphere)
        # at least the box is partially overlapping
        # TODO: find out whether 'infinite' refinement will happen,
        # when we specify hyperplanes close to sides of index-boxes
        # TODO: here, we could check also the vertices of the ndbox
        # (but maybe this should be done if and only if, based on an ε-distance, we can deduce that this is necessary, i.e. when hyperplane is very close to the border of the sphere)?
        return INTERACT


# ----------------------------------------------------------------------------#
def relate_box__sweep_based(planes, box):
    """
    spatial relation between hyperplanes and nd-box
    """
    # pre-condition: some planes to test against
    assert len(planes) > 0
    result = -1
    for plane in planes:
        # first hit point sweeping the box with the hyperplane
        enter = sweep_box_with_plane_enter(plane, box)
        enter_dist = signed_distance(enter, plane)
        if enter_dist >= 0:
            # only override result when not yet seen
            if result == -1:
                result = CONTAINED
        else:
            # perform distance comparison first (diagonal length & side of box)
            # by doing that, we defer sweeping for the far point
            abs_enter_dist = abs(enter_dist)
            if abs_enter_dist > diagonal_dist:
                return NO_INTERACT
# ----------------------------------------------------------------------------#
# - alternative:                                                         -----#
# - not considering the exit point, just accepting possibly interact     -----#
# ----------------------------------------------------------------------------#
##            else:
##                result = INTERACT    
# ----------------------------------------------------------------------------#
            elif abs_enter_dist < box_side_dist:
                result = INTERACT
            else:
                # get last hit point sweeping the box with the hyperplane
                exit = sweep_box_with_plane_exit(plane, box)
                exit_dist = signed_distance(exit, plane)
                if exit_dist >= 0:
                    result = INTERACT
                else:
                    return NO_INTERACT
    # post-condition
    # result here should be either INTERACT / CONTAINED
    assert result in (INTERACT, CONTAINED)
    return result


def sweep_box_with_plane_enter(hyperplane, box):
    """ Return enter points when 'sweeping' the box with the hyperplane

    enter: point that is 'hit' first while sweeping the box with the hyperplane
    along the direction of the normal of the hyperplane
    """
    enter = [0.0] * box.dims
    for i, w in enumerate(hyperplane.w):
        if w < 0:
            enter[i] = box.hi[i]
        else:
            enter[i] = box.lo[i]
    return enter


def sweep_box_with_plane_exit(hyperplane, box):
    """ Return exit point when 'sweeping' the box with the hyperplane

    exit: point where hyperplane leaves the box while sweeping the box with the
    hyperplane along the direction of the normal of the hyperplane
    """
    exit = [0.0] * box.dims
    for i, w in enumerate(hyperplane.w):
        if w < 0:
            exit[i] = box.lo[i]
        else:
            exit[i] = box.hi[i]
    return exit


# -- 
# -- Querying ------------------------------------------------
# -- 

# TODO: also make hquery for the Hilbert curve

def nquery(query_planes, query_hi=1023, use_sphere=True):
    # TODO: find a way to easily express query_hi correctly
    # previously we first found query_hi (largest dim value) based on 
    # max of the query box hi coordinates,
    # this is basically the size of the cube
    #
    # we could evaluate the polytope (query_planes) and see if it is finite
    # and from that determine the corners, and hence its extent
    #
    # otherwise, we need to make sure that it is large enough -> known from
    # the metadata for scaling/translating the cube and the resolution

    from pysfc import _determine_bits, _nchunks_coord, _nchunks_key

    if use_sphere:
        relate = relate_box__sphere_based #relate_box__box_based
    else:
        relate = relate_box__box_based #relate_box__box_based

#    maxbits = 63
    ndims = len(query_planes[0].w) #query.dims
    # get how many bits we need
    # to represent the largest number inside the query box
    # --> 2**(mbits_needed) is the maximum size of a side 
    #     of the domain that we need
    mbits_needed = _determine_bits(query_hi, 2) # max(query.hi)
#    mbits = maxbits // ndims
    npath = ()
    # post order tree traversal gives nodes in order we want
    # for this, two stacks are used (stack + paths)
    # -- https://algorithms.tutorialhorizon.com/binary-tree-postorder-traversal-non-recursive-approach/
    paths = []
    stack = [npath]
    while stack:
        npath = stack.pop() 
        cur_level = len(npath)
        lo = _nchunks_coord(npath, ndims)
        # side_at_depth 
        # how large side of the nd-cube is at this depth
        side_at_depth = 2**(mbits_needed - cur_level)
        # print('side len', side_at_depth)
        lo = tuple(map(lambda x: x * side_at_depth, lo))
        hi = tuple(map(lambda x: x + side_at_depth, lo))
        cur_node = ndbox(lo, hi)
#        print(visualize_box_2d(cur_node))

        # [ ] TODO: == optimization ==
        # we can omit testing a plane for kids boxes if
        # parent box is completely at correct side of this hyperplane already
        # -> result of spatial relate test will again be true for splitted boxes
        ndcmp = relate(query_planes, cur_node)
        # -- full overlap, contained
        if ndcmp == CONTAINED:
            paths.append(npath)
        # -- partial overlap
        elif ndcmp == INTERACT:
            # FIXME: re-introduce maxdepth argument !
            if cur_level < mbits_needed:
                # we have not yet reached the lowest level, recurse
                for ncode in range(2**ndims):
                    new_path = npath + (ncode, )
                    stack.append(new_path)
            else:
                # we are not allowed to go further, so use this partial
                # overlapping range as is
                paths.append(npath)
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp == NO_INTERACT:
            continue
    result = []
    while paths:
        npath = paths.pop()
        # if the npath has less than mbits_needed values,
        # we need to add 0's to the end of the list:
        # -> this happens inside _nchunks_key
        start = _nchunks_key(npath, mbits_needed, ndims)
        side_at_depth = 2 ** (mbits_needed - len(npath))
        range_size = side_at_depth**ndims
        end = start + range_size
#        print npath, start, end, range_size
        result.append((start, end))
    return result


# hp = HyperPlane([0.0, 0.0], -3.9) # <<  has to raise an error, only 0.0's in w (no direction defined at all)

#hp0 = HyperPlane([0.0, -1.0], 5.0)
# hp1 = HyperPlane([+0.70, -0.70], 0.0)
# hp2 = HyperPlane([0.0, 1.0], -5.0)

# --
# -- Unit tests ----------------------------------------------------------------
# --

def test_normalize():
    hp = HyperPlane([10.0], 8.0)

    # before normalizing
    assert norm(hp.w) == 10
    assert hp._b == 8.0
    assert norm(hp.w) != 1 # redundant test
    assert hp.offset == 0.8
    assert hp.through == (0.8,)

    hp.normalize()
    # after normalizing
    assert norm(hp.w) == 1
    assert hp.offset == 0.8
    assert hp.through == (0.8,)
    assert hp._b == 0.8


def test_1d_signed_distance_point_hyperplane():
    hp = HyperPlane([1.0], 0.0)
    for p in range(-10, 11):
        assert signed_distance([p], hp) == p

    hp = HyperPlane([-1.0], 0.0)
    for p in range(-10, 11):
        assert signed_distance([p], hp) == -p

    hp = HyperPlane([1.0], 5.0)
    for p in range(-10, 11):
        assert signed_distance([p], hp) == p - 5.0

    hp = HyperPlane([-1.0], -5.0)
    for p in range(-10, 11):
        assert signed_distance([p], hp) == -p + 5.0

    hp = HyperPlane([-1.0], 5.0)
    for p in range(-10, 11):
        assert signed_distance([p], hp) == -p - 5.0



def test_conversion():
    # when we convert n-d box to n-d sphere

    # redundant with n-d tests below (but these are more explicit, 2D/3D)
    s = as_ndsphere(ndbox([0, 0], [2, 2]))
    assert s.center == (1., ) * 2       # (1, 1)
    assert s.radius == 2.0**0.5         # √2

    s = as_ndsphere(ndbox([0, 0, 0], [2, 2, 2]))
    assert s.center == (1., ) * 3       # (1, 1, 1)
    assert s.radius == 3.0**0.5         # √3

    # -- n-d tests

    # box from [0, 2] in every dimension
    # we should end up at 1 in each dimension, with the radius being √n
    for n in range(1, 11):
        s = as_ndsphere(ndbox([0] * n, [2] * n))
        assert s.center == (1., ) * n       # (1, ..., 1)
        assert s.radius == n ** 0.5         # √n

    # box from [-2, 0] in every dimension
    # we should end up at -1 in each dimension, with the radius being √n
    for n in range(1, 11):
        s = as_ndsphere(ndbox([-2] * n, [0] * n))
        assert s.center == (-1., ) * n       # (-1, ..., -1)
        assert s.radius == n ** 0.5         # √n


def test_1d_anyinteract():
    # 1d sphere, centered at 11 radius 1
    s = as_ndsphere(ndbox([10], [12]))
    hp = HyperPlane([1], 0) # positive half, through [0,]
    assert intersects(hp, s)
    hp = HyperPlane([1], 10) # positive half, through [10,]
    assert intersects(hp, s)
    hp = HyperPlane([1], 12) # positive half, through [12,]
    assert intersects(hp, s)
    hp = HyperPlane([1], 12.00001) # positive half, through [12.00001,]
    assert not intersects(hp, s)
    hp = HyperPlane([-1], 0) # negative half, through [0,]
    assert not intersects(hp, s)
    hp = HyperPlane([-1], -10) # negative half, through [10,]
    assert intersects(hp, s)


def test_1d_relate_sphere():
    s = as_ndsphere(ndbox([10], [12]))
    planes = [
        HyperPlane([1], 0), # positive half, through [0,]
        HyperPlane([-1], -100), # negative half, through [100,]
    ]
    assert relate_sphere(planes, s) == True

    s = as_ndsphere(ndbox([10], [12]))
    planes = [
        HyperPlane([1], 0), # positive half, through [0,]
        HyperPlane([-1], -5), # negative half, through [100,]
    ]
    assert relate_sphere(planes, s) == False


def test_1d_relate_box(relate = relate_box__sphere_based):

    # box interval completely inside
    b = ndbox([10], [12])
    planes = [
        HyperPlane([1], 0), # positive half, through [0,]
        HyperPlane([-1], -100), # negative half, through [100,]
    ]
    assert relate(planes, b) == CONTAINED, relate(planes, b)

    # box interval partly inside
    b = ndbox([10], [101])
    planes = [
        HyperPlane([1], 0), # positive half, through [0,]
        HyperPlane([-1], -100), # negative half, through [100,]
    ]
    assert relate(planes, b) == INTERACT

    # box interval lies outside
    b = ndbox([-10], [-1])
    planes = [
        HyperPlane([1], 0), # positive half, through [0,]
        HyperPlane([-1], -100), # negative half, through [100,]
    ]
    assert relate(planes, b) == NO_INTERACT


def test_2d_relate_box(relate = relate_box__sphere_based):
    b = ndbox([1, 1], [2, 2])
    hp = HyperPlane([1, 0], 0) # +x, through [0, 0]
    hp.normalize()
    planes = [
        hp,
        HyperPlane([0, 1], 0), # +y, through [0, 0]
    ]
    assert relate(planes, b) == CONTAINED

    b = ndbox([0,0], [2,2])
    hp = HyperPlane([-1, 1], 0) # -x,+y, through [0, 0]
    hp.normalize()
    planes = [
        hp,
        HyperPlane([0, 1], 0), # +y, through [0, 0]
    ]
    assert relate(planes, b) == INTERACT

    b = ndbox([-12,-2], [-11,-1])
    hp = HyperPlane([-1, 1], 0) # -x,+y, through [0, 0]
    hp.normalize()
    planes = [
        hp,
        HyperPlane([0, 1], 0), # +y, through [0, 0]
    ]
    assert relate(planes, b) == NO_INTERACT



# --
# -- Playground ----------------------------------------------------------------
# --


### 1D case
def case_1d():
    hp1 = HyperPlane([-1.0], -8.0)
    print('offset  {}'.format(hp1.offset), file=sys.stderr)
    print('through {}'.format(hp1.through), file=sys.stderr)
    left, right = -50, 51
    for x in range(left, right):
        print ( (x,), lies_on_closed_upper_halfspace( (x,), hp1 ) )


def case_1d_range():
    ### 1D case - 'range search'
    hp1 = HyperPlane([1.0], 8.0)    # [8,>
    hp2 = HyperPlane([-1.0], -16.0) #       <,16]
    print('offset  {}'.format(hp1.offset), file=sys.stderr)
    print('through {}'.format(hp1.through), file=sys.stderr)
    print('offset  {}'.format(hp2.offset), file=sys.stderr)
    print('through {}'.format(hp2.through), file=sys.stderr)
    left, right = -50, 51
    for x in range(left, right):
        print ( (x,), lies_on_closed_upper_halfspace( (x,), hp1 ) and lies_on_closed_upper_halfspace( (x,), hp2 ))
    #    print ( (x,), hp1.on((x,)) and hp2.on((x,)) )


def case_2d():
    ### 2D case
    hp = HyperPlane([0.0, -1.0], -5.0)
#    hp = HyperPlane([-1.0, 0.0], -5.0)
    print('offset  {}'.format(hp.offset), file=sys.stderr)
    print('through {}'.format(hp.through), file=sys.stderr)
    # test with some points
    left, right = -10, 11
    IN = ([], [])
    OUT = ([], [])
    for x in range(left, right):
        for y in range(left, right):
            pt = tuple(map(float, (x, y)))
            if lies_on_closed_upper_halfspace(pt, hp):
                for i in range(2):
                    IN[i].append(pt[i])
            else:
                for i in range(2):
                    OUT[i].append(pt[i])
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(IN[0], IN[1], c='g', marker='P')
    ax.scatter(OUT[0], OUT[1], c='r', marker='+')
    ax.set_xlim(left, right)
    ax.set_ylim(left, right)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    plt.show()


def case_3d():
    ## 3D case
    hp3 = HyperPlane([0.0, 0.0, -1.0], 2.5)
    print('offset  {}'.format(hp3.offset), file=sys.stderr)
    print('through {}'.format(hp3.through), file=sys.stderr)
    print(cross([1, 0, 0], [0, 1, 0]), file=sys.stderr)
    X = []
    Y = []
    Z = []

    left, right = -10, 11
    for x in range(left, right):
        for y in range(left, right):
            for z in range(left, right):
                pt = (x, y, z)
                if lies_on_closed_upper_halfspace(pt, hp3):
                    X.append(x)
                    Y.append(y)
                    Z.append(z)

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, c='g', marker='+')

    ax.set_xlim3d(left, right)
    ax.set_ylim3d(left, right)
    ax.set_zlim3d(left, right)
    #ax.axis('equal')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


def case_3d_combine():
    ## 3D case
    planes = [
        # x [-1.5 , + 1.5]
        HyperPlane([1.0, 0.0, 0.0], -1.5),
        HyperPlane([-1.0, 0.0, 0.0], -1.5),
        # y [-3.5 , + 3.5]
        HyperPlane([0.0, 1.0, 0.0], -3.5),
        HyperPlane([0.0, -1.0, 0.0], -3.5),
        # z [-4.5 , + 2.5]
        HyperPlane([0.0, 0.0, 1.0], -4.5),
        #HyperPlane([0.0, 0.0, -1.0], -2.5),  # axis aligned plane
        #HyperPlane([-1.0, -1.0, -1.0], -2.5), # tilted plane
    ]

    for hp in planes:
        print('offset  {}'.format(hp.offset), file=sys.stderr)
        print('through {}'.format(hp.through), file=sys.stderr)

    X = []
    Y = []
    Z = []

    left, right = -10, 11
    for x in range(left, right):
        for y in range(left, right):
            for z in range(left, right):
                pt = (x, y, z)
                if all(lies_on_closed_upper_halfspace(pt, hp) for hp in planes):
                    X.append(x)
                    Y.append(y)
                    Z.append(z)

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, c='b', marker='+')

    ax.set_xlim3d(left, right)
    ax.set_ylim3d(left, right)
    ax.set_zlim3d(left, right)
    #ax.axis('equal')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()


def visualize_box_2d(box):
    corners = [pt for pt in box.vertices]
    coords = ", ".join("{:.1f} {:.1f}".format(x, y) for x, y in [corners[0], corners[2], corners[3], corners[1], corners[0]])
    return 'POLYGON(({}))'.format(coords)


def visualize_2d_boxes_against_some_planes():
    hp0 = HyperPlane([-1, 1], 0) # -x,+y, through [0, 0]
    hp0.normalize()

    from math import cos, sin, radians

    angle_deg = 285. #315.
    hp1 = HyperPlane([cos(radians(angle_deg)), sin(radians(angle_deg))], -20) # +x,-y, through [?,?]
    hp1.normalize()

    planes = [
        hp0,
        HyperPlane([0, 1], 0), # +y, through [0, 0]
        HyperPlane([-1, 0], -10), # -x, through [10, 0]
        hp1,
    ]
    # output halfplanes
    filenm = '/tmp/planes.txt'
    with open(filenm, 'w') as fh:
        pass
    for plane in planes:
        visualize_halfplane_2d(plane, filenm)

    outcomes = []
    step = 1

#    start = -100
#    end = 110

    for x in range (-100, 110, step):
        for y in range (-100, 110, step):
            b = ndbox((x,y), (x+step, y+step))
            relate_s = relate_box__sphere_based(planes, b)
            relate_b = relate_box__box_based(planes, b)
            relate_sw = relate_box__sweep_based(planes, b)
            wkt = visualize_box_2d(b)
            outcomes.append((wkt, relate_s, relate_b, relate_sw))
    # write for visualization in qgis to csv file
    import csv
    with open('/tmp/wkt.csv', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['geometry', 'relate_sphere', 'relate_box', 'relate_sweep'])
        writer.writerows(outcomes)


def visualize_halfplane_2d(h, filenm='/tmp/planes.txt'):
    from pysfc.vectorops import rotate90ccw, rotate90cw, add, mul
    # make 2D line of 2000 units long + the normal
    ccw = rotate90ccw(h.w)
    cw = rotate90cw(h.w)
    start = add(mul(cw, 1000.), h.through)
    end = add(mul(ccw, 1000.), h.through)
    with open(filenm, 'a') as fh:
        fh.write("LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tplane".format(start, end))
        fh.write('\n')
        fh.write("LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tnormal".format(h.through, add(h.through, h.w)))
        fh.write('\n')


def measure_performance(which):
    """
    To run a mini-timing benchmark, in the shell run:

    $ python3 -O -mtimeit -n 5 -r 5  -s'from pysfc import ndgeom' 'ndgeom.measure_performance("sphere")'
    $ python3 -O -mtimeit -n 5 -r 5  -s'from pysfc import ndgeom' 'ndgeom.measure_performance("corner")'
    $ python3 -O -mtimeit -n 5 -r 5  -s'from pysfc import ndgeom' 'ndgeom.measure_performance("sweep")'
    """

    # TODO: == introduce higher dimensions ==
    # -> the higher the dimensions, the more vertices per hyperbox,
    # so the larger the difference with the sphere/sweep representation

    hp0 = HyperPlane([-1, 1, 0, 0], 0) # -x,+y, through [0, 0]
    hp0.normalize()
    planes = [
        hp0,
        HyperPlane([0, 1, 0, 0], 0), # +y, through [0, 0]
        HyperPlane([-1, 0, 0, 0], -10), # +x, through [10, 0]
    ]
    step = 1

    lut = {
        'sphere': relate_box__sphere_based,
        'corner': relate_box__box_based,
        'sweep': relate_box__sweep_based,
    }
    r = lut[which]
    lo, hi = -5, 6
    for x in range (lo, hi, step):
        for y in range (lo, hi, step):
            for z in range (lo, hi, step):
                for k in range (lo, hi, step):
                    b = ndbox((x,y,z, k), (x+step, y+step, z+step, k+step))
                    relate = r(planes, b)


def dist_comp():

    # TODO: == transform into test ==

    # positive plane through x = 0
    hp = HyperPlane([1], 0)

    sphere = as_ndsphere(ndbox([1], [3]))
    assert signed_distance(sphere.center, hp) == 2.0

    sphere = as_ndsphere(ndbox([-3], [0]))
#    assert signed_distance(sphere.center, hp) == -2.0
    print(signed_distance(sphere.center, hp))
    print(signed_distance(sphere.center, hp) + sphere.radius)

    sphere = as_ndsphere(ndbox([0], [2]))
    assert signed_distance(sphere.center, hp) == +1.0

    # partly overlaps positive halfspace
    # why?
    sphere = as_ndsphere(ndbox([-1], [1]))
    assert signed_distance(sphere.center, hp) == 0.0

    sphere = as_ndsphere(ndbox([-1.5], [0.5]))
    assert signed_distance(sphere.center, hp) == -0.5
    print(signed_distance(sphere.center, hp) + sphere.radius >= sphere.radius)

    for i in range(-10, 10):
        i *= 0.5
        s = ndsphere((i,),1)
        print(i, s, '\t', signed_distance(s.center, hp), '\t', 'covers:', covers(hp, s), '\t', 'overlaps:', overlaps(hp, s), '\t', signed_distance(s.center, hp) + s.radius, signed_distance(s.center, hp) - s.radius)


def try_orig_nquery():
    from pysfc.pysfc import nquery as orig_nquery
    from pysfc.pysfc import hquery as orig_hquery
    from pprint import pprint
    #b = ndbox([2.0,2.0],[30.0,30.0])
    b = ndbox([15.4, 15.4], [15.6, 15.6])
    #b = ndbox([0.0, 0.0], [2.0, 2.0])
    #b = ndbox([0.0, 0.0], [4.0, 4.0])
    #b = ndbox([1.0, 1.0], [3.0, 3.0])
    pprint(b)
    pprint(orig_nquery(b))
    pprint(nquery(as_hyperplanes(b)))
    pprint(orig_hquery(b))


def tests():
    test_normalize()
    test_1d_signed_distance_point_hyperplane()
    test_1d_anyinteract()
    test_1d_relate_sphere()
    test_1d_relate_box()
    test_2d_relate_box()
    test_conversion()
    print('tests done.')


def playground():
#    case_1d()
#    case_1d_range()
#    case_2d()
    case_3d()
#    case_3d_combine()
    pass


def main():
    tests()
    # playground()
    # visualize_2d_boxes_against_some_planes()


if __name__ == "__main__":
    main()
