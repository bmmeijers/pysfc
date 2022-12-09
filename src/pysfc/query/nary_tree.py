from __future__ import annotations
from math import ceil    

# from multiprocessing.sharedctypes import Value
from typing import List, Tuple

# from collections import deque

from pysfc.ndgeom.ndbox import (
    ndbox,
)  # , as_hyperplanes, relate_box__sweep_based__with_mask

# from pysfc.ndgeom. import as_hyperplanes, relate_box__sweep_based__with_mask, Interaction, ConvexPolytope, HyperPlane
from pysfc.ndgeom.convex_polytope import ConvexPolytope
from pysfc.ndgeom.hyper_plane import HyperPlane, _visualize_hyperplane_2d
from pysfc.query.relate import relate_planes_box__sweep_based__with_mask, Interaction, relate_box_box

# from pysfc.ndgeom. import visualize_halfplane_2d
# from pysfc.ndgeom import Interaction
# from pysfc.encode_decode import nenc, henc
from pysfc.encode_decode import _transpose_bits, _determine_bits
from pysfc.encode_decode import (
    _hinitial_start_end,
    _hgray_decode_travel,
    _hchild_start_end,
)
from pysfc.encode_decode import _nchunks_to_hchunks, _chunks_key

# hchunks[j] = _hgray_decode_travel(start, end, mask, nchunk)
# start, end = _hchild_start_end(start, end, mask, hchunks[j])


def bitreverse(n, bits):
    N = 1 << bits  # find N: shift left 1 by the number of bits
    nrev = n  # nrev will store the bit-reversed pattern
    for i in range(1, bits):
        n >>= 1
        nrev <<= 1
        nrev |= n & 1  # give LSB of n to nrev
    nrev &= N - 1  # clear all bits more significant than N-1
    return nrev


def format4(s, pad="0"):
    return "{s:{spec}}".format(spec=f"{pad}>4", s=s)


class ndvoxel:
    """Addressing a hypervoxel in nD hypercube spanned by SFC
    (geometrically interpreted as a k-way tree, where k=2**ndims)
    """

    def __init__(self, address, level, mbits):
        # FIXME: does this address follow xyzi..., or ...izyx order???
        # i.e. what is the ordering of the different dimensions?
        #
        # - deal with Hilbert means that we need the path to the root
        #   (as hilbert enc/dec depends on level above recursively)
        #
        self.address = tuple(address)
        self.level = level
        # FIXME: is it necessary to be part of the properties here?
        # maybe better to have a ndSpecification-class that supplies nDims, mBits,
        # and those methods that need it, you pass an instance of this class as argument
        self.mbits = mbits
        self.ndims = len(self.address)

    def __str__(self):
        return f"ndvoxel(address={self.address}, level={self.level}, mbits={self.mbits}, ndims={self.ndims})"

    def zcode(self):
        """The Z-code index of this voxel in this subtree"""
        x = 0
        # the last bits in the address co-ordinates are the zcode of this ndvoxel
        for c in self.address:
            x <<= 1
            x += c & 1
        return x

    def ncode(self):
        """The N-code index of this voxel in this subtree"""
        # the ncode of this ndvoxel are the last bits of this co-ordinates,
        # where bits are traversed in reversed order
        # equivalent to bit_reverse(self.zcode(), self.ndims)
        x = 0
        for c in reversed(self.address):
            x <<= 1
            x += c & 1
        return x

    def ncodes(self):
        """the ncode list from the top of the tree to this subvoxel"""
        return _transpose_bits(self.address, self.level)

    # def hcode(self):
    #     """ using the ncode of this node, we can find the hcode
    #     """
    #     # FIXME:
    #     # -> this needs to traverse the tree up to the root
    #     # as the rotation/mirroring is dependent on parent nodes in the k-way tree
    #     # this traversal upwards is inefficient compared to performing same
    #     # procedure top-down (duplicate operations performed for every level
    #     # in the tree, so it is best to avoid using this method)
    #     # -> It may be more efficient if we keep track of this information, while traversing the
    #     # tree, but we then better stack also the enter, exit and parent H-code during tree traversal
    #     # it seems
    #     # -> it may be necessary to do this in case we want to guarantee the sort order
    #     # that comes out of the tree traversal of the ranges
    #     mask = 2 ** self.ndims - 1
    #     enter, exit = self._hilbert_enter_exit()
    #     return _hgray_decode_travel(enter, exit, mask, self.ncode())

    # def _hilbert_enter_exit(self):
    #     """ look upwards and find rotation / mirror of above lying parent ndvoxel's"""
    #     mask = 2 ** self.ndims - 1
    #     if self.level == 0:
    #         start, end = _hinitial_start_end(self.mbits, self.ndims)
    #     else:
    #         parent = self.get_parent()
    #         start, end = parent._hilbert_enter_exit()
    #         start, end = _hchild_start_end(start, end, mask, parent.hcode())
    #     return start, end

    def as_ndbox(self) -> ndbox:
        side_at_level = 2 ** (self.mbits - self.level)
        lo = [c * side_at_level for c in self.address]
        # FIXME: should we in- or exclude the upper bound inside the box hi-bounds?
        # Maybe inside, if we query with floating point halfplane representation?
        # (or we have a small gap of 1 unit between the boxes)
        # if we subtract 1 here, as_morton_range and as_hilbert_range work as expected
        # if we do not do that here, we have to compensate there...
        hi = [item + side_at_level for item in lo]
        return ndbox(lo, hi)

    def as_morton_range(self) -> Tuple(int, int):
        # it should be possible to get a sfc key for the start and end of this voxel
        # box = self.as_ndbox()
        # return nenc(box.lo), nenc(box.hi)
        side_at_level = 2 ** (self.mbits - self.level)
        steps_in_subcube = side_at_level ** self.ndims - 1
        sfc_key = _chunks_key(self.ncodes(), self.mbits, self.ndims)
        return sfc_key, sfc_key + steps_in_subcube + 1

    def as_hilbert_range(self):
        side_at_level = 2 ** (self.mbits - self.level)
        steps_in_subcube = side_at_level ** self.ndims - 1
        hchunks = _nchunks_to_hchunks(self.ncodes(), self.mbits, self.ndims)
        sfc_key = _chunks_key(hchunks, self.mbits, self.ndims)
        return sfc_key, sfc_key + steps_in_subcube + 1

        # # Why is hilbert so strange? ==>
        # # As the entrance point is elsewhere in the voxel, *stupid*
        # # so the range does not make sense by just looking at the lower left / upper right coordinate...
        # # which 2 coordinates of the nd-box we need to consider depends on the hilbert pattern
        # # (these depends on mirror / rotate depending on what happened for the parent)
        # box = self.as_ndbox()
        # # result = henc(box.lo), henc(box.hi)
        # result = []
        # # an inefficient solution here:
        # # we get all the corner points on the box, and feed these *all* to henc
        # # after which the min and max of these forms the sfc range
        # for corner in box.vertices:
        #     result.append(henc(corner))
        # # raise NotImplementedError('as_hilbert method is not yet correct (we need to twist and turn to make this right)')
        # # FIXME: this is correct, but inefficient; the higher the dimension, the more inefficient it gets
        # return min(result), max(result)

    # get_childs(self):
    # entwine does the following in 3D:
    # Each node at depth D may have up to 8 child nodes at depth D + 1 which represent bisected sub-volumes.
    # To zoom in, depth D is incremented by one, and each of X, Y, and Z are doubled and then possibly incremented by one.
    # Coordinates are treated as Cartesian during this bisection,
    # so X -> 2 * X represents the sub-volume where the new maximum X value is the midpoint from its parent (go left in binary tree),
    # and X -> 2 * X + 1 represents the sub-volume where the new minimum X value is the midpoint from its parent (go right in binary tree).
    def get_childs_along_norder(self) -> List[ndvoxel]:
        # FIXME: can we find a visiting order of child nodes here that corresponds with how the SFC
        # curve runs through the voxels? ---
        # this way we get ranges automatically in sorted order
        #
        # YES, we can, we just need to bit_reverse the child_index variable
        # (maybe making a look-up table for this would be a few magnitudes faster)
        new_level = self.level + 1
        if new_level > self.mbits:
            return []
        childs = [None] * 2 ** self.ndims
        new_mbits = self.mbits
        for child_index in range(2 ** self.ndims):  ## zcode / ncode
            # for each child
            # e.g. 0, 1, 2, ..., 7 for ndims = 3
            # To get the ncode of this child, we reverse the ndims bits inside the child_index
            ncode = bitreverse(child_index, self.ndims)
            # get the address of one of the new childs, depending on how many steps we make in the subcube
            new_address = [None] * self.ndims
            for step, c in enumerate(self.address):
                # tmp = (c << 1) | ((2**self.ndims - child_index) >> dim & 0x1)
                # dim = reverse_bits(dim)
                tmp = (c << 1) | (ncode >> step & 1)
                new_address[step] = tmp
            childs[child_index] = ndvoxel(new_address, new_level, new_mbits)
        return childs

    def get_parent(self):
        if self.level > 0:
            # use right shift ? >> 1
            return ndvoxel([c // 2 for c in self.address], self.level - 1, self.mbits)
        else:
            return None


class ndvoxel_hilbert(ndvoxel):
    """hilbert specifics (most notable entry/exit point at (sub-hypercube) for ndvoxel"""

    def __init__(self, address, level, mbits, enter, exit):
        super().__init__(address, level, mbits)
        self.enter = enter  ## FIXME: enter/exit may be confusing, these are *not* the entry/exit point for the sweep based test
        self.exit = exit
        self.mask = 2 ** self.ndims - 1

    def get_childs_along_horder(self, enter, exit) -> List[ndvoxel_hilbert]:
        new_level = self.level + 1
        if new_level > self.mbits:
            return []
        childs = [None] * 2 ** self.ndims
        for child_index in range(2 ** self.ndims):  ## zcode / ncode
            ncode = bitreverse(
                child_index, self.ndims
            )  # To get the N-order, we reverse the dimensions inside the child_index
            new_hcode = _hgray_decode_travel(enter, exit, self.mask, ncode)
            new_enter, new_exit = _hchild_start_end(enter, exit, self.mask, new_hcode)
            # get the address of one of the new childs, depending on how many steps we make in the subcube
            new_address = [None] * self.ndims
            for step, c in enumerate(self.address):
                # tmp = (c << 1) | ((2**self.ndims - child_index) >> dim & 0x1)
                # dim = reverse_bits(dim)
                tmp = (c << 1) | (child_index >> step & 1)
                new_address[step] = tmp
            childs[new_hcode] = ndvoxel_hilbert(
                new_address, new_level, self.mbits, new_enter, new_exit
            )
        return childs

    def get_parent(self):
        raise NotImplementedError(
            "traversing up, not supported for hilbert due to enter / exit properties unknown in stateless manner -- only by traversing tree top-down"
        )


# def relate_ndbox(rect, qrt, _):
#     """
#     Spatial relationship between two nd-boxes.

#     Outcome can be:
#         0: equal
#         1: contains
#         2: intersects
#         -1: no overlap
#     """
#     # how many dimensions to check?
#     dims = rect.dims

#     # equal, all coordinates are equal
#     # FIXME: is this check necessary?
#     ncmp = 1
#     for d in range(dims):
#         ncmp &= rect.lo[d] == qrt.lo[d] and rect.hi[d] == qrt.hi[d]
#     if ncmp:
#         return Interaction.CONTAINED, _

#     # FIXME: would one loop in which we check for both contained + intersects faster than
#     #        doing the checks one after the other?

#     # fully contains, rect fully contains qrt
#     ncmp = 1
#     for d in range(dims):
#         ncmp &= rect.lo[d] <= qrt.lo[d] and rect.hi[d] >= qrt.hi[d]
#     if ncmp:
#         return Interaction.CONTAINED, _

#     # intersects, the two nd-boxes interact
#     # (either on the boundary or internally)
#     ncmp = 1
#     for d in range(dims):
#         ncmp &= rect.lo[d] < qrt.hi[d] and rect.hi[d] > qrt.lo[d]
#     if ncmp:
#         return Interaction.INTERACT, _

#     # no overlap
#     return Interaction.NO_INTERACT, _


# def playground(is_hilbert):


#     # from polytope import Convex
#     # from polytope import halfSpace
#     # from polytope import point

#     # conv = Convex.convex(2)

#     # if False:
#     #     from math import pi, cos, sin
#     #     circle = 2*pi
#     #     how_many = 50
#     #     increments = circle / how_many
#     #     print(increments)

#     #     normals = []
#     #     for i in range(how_many):
#     #         x, y = cos(i * increments), sin(i * increments)
#     #         normals.append([x, y])

#     #     radius = 100
#     #     for i, normal in enumerate(normals):
#     #         hs = halfSpace.HalfSpace(dim=2, params=[radius] + normal)
#     #         conv.addHalfSpace(hs)
#     # else:
#     #     hs = halfSpace.HalfSpace(dim=2, points=[point.Point(pointCoords=[0., 10.]), point.Point(pointCoords=[1000., 25.])])
#     #     conv.addHalfSpace(hs)
#     #     hs = halfSpace.HalfSpace(dim=2, points=[point.Point(pointCoords=[1000., 0.]), point.Point(pointCoords=[0., 25.])])
#     #     conv.addHalfSpace(hs)


#     # from hplane import HyperPlane
#     # print(conv)
#     # query_planes = []
#     # for hs in conv.halfSpaces:
#     #     query_planes.append(HyperPlane(hs.aVal[1:], hs.aVal[0]))
#     # for plane in query_planes:
#     #     print(plane)
#     # nodes, edges = Convex.as_nodes_and_edges(conv)
#     # print("**** ", len(nodes), len(edges), len(query_planes))
#     # for node in nodes:
#     #     print("node", node)

#     # def as_linestring(start, end):
#     #     return f"LINESTRING({start[0]} {start[1]}, {end[0]} {end[1]})"

#     # # translate the convex polytope
#     # with open("/tmp/ch_edges_2d.csv", "w") as fh:
#     #     print("wkt", file=fh)
#     #     from pysfc.vectorops import add
#     #     for edge in edges:
#     #         # print("edge", edge)
#     #         print(as_linestring(add(nodes[edge[0]], (5000, 5000)), add(nodes[edge[1]], (5000, 5000))), file=fh)# , edge)

#     # translated_planes = []
#     # for query_plane in query_planes:
#     #     new_pt = add(query_plane.through, (5000, 5000))
#     #     print(new_pt)
#     #     # hp = HyperPlane.from_normal_and_distance(query_plane.n, norm(new_pt))
#     #     hp = HyperPlane.from_normal_and_point(query_plane.n, new_pt)
#     #     translated_planes.append(hp)
#     #     print(hp)
#     # query_planes = translated_planes[:]
#     # conv = Convex.convex(2)
#     # with open("/tmp/planes_2d.csv", "w") as fh:
#     #     for hs in query_planes:
#     #         print(hs)
#     #         from hplane import visualize_halfplane_2d
#     #         visualize_halfplane_2d(hs, size=80000, filenm="/tmp/planes_2d.csv")
#     #         hs_ = halfSpace.HalfSpace(dim=2, params=[hs.d] + list(hs.n))
#     #         conv.addHalfSpace(hs_)


#     # return

#     # for i in range(16):
#     #     print(format4(i,""),  format4(bin(i)[2:]), "<==>", format4(bitreverse(i, 4),""), format4(bin(bitreverse(i, 4))[2:]))
#     mbits = 25
#     # query_box = ndbox([1, 1, 1, 1], [5, 3301, 5, 5])
#     #query_box = ndbox([1, 1], [4, 4])

#     # query_box = ndbox([100, 100], [104, 104])
#     # ndims = len(query_box.lo)
#     # query_box = ndbox([0,0,3, 3], [0,0,4, 4])
#     ndims = 2

#     # query_planes = as_hyperplanes(query_box)
#     # FIXME: rename/refactor to ndbox.as_convex_polytope() (and have type that represents ndConvexPolytope)
#     # FIXME: maybe even better is to have also HalfSpace class, as HyperPlane is just the plane,
#     #        a HalfSpace makes a distinction which side of the HyperPlane is considered Interior (part of HalfSpace) and which not
#     #        the convexPolytope is then a list of HalfSpaces (for which we can have .on_interior(ndpoint) and .on_exterior(ndpoint) methods)
#     # for plane in query_planes:
#     #     plane.normalize()

#     if is_hilbert:
#         enter, exit = _hinitial_start_end(mbits, ndims)
#         root = ndvoxel_hilbert([0] * ndims, level=0, mbits=mbits, enter=enter, exit=exit)
#     else:
#         root = ndvoxel([0] * ndims, level=0, mbits=mbits)

#     max_level = mbits # root.mbits # - 2
#     max_level = 20

#     # root = ndvoxel([0, 0], level=0, mbits=8)
#     # FIXME: which dim of address corresponds with which dimension?
#     # is it:
#     # x, y, [z, ...]            <=> Z-order morton curve
#     # or is it:
#     # [... , z], y, x           <=> N-order morton

#     # convex = ConvexPolytope.from_ndbox(ndbox([-1]*ndims, [2**mbits + 1]*ndims))
#     # hs = HyperPlane.from_points([[0., 0.], [0., 1000.]])
#     # convex.add_halfplane_to_convex(hs)
#     # hs = HyperPlane.from_points([[0.,1000.],[1000., 25.]])
#     # convex.add_halfplane_to_convex(hs)
#     # print(convex.show()
#     # )

#     # convex = ConvexPolytope.from_ndbox(ndbox([-10, -10], [12500, 7]))
#     # hs = HyperPlane.from_points([[24., 6.8], [125.,0.]])
#     # convex.add_halfplane_to_convex(hs)

#     convex = ConvexPolytope.from_ndbox(ndbox([-1]*ndims, [2**mbits + 1]*ndims))

#     pts = []
#     ## -- circle
#     center = 2**14
#     r = 2**8
#     if False:
#         from math import pi, cos, sin
#         circle = 2*pi
#         how_many = 40
#         increments = circle / how_many

#         print(increments)
#         normals = []
#         for i in range(how_many + 1):
#             x, y = cos(i * increments), sin(i * increments)
#             normals.append([x, y])
#             x, y = center + r * x, center + r * y
#             pts.append([x,y])
#     else:
#         # # # diamond
#         x, y = center + r, center
#         pts.append([x,y])
#         x, y = center, center + 0.1 * r
#         pts.append([x,y])
#         x, y = center - r, center
#         pts.append([x,y])
#         x, y = center, center - 0.1 * r
#         pts.append([x,y])
#         # make ring
#         pts = pts + [pts[0]]
#     # pts in cw order
#     pts.reverse()
#     for start, end in zip(pts, pts[1:]):
#         hs = HyperPlane.from_points([start, end])

#     # interesting link
#     # https://github.com/kovacsv/JSModeler/blob/master/src/extras/solidgenerator.js


#     #     convex.add_halfplane_to_convex(hs)
#     # # radius = 10000
#     # for i, (normal, pt) in enumerate(zip(normals, pts)):
#     #     #hs = HyperPlane(normal, -radius)
#     #     hs = HyperPlane.from_normal_and_point(normal, pt)

#         convex.add_halfplane_to_convex(hs)
#     print(convex.show())


def _make_convex_test_object(ndims, mbits):
    convex = ConvexPolytope.from_ndbox(ndbox([-1] * ndims, [2 ** mbits + 1] * ndims))

    pts = []
    ## -- circle
    center = 2 ** 14
    r = 2 ** 8
    if True:
        # # # circle object
        from math import pi, cos, sin

        circle = 2 * pi
        how_many = 40
        increments = circle / how_many

        print(increments)
        normals = []
        for i in range(how_many + 1):
            x, y = cos(i * increments), sin(i * increments)
            normals.append([x, y])
            x, y = center + r * x, center + r * y
            pts.append([x, y])
    else:
        # # # diamond
        x, y = center + r, center
        pts.append([x, y])
        x, y = center, center + 0.1 * r
        pts.append([x, y])
        x, y = center - r, center
        pts.append([x, y])
        x, y = center, center - 0.1 * r
        pts.append([x, y])
        # make ring
        pts = pts + [pts[0]]
    # pts in cw order
    pts.reverse()
    for start, end in zip(pts, pts[1:]):
        hp = HyperPlane.from_points([start, end])
        convex.add_halfplane_to_convex(hp)

    return convex


def range_generation(query_planes_or_box, ndims, mbits, max_level=10, is_hilbert=False, verbose=False):
    #     query_planes = convex.hyperplanes()
    """Perform a traversal of the n-ary tree
    checking whether a index node is inside the query_planes
    """
    if is_hilbert:
        enter, exit = _hinitial_start_end(mbits, ndims)
        root = ndvoxel_hilbert(
            [0] * ndims, level=0, mbits=mbits, enter=enter, exit=exit
        )
    else:
        root = ndvoxel([0] * ndims, level=0, mbits=mbits)
    
    if verbose:
        print("+-- root ", root, "--+")
        if is_hilbert:
            print("HILBERT curve")
        else:
            print("MORTON N-curve")
        print(" ndims", root.ndims)
        print(" mbits", root.mbits)
        print(" branching factor", 2 ** root.ndims)
        print(" root", root.as_ndbox())
        print(root)
        print(" querying with:")
        for hp in query_planes_or_box:
            print("  ", hp)

    # the fastest nd-test we came up with, remembering if hyperplane contains ndbox of parent already, then no need to test again
    if isinstance(query_planes_or_box, ndbox):
        relate=relate_box_box
    else:
        relate = relate_planes_box__sweep_based__with_mask
    # relate=relate_ndbox

    result: list[tuple(ndvoxel | ndvoxel_hilbert, int)] = []
    stack = [(root, 0)]
    while stack:
        node, remembering_mask = stack.pop()
        ndcmp, remembering_mask = relate(
            query_planes_or_box, node.as_ndbox(), remembering_mask
        )
        # print(ndcmp, node.as_ndbox())
        # ndcmp, remembering_mask = relate(query_box, node.as_ndbox(), remembering_mask)
        # ndcmp = relate_ndbox(query_box, node.as_ndbox())
        # print(ndcmp, remembering_mask, query_box, node.as_ndbox())
        # print(ndcmp)
        # print(f"{node.level * ' '}at: {node} {ndcmp}")
        # -- full overlap := accept
        if ndcmp == Interaction.CONTAINED:
            # pass
            # print(f"{node.level * ' '}full @ level = {node.level}")
            result.append((node, ndcmp, remembering_mask))
            # print("")
        # -- partial overlap := refine or accept
        elif ndcmp == Interaction.INTERACT:
            # print(f"{node.level * ' '}partial @ level = {node.level}")
            # print("")
            if node.level < max_level:  # refine
                # get child nodes in correct visiting order
                if is_hilbert:
                    childs = node.get_childs_along_horder(node.enter, node.exit)
                else:
                    childs = node.get_childs_along_norder()
                # schedule childs for visiting them
                for child in childs:
                    stack.append((child, remembering_mask))
            else:  # accept
                # print('we accept due to max_level')
                result.append((node, ndcmp, remembering_mask))
                # FIXME: as this can be partial intersecting node, maybe we should also store the
                # 'remembering_mask' value here
                # --> partial overlapping needs additional testing of points
                # --> this way we do not need to address all hyperplanes,
                #     but only relevant hyperplanes of node

        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp == Interaction.NO_INTERACT:
            # print(f"{node.level * ' '}disjoint @ level = {node.level}")
            # print("")
            continue
        # -- error (unhandled relate case)
        else:
            raise ValueError("unknown relate outcome")

    if verbose:
        print("--- result of range generation ---")
        print(len(result))
    # node: ndvoxel | ndvoxel_hilbert
    # interaction: int
    # for node, interaction in reversed(result):
    #     if is_hilbert:
    #         print(interaction, node.as_ndbox(), node.as_hilbert_range())
    #     else:
    #         print(interaction, node.as_ndbox(), node.as_morton_range())
    
    # with open("/tmp/query_planes_2d.csv", "w") as fh:
    #     print("geometry;what", file=fh)
    #     for hs in query_planes:
    #         visualize_halfplane_2d(hs, fh, size=2**mbits)

    # with open("/tmp/query_result.csv", "w") as fh:
    #     print("interaction;geometry;sfc_start;sfc_end;mask", file=fh)
    #     for node, interaction, mask in reversed(result):
    #         if is_hilbert:
    #             range_gen = node.as_hilbert_range
    #         else:
    #             range_gen = node.as_morton_range
    #         line = "{0};{1};{2[0]};{2[1]};{3}".format(interaction, node.as_ndbox().as_wkt_2d(), range_gen(), mask)
    #         print(line, file=fh)
    return result


# def play_intersection():
#     import numpy as np
#     import polyscope as ps

#     mbits = 25
#     ndims = 2

#     convex = ConvexPolytope.from_ndbox(ndbox([-(2**mbits +1)]*ndims, [2**(mbits + 1)]*ndims))
#     # hs = HyperPlane.from_points([[0., 10.], [1000., 25.]])
#     # convex.add_halfplane_to_convex(hs)
#     # hs = HyperPlane.from_points([[1000., 0.],[0., 25.]])
#     # convex.add_halfplane_to_convex(hs)
#     nodes, edges = convex.as_nodes_and_edges()
#     ps.register_curve_network("the cube", np.array(nodes), np.array(edges))


#     root = ndvoxel([0] * ndims, level=0, mbits=mbits)
#     index_box = root.as_ndbox()
#     index = ConvexPolytope.from_ndbox(index_box)
#     nodes, edges = index.as_nodes_and_edges()
#     ps.register_curve_network("the index root", np.array(nodes), np.array(edges))


#     result = relate_box__sweep_based__with_mask(convex.hyperplanes(), index_box, 0)
#     print(result)


#     ps.set_up_dir("z_up")
#     ps.init()
#     ps.show()


# def visualize():
#     convex, boxes = playground(is_hilbert=True)
#     import numpy as np
#     import polyscope as ps
#     nodes, edges = convex.as_nodes_and_edges()
#     ps.register_curve_network("the query planes", np.array(nodes), np.array(edges))
#     for box in boxes:
#         print(box)
#         indexnode, interact, mask = box
#         print(indexnode)
#         print(indexnode.as_ndbox())
#     ps.set_up_dir("z_up")
#     ps.init()
#     ps.show()


def _playground(is_hilbert):
    import time
    ndims = 2
    mbits = 25
    max_level = 25

    timings = [(time.perf_counter_ns(), 'start')]
    convex = _make_convex_test_object(ndims, mbits)
    query_planes = convex.hyperplanes()
    timings.append((time.perf_counter_ns(), 'convex'))
    result = range_generation(query_planes, ndims, mbits, max_level, is_hilbert)
    timings.append((time.perf_counter_ns(), 'range_gen'))

    for first, second in zip(timings, timings[1:]):
        print("{:10.3f}".format((second[0] - first[0]) / 1e6), "ms  ", second[1])

    with open("/tmp/query_planes_2d.csv", "w") as fh:
        print("geometry;what", file=fh)
        for hs in query_planes:
            _visualize_hyperplane_2d(hs, fh, size=2 ** mbits)

    with open("/tmp/query_result.csv", "w") as fh:
        print("interaction;geometry;sfc_start;sfc_end;mask", file=fh)
        for node, interaction, mask in reversed(result):
            if is_hilbert:
                range_gen = node.as_hilbert_range
            else:
                range_gen = node.as_morton_range
            line = "{0};{1};{2[0]};{2[1]};{3}".format(
                interaction, node.as_ndbox().as_wkt_2d(), range_gen(), mask
            )
            print(line, file=fh)

def nquery(box):
    return _query(box, False)

def hquery(box):
    return _query(box, True)

def _query(mybox, is_hilbert):
    # mybox = ndbox([1, 0, 0, 0], [2, 1, 1, 1])
    if isinstance(mybox, ndbox):
        ndims = mybox.dims
    elif isinstance(mybox, list):
        # assume it is a list of hyperplanes
        ndims = len(mybox[0].n)
    # get the number of needed bits, by looking at the maximum value of the coordinate
    # this should be done differently when we have a dataset fixed, then we should
    # look at the metadata, and use the number of bits,
    # also because we then have a nd histogram, which records at which level how many
    # points are present, and this could otherwise go 'skewed'
    # mbits = _determine_bits(int(ceil(max(mybox.hi))), 1) 
    # mybox=ndbox(mybox.lo, [h - 1 for h in mybox.hi])
    mbits = 25
    # mbits = 25 ## could generate maybe this from the dimensions of the box...
    max_level = mbits
    # query_planes = mybox.as_hyperplanes()
    #query_planes = convex.hyperplanes()
    #timings.append((time.perf_counter_ns(), 'convex'))
    result = range_generation(mybox, ndims, mbits, max_level, is_hilbert)
    ranges = []
    for node, interaction, mask in reversed(result):
        if is_hilbert:
            range_gen = node.as_hilbert_range
        else:
            range_gen = node.as_morton_range
        ranges.append(range_gen())
    return ranges

if __name__ == "__main__":
    # _playground(is_hilbert=True)
    mybox = ndbox([1, 0, 0, 0], [2, 1, 1, 1])
    ranges = nquery(mybox) 
    for rng in ranges:
        print(rng)
    ranges = hquery(mybox)
    for rng in ranges:
        print(rng)
#     # visualize()
#     # playground(is_hilbert=False)

#     # play_intersection()
