from math import isclose
from enum import Enum

from pysfc.linalg import mul, dot, unit, sub, add
from pysfc.linalg import determinant

import math


class HyperPlane:
    def __init__(self, n, d):
        """
        This class represents an hyperplane as the zero set of the implicit equation:
            n⋅x + d = 0
        where:
            n   is a unit normal vector of the plane (linear part), as list of floats
            d   is the distance (offset) to the origin, as float

        Any hyperplane can be written as the set of points x satisfying
            n⋅x = -d
            =>
            n⋅x + d = 0,

        We keep it in Hesse normal form, so:
            . n is a unit vector
            . the sign of d determines the side of the plane on which the origin is located
                when d > 0 it is in the half space determined by the normal n
                when d < 0 it is in the other half space
            . d is the distance from the origin

        This way the two halfspaces that go with this hyperplane are also defined:
            . normal points towards = exterior (upper halfspace)
            . normal points away from = interior (lower halfspace)

        """
        assert isclose(dot(n, n), 1.0)  # vec is normalized / unit vec
        self.n = tuple(map(float, n))  # normal vector, should be normalized!
        self.d = float(d)

    def __repr__(self):
        return "HyperPlane({}, {})".format(repr(self.n), repr(self.d))

    @staticmethod
    def from_normal_and_point(normal, point):
        """create a hyperplane from unit vector and point through which the plane goes"""
        return HyperPlane(normal, -dot(normal, point))

    # @staticmethod
    # def from_normal_and_distance(normal, distance):
    #     """constructs a plane from its normal n and distance to the origin d"""
    #     return HyperPlane(normal, distance)

    # @property
    # def offset(self):
    #     return self.d

    # @property
    # def normal(self):
    #     return self.n

    def signed_distance(self, point):
        """Signed distance of the point to the plane"""
        # we can save dividing by the length of the vector n, as n is unit vec and always has length == 1.0
        return dot(self.n, point) + self.d  # / norm(self.n)

    def projection(self, point):
        """Get projection of point onto plane"""
        return sub(point, mul(self.n, self.signed_distance(point)))

    @property
    def through(self):
        """
        Returns a nD-point through which this hyperplane passes
        """
        return mul(self.n, -self.d)

    @staticmethod
    def from_points(points):
        """
        construct a hyperplane from nd-points, the order in which the points are given matters
        (dependent on the order, the distance to the origin will be different in sign and thus the orientation in space)
        """
        # FIXME:
        # - rename aVal
        # - do not insert(0, dist), but keep separate
        dim = len(points[0])
        # self.dim = points[0].dim
        # if dim != None:
        #     self.dim = dim
        # assert self.dim >= points[0].dim, (
        #     "Points " + str(points[0].dim) + " dim " + str(self.dim)
        # )
        normal = []
        for i in range(dim):  # for each column
            mtx = []
            for pt in points:  # for each row(point) leave out column i
                row = pt[:i]
                row.append(1)
                row += pt[i + 1 :]
                mtx.append(row)
            val = determinant(mtx)
            # print(mtx, val)
            normal.append(-val)
        for i in range(dim):  # for each column
            mtx = []
            for pt in points:  # for each row(point)
                row = pt[:]
                mtx.append(row)
        dist = determinant(mtx)
        dd = 0
        point0 = points[0]
        for i in range(dim):
            dd += -normal[i] * point0[i]
        ## compare d and dd
        if not isclose(dist, dd):
            raise ValueError("dist " + str(dist) + " <> dd " + str(dd))
        ## normalize, so that we obtain a unit vector
        sq = 0
        for aa in normal:
            sq += aa ** 2
        # print(sq)
        # print(dot(normal, normal))
        sq = math.sqrt(sq)
        if sq > 0.0:
            dist /= sq
            for i in range(dim):
                normal[i] /= sq

        # FIXME: add extra argument here for the ambient_dimension?
        # deal with the case when we have less dimensions in point than in (ambient) halfspace
        # but then, it would be nice to specify which of the dimensions is 0...
        # for i in range(dim, dim):
        #     normal.append(0.0)
        if not any(normal):
            raise ValueError("all zeros for the normal")
        return HyperPlane(unit(normal[:]), dist)


def _visualize_hyperplane_2d(h, fh, size=1000.0):
    from pysfc.linalg import rotate90ccw, add, mul

    # make 2D line of 2000 units long + the normal
    ccw = rotate90ccw(h.n)
    # cw = rotate90cw(h.n) # direction vector on the line

    start = add(mul(ccw, size), h.through)
    end = add(mul(ccw, -size), h.through)

    fh.write("LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tplane".format(start, end))
    fh.write("\n")
    fh.write(
        "LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tnormal".format(
            h.through, add(h.through, h.n)
        )
    )
    fh.write("\n")

# def test_from_points():
#     print("*")
#     points = [
#         [1000.0, -1.0],
#         [1000.0, +1.0],
#     ]
#     hp = HyperPlane.from_points(points)
#     for pt in points:
#         print(hp.signed_distance(pt))
#         print(hp.through)

#     print("*")
#     points = [
#         [1000.0, +1.0],
#         [1000.0, -1.0],
#     ]
#     hp = HyperPlane.from_points(points)
#     for pt in points:
#         print(hp.signed_distance(pt))
#         print(hp.through)

#     print("*")
#     points = [
#         [1.0, 100.0],
#         [-1.0, 100.0],
#     ]
#     hp = HyperPlane.from_points(points)
#     for pt in points:
#         print(hp.signed_distance(pt))
#         print(hp.through)

#     print("*")
#     points = [
#         [-1.0, 100.0],
#         [1.0, 100.0],
#     ]
#     hp = HyperPlane.from_points(points)
#     for pt in points:
#         print(hp.signed_distance(pt))
#         print(hp.through)

#     points = [
#         [4.0, 1.0],
#         [2.0, 2.0],
#     ]
#     hp = HyperPlane.from_points(points)
#     for pt in points:
#         print(hp.signed_distance(pt))
#         print(hp.through)
#         print(hp.signed_distance(hp.through))


# # def normalize(n, d):
# #     """normalize plane equation, after this n is a unit vector and d represents the correct distance from the origin"""
# #     normalized_d = d / norm(n)
# #     normalized_n = unit(n)
# #     return normalized_n, normalized_d


# def playground():
#     plane = HyperPlane.from_normal_and_point(normal=(1, 0), point=(5, 0))
#     print(plane)
#     print(plane.signed_distance((0, 0)))
#     print(plane.through)

#     print("***")

#     plane = HyperPlane.from_normal_and_point(normal=(1, 0), point=(-5, 0))
#     print(plane)
#     print(plane.signed_distance((0, 0)))
#     print(plane.through)

#     print("***")

#     plane = HyperPlane.from_normal_and_point((-1, 0), (5, 0))
#     print(plane)
#     print(plane.through)
#     print(plane.signed_distance((0, 0)))
#     print(plane.projection((0, 7)))
#     print(plane.projection((15, 7)))
#     print(plane.projection((5, 0)))

#     print("***")
#     HyperPlane.from_normal_and_distance(plane.n, plane.d)
#     print(plane)
#     print(plane.through)

#     from random import random

#     vec = unit([random() for _ in range(4)])
#     pt = [random() for _ in range(4)]
#     print(pt)
#     plane = HyperPlane.from_normal_and_point(vec, pt)
#     print(plane.signed_distance(pt))
#     assert isclose(plane.signed_distance(plane.through), 0.0, abs_tol=1e-15)

#     print(plane.through)
#     print(sub(pt, plane.through))


# class Interaction(Enum):
#     NO_INTERACT = 0  # no interaction
#     INTERACT = 1  # possibly some interaction
#     CONTAINED = 2  # interaction, for sure


# def relate_box__sweep_based__with_mask(planes, box, mask=0, stats=None):
#     """
#     spatial relation between hyperplanes and nd-box
#     remembering outcome of relation test in mask so that
#     some hyperplanes might be skipped for testing
#     when nd-box is split in 2^n smaller boxes
#     """
#     # global TESTED_CT, CONTAINED_CT
#     from math import sqrt

#     # pre-condition: some planes to test against
#     assert len(planes) > 0
#     result = -1
#     # another pre-condition:
#     # box has same edge lengths in all dimensions
#     # (not tested)
#     # side length of box: E
#     box_side_dist = box.hi[0] - box.lo[0]
#     # length of diagonal: E * √n (:= √(n * E * E), with n=nDims
#     diagonal_dist = box_side_dist * sqrt(box.dims)
#     # mask_out will be changed, if we detect that box now is contained
#     mask_out = mask
#     for index, plane in enumerate(planes):
#         # print(f'testing {index} {plane} result so far: {result}')
#         if stats is not None:
#             stats.tested += 1
#         # TESTED_CT += 1
#         # we know that the plane is contained from parent index boxes,
#         # skip checking this particular hyperplane again
#         is_contained = bool(mask & (1 << index))
#         # print(f'{is_contained}')
#         if is_contained:
#             if stats is not None:
#                 stats.contained += 1
#             # CONTAINED_CT += 1
#             if result == -1:
#                 # print('l280')
#                 result = Interaction.CONTAINED
#             continue
#         # first hit point sweeping the box with the hyperplane
#         # FIXME: for this logic to be correct, the
#         # point that is tested first
#         # needs to be the corner point that is
#         # 'least-likely-to-be-inside'
#         # along the direction of the normal
#         #
#         # then we can add the two distances (box size / diagonal size)
#         # and eventually test the point that is 'most-likely-to-be-inside'
#         # enter = sweep_box_with_plane_enter(plane, box)
#         # enter_dist = plane.signed_distance(enter)

#         exit = sweep_box_with_plane_exit(plane, box)
#         exit_dist = plane.signed_distance(exit)


#         # print(plane, 'to', enter, 'with', enter_dist)
#         if exit_dist <= 0:
#             mask_out += pow(2, index)
#             # only override result when not yet seen
#             if result == -1:
#                 # print('l289')
#                 result = Interaction.CONTAINED
#         else:
#             # perform distance comparison first (diagonal length & side of box)
#             # by doing that, we defer sweeping for the far point
#             abs_enter_dist = abs(exit_dist)
#             if abs_enter_dist > diagonal_dist:
#                 # print('enter > diag')
#                 return Interaction.NO_INTERACT, mask_out
#             elif abs_enter_dist < box_side_dist:
#                 # print('enter < side')
#                 result = Interaction.INTERACT
#             else:
#                 # get last hit point sweeping the box with the hyperplane
#                 enter = sweep_box_with_plane_enter(plane, box)
#                 enter_dist = plane.signed_distance(enter)
#                 if enter_dist <= 0:
#                     result = Interaction.INTERACT
#                 else:
#                     return Interaction.NO_INTERACT, mask_out
#     # post-condition
#     # result here should be either INTERACT / CONTAINED
#     assert result in (Interaction.INTERACT, Interaction.CONTAINED)
#     return result, mask_out


# # fixme: rename to caliper touch points?
# def sweep_box_with_plane_enter(hyperplane, box):
#     """Return enter point when 'sweeping' the box with the hyperplane

#     enter: point that is 'hit' first while sweeping the box with the hyperplane
#     along the direction of the normal of the hyperplane
#     (i.e. *first* encounter point along ray swept by hyperplane)
#     """
#     enter = [0.0] * box.dims
#     for i, n in enumerate(hyperplane.n):
#         # if n >= 0:
#         if n < 0:
#             enter[i] = box.hi[i]
#         else:
#             enter[i] = box.lo[i]
#     return enter


# def sweep_box_with_plane_exit(hyperplane, box):
#     """Return exit point 'sweeping' the box with the hyperplane

#     exit: point where hyperplane leaves the box while sweeping the box with the
#     hyperplane along the direction of the normal of the hyperplane
#     (i.e. *second* encounter point along ray swept by hyperplane)
#     """
#     exit = [0.0] * box.dims
#     for i, n in enumerate(hyperplane.n):
#         # if n >= 0:
#         if n < 0:
#             exit[i] = box.lo[i]
#         else:
#             exit[i] = box.hi[i]
#     return exit


# # from pysfc.ndgeom import ndbox


# def visualize_halfplane_2d(h, fh, size=1000.0):
#     from pysfc.vectorops import rotate90ccw, add, mul

#     # make 2D line of 2000 units long + the normal
#     ccw = rotate90ccw(h.n)
#     # cw = rotate90cw(h.n) # direction vector on the line

#     start = add(mul(ccw, size), h.through)
#     end = add(mul(ccw, -size), h.through)

#     fh.write("LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tplane".format(start, end))
#     fh.write("\n")
#     fh.write(
#         "LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tnormal".format(
#             h.through, add(h.through, h.n)
#         )
#     )
#     fh.write("\n")


# def as_hyperplanes(ndbox):
#     """
#     Given a ndbox, return a list of HyperPlanes that bound the same space
#     as the box
#     """
#     planes = []
#     dims = [0.0] * len(ndbox.lo)
#     # lo dimensions, normal negative, towards origin (assuming lo > 0)

#     # FIXME: does this also work with the box on both sides of the origin
#     # (i.e. having negative coordinates)
#     for dim, _ in enumerate(ndbox.lo):
#         w = dims[:]
#         w[dim] = +1.0
#         planes.append(HyperPlane.from_normal_and_point(w, ndbox.lo))
#     # hi dimensions, negative normal
#     for dim, _ in enumerate(ndbox.hi):
#         w = dims[:]
#         w[dim] = -1.0
#         planes.append(HyperPlane.from_normal_and_point(w, ndbox.hi))
#     return planes


# def project_3d_to_2d():
#     # planes = [
#     #     # x [-1.5 , + 1.5]
#     #     HyperPlane([1.0, 0.0, 0.0], -1.5),
#     #     HyperPlane([-1.0, 0.0, 0.0], -1.5),
#     #     # y [-3.5 , + 3.5]
#     #     HyperPlane([0.0, 1.0, 0.0], -3.5),
#     #     HyperPlane([0.0, -1.0, 0.0], -3.5),
#     #     # z [-4.5 , + 2.5]
#     #     HyperPlane([0.0, 0.0, 1.0], -4.5),
#     #     #HyperPlane([0.0, 0.0, -1.0], -2.5),  # axis aligned plane
#     #     #HyperPlane([-1.0, -1.0, -1.0], -2.5), # tilted plane
#     # ]
#     planes = as_hyperplanes(ndbox([-5, -20, 30], [15, 25, 35]))
#     from itertools import combinations

#     ndims = len(planes[0].n)
#     combis = tuple(combinations(range(ndims), 2))
#     names = "xyzuvwklmnop"
#     filenm = "/tmp/wopwop"
#     origin = (0, 0)
#     for combi in combis:
#         side = names[combi[0]], names[combi[1]]
#         with open(filenm + str(side) + ".txt", "w") as fh:
#             pass
#         print("###", side)

#         for plane in planes:
#             print("")
#             print("original", plane)
#             through = plane.through
#             print(" through", through)
#             new_through = []
#             new_w = []
#             for dim in combi:
#                 new_w.append(plane.n[dim])
#                 new_through.append(through[dim])
#             print(new_w, new_through)
#             # print(dot(new_w, new_w))
#             # FIXME: apparently I am doing something wrong...
#             # as the lines do not end up with the original sides of the ndbox :|
#             if dot(new_w, new_w) > 0.0:
#                 # v = make_vector(origin, new_through)
#                 # new_d = norm(proj(v, new_w))
#                 # norm( proj new_n onto -new_through )
#                 new_plane = HyperPlane.from_normal_and_point(new_w, new_through)
#                 print("new equation", new_plane, new_plane.through)
#                 # print(combi, new_plane)
#                 visualize_halfplane_2d(new_plane, 100.0, filenm + str(side) + ".txt")
#         print("")


# def play__construct_from_points():
#     hplanes = as_hyperplanes(ndbox([-10, -10], [10, 10]))
#     print(hplanes)

#     hplanes = as_hyperplanes(ndbox([5, 5], [10, 10]))  # TODO: as_convex()
#     print(hplanes)
#     for hp in hplanes:
#         print(hp.through, "->", hp.signed_distance(hp.through))

#     # print(list(ndbox([-10,-10], [10, 10]).vertices))

#     print("--")
#     c = [(-10, -10), (10, -10), (10, 10), (-10, 10), (-10, -10)]
#     c = list(map(list, c))
#     for start, end, axis in zip(c, c[1:], ["-y", "+x", "+y", "-x"]):
#         hp = HyperPlane.from_points([start, end])
#         print(axis, hp)
#         print(hp.through, "  ↹ dist:", hp.signed_distance(hp.through))
#         print("                ->", start, "↹ dist:", hp.signed_distance(start))
#         print("                ->", end, "↹ dist:", hp.signed_distance(end))

#     print("--")
#     c = [(5, 5), (10, 5), (10, 10), (5, 10), (5, 5)]
#     c = list(map(list, c))
#     for start, end, axis in zip(c, c[1:], ["-y", "+x", "+y", "-x"]):
#         hp = HyperPlane.from_points([start, end])
#         print(axis, hp)


# def test_diamant():
#     print("--")
#     c = [(0, -10), (10, 0), (0, 10), (-10, 0), (0, -10)]
#     c = list(map(list, c))

#     ll = ndbox((-20, -20), (-10, -10))
#     for start, end, axis in zip(c, c[1:], ["-y", "+x", "+y", "-x"]):
#         hp = HyperPlane.from_points([start, end])
#         print(axis, hp)
#         print(hp.through, "  ↹ dist:", hp.signed_distance(hp.through))
#         print("                ->", start, "↹ dist:", hp.signed_distance(start))
#         print("                ->", end, "↹ dist:", hp.signed_distance(end))

#         print(sweep_box_with_plane_enter(hp, ll))
#         print(sweep_box_with_plane_exit(hp, ll))


# def test_what_is_inside():
#     # d, n = 10.0, [-1.0, 0.0]
#     # 1000.0, [1.0, 0.0]
#     # 1000.0, [0.0, -1.0]
#     # 1000.0, [0.0, 1.0]
#     # hp = HyperPlane.from_normal_and_distance(n,d)

#     # hp = HyperPlane.from_normal_and_point((0.5*math.sqrt(2), 0.5*math.sqrt(2)), (5,5))

#     hp = HyperPlane.from_points(([0.0, 10.0 - 50], [1000.0, 25.0 - 50]))

#     # print()
#     with open("/tmp/points.txt", "w") as fh:
#         print("x; y; outcome_nx-p ; outcome_nx_c", file=fh)
#         n = hp.n
#         p = hp.through
#         for x in range(-50, 50):
#             for y in range(-50, 50):
#                 pt = (x, y)
#                 result1 = dot(n, add(pt, mul(p, -1))) <= 0.0
#                 # result2 = dot(n, p) <= hp.d
#                 result2 = hp.signed_distance(pt) <= 0.0
#                 print(x, ";", y, ";", result1, ";", result2, file=fh)

#     with open("/tmp/hp.txt", "w") as fh:
#         print("geometry\twhat", file=fh)
#     visualize_halfplane_2d(hp, size=1000.0, filenm="/tmp/hp.txt")


# def vis_halfplane_2d(fh, hp, size=1000.0, normal_factor=1.0):
#     from pysfc.vectorops import rotate90ccw, add, mul

#     # make 2D line of 2000 units long + the normal
#     ccw = rotate90ccw(hp.n)
#     # cw = rotate90cw(h.n) # direction vector on the line

#     start = add(mul(ccw, size), hp.through)
#     end = add(mul(ccw, -size), hp.through)
#     print(
#         "LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tplane".format(start, end), file=fh
#     )
#     print(
#         "LINESTRING({0[0]} {0[1]} , {1[0]} {1[1]})\tnormal".format(
#             hp.through, add(hp.through, mul(hp.n, normal_factor))
#         ),
#         file=fh,
#     )


# # class Convex:
# #     def __init__(self, v, p, v2e, e2v, combis):
# #         self.vertices_lut = v
# #         self.planes_lut = p
# #         # the topology structure
# #         # vertex -> edge
# #         self.vertex2edge_lut = v2e
# #         # edge -> 2 vertices
# #         self.edge2vertex_lut = e2v
# #         # the edges (hyperplane combis that, when intersected form the edge)
# #         self.edges = combis



# # @profile
# def profile():
#     boxes = []
#     for d in range(2, 10):
#         box = ndbox([-5] * d, [+5] * d)
#         boxes.append(box)
#         box = ndbox([5] * d, [10] * d)
#         boxes.append(box)
#         box = ndbox([-10] * d, [-5] * d)
#         boxes.append(box)
#     for box in boxes:
#         # polytope =
#         ConvexPolytope.from_ndbox(box)
#         # polytope.is_consistent()


# def around_origin():
#     box = ndbox([-5, -5, -5], [+5, +5, +5])
#     hp0 = HyperPlane.from_normal_and_point([0, 0, -1], [0, 0, 10])  # make convex empty
#     hp1 = HyperPlane.from_normal_and_point([0, 0, -1], [0, 0, -10])  # ignore hyperplane
#     hp2 = HyperPlane.from_normal_and_point(
#         [0, 0, -1], [0, 0, 0]
#     )  # embed hyperplane in structure

#     for hp in [hp0, hp1, hp2]:
#         polytope = ConvexPolytope.from_ndbox(box)
#         polytope.add_halfplane_to_convex(hp)
#         assert polytope.is_consistent()


# def around_origin_2d():
#     box = ndbox([-5, -5], [+5, +5])
#     polytope = ConvexPolytope.from_ndbox(box)
#     assert polytope.is_consistent()
#     nodes, edges = polytope.as_nodes_and_edges()
#     assert nodes == [(-5.0, -5.0), (-5.0, 5.0), (5.0, -5.0), (5.0, 5.0)]
#     assert edges == [(0, 2), (0, 1), (1, 3), (2, 3)]


# def play_with_cone():
#     import numpy as np
#     import polyscope as ps

#     ps.set_up_dir("z_up")
#     ps.init()
#     conv = ConvexPolytope.from_ndbox(ndbox([-10000] * 3, [10000] * 3))
#     # output_conv(conv, "/tmp/halfedges_cube.csv")
#     nodes, edges = conv.as_nodes_and_edges()
#     ps.register_curve_network("my cube", np.array(nodes), np.array(edges))

#     # hp = HyperPlane.from_normal_and_point([0,0,-1],[0,0,0]) # embed hyperplane in structure
#     # conv.add_halfplane_to_convex(hp)
#     # nodes, edges = conv.as_nodes_and_edges()
#     # ps.register_curve_network("bottom added", np.array(nodes), np.array(edges))

#     from math import pi, cos, sin

#     circle = 2 * pi
#     radius = 50.0
#     how_many = 25
#     increments = circle / how_many
#     # print(increments)
#     normals = []
#     for i in range(how_many):
#         x, y = cos(i * increments), sin(i * increments)
#         normals.append(unit([-x, -y, 0.5]))
#     print(normals)
#     for n in normals:
#         hp = HyperPlane(n, radius)
#         # we can make a copy, to skip insertion of a hyperplane
#         # and rollback to before state
#         copy = conv.copy()
#         # try:
#         conv.add_halfplane_to_convex(hp)
#         # except KeyError:
#         #     conv = copy
#         # for i, pt in enumerate(pts):
#         #     hs = HalfSpace(dim=3, params=[radius] + pt)
#         #     hs.normalise()
#         #     conv.addHalfSpace(hs)
#         if not conv.is_consistent():
#             print("rolling back inserting of hp, as it leads to inconsistent polytope")
#             conv = copy
#     nodes, edges = conv.as_nodes_and_edges()
#     # conv.print()

#     conv.show()
#     assert conv.is_consistent()

#     ps.register_curve_network("circle added", np.array(nodes), np.array(edges))
#     ps.show()


# def play_with_frustum():
#     import numpy as np
#     import polyscope as ps

#     from random import randint

#     ps.set_up_dir("z_up")
#     ps.init()
#     conv = ConvexPolytope.from_ndbox(ndbox([-10000] * 3, [10000] * 3))
#     nodes, edges = conv.as_nodes_and_edges()
#     ps.register_curve_network("the cube", np.array(nodes), np.array(edges))
#     hp = HyperPlane.from_points([[1, 0, 0], [0, 10, 0], [-13.5, 13.5, 13.5]])  # RHR
#     conv.add_halfplane_to_convex(hp)

#     for i in range(5):
#         try:
#             hp = HyperPlane.from_points(
#                 [
#                     [-randint(0, 10), 0, 0],
#                     [0, randint(0, 5), 0],
#                     [randint(0, 5), randint(0, 5), randint(0, 5)],
#                 ]
#             )  # RHR
#             conv.add_halfplane_to_convex(hp)
#         except ValueError:
#             pass

#     nodes, edges = conv.as_nodes_and_edges()
#     ps.register_curve_network("cut", np.array(nodes), np.array(edges))

#     ps.show()


# def play_with_projection_to_lesser_dimension():
#     import numpy as np
#     import polyscope as ps

#     dims = 4
#     ps.set_up_dir("z_up")
#     ps.init()
#     conv = ConvexPolytope.from_ndbox(ndbox([-10000] * dims, [10000] * dims))
#     nodes, edges = conv.as_nodes_and_edges()
#     # which cols to delete
#     cols = [0, 1]
#     ps.register_curve_network(
#         "the cube", np.delete(np.array(nodes), cols, axis=1), np.array(edges)
#     )
#     ps.show()


# def play_with_deca():
#     import numpy as np
#     import polyscope as ps

#     dims = 3
#     conv = ConvexPolytope.from_ndbox(ndbox([-100] * dims, [100] * dims))
#     nodes, edges = conv.as_nodes_and_edges()
#     ps.register_curve_network("the cube", np.array(nodes), np.array(edges))

#     ps.register_curve_network(
#         "the cube 2d", np.delete(np.array(nodes), [2], axis=1), np.array(edges)
#     )

#     from math import pi, cos, sin

#     circle = 2 * pi
#     radius = 50.0
#     how_many = 5
#     deg30 = 2 * pi / 12
#     increments = circle / how_many
#     # print(increments)

#     normals_top = []
#     for i in range(how_many):
#         x, y = cos(i * increments + deg30), sin(i * increments + deg30)
#         normals_top.append(unit([x, y, +0.5]))

#     normals_bottom = []
#     for i in range(how_many):
#         x, y = cos(i * increments), sin(i * increments)
#         normals_bottom.append(unit([x, y, -0.5]))

#     pts = (mul(n, radius) for n in normals_top)
#     for n, pt in zip(normals_top, pts):
#         # hp = HyperPlane(n, radius)
#         hp = HyperPlane.from_normal_and_point(n, pt)
#         conv.add_halfplane_to_convex(hp)

#     pts = (mul(n, radius) for n in normals_bottom)
#     for n, pt in zip(normals_bottom, pts):
#         # hp = HyperPlane(n, radius)
#         hp = HyperPlane.from_normal_and_point(n, pt)
#         conv.add_halfplane_to_convex(hp)

#     hp = HyperPlane.from_normal_and_point((0, 0, -1), (0, 0, -40))
#     conv.add_halfplane_to_convex(hp)

#     hp = HyperPlane.from_normal_and_point((0, 0, +1), (0, 0, 40))
#     conv.add_halfplane_to_convex(hp)

#     nodes, edges = conv.as_nodes_and_edges()

#     ps.register_curve_network("pentagon", np.array(nodes), np.array(edges))
#     ps.register_curve_network(
#         "penta 2d xy", np.delete(np.array(nodes), [2], axis=1), np.array(edges)
#     )
#     ps.register_curve_network(
#         "penta 2d xz", np.delete(np.array(nodes), [1], axis=1), np.array(edges)
#     )
#     ps.register_curve_network(
#         "penta 2d yz", np.delete(np.array(nodes), [0], axis=1), np.array(edges)
#     )

#     ps.set_up_dir("z_up")
#     ps.init()
#     ps.show()


# def test_relate():
#     # this is clearly wrong, is there a missing 1 in the matrix of point coordinates?
#     # as Rod did add 1 into the matrix...
#     # maybe good to check the original code of Rod, and see what it produces
#     hp = HyperPlane.from_points([[0, 10], [10, 0]])
#     print(hp.through)

#     # this is fine, if we construct the plane, then we get a point through it
#     pt = [-5, -5]
#     hp = HyperPlane.from_normal_and_point(unit([1, -1]), pt)
#     for dim in range(len(pt)):
#         assert isclose(hp.signed_distance(hp.through), 0)
#     box1 = ndbox([0, 0], [1, 1])
#     box2 = ndbox([4.5, 4.5], [5.5, 5.5])
#     box3 = ndbox([10.0, 10.0], [11, 11])

#     print(hp)
#     print(hp.through)
#     print(box1)
#     # distance calculation to all vertices
#     # all distances <= 0: inside
#     # mix of signs: partly in/partly out
#     # all distances >0 for sure outside
#     for bid, box in enumerate([box1, box2, box3]):
#         for vertex in box.vertices:
#             print(bid, hp.signed_distance(vertex))
#         print("")

#     # for bid, box in enumerate([box1, box2, box3]):

#     # the entrance point is the origin (even though the plane is at 5,5, we would
#     # hit the origin first when we move over the domain from -inf to +inf with the
#     # plane)
#     enter_point = sweep_box_with_plane_enter(hp, box1)
#     print(enter_point)
#     assert enter_point == [1.0, 1.0]
#     exit_point = sweep_box_with_plane_exit(hp, box1)
#     assert exit_point == [0.0, 0.0]
#     interaction, mask = relate_box__sweep_based__with_mask([hp], box1, 0)
#     assert interaction == Interaction.CONTAINED
#     assert mask == 1
#     # unittest.main()


# def relate_box__box_based(planes, box):
#     """
#     Spatial relation between hyperplanes and nd-box

#     note, hyperplanes should form a minimal description,
#     i.e. no redundant planes allowed
#     """

#     #       [ ] small optimization: 3^d - 2^d (may not be needed when using sphere around box?)
#     #            -- test only new corners for 'kids' boxes in traversal algo of n-ary tree
#     #           2D =  9 -  4 =  5 refined
#     #           3D = 27 -  8 = 19 refined
#     #           4D = 81 - 16 = 65 refined

#     all_dists = []
#     # -- test if the vertices of the box are outside of one of the planes
#     for plane in planes:
#         dists = []  # distances for this particular plane
#         for pt in box.vertices:
#             dist = plane.signed_distance(pt)
#             dists.append(dist)
#         # if so, for sure no interaction
#         if all(map(lambda x: x > 0, dists)):
#             return Interaction.NO_INTERACT
#         all_dists.extend(dists)
#     # -- see what the type of interaction is,
#     # if all vertices all inside:
#     #   then
#     #       fully contained
#     #   otherwise
#     #       some interaction
#     if all(map(lambda x: x <= 0, all_dists)):
#         # all points on correct side of all planes
#         # -> box is fullly contained by planes
#         return Interaction.CONTAINED
#     else:
#         # some interaction
#         return Interaction.INTERACT


# def test_relate2():
#     box = ndbox((0.0, 0.0), (33554432.0, 33554432.0))
#     planes = [
#         HyperPlane((-1.0, 0.0), -1.0),
#         HyperPlane((0.0, -1.0), -1.0),
#         HyperPlane((1.0, 0.0), -9.0),
#         HyperPlane((0.0, 1.0), -9.0),
#     ]
#     interaction = relate_box__box_based(planes, box)
#     print(interaction)

#     interaction, mask = relate_box__sweep_based__with_mask(planes, box, 0)
#     print(interaction)


# def test_hyperplane_construction():
#     pt = [5, 5]
#     hp = HyperPlane.from_normal_and_point(unit([1, 1]), pt)
#     for dim in range(len(pt)):
#         assert isclose(hp.through[dim], pt[dim])
#     # a hyperplane that runs through 5,5 as well, but constructed based on determinant
#     hp = HyperPlane.from_points([[0, 10], [10, 0]])
#     for dim in range(len(pt)):
#         assert isclose(hp.through[dim], pt[dim])


# def test_enter_exit():
#     hp = HyperPlane.from_normal_and_point([1], [0])
#     enter = sweep_box_with_plane_enter(hp, ndbox([0], [1]))
#     exit = sweep_box_with_plane_exit(hp, ndbox([0], [1]))
#     print(enter)
#     print(exit)

#     hp = HyperPlane.from_normal_and_point([-1], [0])
#     enter = sweep_box_with_plane_enter(hp, ndbox([0], [1]))
#     exit = sweep_box_with_plane_exit(hp, ndbox([0], [1]))
#     print(enter)
#     print(exit)


# def as_signs(v):
#     return "".join(["+" if component >= 0 else "-" for component in v])


# def relate__plane_box(plane, box):
#     # side length of box: E
#     box_side_dist = box.hi[0] - box.lo[0]
#     # length of diagonal: E * √n (:= √(n * E * E), with n=nDims
#     diagonal_dist = box_side_dist * box.dims ** 0.5

#     exit = sweep_box_with_plane_exit(plane, box)
#     exit_dist = plane.signed_distance(exit)

#     # print(plane, 'running through', plane.through, 'to', exit, 'with', exit_dist)
#     if exit_dist <= 0:
#         # mask_out += pow(2, index)
#         # only override result when not yet seen
#         print("exit l1506")  # √ tested
#         return Interaction.CONTAINED
#     else:
#         # perform distance comparison first (diagonal length & side of box)
#         # by doing that, we defer sweeping for the far point
#         abs_exit_dist = abs(exit_dist)
#         if abs_exit_dist > diagonal_dist:
#             # print('enter > diag')
#             print("exit l1514")  # √ tested
#             return Interaction.NO_INTERACT  # , mask_out
#         elif abs_exit_dist < box_side_dist:
#             # print('enter < side')
#             print("exit l1518")  # √ tested
#             return Interaction.INTERACT
#         else:
#             # get first hit point sweeping the box with the hyperplane
#             enter = sweep_box_with_plane_enter(plane, box)
#             enter_dist = plane.signed_distance(enter)
#             if enter_dist <= 0:
#                 print("exit l1525")
#                 return Interaction.INTERACT
#             else:
#                 print("exit l1528")
#                 return Interaction.NO_INTERACT  # , mask_out


# def test_directions():
#     for dims in range(1, 6):
#         # dims = 2
#         print("")
#         print("*" * dims, dims)
#         print("=" * 75)
#         negdir = -1
#         posdir = +1
#         directions = [[posdir], [negdir]]
#         for _ in range(1, dims):
#             new_directions = []
#             for v in directions:
#                 v1 = v[:]
#                 v1.append(negdir)
#                 v2 = v[:]
#                 v2.append(posdir)
#                 new_directions.extend([v1, v2])
#             directions = new_directions[:]
#         directions = [unit(v) for v in directions]
#         # print(directions)
#         for n in directions:
#             lo = [-1] * dims
#             hi = [1] * dims
#             box = ndbox(lo, hi)
#             support_points = [
#                 mul(n, central) for central in [-(10.0)**dims, 0.0, (10.0)**dims]
#             ]
#             expected_outcome = [
#                 Interaction.NO_INTERACT, Interaction.INTERACT, Interaction.CONTAINED
#             ]
#             for support_point, expectation in zip(support_points, expected_outcome):

#                 hp = HyperPlane.from_normal_and_point(n, support_point)
#                 # print('dist to ', support_point, 'is :=', hp.signed_distance(support_point), hp.signed_distance(hp.through))
#                 # print(hp.through)
#                 enter = sweep_box_with_plane_enter(hp, box)
#                 exit = sweep_box_with_plane_exit(hp, box)
#                 print(
#                     "sign of normal:", as_signs(n), " __ near:", enter, "far:", exit
#                 )  # ,
#                 print("plane at point", support_point)
#                 print(relate__plane_box(hp, box))
#                 assert relate__plane_box(hp, box) == expectation
#             print("")


# def test_non_obvious_cases():
#     """
#     Test the not obvious sweep cases, where we cannot just
#     based on second encounter and distances get what we are after

#     """
#     dims = 2
#     lo = [0] * dims
#     hi = [1] * dims
#     box = ndbox(lo, hi)
#     print(box)

#     # n = [0,+1]
#     # support_point = [0, 1.25]
#     # hp = HyperPlane.from_normal_and_point(n, support_point)
#     # hp = HyperPlane.from_points([ [0.0, -1.0], [1.25, 0] ])
#     # print(hp)
#     # print(relate__plane_box(hp, box))

#     hp = HyperPlane.from_points([ [0, 0], [1.25, 0] ])
#     print(hp)
#     print(relate__plane_box(hp, box))

#     hp = HyperPlane.from_points([ [0, -0.1], [1.25, 0] ])
#     print(hp)
#     print(relate__plane_box(hp, box))

#     box = ndbox([0,0], [2,2])
#     hp = HyperPlane.from_points([ [0, 2.05], [1, 2.1] ])
#     print(hp)
#     print(relate__plane_box(hp, box))

#     box = ndbox([0,0], [2,2])
#     hp = HyperPlane.from_points([ [0, 2.05 -2.15], [2, -0.0001] ])
#     print(hp)
#     print(relate__plane_box(hp, box))
#     with open("/tmp/box.txt", "w") as fh:
#         print(box.as_wkt_2d(), file=fh)
#     with open("/tmp/hp.txt", "w") as fh:
#         visualize_halfplane_2d(hp, fh)


#     box = ndbox([0,0], [2,2])
#     hp = HyperPlane.from_points([ [0, 2.05 -2.15], [2, 1] ])
#     print(hp)
#     print(relate__plane_box(hp, box))
#     with open("/tmp/box.txt", "w") as fh:
#         print(box.as_wkt_2d(), file=fh)
#     with open("/tmp/hp.txt", "w") as fh:
#         visualize_halfplane_2d(hp, fh)


#     box = ndbox([0,0], [2,2])
#     hp = HyperPlane.from_points([ [0, -0.001], [2, -0.001] ])
#     print(hp)
#     print(relate__plane_box(hp, box))
#     with open("/tmp/box.txt", "w") as fh:
#         print(box.as_wkt_2d(), file=fh)
#     with open("/tmp/hp.txt", "w") as fh:
#         visualize_halfplane_2d(hp, fh)

#     # for y in range(10, 100, 10):
#     #     n = unit([y,1])
#     #     print(n)
#     #     support_point = [0, 1.25]

#     #     hp = HyperPlane.from_normal_and_point(n, support_point)
#     #     # print('dist to ', support_point, 'is :=', hp.signed_distance(support_point), hp.signed_distance(hp.through))
#     #     # print(hp.through)
#     #     enter = sweep_box_with_plane_enter(hp, box)
#     #     exit = sweep_box_with_plane_exit(hp, box)
#     #     print(enter, '->', exit)
#     #     print(relate__plane_box(hp, box))

# def test_convex():
#     convex = ConvexPolytope.from_ndbox(ndbox([0, 0, 0], [10, 10, 10]))
#     print(convex)



# if __name__ == "__main__":
#     # test_directions()
    
#     test_non_obvious_cases()

#     # test_convex()
#     # test_enter_exit()
#     # play_with_deca()
#     # test_relate()
#     # test_hyperplane_construction()

#     # around_origin_2d()
#     # test_convex_polytope()

#     # play_with_cone()

#     # play_with_frustum()

#     # playground()
#     # project_3d_to_2d()
#     # test_from_points()

#     # test_diamant()
#     # test_what_is_inside()

#     # box = ndbox([5, 5], [15, 15])
#     # polytope = ConvexPolytope.from_ndbox(box)
#     # hp = HyperPlane.from_normal_and_point([0,1], [10,10])
#     # polytope.add_halfplane_to_convex(hp)
#     # hp = HyperPlane.from_normal_and_point([1,0], [10,10])
#     # polytope.add_halfplane_to_convex(hp)

#     # # FIXME: hid__vids
#     # # how to update it properly?
#     # # does it work when we have a tilted plane?

#     # # box = ndbox(range(6), range(10,16))
#     # # polytope = ConvexPolytope.from_ndbox(box)
