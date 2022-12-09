import unittest
from math import isclose

from pysfc.ndgeom.ndbox import ndbox
from pysfc.ndgeom.convex_polytope import ConvexPolytope
from pysfc.ndgeom.hyper_plane import HyperPlane

class ConvexPolytopeConstruction(unittest.TestCase):
    def setUp(self):
        self.boxes = []
        for d in range(2, 10):
            box = ndbox([-5] * d, [+5] * d)
            self.boxes.append(box)
            box = ndbox([5] * d, [10] * d)
            self.boxes.append(box)
            box = ndbox([-10] * d, [-5] * d)
            self.boxes.append(box)

    def test_construction(self):
        for box in self.boxes:
            polytope = ConvexPolytope.from_ndbox(box)
            self.assertTrue(polytope.is_consistent())


class ConvexPolytopeNodesAndEdges(unittest.TestCase):
    def test_nodes_and_edges(self):
        box = ndbox([-5, -5], [+5, +5])
        polytope = ConvexPolytope.from_ndbox(box)
        self.assertTrue(polytope.is_consistent())
        nodes, edges = polytope.as_nodes_and_edges()
        self.assertEqual(nodes, [(-5.0, -5.0), (-5.0, 5.0), (5.0, -5.0), (5.0, 5.0)])
        self.assertEqual(edges, [(0, 2), (0, 1), (1, 3), (2, 3)])


class ConvexPolytopeInsertion(unittest.TestCase):
    def setUp(self):
        self.boxes = []
        for d in range(2, 10):
            box = ndbox([-5] * d, [+5] * d)
            self.boxes.append(box)
            box = ndbox([5] * d, [10] * d)
            self.boxes.append(box)
            box = ndbox([-10] * d, [-5] * d)
            self.boxes.append(box)

    def test_insertion_cut_off_bottom(self):
        for box in self.boxes:
            polytope = ConvexPolytope.from_ndbox(box)
            self.assertTrue(polytope.is_consistent())
            hp = HyperPlane.from_normal_and_point(
                [0] * (box.dims - 1) + [-1], [0] * box.dims
            )  # embed hyperplane in structure
            polytope.add_halfplane_to_convex(hp)
            self.assertTrue(polytope.is_consistent())

    def test_insertion_cut_off_top(self):
        for box in self.boxes:
            polytope = ConvexPolytope.from_ndbox(box)
            self.assertTrue(polytope.is_consistent())
            hp = HyperPlane.from_normal_and_point(
                [0] * (box.dims - 1) + [1], [0] * box.dims
            )  # embed hyperplane in structure
            polytope.add_halfplane_to_convex(hp)
            self.assertTrue(polytope.is_consistent())

    def test_insertion_cut_off_positive_corner_top(self):
        for box in self.boxes:
            polytope = ConvexPolytope.from_ndbox(box)
            self.assertTrue(polytope.is_consistent())
            hp = HyperPlane.from_normal_and_point(
                [1 / box.dims ** 0.5] * box.dims, [0] * box.dims
            )  # embed hyperplane in structure
            polytope.add_halfplane_to_convex(hp)
            self.assertTrue(polytope.is_consistent())

    def test_insertion_cut_off_negative_corner_bottom(self):
        for box in self.boxes:
            polytope = ConvexPolytope.from_ndbox(box)
            self.assertTrue(polytope.is_consistent())
            hp = HyperPlane.from_normal_and_point(
                [-1 / box.dims ** 0.5] * box.dims, [0] * box.dims
            )  # embed hyperplane in structure
            polytope.add_halfplane_to_convex(hp)
            self.assertTrue(polytope.is_consistent())


if __name__ == "__main__":
    unittest.main()
