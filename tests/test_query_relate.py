import unittest
from pysfc.ndgeom.ndbox import ndbox
from pysfc.ndgeom.hyper_plane import HyperPlane
from pysfc.query.relate import Interaction
from pysfc.query.relate import relate_planes_box__sweep_based__with_mask
from pysfc.query.relate import relate_box_box
from pysfc.linalg import unit, mul

class RelateBoxBoxTest(unittest.TestCase):
    def test_relate_box_box(self):
        #assert False
        space = ndbox([0], [10])
        
        for query, expected_outcome in zip([
            ndbox([0], [10]),
            ndbox([5], [6]),
            ndbox([-2], [-1]),
            ndbox([-1], [0]),
            ndbox([-1], [1]),
            ndbox([9], [11]),
            ndbox([10], [11]),
            ndbox([11], [12]),
            ], [
                Interaction.CONTAINED,
                Interaction.CONTAINED,
                Interaction.NO_INTERACT,
                Interaction.NO_INTERACT, # TOUCHING BOX 
                Interaction.INTERACT,
                Interaction.INTERACT,
                Interaction.NO_INTERACT, # TOUCHING BOX 
                Interaction.NO_INTERACT,
            ]):
            relation = relate_box_box(space, query)
            self.assertEqual(relation[0], expected_outcome)

class RelatePlanesBoxTest(unittest.TestCase):
    def test_relate_in_all_directions(self):
        for dims in range(1, 10):
            # dims = 2
            # print("")
            # print("*" * dims, dims)
            # print("=" * 75)
            negdir = -1
            posdir = +1
            directions = [[posdir], [negdir]]
            for _ in range(1, dims):
                new_directions = []
                for v in directions:
                    v1 = v[:]
                    v1.append(negdir)
                    v2 = v[:]
                    v2.append(posdir)
                    new_directions.extend([v1, v2])
                directions = new_directions[:]
            directions = [unit(v) for v in directions]
            # print(directions)
            for n in directions:
                lo = [-1] * dims
                hi = [1] * dims
                box = ndbox(lo, hi)
                support_points = [
                    mul(n, central) for central in [-(10.0)**dims, 0.0, (10.0)**dims]
                ]
                expected_outcome = [
                    Interaction.NO_INTERACT, Interaction.INTERACT, Interaction.CONTAINED
                ]
                for support_point, expectation in zip(support_points, expected_outcome):

                    hp = HyperPlane.from_normal_and_point(n, support_point)
                    # print('dist to ', support_point, 'is :=', hp.signed_distance(support_point), hp.signed_distance(hp.through))
                    # print(hp.through)
                    # enter = sweep_box_with_plane_first_encounter(hp, box)
                    # exit = sweep_box_with_plane_second_encounter(hp, box)
                    # print(
                    #     "sign of normal:", as_signs(n), " __ near:", enter, "far:", exit
                    # )  # ,
                    # print(hp)
                    # print(box)
                    # print("plane at point", support_point)
                    # print(relate_box__sweep_based__with_mask([hp], box))
                    self.assertEqual(relate_planes_box__sweep_based__with_mask([hp], box)[0], expectation)
                # print("")

    # def test_relate(self):
    #     hp = HyperPlane.from_points([[0, 10], [10, 0]])
    #     print(hp.through)

    #     # this is fine, if we construct the plane, then we get a point through it
    #     pt = [-5, -5]
    #     hp = HyperPlane.from_normal_and_point(unit([1, -1]), pt)
    #     for dim in range(len(pt)):
    #         self.assertAlmostEqual(hp.signed_distance(hp.through), 0)
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
    #     enter_point = sweep_box_with_plane_first_encounter(hp, box1)
    #     print(enter_point)
    #     assert enter_point == [1.0, 1.0]
    #     exit_point = sweep_box_with_plane_second_encounter(hp, box1)
    #     assert exit_point == [0.0, 0.0]
    #     interaction, mask = relate_box__sweep_based__with_mask([hp], box1, 0)
    #     assert interaction == Interaction.CONTAINED
    #     assert mask == 1
    #     # unittest.main()


    def test_with_mask(self):
        hp = HyperPlane((1.0,), -10.0)
        box = ndbox((-1.0,), (1.0,))
        # with mask
        outcome = relate_planes_box__sweep_based__with_mask([hp], box, mask=1)
        assert outcome[0] == Interaction.CONTAINED

    def test_first_encounter1(self):
        box = ndbox([0,0], [2,2])
        hp = HyperPlane.from_points([ [0, 2.05 -2.15], [2, 1] ])
        self.assertEqual(relate_planes_box__sweep_based__with_mask([hp], box)[0], Interaction.INTERACT)

    def test_first_encounter2(self):
        box = ndbox([0,0], [2,2])
        hp = HyperPlane.from_points([ [0, -0.001], [2, -0.001] ])
        self.assertEqual(relate_planes_box__sweep_based__with_mask([hp], box)[0], Interaction.NO_INTERACT)
        # print(hp)
        # print(relate__plane_box(hp, box))
        # with open("/tmp/box.txt", "w") as fh:
        #     print(box.as_wkt_2d(), file=fh)
        # with open("/tmp/hp.txt", "w") as fh:
        #     visualize_halfplane_2d(hp, fh)


    # box = ndbox([0,0], [2,2])
    # hp = HyperPlane.from_points([ [0, -0.001], [2, -0.001] ])
    # print(hp)
    # print(relate__plane_box(hp, box))
    # with open("/tmp/box.txt", "w") as fh:
    #     print(box.as_wkt_2d(), file=fh)
    # with open("/tmp/hp.txt", "w") as fh:
    #     visualize_halfplane_2d(hp, fh)


if __name__ == "__main__":
    unittest.main()
    # small_test()