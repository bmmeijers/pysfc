from math import isclose
from pysfc.ndgeom.hyper_plane import HyperPlane
from pysfc.linalg import add, mul, sub


class ConvexPolytope:
    """A convex polytope, represented by hyperplanes

    During incremental construction the 0D (vertices) and 1D (edges) primitives
    are derived from the hyperplane equations
    """

    def __init__(self):
        self.vid__vertex = {}  # vertices, indexed by id
        self.hid__hyperplane = {}  # hyperplanes, indexed by id
        self.vid__stars = {}  # stars (halfedges), indexed by vertex id
        self.vid_star__twin_vid = {}  # vid+star (halfedge) to twin vertex
        self.hid__star_count = {}  # hyperplane to vertex references
        self.dim = 0  # dimension

    def hyperplanes(self):
        return list(self.hid__hyperplane.values())

    def show(self):
        """print a representation of this polytope"""
        print("Polytope of dimension", self.dim)
        names = [
            "vid__vertex",
            "hid__hyperplane",
            "vid__stars",
            "vid_star__twin_vid",
            "hid__star_count",
        ]
        to_show = [
            self.vid__vertex,
            self.hid__hyperplane,
            self.vid__stars,
            self.vid_star__twin_vid,
            self.hid__star_count,
        ]
        for name, show in zip(names, to_show):
            print(name)
            print("-" * 80)
            for key, val in sorted(show.items()):
                print(key, "  |  ", val)
            print()

    def copy(self):
        from copy import deepcopy

        copy = ConvexPolytope()
        copy.vid__vertex = deepcopy(self.vid__vertex)
        copy.hid__hyperplane = deepcopy(self.hid__hyperplane)
        copy.vid__stars = deepcopy(self.vid__stars)
        copy.vid_star__twin_vid = deepcopy(self.vid_star__twin_vid)
        copy.hid__star_count = deepcopy(self.hid__star_count)
        copy.dim = self.dim
        return copy

    def as_nodes_and_edges(self):
        """Getting a more topology oriented view from this convex polytope"""
        # construct nodes
        vertices_lut = {}
        nodes = []
        for vid, vertex in self.vid__vertex.items():
            vertices_lut[vid] = len(nodes)
            nodes.append(vertex)
        # an edge has 2 vertices, start and end
        # only output an edge if the start index < end index
        edges = []
        for vid in self.vid__vertex:
            stars = self.vid__stars[vid]
            for star in stars:
                start = vertices_lut[vid]
                try:
                    end = vertices_lut[self.vid_star__twin_vid[(vid, star)]]
                except KeyError:
                    print("error with (vid, star)", (vid, star))
                    continue
                if start < end:
                    edges.append((start, end))
        return nodes, edges

    @staticmethod
    def from_ndbox(box):
        """constructor - from an ndbox get a convex polytope"""
        conv = ConvexPolytope()
        # create a convex polytope, initialized as hypercube with the size of the nd-box
        conv.dim = box.dims

        # create vertices, 2**d vertices (3d = 8 vertices)
        # using Π_{0..d} of (min, max)
        vertices = [[box.lo[0]], [box.hi[0]]]
        for d in reversed(range(1, conv.dim)):
            new_vertices = []
            for v in vertices:
                v1 = v[:]
                v1.append(box.lo[d])
                v2 = v[:]
                v2.append(box.hi[d])
                new_vertices.extend([v1, v2])
            vertices = new_vertices[:]
        conv.vid__vertex = dict(
            list(enumerate((tuple(vertex) for vertex in vertices), start=0))
        )
        assert len(conv.vid__vertex) == 2 ** conv.dim
        # in 3D this creates vertices as follows
        # bit pattern (of the vertex identifier) goes slowest for x;
        # z goes fastest, i.e.:
        # xyz
        # ---
        # 000 0 (0.0, 0.0, 0.0)
        # 001 1 (0.0, 0.0, 10.0)
        # 010 2 (0.0, 10.0, 0.0)
        # 011 3 (0.0, 10.0, 10.0)
        #
        # 100 4 (10.0, 0.0, 0.0)
        # 101 5 (10.0, 0.0, 10.0)
        # 110 6 (10.0, 10.0, 0.0)
        # 111 7 (10.0, 10.0, 10.0)
        #
        # for vid, vertex in conv.vid__vertex.items():
        #     print(bin(vid)[2:].zfill(conv.dim), vid, vertex)

        # create hyperplanes, 2*d (3d = 6 planes)
        # indices = []
        # hplanes = []
        hplane_at_side = []
        for dim in range(conv.dim):
            tmp = []
            for _ in range(2):
                tmp.append(None)
            hplane_at_side.append(tmp)
        # vector with all 0.0 of dims length
        null = [0.0] * conv.dim
        conv.hid__hyperplane = {}
        # lo dimension, normal negative
        for hid, dim in enumerate(range(conv.dim), start=1024):
            w = null[:]
            w[dim] = -1.0
            # hplanes.append()
            hplane = HyperPlane.from_normal_and_point(w, box.lo)
            # conv.hid__hyperplane[hid]
            # hplanes.append(hplane)
            conv.hid__hyperplane[hid] = hplane
            # indices.append(hid)
            hplane_at_side[dim][0] = hid
        # hi dimension side = positive normal for outwards
        for hid, dim in enumerate(range(conv.dim), start=hid + 1):
            w = null[:]
            w[dim] = +1.0
            hplane = HyperPlane.from_normal_and_point(w, box.hi)
            conv.hid__hyperplane[hid] = hplane
            # hplanes_.append(hplane)
            # hplanes.append(hplane)
            # indices.append(hid)
            hplane_at_side[dim][1] = hid
        # conv.hid__hyperplane = dict(enumerate(hplanes, start=1024))
        assert len(conv.hid__hyperplane) == 2 * conv.dim

        # dims = "xyz"
        for vid, vertex in conv.vid__vertex.items():
            # print(vid)
            match = [None] * conv.dim
            for dim in range(conv.dim):
                # print(bin(vid))
                # this comes from the pattern in which dimension order the vertices are constructed
                to_shift = conv.dim - 1 - dim
                # print(to_shift)
                is_hi = int(bool(vid & (1 << to_shift)))
                # print (dims[dim], is_hi)
                hid = hplane_at_side[dim][is_hi]
                hplane = conv.hid__hyperplane[hid]
                # print("getting", hid,"at dim", dim, "at side", is_hi)
                # print("the distance", vid, "vtx", vertex, hplane.signed_distance(vertex))
                assert hplane.signed_distance(vertex) == 0

                # print(vid, hid)
                match[dim] = hid
            # print(vid, match)
            assert len(match) == conv.dim
            # for the hyperplanes that run through this vertex,
            # group them into a series of stars

            # a star is a sorted tuple of hyper plane ids (ordered by id asc)
            # these planes, when intersected geometrically
            # should give the edge support line geometry
            # (:= infinite 1D line in nD space, or 1D hyperplane)
            for h1 in match:
                star = [None] * (conv.dim - 1)
                sid = 0
                for h2 in match:
                    if h1 != h2:
                        star[sid] = h2
                        sid += 1
                star.sort()
                star = tuple(star)
                for hid in star:
                    if hid not in conv.hid__star_count:
                        # conv.hid__stars[hid] = []
                        conv.hid__star_count[hid] = 0
                    # conv.hid__stars[hid].append(star)
                    conv.hid__star_count[hid] += 1
                assert len(star) == conv.dim - 1
                if vid not in conv.vid__stars:
                    conv.vid__stars[vid] = []
                conv.vid__stars[vid].append(star)
        # we now have to have dim stars per vertex
        for vid in conv.vid__stars:
            assert len(conv.vid__stars[vid]) == conv.dim

        ## alternative implementation, brute force matching
        ## of vertices against a list of hyperplanes
        # # find all the hyperplanes on which this vertex lies
        # # and make stars of each vertex
        # for vid, vertex in conv.vid__vertex.items():
        #     # by brute force geometry matching
        #     # (a vertex lies on a hyperplane if its signed distance is 0.0)
        #     # (note, order: O(v h d))
        #     #
        #     # this could be done a bit smarter, if we make a temporary structure
        #     # of hyperplanes, organized by dimension number, min max
        #     # [(min-x, max-x), (min-y, max-y), ... ]
        #     # then we can check this list for the normal vector, if it is
        #     # not 0.0 it should be a match
        #     match = []
        #     for hid, hplane in conv.hid__hyperplane.items():
        #         if hplane.signed_distance(vertex) == 0:
        #             match.append(hid)
        #     assert len(match) == conv.dim
        #     # for the hyperplanes that run through this vertex,
        #     # group them into a series of stars

        #     # a star is a sorted tuple of hyper plane ids (ordered by id asc)
        #     # these planes, when intersected geometrically
        #     # should give the edge support line geometry
        #     # (:= infinite 1D line in nD space, or 1D hyperplane)
        #     for h1 in match:
        #         star = [None] * (conv.dim - 1)
        #         sid = 0
        #         for h2 in match:
        #             if h1 != h2:
        #                 star[sid] = h2
        #                 sid += 1
        #         star.sort()
        #         star = tuple(star)
        #         for hid in star:
        #             if hid not in conv.hid__star_count:
        #                 # conv.hid__stars[hid] = []
        #                 conv.hid__star_count[hid] = 0
        #             # conv.hid__stars[hid].append(star)
        #             conv.hid__star_count[hid] += 1
        #         assert len(star) == conv.dim - 1
        #         if vid not in conv.vid__stars:
        #             conv.vid__stars[vid] = []
        #         conv.vid__stars[vid].append(star)
        # # we now have to have dim stars per vertex
        # for vid in conv.vid__stars:
        #     assert len(conv.vid__stars[vid]) == conv.dim

        # invert the vid__stars dict (this is a temporary dict)
        stars__vids = {}
        for vid, stars in conv.vid__stars.items():
            for star in stars:
                if star not in stars__vids:
                    stars__vids[star] = []
                stars__vids[star].append(vid)
        # and with this temp dict, we build the vid_star__twin_vid
        # which contains a vid + star to point to the twin vid
        for (star, vids) in stars__vids.items():
            for vid in vids:
                for twin_vid in vids:
                    if vid != twin_vid:
                        conv.vid_star__twin_vid[(vid, star)] = twin_vid
        # all edges should be represented twice, number of edges in hypercube is: 2 ** (d-1) * d
        assert len(conv.vid_star__twin_vid) / 2 == 2 ** (conv.dim - 1) * conv.dim

        # can we state something about the number of edges per hyperplane for a hypercube
        # I guess we can:
        # number of edges = 2 ** (conv.dim-1) * conv.dim
        # number of planes = 2 * dim
        #
        # so edges / planes :
        # (2 ** (conv.dim-1) * conv.dim) / (2 * dim)

        # can we state something about the number of vertices per plane (in 2D: 2, in 3D: 4, in 4D: ?... what is this series?)
        # print(conv.hid__vids)
        # conv.show()
        # assert conv.is_consistent()
        return conv

    def is_consistent(self):
        """
        consistency check of combinatorial structure (topology) of the
        convex polytope represented
        """
        # consistency check
        # vertices represented
        vids1 = set(self.vid__vertex.keys())
        vids2 = set(self.vid__stars.keys())
        vids3 = set([vid_star[0] for vid_star in self.vid_star__twin_vid.keys()])
        if vids1 != vids2:
            return False
        if vids1 != vids3:
            return False
        # hyperplanes represented
        hids1 = set(self.hid__hyperplane.keys())
        # hids2 = set(self.hid__stars.keys())
        # if hids1 != hids2:
        #     return False
        # get all stars and do check if they are referring to existing hyperplanes
        hids3 = set()
        for stars in self.vid__stars.values():
            for star in stars:
                for hid in star:
                    hids3.add(hid)
        if hids1 != hids3:
            return False
        tmp_count = {}
        for stars in self.vid__stars.values():
            for star in stars:
                for hid in star:
                    if hid not in tmp_count:
                        tmp_count[hid] = 0
                    tmp_count[hid] += 1
        if tmp_count != self.hid__star_count:
            return False

        # twin vertices should point to each other
        for vid_star in self.vid_star__twin_vid:
            vid1, star = vid_star
            twin_vid = self.vid_star__twin_vid[vid_star]
            vid2 = self.vid_star__twin_vid[(twin_vid, star)]
            if vid1 != vid2:
                return False

        return True

    def clear(self):
        """Makes the convex polytope empty (ø)"""
        # clear all internal dict's
        self.vid__vertex.clear()
        self.hid__hyperplane.clear()
        self.vid__stars.clear()
        self.vid_star__twin_vid.clear()
        self.hid__star_count.clear()

    def add_halfplane_to_convex(self, hp: HyperPlane):
        """add a halfplane to an existing convex polytope"""
        some_included = False
        some_excluded = False

        vids_to_delete = set([])

        distances = {}  # record the distance
        # we need this when we have a proper plane to insert
        for vid, vertex in self.vid__vertex.items():
            dist = hp.signed_distance(vertex)
            # print(vid, vertex, dist)
            distances[vid] = dist
            if dist <= 0:
                some_included = True
            else:
                some_excluded = True
                vids_to_delete.add(vid)
            # remove name bindings
            del vid
            del vertex

        # deal with the plane now that we have all distances
        if not some_included:
            # hyperplane makes convex empty
            self.clear()
            return

        elif not some_excluded:
            # hyperplane does not add anything, no changes needed
            return

        else:
            # hyperplane is intersecting edges, let us add it
            new_hid = max(self.hid__hyperplane) + 1
            self.hid__hyperplane[new_hid] = hp
            # try:
            new_vid = max(self.vid__vertex) + 9

            # find which hyperplanes will not be needed after insertion of new hyperplane
            # FIXME: could it be sufficient to just use a counter for how many vertices
            # are associated with the hyperplane?
            # temp dict to match end points of edges, after all new vertices have been put
            star__vids = {}
            stars_to_delete = []
            # /--+ make a pass to shrink all edges related to the evicted vertices
            for evid in vids_to_delete:
                stars = self.vid__stars[evid]
                for star in stars:
                    # print(star)
                    twin_vid = self.vid_star__twin_vid[(evid, star)]
                    # print(twin_vid) # other end of edge

                    if twin_vid in vids_to_delete:

                        stars_to_delete.append(((evid, star), True))
                        # True because we need to remove associated hyperplanes

                    else:

                        stars_to_delete.append(((evid, star), False))
                        # False because we need to keep associated hyperplanes

                        # twin_vid will be kept, interpolate the position for new vertex geometry
                        # print("shrink this edge", evid, twin_vid, star)
                        assert (
                            distances[evid] > 0
                        )  # this vertex is to be deleted (outside of hyperplane), so + distance (it can not be 0)
                        assert (
                            distances[twin_vid] <= 0
                        )  # other vertex is to be kept, so - distance
                        #
                        total_length = distances[evid] - distances[twin_vid]
                        # print(distances[evid], distances[twin_vid], total_length)
                        #
                        if total_length <= 0.0:
                            new_vertex_geom = self.vid__vertex[evid]
                            raise NotImplementedError("plane through 0 length edge?")
                        else:
                            outside_frac = distances[evid] / total_length
                            inside_frac = -distances[twin_vid] / total_length
                            assert isclose(outside_frac + inside_frac, 1.0)
                            new_vertex_geom = interpolate(
                                self.vid__vertex[twin_vid],
                                self.vid__vertex[evid],
                                inside_frac,
                            )

                        # put the new vertex in relevant dict's
                        new_vid += 1
                        self.vid__vertex[new_vid] = new_vertex_geom

                        # set up the vertex / other edge relation ship
                        # this is on the part of the edge we know
                        # the edges that lie in the hyperplane to other new vertices will follow later
                        # assert self.vid_star__twin_vid[(twin_vid, star)] == evid # does not always hold...
                        self.vid_star__twin_vid[(new_vid, star)] = twin_vid
                        self.vid_star__twin_vid[(twin_vid, star)] = new_vid
                        assert self.vid_star__twin_vid[(twin_vid, star)] == new_vid

                        # check distances to hyperplanes in the current star (maybe use as robustness measure)
                        # hp_dists = []
                        # for hid in star:
                        #     hp = self.hid__hyperplane[hid]
                        #     d = abs(hp.signed_distance(new_vertex_geom))
                        #     hp_dists.append(d)
                        # print("dists to new hp", hp_dists, "eps:", max(hp_dists))

                        # we know that the plane we add and the hyperplanes
                        # on the edge we are intersecting are
                        # the matching hyperplanes to this new vertex
                        #
                        # let us set up the appropriate stars
                        assert new_vid not in self.vid__stars
                        # print("current star", star)
                        # print("     new hid", new_hid)
                        match = list(star) + [new_hid]

                        def process_match(self, match, current_star, star__vids):
                            # in a inner function, not to bother with 'star' variable in this loop
                            for h1 in match:
                                star = [None] * (self.dim - 1)
                                sid = 0
                                for h2 in match:
                                    if h1 != h2:
                                        star[sid] = h2
                                        sid += 1
                                star.sort()
                                star = tuple(star)
                                assert len(star) == self.dim - 1
                                if star != current_star:
                                    # the shrinking edge is already there,
                                    # thus do not add it again in the hid__star_count dict
                                    for hid in star:
                                        if hid not in self.hid__star_count:
                                            self.hid__star_count[hid] = 0
                                        self.hid__star_count[hid] += 1
                                if new_vid not in self.vid__stars:
                                    self.vid__stars[new_vid] = []
                                self.vid__stars[new_vid].append(star)
                                # record the star / new vertex combi, so we
                                # can update the twin vertex at the other end
                                if star not in star__vids:
                                    star__vids[star] = []
                                star__vids[star].append(new_vid)

                        process_match(self, match, star, star__vids)
            # match the edges based on the star__vids dict made earlier
            # (these are the edges that run 'over' the inserted hyperplane)
            for (star, vids) in star__vids.items():
                for vid in vids:
                    for twin_vid in vids:
                        if vid != twin_vid:
                            self.vid_star__twin_vid[(vid, star)] = twin_vid
            # \--+ end of shrinking edges

            # cleaning up vertex data not needed any more
            for evid in vids_to_delete:
                self.vid__vertex.pop(evid)
                self.vid__stars.pop(evid)

            # remove hyperplanes and edges that we do not need any more
            for vid_star, remove_hyperplane_reference in stars_to_delete:
                vid, star = vid_star
                # print("v", vid, "s", star)
                if remove_hyperplane_reference:
                    for hid in star:
                        self.hid__star_count[hid] -= 1
                        if self.hid__star_count[hid] == 0:
                            self.hid__hyperplane.pop(hid)
                            self.hid__star_count.pop(hid)
                self.vid_star__twin_vid.pop(vid_star)

            # after all modifications, check if we are still consistent
            # self.is_consistent()


def interpolate(one, other, fraction):
    """interpolate vertex coordinates, starting at one, going fraction of vec towards other"""
    return add(one, mul(sub(other, one), fraction))
