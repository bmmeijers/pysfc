from enum import Enum


class Interaction(Enum):
    UNINITIALIZED = -1  # interaction not yet determined
    NO_INTERACT = 0  # no interaction
    INTERACT = 1  # possibly some interaction
    CONTAINED = 2  # interaction, for sure



def relate_box_box(rect, qrt, _ = None):
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
    # ncmp = 1
    # for d in range(dims):
    #     ncmp &= rect.lo[d] == qrt.lo[d] and rect.hi[d] == qrt.hi[d]
    # if ncmp:
    #     return Interaction.CONTAINED

    # fully contains (or equal), rect fully contains qrt
    ncmp = 1
    for d in range(dims):
        ncmp &= rect.lo[d] <= qrt.lo[d] and rect.hi[d] >= qrt.hi[d]
    if ncmp:
        return Interaction.CONTAINED, None

    # intersects, the two nd-boxes interact 
    # (either on the boundary or internally)
    ncmp = 1
    for d in range(dims):
        ncmp &= rect.lo[d] < qrt.hi[d] and rect.hi[d] > qrt.lo[d]
    if ncmp:
        return Interaction.INTERACT, None

    # no overlap
    return Interaction.NO_INTERACT, None


def relate_planes_box__box_based(planes, box):
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

    all_dists = []
    # -- test if the vertices of the box are outside of one of the planes
    for plane in planes:
        dists = []  # distances for this particular plane
        for pt in box.vertices:
            dist = plane.signed_distance(pt)
            dists.append(dist)
        # if so, for sure no interaction
        if all(map(lambda x: x > 0, dists)):
            return Interaction.NO_INTERACT
        all_dists.extend(dists)
    # -- see what the type of interaction is,
    # if all vertices all inside:
    #   then
    #       fully contained
    #   otherwise
    #       some interaction
    if all(map(lambda x: x <= 0, all_dists)):
        # all points on correct side of all planes
        # -> box is fullly contained by planes
        return Interaction.CONTAINED
    else:
        # some interaction
        return Interaction.INTERACT




def relate_planes_box__sweep_based__with_mask(planes, box, mask=0, stats=None):
    """
    spatial relation between hyperplanes and nd-box
    remembering outcome of relation test in mask so that
    some hyperplanes might be skipped for testing
    when nd-box is split in 2^n smaller boxes
    """
    # global TESTED_CT, CONTAINED_CT
    from math import sqrt

    # pre-condition: some planes to test against
    assert len(planes) > 0
    result = Interaction.UNINITIALIZED
    # another pre-condition:
    # box has same edge lengths in all dimensions
    # (not tested)
    # side length of box: E
    box_side_dist = box.hi[0] - box.lo[0]
    # length of diagonal: E * √n (:= √(n * E * E), with n=nDims
    diagonal_dist = box_side_dist * sqrt(box.dims)
    # mask_out will be changed, if we detect that box now is contained
    mask_out = mask
    for index, plane in enumerate(planes):
        # print(f'testing {index} {plane} result so far: {result}')
        if stats is not None:
            stats.tested += 1
        # TESTED_CT += 1
        # we know that the plane is contained from parent index boxes,
        # skip checking this particular hyperplane again
        is_contained = bool(mask & (1 << index))
        # print(f'{is_contained}')
        if is_contained:
            if stats is not None:
                stats.contained += 1
            # CONTAINED_CT += 1
            if result == Interaction.UNINITIALIZED:
                # print('l280')
                result = Interaction.CONTAINED
            continue
        # first hit point sweeping the box with the hyperplane
        # FIXME: for this logic to be correct, the
        # point that is tested first
        # needs to be the corner point that is
        # 'least-likely-to-be-inside'
        # along the direction of the normal
        #
        # then we can add the two distances (box size / diagonal size)
        # and eventually test the point that is 'most-likely-to-be-inside'
        # enter = sweep_box_with_plane_enter(plane, box)
        # enter_dist = plane.signed_distance(enter)

        exit = sweep_box_with_plane_second_encounter(plane, box)
        exit_dist = plane.signed_distance(exit)


        # print(plane, 'to', enter, 'with', enter_dist)
        if exit_dist < 0:
            mask_out += pow(2, index)
            # only override result when not yet seen
            if result == Interaction.UNINITIALIZED:
                # print('l289')
                result = Interaction.CONTAINED
        else:
            # perform distance comparison first (diagonal length & side of box)
            # by doing that, we defer sweeping for the far point
            abs_enter_dist = abs(exit_dist)
            if abs_enter_dist > diagonal_dist:
                # print('enter > diag')
                return Interaction.NO_INTERACT, mask_out
            elif abs_enter_dist < box_side_dist:
                # print('enter < side')
                result = Interaction.INTERACT
            else:
                # get last hit point sweeping the box with the hyperplane
                enter = sweep_box_with_plane_first_encounter(plane, box)
                enter_dist = plane.signed_distance(enter)
                if enter_dist < 0:
                    result = Interaction.INTERACT
                else:
                    return Interaction.NO_INTERACT, mask_out
    # post-condition
    # result here should be either INTERACT / CONTAINED
    # (no longer uninitialized)
    # print('current result', result)
    assert result in (Interaction.CONTAINED, Interaction.INTERACT)
    return result, mask_out


def sweep_box_with_plane_first_encounter(hyperplane, box):
    """Return enter point when 'sweeping' the box with the hyperplane

    enter: point that is 'hit' first while sweeping the box with the hyperplane
    along the direction of the normal of the hyperplane
    (i.e. *first* encounter point along ray swept by hyperplane)
    """
    enter = [0.0] * box.dims
    for i, n in enumerate(hyperplane.n):
        if n < 0:
            enter[i] = box.hi[i]
        else:
            enter[i] = box.lo[i]
    return enter


def sweep_box_with_plane_second_encounter(hyperplane, box):
    """Return exit point 'sweeping' the box with the hyperplane

    exit: point where hyperplane leaves the box while sweeping the box with the
    hyperplane along the direction of the normal of the hyperplane
    (i.e. *second* encounter point along ray swept by hyperplane)
    """
    exit = [0.0] * box.dims
    for i, n in enumerate(hyperplane.n):
        if n < 0:
            exit[i] = box.lo[i]
        else:
            exit[i] = box.hi[i]
    return exit
