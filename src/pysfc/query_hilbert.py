from pysfc.speedups.relate import ndbox, relate
from collections import deque

from traversal import morton_zcode_ncode
from traversal import path_to_hilbert

from hilbert import initial_start_end, gray_decode_travel, child_start_end

from operator import itemgetter


def hquery(query, maxdepth=None, maxbits=63):
    """
    Traverse n-ary tree, each node its childs are visited in H-order
    """
    dim = query.dims            # the number of dimensions of the top box
    hmask = 2 ** dim - 1        # a number with all 1's for every bit

    # >> # FIXME: some options:
    # - from metadata table, or
    # - given by user (specifying max depth to visit), or
    # - inferred from largest coordinate inside the query box?
    # how many levels we want max to visit in the tree
    order = maxbits / dim       # -- influences maximum size of the domain
    # <<

    # -- determine how deep to go at maximum
    if maxdepth is None:
        maxdepth = order + 1
    maxdepth = min(maxdepth + 1, order + 1)

    # min and max, min and max coordinate values as list
    # each has a length equal to the dimension
    l_root = [0] * dim
    h_root = [2 ** order] * dim
    root_box = ndbox(l_root, h_root)  # Here we determine the size of the cube

    # >> dependent on order (depth of the nary-tree) !!!
    # (Note, otherwise curve will be running in wrong direction)
    hstart, hend = initial_start_end(order, dim)
    # <<

    # >> # FIXME: instead of having the path as outcome, we probably could
    # already obtain the hilbert code so far obtained? ->
    # will reduce the computation that is the same while going down the path
    result = []
    # <<

    # What we place in the stack:
    # (Hilbert-code of the node at this level -- note: not used,
    #  node its box,
    #  (start, end) -- how does the curve enter/exit at this node),
    #  level of the node in the tree,
    #  path traversed -- H-codes of the parents)
    traverse_childs_order = [ncode for _, ncode in morton_zcode_ncode(dim, False)]
    root = (None, root_box, (hstart, hend), 1, ())
    deq = deque([root])
    while deq:
        node = deq.popleft()
        cur_hcode, cur_box, (cur_hstart, cur_hend), cur_level, cur_hpath = node
        # determine the relation of the query box and the current node box
        ndcmp = relate(query, cur_box)
        # -- equal or fully contains, done: add to result
        if ndcmp in (0, 1,):
            result.append(cur_hpath)
        # -- partial overlap
        elif ndcmp in (2,):
            if cur_level < maxdepth:
                # we have not yet reached the lowest level, recurse
                childs = []
                next_level = cur_level + 1
                for ncode in traverse_childs_order:
                # for _, ncode in morton_zcode_ncode(dim, False):
                    hcode = gray_decode_travel(
                        cur_hstart, cur_hend, hmask, ncode)
                    next_hpath = cur_hpath + (hcode, )
                    next_box = cur_box.nth_child_norder(ncode)
                    next_start, next_end = child_start_end(
                        cur_hstart, cur_hend, hmask, hcode)
                    childs.append(
                        (hcode, next_box, (next_start, next_end), next_level, next_hpath))
                # -- visit the childs in order of the hilbert curve
                # (we need to sort the childs list for this)
                # !! in the tuple 0 = hcode, sort on local hilbert index (hcode)
                childs.sort(key=itemgetter(0))
                for child in childs:
                    deq.append(child)
            else:
                result.append(cur_hpath)
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp in (-1,):
            continue

    # Transform paths into ranges
    ranges = []
    for path in result:
        start, end = path_to_hilbert(
            path, dim, order)  # << IMPORTANT: order = total depth of the tree
        ranges.append((start, end))

    # Return found ranges
    return ranges


def _test_huge():
    qrybox = ndbox([1066052,1642769,1899], [1083529,1677722,2057])
    hquery(qrybox, maxdepth=16)


if __name__ == "__main__":
    _test_huge()
#    from pprint import pprint
#    pprint(hquery(query=ndbox([1, 1], [2, 2])))
