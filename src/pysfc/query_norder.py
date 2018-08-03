from relate import ndbox, relate

from collections import deque

from traversal import path_to_morton
from traversal import nary

from morton_norder import transpose_bits


def nquery(query, maxdepth=None, maxbits=63):
    """Query with N-order curve

    Arguments:

        query:
        maxdepth:
        maxbits:

    """

    dim = query.dims        # the number of dimensions of the top box
    order = maxbits / dim   # order := how many levels we want in the full tree
                            # (this depends on max number of bits that fit in the type in-use for representing the ranges + the number of dimensions used)
                            # FIXME: the numbering of levels starts at 1, however,
                            # for searching that leads to a strange condition that
                            # needs order + 1
    mbits = maxbits / dim

    # -- determine how deep to go at maximum
    if maxdepth is None:
        maxdepth = order + 1
    maxdepth = min(maxdepth + 1, order + 1)

    # min and max, min and max coordinate values as list
    # each with length of the dimension
    l_root = [0] * dim
    h_root = [2**order] * dim
    # What we place in the stack:
    # (N-order code of the node at this level -- note: not used,
    #  node its box,
    #  level of the node in the tree,
    #  path traversed -- N-codes of the parents)
    traverse_child_order = range(nary(dim))
    #root = (None, 1, ())
    root = (None, ndbox(l_root, h_root), 1, ())
    deq = deque([root])
    result = []
    while deq:
        #
        # _, cur_level, cur_path = deq.popleft()
        _, cur_node, cur_level, cur_path = deq.popleft()



        # derive cur_node from cur_path (path with nchunk's)
        # determine the size of a cell at this depth
        _path = [0] * mbits
        i = 0
        v = 0
        for i, v in enumerate(cur_path):
            _path[i] = v
        lo = transpose_bits(_path, dim)
        #
        range_width_at_depth = 2**(mbits - len(cur_path))
        hi = [_ + range_width_at_depth for _ in lo]
        new_cur_node = ndbox(lo, hi)
        del lo
        del hi
        del i
        del v
        del _path
        # 

#        print cur_node
#        print new_cur_node, cur_path
#        print ""

        # compare the query geometry to the node
        # derive the relationship based on the geometry
        # comparison can be: -1, 0, 1, 2
        ndcmp = relate(query, cur_node)

        # handle the relation between the query geometry and the node geometry
        # -- equal or fully contains, done: add to result
        if ndcmp in (0, 1,):
            result.append(cur_path)
                          # FIXME: could use 'yield' here
                          # -> generator would be nice to have as downstream caller
        # -- partial overlap
        elif ndcmp in (2,):
            if cur_level < maxdepth:
                # we have not yet reached the lowest level, recurse
                for ncode in traverse_child_order:
                    new_path = cur_path + (ncode, )
                    deq.append(
                        (ncode, cur_node.nth_child_norder(ncode), cur_level + 1, new_path))
            else:
                # we are not allowed to go further, so use this partial
                # overlapping range as is
                result.append(cur_path)
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp in (-1,):
            continue

    # Transform paths into ranges
    ranges = []
    for path in result:
        rng = path_to_morton(path, dim, order)
        ranges.append(rng)

    return ranges


if __name__ == "__main__":
    print nquery(query=ndbox([0, 0], [8, 8]))
    print nquery(query=ndbox([1, 1], [4, 4]))
    print nquery(query=ndbox([1, 1], [3, 6]))
    print nquery(query=ndbox([1, 1], [3, 3]))
