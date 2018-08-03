#from pysfc.speedups.relate import ndbox, relate
from relate import ndbox, relate

#from traversal import morton_zcode_ncode
from traversal import path_to_hilbert

from hilbert import initial_start_end, gray_decode_travel, gray_encode_travel, child_start_end

from morton_norder import transpose_bits

from operator import itemgetter

from math import ceil, log

def hchunks_to_nchunks(hchunks, nD):
    hchunks_len = len(hchunks)
    mask = 2 ** nD - 1
    start, end = initial_start_end(hchunks_len, nD)
    nchunks = [0] * hchunks_len
    for j in range(hchunks_len):
        i = hchunks[j]
        nchunks[j] = gray_encode_travel(start, end, mask, i)
        start, end = child_start_end(start, end, mask, i)
    return nchunks


def nchunks_to_ndbox(nchunks, mbits, ndims):
    # to low coordinate
    _path = [0] * mbits
    i = 0
    v = 0
    for i, v in enumerate(nchunks):
        _path[i] = v
    lo = transpose_bits(_path, ndims)

    # to ndbox from low coordinate + width
    range_width_at_depth = 2**(mbits - len(nchunks))
    hi = [_ + range_width_at_depth for _ in lo]
    cur_node = ndbox(lo, hi)
    # end transfrom
    return cur_node


def hquery(query, maxdepth=None, maxbits=63):
    """
    Traverse n-ary tree, each node its childs are visited in H-order
    """

#    print query
#    print "*" * 88

    dim = query.dims            # the number of dimensions of the top box
    hmask = 2 ** dim - 1        # a number of dim bits, with every bit set to 1

    # >> # FIXME: some options:
    # - from metadata table, or
    # - given by user (specifying max depth to visit), or
    # - inferred from largest coordinate inside the query box?
    # how many levels we want max to visit in the tree
    order = maxbits / dim       # -- influences maximum size of the domain
    mbits = maxbits / dim
#    print "mbits mod 2**dim", mbits % 2**dim
    # <<

    # -- determine how deep to go at maximum
    if maxdepth is None:
        maxdepth = order + 1
    maxdepth = min(maxdepth + 1, order + 1)

#    order = max(1, int(ceil(log(max(query.hi)+1, 2**dim))))

    # min and max, min and max coordinate values as list
    # each has a length equal to the dimension
    l_root = [0] * dim
    h_root = [2 ** order] * dim
    root_box = ndbox(l_root, h_root)  # Here we determine the size of the cube

    # >> dependent on order (depth of the nary-tree) !!!
    # (Note, otherwise curve will be running in wrong direction)
    hstart, hend = initial_start_end(1, dim) # order = 1, 5, 9, 13, ... (these all work!)
#    print hstart, hend
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
#    traverse_childs_order = [ncode for _, ncode in morton_zcode_ncode(dim, False)]
#    print "child order", traverse_childs_order
    root = (root_box, (hstart, hend), 1, None, (), None, ())
    stack = [root]
    while stack:
        node = stack.pop()
        cur_box, (cur_hstart, cur_hend), cur_level, cur_hcode, cur_hpath, cur_ncode, cur_npath = node

#        lo = _nchunks_coord(cur_npath, ndims)

        # transform hpath into path with nchunks
        nchunks = hchunks_to_nchunks(cur_hpath, dim)
#        print nchunks
        new_cur_box = nchunks_to_ndbox(cur_npath, mbits, dim)
#        assert cur_box.lo == new_cur_box.lo
#        assert cur_box.hi == new_cur_box.hi

#        print new_cur_box, cur_box, cur_hcode, cur_hpath, cur_ncode, cur_npath

#        print cur_box, nchunks, cur_hpath, new_cur_box
#        print ""
        # from cur_hpath obtain cur_box <--> to remove cur_box from stackue

        # determine the relation of the query box and the current node box
        ndcmp = relate(query, new_cur_box)
        # -- equal or fully contains, done: add to result
        if ndcmp in (0, 1,):
#            print new_cur_box, cur_npath, cur_hpath
            result.append((cur_npath, cur_hpath))
        # -- partial overlap
        elif ndcmp in (2,):
            if cur_level <= maxdepth:
                # we have not yet reached the lowest level, recurse
                childs = []
                next_level = cur_level + 1
#                for ncode in traverse_childs_order:
                for ncode in range(dim**2):
                # for _, ncode in morton_zcode_ncode(dim, False):
                    hcode = gray_decode_travel(
                        cur_hstart, cur_hend, hmask, ncode)
#                    print new_cur_box, 'n:', ncode, (cur_hstart, cur_hend), 'h:', hcode, "[", hmask, "]", cur_hpath, cur_npath
                    next_hpath = cur_hpath + (hcode, )
                    next_npath = cur_npath + (ncode, )
#                    next_box = cur_box.nth_child_norder(ncode)
                    next_start, next_end = child_start_end(
                        cur_hstart, cur_hend, hmask, hcode)
                    childs.append(
                        (None, (next_start, next_end), next_level, hcode, next_hpath, ncode, next_npath))
                # -- visit the childs in order of the hilbert curve
                # (we need to sort the childs list for this)
                # !! in the tuple 0 = hcode, sort on local hilbert index (hcode)
                childs.sort(key=itemgetter(3))
                for child in childs:
#                    print child
                    stack.append(child)
            else:
                print new_cur_box, cur_npath, cur_hpath
                result.append((cur_npath, cur_hpath))
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp in (-1,):
            continue

    # Transform paths into ranges
    ranges = []
    while result:
        npath, hpath = result.pop()
        start, end = path_to_hilbert(hpath, dim, order)
        print npath, hpath, "->", start, "order:", order
            # << IMPORTANT: order = total depth of the tree
        ranges.append((start, end))

    # Return found ranges
    return ranges


def _test_huge():
    qrybox = ndbox([1066052, 1642769, 1899], [1083529, 1677722, 2057])
    hquery(qrybox, maxdepth=16)


if __name__ == "__main__":

#    print hquery(query=ndbox([1, 1],
#                             [4, 4]), maxdepth=None, maxbits = 63)


#    print hquery(query=ndbox([0, 1], [2, 4]), maxdepth=None, maxbits = 63)

#    print hquery(query=ndbox([2, 0], [4, 3]), maxdepth=None, maxbits = 4)

    print "range to obtain := 4,8"
    for i in range(7):
        print ""
        print "maxbits := ", i
        print hquery(query=ndbox([0, 2], [2, 4]), maxdepth=None, maxbits = i)

    print "-----"

    print "range to obtain := 12,16"
    for i in range(7):
        print ""
        print "maxbits := ", i
        print hquery(query=ndbox([2,0], [4, 2]), maxdepth=None, maxbits = i)


#    print hquery(query=ndbox([7,9,15], [8,10,16]), maxdepth=None)
#    _test_huge()
#    from pprint import pprint
#    pprint(hquery(query=ndbox([1, 1], [2, 2])))
