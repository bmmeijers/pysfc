#cython: infer_types=True, boundscheck=False, profile=True

from pysfc.speedups.relate import ndbox, relate
from collections import deque

#from pysfc.traversal import nth_child_coordinates_norder
from pysfc.traversal import morton_zcode_ncode
#from pysfc.traversal import path_to_hilbert

from operator import itemgetter

import cython

cdef pack_index(tuple chunks, int nD):
    if chunks:
        p = 2 ** nD
        tally = chunks[0]
        for nxt in chunks[1:]:
            tally = tally * p + nxt
        return tally
    else:
        return 0

cdef path_to_hilbert(tuple path, int nD, int order):
    """Given a traversed path in the tree (with integers that represent h-codes),
    determine the hilbert code that starts in this box of the path
    """
    level = len(path)
    # width: the size of the range belonging to this node at this level
    width = 2**(nD*(order-level))
    hcode = pack_index(path, nD) * width
    return hcode, hcode + width 

cdef int initial_end(int nChunks, int nD):
    # This orients the largest cube so that
    # its start is the origin (0 corner), and
    # the first step is along the x axis, regardless of nD and nChunks:
    return 2**((-nChunks - 1) % nD)  # in Python 0 <=   a % b   < b.

cdef inline int gray_decode(int n):
    sh = 1
    while True:
        div = n >> sh
        n ^= div
        if div <= 1:
            return n
        sh <<= 1

@cython.cdivision(True)
cdef inline int gray_decode_travel(int start, int end, int mask, int g):
    cdef int travel_bit = start ^ end
    cdef int modulus = mask + 1          # == 2**nBits
    cdef int rg = (g ^ start) * (modulus / (travel_bit * 2))
#    print 'rg', rg
    cdef int result = gray_decode((rg | (rg / modulus)) & mask)
#    print 'rgresult', result
    return result


cdef inline int gray_encode(int bn):
    # Gray encoder and decoder from http://en.wikipedia.org/wiki/Gray_code
#    assert bn >= 0
#    assert type(bn) in [int, long]
    return bn ^ (bn / 2)


@cython.cdivision(True)
cdef inline int gray_encode_travel(int start, int  end, int mask, int i):
    cdef int travel_bit = start ^ end
    cdef int modulus = mask + 1          # == 2**nBits
    # travel_bit = 2**p, the bit we want to travel.
    # Canonical Gray code travels the top bit, 2**(nBits-1).
    # So we need to rotate by ( p - (nBits-1) ) == (p + 1) mod nBits.
    # We rotate by multiplying and dividing by powers of two:
    cdef int g = gray_encode(i) * (travel_bit * 2)
    return ((g | (g / modulus)) & mask) ^ start


cdef int  child_start(int parent_start, int parent_end, int mask, int i):
    cdef int start_i = max(0, (i - 1) & ~1)  # next lower even number, or 0
    cdef int child_start = gray_encode_travel(parent_start, parent_end, mask, start_i)
    return child_start


cdef int child_end(int parent_start, int parent_end, int mask, int i):
    cdef int end_i = min(mask, (i + 1) | 1)  # next higher odd number, or mask
    cdef int child_end = gray_encode_travel(parent_start, parent_end, mask, end_i)
    return child_end




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
    hstart = 0
    hend = initial_end(order, dim)
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
                    hcode = gray_decode_travel(
                        cur_hstart, cur_hend, hmask, ncode)
                    next_hpath = cur_hpath + (hcode, )
                    next_box = cur_box.nth_child_norder(ncode)
                    next_start = child_start(cur_hstart, cur_hend, hmask, hcode)
                    next_end = child_end(cur_hstart, cur_hend, hmask, hcode)
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
            pass

    # Transform paths into ranges
    ranges = []
    for path in result:
        start, end = path_to_hilbert(
            path, dim, order)  # << IMPORTANT: order = total depth of the tree
        ranges.append((start, end))

    return ranges


def _test_huge():
    qrybox = ndbox([1066052,1642769,1899], [1083529,1677722,2057])
    hquery(qrybox, maxdepth=16)


if __name__ == "__main__":
    _test_huge()
#    from pprint import pprint
#    pprint(hquery(query=ndbox([1, 1], [2, 2])))
