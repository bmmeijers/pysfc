from math import ceil, log
from collections import deque
from operator import itemgetter
from pysfc.relate import ndbox, relate

#
# Conversion of nD coordinate to SFC key and vice versa
#
# * Route for Morton:
#
# coord <> nchunks <> key
#
#
# * Route for Hilbert:
#
# coord <> nchunks <> hchunks <> key
#
#
# where:
#
# coord -- nD coordinate with len() == ndims 
# nchunks -- list of integers, where each int is location of cell to descend 
#            (following the Morton N-order)
# hchunks -- list of integers, where each int is location of cell to descend 
#            (following the Hilbert order) 
#            *Note*, this depends on depth of tree (mbits)
# key -- Space Filling Curve key 
#        (1D value, steps on the SFC from the start)
#

def _determine_bits(val, base):
    # gives exponent for next power of 2
    # next = pow(2, ceil(log(x)/log(2)));
    return max(1, int(ceil(log(val+1, base))))


def _transpose_bits(srcs, nDests):
    srcs = list(srcs)                   # Make a copy we can modify safely.
    nSrcs = len(srcs)
    dests = [0] * nDests
    for j in range(nDests - 1, -1, -1): # (to first deal with y, then with x)
        dest = 0
        for k in range(nSrcs):
            dest = dest * 2 + srcs[k] % 2
            srcs[k] //= 2 # divide by two, i.e. shift right (>>) with 1
        dests[j] = dest
    return tuple(dests)


def _chunks_key(chunks, mbits, ndims):
    key = 0
    nd = 2**ndims
    for m in range(mbits):
        if m >= len(chunks):
            chunk = 0
        else:
            chunk = chunks[m]
        key = key * nd + chunk
    return key

# equivalent code:
#def pack_index(chunks, nD):
#    p = 2**nD  # Turn digits mod 2**nD back into a single number:
#    return reduce(lambda n, chunk: n * p + chunk, chunks)


_coord_nchunks = _transpose_bits # args: nDests = mbits


_nchunks_coord = _transpose_bits # args: nDests = ndims


_nchunks_key = _chunks_key
# equivalent code:
#def _nchunks_key(nchunks, mbits, ndims):
#    diff = mbits - 1 # mbits gives depth of the tree, starting from 1 (with the level that has nary boxes)
#    key = 0
#    for nchunk in nchunks:
#        width = (diff * ndims)
#        key += nchunk << width
#        diff -= 1    # we move 1 closer to the bottom of the tree
#    return key


def _key_nchunks(key, mbits, ndims):
    """Convert Morton key into nchunks list

    For example:

        >>> _key_nchunks(45, 3, 2)
        [2, 3, 1]

    """
    nchunks = [0] * mbits
    mask = (1 << ndims) - 1 # number with all bits set to 1, with ndims bits, e.g. ndims = 3 == 0b111
    for i in range(mbits):
        nchunks[mbits - i - 1] = ((key >> (i * ndims)) & mask)
    return nchunks


# -- Hilbert specifics: for rotation and mirroring of the curve pattern --------


def _hgray_encode(bn):
    # Gray encoder and decoder from http://en.wikipedia.org/wiki/Gray_code
    assert bn >= 0
    assert type(bn) in [int], "Expected int, but found: {}".format(type(bn))
    return bn ^ (bn // 2)


def _hgray_decode(n):
    assert type(n) in [int], "Expected int, but found: {}".format(type(n))
    sh = 1
    while True:
        div = n >> sh
        n ^= div
        if div <= 1:
            return n
        sh <<= 1


def _hgray_encode_travel(start, end, mask, i):
    assert type(start) in [int], "Expected int, but found: {}".format(type(start))
    assert type(end) in [int], "Expected int, but found: {}".format(type(end))
    assert type(mask) in [int], "Expected int, but found: {}".format(type(mask))
    assert type(i) in [int], "Expected int, but found: {}".format(type(i))

    # gray_encode_travel -- gray_encode given start and end using bit rotation.
    #    Modified Gray code.  mask is 2**nbits - 1, the highest i value, so
    #        gray_encode_travel( start, end, mask, 0 )    == start
    #        gray_encode_travel( start, end, mask, mask ) == end
    #        with a Gray-code-like walk in between.
    #    This method takes the canonical Gray code,
    #    rotates the output word bits,
    #    then xors ("^" in Python) with the start value.
    travel_bit = start ^ end
    modulus = mask + 1          # == 2**nBits
    # travel_bit = 2**p, the bit we want to travel.
    # Canonical Gray code travels the top bit, 2**(nBits-1).
    # So we need to rotate by ( p - (nBits-1) ) == (p + 1) mod nBits.
    # We rotate by multiplying and dividing by powers of two:
    g = _hgray_encode(i) * (travel_bit * 2)
    return ((g | (g // modulus)) & mask) ^ start


def _hgray_decode_travel(start, end, mask, g):
    travel_bit = start ^ end
    modulus = mask + 1          # == 2**nBits
    rg = (g ^ start) * (modulus // (travel_bit * 2))
    result = _hgray_decode((rg | (rg // modulus)) & mask)
    return result


def _hchild_start_end(parent_start, parent_end, mask, i):
    # child_start_end( parent_start, parent_end, mask, i ) 
    # -- Get start & end for child.
    #    i is the parent's step number, between 0 and mask.
    #    Say that parent( i ) =
    #           gray_encode_travel( parent_start, parent_end, mask, i ).
    #    And child_start(i) and child_end(i) are what child_start_end()
    #    should return -- the corners the child should travel between
    #    while the parent is in this quadrant or child cube.
    #      o  child_start( 0 ) == parent( 0 )       (start in a corner)
    #      o  child_end( mask ) == parent( mask )   (end in a corner)
    #      o  child_end(i) - child_start(i+1) == parent(i+1) - parent(i)
    #         (when parent bit flips, same bit of child flips the opposite way)
    # Those constraints still leave choices when nD (# of bits in mask) > 2.
    #    Here is how we resolve them when nD == 3 (mask == 111 binary),
    #    for parent_start = 000 and parent_end = 100 (canonical Gray code):
    #         i   parent(i)    child_
    #         0     000        000   start(0)    = parent(0)
    #                          001   end(0)                   = parent(1)
    #                 ^ (flip)   v
    #         1     001        000   start(1)    = parent(0)
    #                          010   end(1)                   = parent(3)
    #                ^          v
    #         2     011        000   start(2)    = parent(0)
    #                          010   end(2)                   = parent(3)
    #                 v          ^
    #         3     010        011   start(3)    = parent(2)
    #                          111   end(3)                   = parent(5)
    #               ^          v
    #         4     110        011   start(4)    = parent(2)
    #                          111   end(4)                   = parent(5)
    #                 ^          v
    #         5     111        110   start(5)    = parent(4)
    #                          100   end(5)                   = parent(7)
    #                v          ^
    #         6     101        110   start(6)    = parent(4)
    #                          100   end(6)                   = parent(7)
    #                 v          ^
    #         7     100        101   start(7)    = parent(6)
    #                          100   end(7)                   = parent(7)
    #    This pattern relies on the fact that gray_encode_travel()
    #    always flips the same bit on the first, third, fifth, ...
    #    and last flip.
    #    The pattern works for any nD >= 1.
    start_i = max(0, (i - 1) & ~1)  # next lower even number, or 0
    end_i = min(mask, (i + 1) | 1)  # next higher odd number, or mask
    child_start = _hgray_encode_travel(parent_start, parent_end, mask, start_i)
    child_end = _hgray_encode_travel(parent_start, parent_end, mask, end_i)
    return child_start, child_end


def _hinitial_start_end(nChunks, nD):
    # This orients the largest cube so that
    # its start is the origin (0 corner), and
    # the first step is along the x axis, regardless of nD and nChunks:
    return 0, 2**((-nChunks - 1) % nD)  # in Python 0 <=   a % b   < b.


def _hchunks_to_nchunks(hchunks, mbits, ndims):
    hchunks_len = len(hchunks)
    # number with all bits set to 1, with ndims bits, 
    # e.g. ndims:
    # 2 == 0b11
    # 3 == 0b111
    # 4 == ...
    mask = 2 ** ndims - 1 
    start, end = _hinitial_start_end(mbits, ndims)
    nchunks = [0] * hchunks_len
    for j, hchunk in enumerate(hchunks):
        nchunks[j] = _hgray_encode_travel(start, end, mask, hchunk)
        start, end = _hchild_start_end(start, end, mask, hchunk)
    return tuple(nchunks)


def _nchunks_to_hchunks(nchunks, mbits, ndims):
    nchunks_len = len(nchunks)
    mask = 2 ** ndims - 1
    start, end = _hinitial_start_end(mbits, ndims)
    hchunks = [0] * nchunks_len
    for j, nchunk in enumerate(nchunks):
        hchunks[j] = _hgray_decode_travel(start, end, mask, nchunk)
        start, end = _hchild_start_end(start, end, mask, hchunks[j])
    return tuple(hchunks)


_hchunks_key = _chunks_key


def _key_hchunks(key, mbits, ndims):
    # pre conditions
    assert type(key) in [int], "Expected int, but found: {}".format(type(key))
    assert type(mbits) in [int], "Expected int, but found: {}".format(type(mbits))
    assert type(ndims) in [int], "Expected int, but found: {}".format(type(ndims))

    p = 2**ndims
    hchunks = [0] * mbits
    for j in range(mbits - 1, -1, -1):
        hchunks[j] = key % p
        key //= p

    # post condition: all int's
    for h in hchunks:
        assert type(h) in [int], "Expected int, but found: {}".format(type(p))

    return hchunks


# Public api, encode and decode, n-order (nenc+ndec) and hilbert (henc + hdec)


# -- Morton (N-order)
def nenc(coord):
    mbits = _determine_bits(max(coord), 2)
    ndims = len(coord)
    #
    nchunks = _coord_nchunks(coord, mbits)
    key = _hchunks_key(nchunks, mbits, ndims)
    return key


def ndec(key, ndims):
    mbits = _determine_bits(key, 2**ndims)
    #
    nchunks = _key_nchunks(key, mbits, ndims)
    coord = _nchunks_coord(nchunks, ndims)
    return tuple(coord)


# -- Hilbert
def henc(coord):
    mbits = _determine_bits(max(coord), 2)
    ndims = len(coord)
    #
    nchunks = _coord_nchunks(coord, mbits)
    hchunks = _nchunks_to_hchunks(nchunks, mbits, ndims)
    key = _hchunks_key(hchunks, mbits, ndims)
    return key


def hdec(key, ndims):
    mbits = _determine_bits(key, 2**ndims)
    #
    assert type(mbits) in [int]
    assert type(key) in [int]
    assert type(ndims) in [int]
    hchunks = _key_hchunks(key, mbits, ndims)

    nchunks = _hchunks_to_nchunks(hchunks, mbits, ndims)
    coord = _nchunks_coord(nchunks, ndims)
    return tuple(coord)


def nquery(query):
#    maxbits = 63
    ndims = query.dims
    # get how many bits we need
    # to represent the largest number inside the query box
    # --> 2**(mbits_needed) is the maximum size of a side 
    #     of the domain that we need
    mbits_needed = _determine_bits(max(query.hi), 2)
#    mbits = maxbits // ndims
    npath = ()
    # post order tree traversal gives nodes in order we want
    # for this, two stacks are used (stack + paths)
    # -- https://algorithms.tutorialhorizon.com/binary-tree-postorder-traversal-non-recursive-approach/
    paths = []
    stack = [npath]
    while stack:
        npath = stack.pop() 
        cur_level = len(npath)
        lo = _nchunks_coord(npath, ndims)
        # side_at_depth 
        # how large side of the nd-cube is at this depth
        side_at_depth = 2**(mbits_needed - cur_level)
        lo = tuple(map(lambda x: int(x * side_at_depth), lo))
        hi = tuple(map(lambda x: int(x + side_at_depth), lo))
        cur_node = ndbox(lo, hi)
        ndcmp = relate(query, cur_node)
        if ndcmp in (0, 1,):
            paths.append(npath)
        # -- partial overlap
        elif ndcmp in (2,):
            # FIXME: re-introduce maxdepth argument !
            if cur_level < mbits_needed:
                # we have not yet reached the lowest level, recurse
                for ncode in range(2**ndims):
                    new_path = npath + (ncode, )
                    stack.append(new_path)
            else:
                # we are not allowed to go further, so use this partial
                # overlapping range as is
                paths.append(npath)
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp in (-1,):
            continue
    result = []
    while paths:
        npath = paths.pop()
        # if the npath has less than mbits_needed values,
        # we need to add 0's to the end of the list:
        # -> this happens inside _nchunks_key
        start = _nchunks_key(npath, mbits_needed, ndims)
        side_at_depth = 2 ** (mbits_needed - len(npath))
        range_size = side_at_depth**ndims
        end = start + range_size
#        print npath, start, end, range_size
        result.append((start, end))
    return result


def hquery(query):
#    maxbits = 63
    ndims = query.dims
    # get how many bits we need
    # to represent the largest number inside the query box
    # --> 2**(mbits_needed) is the maximum size of a side 
    #     of the domain that we need
    mbits_needed = _determine_bits(max(query.hi), 2)
#    mbits = maxbits // ndims
    npath = ()
    # post order tree traversal gives nodes in order we want
    # for this, two stacks are used (stack + paths)
    stack = [npath]
    paths = []
    while stack:
        npath = stack.pop() 
#        npath = stack.popleft()
        cur_level = len(npath)
        lo = _nchunks_coord(npath, ndims)
        # side_size_at_depth how large is a side of the cube at this depth?
        side_at_depth = 2 ** (mbits_needed - cur_level)
        lo = tuple(map(lambda x: int(x * side_at_depth), lo))
        hi = tuple(map(lambda x: int(x + side_at_depth), lo))
        cur_node = ndbox(lo, hi)
        ndcmp = relate(query, cur_node)
        if ndcmp in (0, 1,):
            paths.append(npath)
        # -- partial overlap
        elif ndcmp in (2,):
            # FIXME: re-introduce maxdepth argument !
            if cur_level < mbits_needed:
                # we have not yet reached the lowest level, recurse
                childs = []
                for ncode in range(2**ndims):
                    new_path = npath + (ncode, )
                    hpath = _nchunks_to_hchunks(new_path, mbits_needed, ndims)
                    hcode = _hchunks_key(hpath, mbits_needed, ndims)
                    childs.append((hcode, new_path))
                childs.sort(key=itemgetter(0))
                for child in childs:
                    stack.append(child[1])
            else:
                # we are not allowed to go further, so use this partial
                # overlapping range as is
                paths.append(npath)
        # -- disjoint, we do not need to process further this path down the tree
        elif ndcmp in (-1,):
            continue
    #
    result = []
    while paths:
        npath = paths.pop()
        hpath = _nchunks_to_hchunks(npath, mbits_needed, ndims)
        start = _hchunks_key(hpath, mbits_needed, ndims)
        side_at_depth = 2**(mbits_needed - len(npath))
        range_size = side_at_depth**ndims
        end = start + range_size
        result.append((start, end))
    #
    return result



#if __name__ == "__main__":
#    pass
#    ndims = 2
#    for i in range(4):
#        for j in range(4):
#            print [i, j], _nchunks_to_hchunks([i, j], ndims)

#    print hquery(ndbox((1, 1, 1, 1), (3, 3, 3, 3)))

#    print hquery(query=ndbox([0, 2], [2, 4]))
#    print nchunk_table()

#    print _nchunks_key((1,3,0), 3, 2)

#    import doctest
#    doctest.testmod()
#    test()

#    hchunk_table()
#    nchunk_table()

#    nquery(query=ndbox([0, 0], [8, 8]))
#    nquery(query=ndbox([1, 1], [4, 4]))
#    nquery(query=ndbox([1, 1], [3, 6]))
#    nquery(query=ndbox([1, 1], [3, 3]))

#    print "***", _nchunks_hchunks((1,0,), 2)
#    print "***", _hchunks_key((1,0,), 2, 2)

##    for i in range(64):
##        print i, hdec(i, 2)
#    hquery_new(query=ndbox([0, 0], [1, 1]))
#    hquery_new(query=ndbox([1, 1], [4, 4]))
#    hquery_new(query=ndbox([0, 1], [2, 4]))
#    hquery(query=ndbox([1, 1], [4, 4]))
#    hquery(query=ndbox([1, 1], [3, 6]))
#    hquery(query=ndbox([1, 1], [3, 3]))
