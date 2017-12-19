from math import ceil, log



def decode(val, ndims=2):
    """Decode a N-order Morton code as a nD coordinate"""
    # make from the Morton code a list that is the path in the nary tree
    # with at each level the index to which val overlaps
    p = 2**ndims     # Chunks are like digits in base 2**nD.
    mbits = max(1, int(ceil(log(val + 1, p))))  # nr of bits we need
    path = mortoncode_to_path(val, mbits, ndims)
#    print 'd:', path
    # convert the path into nd Coordinate
    result = tuple(transpose_bits(path, ndims))
    return result


def encode(coords):
    """Encode a nD coordinate as N-order Morton code"""
    nD = len(coords)
    # biggest = reduce( max, coords )  # the max of all coords
    biggest = max(coords)
    nChunks = max(1, int(ceil(log(biggest + 1, 2))))  # max nr of bits needed
    path = transpose_bits(coords, nChunks)
#    print ""
#    print 'en:', path
    result = path_to_mortoncode(path, nChunks, nD)
#    print 'result', result
    return result


def path_to_mortoncode(path, mbits, ndims):
    """Convert path in nary-tree, top down, each with morton index
    for that level to a complete Morton code.

    The max index in each step of the path depends on the number of dims:

        2D: (2 ** 2) - 1 = 3
        3D: (2 ** 3) - 1 = 7
        4D: (2 ** 4) - 1 = 15
        etc.

    In 2D indexing the quads like this:

        1   3
        | \ |
        0   2

    For example:

        >>> path_to_mortoncode([2, 3, 1], 3, 2)
        45

    """
    # transform the path with indexes to a number
    # FIXME: is mbits argument needed ??? Normally it is the length of the path (so known)...
    diff = mbits - 1 # mbits gives depth of the tree, starting from 1 (with the level that has nary boxes)
    morton = 0
    for zorder in path:
        width = (diff * ndims)
        morton += zorder << width
        diff -= 1    # we move 1 closer to the bottom of the tree
    return morton


def mortoncode_to_path(code, mbits, ndims): # FIXME: do we need mbits as argument?
    """Convert complete Morton code into path (list with each morton index)

    For example:

        >>> mortoncode_to_path(45, 3, 2)
        [2, 3, 1]

    """
    path = [0] * mbits
    mask = (1 << ndims) - 1
    for i in range(mbits):
        path[mbits - i - 1] = ((code >> (i*ndims)) & mask)
    return path


def transpose_bits(srcs, nDests):
    srcs = list(srcs)  # Make a copy we can modify safely.
    nSrcs = len(srcs)
    dests = [0] * nDests
    # Break srcs down most-significant value first, shifting up:
    for j in range(nDests - 1, -1, -1): # (to first deal with y, then with x)
        dest = 0
        for k in xrange(nSrcs): # FIXME: Check is this same code as in HILBERT.PY 
            dest = dest * 2 + srcs[k] % 2
            srcs[k] /= 2 # divide by two, i.e. shift right (>>) with 1
        dests[j] = dest
    return dests


def _test_lengthy(maxside):
    """Test that encode and decode functions agree with each other,
    also in higher dims"""
    # 2D
    for x in xrange(maxside):
        for y in xrange(maxside):
            c = (x, y)
            e = encode(c)
            d = decode(e, 2)
            assert c == d
    # 3D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                c = (x, y, z)
                e = encode(c)
                d = decode(e, 3)
                assert c == d
    # 4D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                for t in xrange(maxside):
                    c = (x, y, z, t)
                    e = encode(c)
                    d = decode(e, 4)
                assert c == d
    # 5D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                for t in xrange(maxside):
                    for s in xrange(maxside):
                        c = (x, y, z, t, s)
                        e = encode(c)
                        d = decode(e, 5)
                    assert c == d

def _test_small():
    # test first 9 steps of the n-order curve in 2D
    tests = [
        ((0, 0), 0),
        ((0, 1), 1),
        ((1, 0), 2),
        ((1, 1), 3), 
        ((0, 2), 4),
        ((0, 3), 5),
        ((1, 2), 6),
        ((1, 3), 7),
        ((2, 0), 8)]
    for c, outcome in tests:
        assert encode(c) == outcome
        assert c == decode(outcome, 2)

if __name__ == "__main__":
#    print transpose_bits((7, 2), 100)
#    print transpose_bits((2, 7), 100)
#    print mortoncode_to_path(45, 3, 2)

# n order ==
#(0, 0) 	0
#(0, 1) 	1
#(1, 0) 	2
#(1, 1) 	3
#
# 5  7
# |\ |\
# | \| 
# 4  6 \
#  \    
#   \   \
# 1  3    
# |\ |   \ 
# | \|    
# 0  2    8


    import doctest
    doctest.testmod()

    _test_small()
    _test_lengthy(maxside = 2**4)
