from math import ceil, log


#def encode_2d_reference(ptcoord):
#    """Encode a 2D coordinate as Z-order Morton code"""
#    ndims = len(ptcoord)
#    biggest = max(ptcoord)
#    mbits = max( 1, int( ceil( log( biggest + 1, 2 ) ) ) )
#    result = 0
#    for i in range(ndims*mbits):
#        result |= (ptcoord[0] & 1L << i) << i | (ptcoord[1] & 1L << i) << (i + 1)
#    return result


# FIXME: Do we need the two versions of transpose_bits / _rev ???


def decode(val, ndims=2):
    """Decode a Z-order Morton code as a nD coordinate"""
    # make from the Morton code a list that is the path in the nary tree
    # with at each level the index to which val overlaps
    p = 2**ndims     # Chunks are like digits in base 2**nD.
    mbits = max(1, int(ceil(log(val + 1, p))))  # nr of bits we need
    path = mortoncode_to_path(val, mbits, ndims)
    # convert the path into nd Coordinate
    result = tuple(transpose_bits(path, ndims))
    return result


def encode(coords):
    """Encode a nD coordinate as Z-order Morton code"""
    nD = len(coords)
    # biggest = reduce( max, coords )  # the max of all coords
    biggest = max(coords)
    nChunks = max(1, int(ceil(log(biggest + 1, 2))))  # max # of bits needed
    path = transpose_bits_rev(coords, nChunks)
    result = path_to_mortoncode(path, nChunks, nD)
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

        2 - 3
          \  
        0 - 1

    E.g. a path [2, 3, 1] means that the path indicates:

        - top level, overlap with top left
        - middle level, overlap with top right
        - lowest level, overlap with bottom right

    For example:

        >>> path_to_mortoncode([2, 3, 1], 3, 2)
        45

    """
    # transform the path with indexes to a number
    # FIXME: is mbits argument needed ??? Normally it is the length of the path...
    diff = mbits - 1 # mbits gives depth of the tree, starting from 1 (with the level that has nary boxes)
    morton = 0
    for zorder in path:
        width = (diff * ndims)
        morton += zorder << width
        diff -= 1    # we move 1 closer to the bottom of the tree
    return morton


def mortoncode_to_path(code, mbits, ndims):
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


def transpose_bits_rev(srcs, nDests):
    srcs = list(srcs)  # Make a copy we can modify safely.
    nSrcs = len(srcs)
    dests = [0] * nDests
    # Break srcs down most-significant value first, shifting up:
    for j in range(nDests - 1, -1, -1):
        dest = 0
        for k in reversed(xrange(nSrcs)): # <<< THIS LOOP IS REVERSED FROM CODE IN HILBERT.PY (to first deal with y, then with x)
            dest = dest * 2 + srcs[k] % 2
            srcs[k] /= 2 # divide by two, i.e. shift right (>>) with 1
        dests[j] = dest
    return dests


def transpose_bits(srcs, nDests):
    srcs = list(srcs)  # Make a copy we can modify safely.
    nSrcs = len(srcs)
    dests = [0] * nDests
    for j in xrange(nDests): # <<< THIS LOOP IS REVERSED FROM CODE IN HILBERT.PY (to first deal with x, then with y)
        # Put dests together most-significant first, shifting up:
        dest = 0
        for k in xrange(nSrcs):
            dest = dest * 2 + srcs[k] % 2
            srcs[k] /= 2 # divide by two, i.e. shift right (>>) with 1
        dests[j] = dest
    return dests


def _test_encode():
    assert encode((8, 8, 8)) == 3584
    assert encode((0, 0)) == 0
    assert encode((2, 2)) == 12
    assert encode((3, 3)) == 15
    assert encode((2, 7)) == 46
    assert encode((7, 7)) == 63


def _test_decode():
    assert decode(3584, 3) == (8,8,8)
    assert decode(0 , 2) == (0, 0)
    assert decode(12, 2) == (2, 2)
    assert decode(15, 2) == (3, 3)
    assert decode(46, 2) == (2, 7)
    assert decode(63, 2) == (7, 7)


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


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    _test_encode()
    _test_decode()
    _test_lengthy(maxside = 2**4)
