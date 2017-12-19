#!/usr/bin/env python
#cython: infer_types=True, boundscheck=False
# hilbert.py -- Hilbert walk coordinate in multiple dimensions.
#    decode( i, 3 ) ==> ( x, y, z )  ### decode
#        decode( 0, nD ) ==> ( 0, 0, 0, ... 0 ) Start at origin.
#        decode( 1, nD ) ==> ( 1, 0, 0, ... 0 ) 1st step is along x.
#
#    encode( ( x, y, z ) ) ==> i     ### encode
#
# Steve Witham ess doubleyou at tiac remove-this dot net.
# http://www.tiac.net/~sw/2008/10/Hilbert


from libc.math cimport log, ceil
#import array
#from cpython cimport array
import cython


def decode(object i, int nD=2):  # Default is the 2D Hilbert walk.
    """Decode a Hilbert key as nD coordinate"""
    cdef int j, nChunks, mask, start, end, parent_start, parent_end
    index_chunks = unpack_index(i, nD)
    nChunks = len(index_chunks)
    mask = 2 ** nD - 1
    start = 0
    end = initial_end(nChunks, nD)
    coord_chunks = [0] * nChunks
    for j in range(nChunks):
        i = index_chunks[j]
        coord_chunks[j] = gray_encode_travel(start, end, mask, i)
        parent_start = start
        parent_end = end
        start = child_start(parent_start, parent_end, mask, i)
        end = child_end(parent_start, parent_end, mask, i)
#    print "d: index chunks", index_chunks
#    print "d: coord chunks", coord_chunks
#    return tuple(pack_coords(coord_chunks, nD)), coord_chunks, index_chunks, val
    return tuple(pack_coords(coord_chunks, nD))


def encode(coords):
    """Encode a nD coordinate as Hilbert key"""
    cdef int i, j, nChunks, mask, start, end, parent_start, parent_end
    nD = len(coords)
    coord_chunks = unpack_coords(coords)
    nChunks = len(coord_chunks)
    mask = 2 ** nD - 1 # a number with all bits set to 1 and has as many bits as there are dimensions, so for 2D the mask is 11, for 3D: 111, for 4D: 1111, etc.
    start = 0 
    end = initial_end(nChunks, nD)
    index_chunks = [0] * nChunks
    for j in range(nChunks):
        i = gray_decode_travel(start, end, mask, coord_chunks[j]) # << coord_chunks[j] :=> NORDER VALUE
        index_chunks[j] = i
        parent_start = start
        parent_end = end
        start = child_start(parent_start, parent_end, mask, i)
        end = child_end(parent_start, parent_end, mask, i)
#        start, end = child_start_end(start, end, mask, i) # << i :=> CURRENT HINDEX 
#    print ""
#    print "e: coord chunks", coord_chunks # << in the n-ary tree this is the local N-ORDER identifier
#    print "e: index chunks", index_chunks # << in the n-ary tree this is the local Hilbert identifier 
    # (Note, it is *not* possible to us this index value to know which field you are in (as is the case with morton) how the curve was oriented previously has an influence on this...)!
#    return coords, coord_chunks, index_chunks, pack_index(index_chunks, nD)
    return pack_index(index_chunks, nD)

cdef int initial_end(int nChunks, int nD):
    # This orients the largest cube so that
    # its start is the origin (0 corner), and
    # the first step is along the x axis, regardless of nD and nChunks:
    return 2**((-nChunks - 1) % nD)  # in Python 0 <=   a % b   < b.


# Unpacking arguments and packing results of int <-> Hilbert functions.
# nD == # of dimensions.
# A "chunk" is an nD-bit int (or Python long, aka bignum).
# Lists of chunks are highest-order first.
# Bits within "coord chunks" are x highest-order, y next, etc.,
# i.e., the same order as coordinates input to Hilbert_to_int()
# and output from int_to_Hilbert().


# unpack_index( int index, nD ) --> list of index chunks.
cdef inline unpack_index(object i, int nD):
    cdef int p, j, nChunks
    p = 2**nD     # Chunks are like digits in base 2**nD.
    nChunks = max(1, int(ceil(log(i + 1)/log(p))))  # of digits
    chunks = [0] * nChunks
    for j in range(nChunks - 1, -1, -1):
        chunks[j] = i % p
        i = i / p
    return chunks


cdef inline pack_index(chunks, nD):
    p = 2**nD  # Turn digits mod 2**nD back into a single number:
    return reduce(lambda n, chunk: n * p + chunk, chunks)


# unpack_coords( list of nD coords ) --> list of coord chunks each nD bits.
cdef inline unpack_coords(coords):
    cdef int nChunks, nD
    nD = len(coords)
    # biggest = reduce( max, coords )  # the max of all coords
    biggest = max(coords)
    nChunks = max(1, int(ceil(log(biggest + 1)/ log(2))))  # max # of bits
    return transpose_bits(coords, nChunks)


cdef inline pack_coords(chunks, nD):
    return transpose_bits(chunks, nD)


# transpose_bits --
#    Given nSrcs source ints each nDests bits long,
#    return nDests ints each nSrcs bits long.
#    Like a matrix transpose where ints are rows and bits are columns.
#    Earlier srcs become higher bits in dests;
#    earlier dests come from higher bits of srcs.
cdef inline transpose_bits(object srcarg, int nDests):
    cdef int dest, k, j, nSrcs
    srcs = list(srcarg) # Make a copy we can modify safely.
    nSrcs = len(srcs)
    dests = [0] * nDests
    #dests = [0] * nDests
    # Break srcs down least-significant bit first, shifting down:
    for j in range(nDests - 1, -1, -1):
        # Put dests together most-significant first, shifting up:
        dest = 0
        for k in range(nSrcs):
            dest = dest * 2 + srcs[k] % 2
            srcs[k] = srcs[k] / 2
        dests[j] = dest
    return dests


cdef inline int gray_encode(int bn):
    # Gray encoder and decoder from http://en.wikipedia.org/wiki/Gray_code
#    assert bn >= 0
#    assert type(bn) in [int, long]
    return bn ^ (bn / 2)


cdef inline int gray_decode(int n):
    sh = 1
    while True:
        div = n >> sh
        n ^= div
        if div <= 1:
            return n
        sh <<= 1


# gray_encode_travel -- gray_encode given start and end using bit rotation.
#    Modified Gray code.  mask is 2**nbits - 1, the highest i value, so
#        gray_encode_travel( start, end, mask, 0 )    == start
#        gray_encode_travel( start, end, mask, mask ) == end
#        with a Gray-code-like walk in between.
#    This method takes the canonical Gray code, rotates the output word bits,
#    then xors ("^" in Python) with the start value.

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

@cython.cdivision(True)
cdef inline int gray_decode_travel(int start, int end, int mask, int g):
    cdef int travel_bit = start ^ end
    cdef int modulus = mask + 1          # == 2**nBits
    cdef int rg = (g ^ start) * (modulus / (travel_bit * 2))
#    print 'rg', rg
    cdef int result = gray_decode((rg | (rg / modulus)) & mask)
#    print 'rgresult', result
    return result


# child_start_end( parent_start, parent_end, mask, i ) -- Get start & end for child.
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
#    always flips the same bit on the first, third, fifth, ... and last flip.
#    The pattern works for any nD >= 1.
#
cdef int  child_start(int parent_start, int parent_end, int mask, int i):
    cdef int start_i = max(0, (i - 1) & ~1)  # next lower even number, or 0
    cdef int child_start = gray_encode_travel(parent_start, parent_end, mask, start_i)
    return child_start

cdef int child_end(int parent_start, int parent_end, int mask, int i):
    cdef int end_i = min(mask, (i + 1) | 1)  # next higher odd number, or mask
    cdef int child_end = gray_encode_travel(parent_start, parent_end, mask, end_i)
    return child_end



def _test_lengthy(maxside):
    """
    Test that encode and decode functions agree with each other,
    also in higher dims
    """
    # 2D
    for x in xrange(maxside):
        for y in xrange(maxside):
            c = (x, y)
            e = encode(c)
            d = decode(e, 2)
            assert c == d, (c, d)
    # 3D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                c = (x, y, z)
                e = encode(c)
                d = decode(e, 3)
                assert c == d, (c, d)
    # 4D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                for t in xrange(maxside):
                    c = (x, y, z, t)
                    e = encode(c)
                    d = decode(e, 4)
                assert c == d, (c, d)
    # 5D
    for x in xrange(maxside):
        for y in xrange(maxside):
            for z in xrange(maxside):
                for t in xrange(maxside):
                    for s in xrange(maxside):
                        c = (x, y, z, t, s)
                        e = encode(c)
                        d = decode(e, 5)
                    assert c == d, (c, d)

    print pow(maxside, 2) + pow(maxside, 3)+ pow(maxside, 4)+ pow(maxside, 5)


def _test_first_step_in_x():
    # if we increase the dimension, the first step should be in x-direction
    # the other dimensions should be 0
    for d in range(1, 10):
        result = decode(1, d)
        assert result[0] == 1
        for digit in result[1:]:
            assert digit == 0


if __name__ == "__main__":
    _test_lengthy(maxside = 2**4)
    _test_first_step_in_x()

