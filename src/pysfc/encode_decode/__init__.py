from functools import lru_cache
# from operator import itemgetter

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

def _determine_bits(n, base):
    """Determine how many words we need to represent integer *n*
    *base* encodes how many bits a word contains
    """
    count = 0
    while (n):
        count +=1
        n >>= base
    return count


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


# def _chunks_key(chunks, mbits, ndims):
#     key = 0
#     nd = 2**ndims
#     for m in range(mbits):
#         if m >= len(chunks):
#             chunk = 0
#         else:
#             chunk = chunks[m]
#         key = key * nd + chunk
#     return key

def _chunks_key(chunks, mbits, ndims):
    key = 0
    nd = 2**ndims
    for m in range(mbits):
        chunk = 0
        if m < len(chunks):
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
    """Convert Morton key into nchunks tuple

    For example:

        >>> _key_nchunks(45, 3, 2)
        (2, 3, 1)

    """
    # EXACT SAME behaviour as _key_hchunks ???
    nchunks = [0] * mbits
    # mask is a number with all bits set to 1,
    # with ndims bits, e.g. ndims = 3 -> (1 << 3) - 1 == 8 - 7 == 7 == 0b111
    mask = (1 << ndims) - 1 
    for i in range(mbits):
        nchunks[mbits - i - 1] = ((key >> (i * ndims)) & mask)
    return tuple(nchunks)


# -- Hilbert specifics: for rotation and mirroring of the curve pattern --------

@lru_cache(maxsize=None)
def _hgray_encode(bn):
    # Gray encoder and decoder from http://en.wikipedia.org/wiki/Gray_code
    assert bn >= 0
    assert type(bn) in [int], "Expected int, but found: {}".format(type(bn))
    return bn ^ (bn // 2)


@lru_cache(maxsize=None)
def _hgray_decode(n):
    assert type(n) in [int], "Expected int, but found: {}".format(type(n))
    sh = 1
    while True:
        div = n >> sh
        n ^= div
        if div <= 1:
            return n
        sh <<= 1


@lru_cache(maxsize=None)
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
    #print(f"travel_bit {travel_bit}")
    #print(f"modulus {modulus}")
    # travel_bit = 2**p, the bit we want to travel.
    # Canonical Gray code travels the top bit, 2**(nBits-1).
    # So we need to rotate by ( p - (nBits-1) ) == (p + 1) mod nBits.
    # We rotate by multiplying and dividing by powers of two:
    g = _hgray_encode(i) * (travel_bit * 2)
    #print(f"g {g}")
    result =  ((g | (g // modulus)) & mask) ^ start
    #print(f"result {result}")
    return result


@lru_cache(maxsize=None)
def _hgray_decode_travel(start, end, mask, g):
    travel_bit = start ^ end
    modulus = mask + 1          # == 2**nBits
    rg = (g ^ start) * (modulus // (travel_bit * 2))
    result = _hgray_decode((rg | (rg // modulus)) & mask)
    return result


@lru_cache(maxsize=None)
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
    #print(f"start_i {start_i}")
    end_i = min(mask, (i + 1) | 1)  # next higher odd number, or mask
    #print(f"end_i {end_i}")
    # print(parent_start, parent_end, mask, start_i)
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
    #print(f"nchunks_len {nchunks_len}")
    mask = 2 ** ndims - 1
    #print(f"mask {mask}")
    start, end = _hinitial_start_end(mbits, ndims)
    #print(f"start {start} end {end}")
    hchunks = [0] * nchunks_len
    for j, nchunk in enumerate(nchunks):
        #print(f"nchunk {nchunk}")
        hchunks[j] = _hgray_decode_travel(start, end, mask, nchunk)
        start, end = _hchild_start_end(start, end, mask, hchunks[j])
        #print(f"start {start}, end {end}, mask {mask}, hchunks[j] {hchunks[j]}")
        #print("")
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

    return tuple(hchunks)


# Public api, encode and decode, n-order (nenc+ndec) and hilbert (henc + hdec)

# -- Morton (N-order)
# def nenc_fixed(coord, mbits, ndims): # ndims is also the len(coord) -- should it be in the api?
#     #mbits = _determine_bits(max(coord), 2)
#     # find next power of two that surrounds cube fully, i.e. max(coord)
#     # mbits = _determine_bits(max(coord), 1) 
#     # print(mbits)
#     ndims = len(coord)
#     #
#     nchunks = _coord_nchunks(coord, mbits)
#     #print(nchunks)
    
#     key_arr_u8 = pack(nchunks, mbits, ndims)
#     return key_arr_u8


# def ndec_fixed(key_arr_u8, mbits, ndims):
#     #packer = BitPacker(mbits, ndims)
#     nchunks = unpack(key_arr_u8, mbits, ndims)
#     #mbits = _determine_bits(key, 2**ndims)
#     #mbits = _determine_bits(key, ndims)
#     # print(mbits)
#     #
#     #nchunks = _key_nchunks(key, mbits, ndims)
#     coord = _nchunks_coord(nchunks, ndims)

#     # print("   ndec --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", nchunks, "key:", key, "coord:", coord)
#     # print(nchunks, coord)
#     # return mbits, nchunks, coord
#     return coord


# -- Hilbert
# def henc_fixed(coord, mbits, ndims):
#     #mbits = _determine_bits(max(coord), 2)
#     #mbits = _determine_bits(max(coord), 1)
#     #ndims = len(coord)
#     #
#     #print("coord", coord, "mbits", mbits)
#     nchunks = _coord_nchunks(coord, mbits)
#     #print(nchunks, mbits, ndims)
#     # FIXME: can we replace the following function call
#     # by look ups in dictionaries? size of the lut's does depend on ndims
#     # so memory consumption will grow with higher dims (but runtime will most
#     # likely be some factor faster)
#     hchunks = _nchunks_to_hchunks(nchunks, mbits, ndims)
#     #print(hchunks)
#     #key = _hchunks_key(hchunks, mbits, ndims)
#     # print("   henc --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", hchunks, "key:", key, "coord:", coord)
#     #return key
#     key_arr_u8 = pack(hchunks, mbits, ndims)
#     return key_arr_u8


# def hdec_fixed(key_arr_u8, mbits, ndims):
#     #mbits = _determine_bits(key, 2**ndims)
#     #mbits = _determine_bits(key, 2**ndims)
#     #mbits = _determine_bits(key, ndims)
#     #
#     #assert type(mbits) in [int]
#     #assert type(key) in [int]
#     #assert type(ndims) in [int]
#     hchunks = unpack(key_arr_u8, mbits, ndims)
#     #print(hchunks)
#     nchunks = _hchunks_to_nchunks(hchunks, mbits, ndims)
#     #print(nchunks)
#     coord = _nchunks_coord(nchunks, ndims)
#     # print("   hdec --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", hchunks, "key:", key, "coord:", coord)
#     return tuple(coord)





# -- Morton (N-order)
def nenc(coord):
    #mbits = _determine_bits(max(coord), 2)
    # find next power of two that surrounds cube fully, i.e. max(coord)
    mbits = _determine_bits(max(coord), 1) 
    # print(mbits)
    ndims = len(coord)
    #
    nchunks = _coord_nchunks(coord, mbits)
    key = _hchunks_key(nchunks, mbits, ndims)
    # print("   nenc --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", nchunks, "key:", key, "coord:", coord)
    return key


def ndec(key, ndims):
    #mbits = _determine_bits(key, 2**ndims)
    mbits = _determine_bits(key, ndims)
    # print(mbits)
    #
    nchunks = _key_nchunks(key, mbits, ndims)
    coord = _nchunks_coord(nchunks, ndims)

    # print("   ndec --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", nchunks, "key:", key, "coord:", coord)
    # print(nchunks, coord)
    # return mbits, nchunks, coord
    return coord


# -- Hilbert
def henc(coord):
    #mbits = _determine_bits(max(coord), 2)
    mbits = _determine_bits(max(coord), 1)
    ndims = len(coord)
    #
    #print("coord", coord, "mbits", mbits)
    nchunks = _coord_nchunks(coord, mbits)
    #print(nchunks, mbits, ndims)
    # FIXME: can we replace the following function call
    # by look ups in dictionaries? size of the lut's does depend on ndims
    # so memory consumption will grow with higher dims (but runtime will most
    # likely be some factor faster)
    hchunks = _nchunks_to_hchunks(nchunks, mbits, ndims)
    #print(hchunks)
    key = _hchunks_key(hchunks, mbits, ndims)
    # print("   henc --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", hchunks, "key:", key, "coord:", coord)
    return key


def hdec(key, ndims):
    #mbits = _determine_bits(key, 2**ndims)
    #mbits = _determine_bits(key, 2**ndims)
    mbits = _determine_bits(key, ndims)
    #
    assert type(mbits) in [int]
    assert type(key) in [int]
    assert type(ndims) in [int]
    hchunks = _key_hchunks(key, mbits, ndims)
    #print(hchunks)
    nchunks = _hchunks_to_nchunks(hchunks, mbits, ndims)
    #print(nchunks)
    coord = _nchunks_coord(nchunks, ndims)
    # print("   hdec --", "bitlength(key):", key.bit_length(), "mbits:", mbits, "nchunks:", hchunks, "key:", key, "coord:", coord)
    return tuple(coord)
