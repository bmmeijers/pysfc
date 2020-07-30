from pysfc.pysfc import _determine_bits, _key_nchunks, _key_hchunks, _hchunks_to_nchunks, _nchunks_coord
from pysfc.ndgeom import ndbox

# FIXME: 

# [ ] does the histogram preserve correctly space <-> (partial) key mapping?
#     [*] what does chunks represent: N-chunks / H-chunks? --> N-chunks!
#     [ ] unit tests
#         [ ] 2d
#         [ ] 3d
#         ......
#         [ ] nd
# [ ] take sample of records for histogram -> how large should such a sample be?
#     https://www.surveymonkey.com/mp/sample-size-calculator/ uses z-score
#     https://en.wikipedia.org/wiki/Margin_of_error
#

class Histogram(object):
    """nD-Histogram"""
    def __init__(self, mbits, ndims, max_depth, is_hilbert = True):
        """ """
        self.mbits = mbits
        self.ndims = ndims

        self.max_depth = max_depth

        self.is_hilbert = is_hilbert
        # FIXME: 
        # should we store how large the sample size is, so we 
        # can get to an estimate if we do not count *all* records inside
        # a cell?
        self.counts = {(): 0} # dict: n-chunks (tuple of int's) -> count (int)

    def add_key(self, key):
        """ Add a SFC key to the histogram """

        # FIXME: this unpacks the key completely considering mbits/ndims
        # we could make this more lazily (only unpack while iterating
        # over part of the key we *really* need for indexing in the histogram)
        if self.is_hilbert:
            hchunks = tuple(_key_hchunks(key, self.mbits, self.ndims))
            nchunks = _hchunks_to_nchunks(hchunks, self.mbits, self.ndims)
        else:
            nchunks = tuple(_key_nchunks(key, self.mbits, self.ndims))

        self.counts[()] += 1
        for i in range(min(len(nchunks), self.max_depth+1)):
            partkey = nchunks[:i+1]
            if partkey in self.counts:
                self.counts[partkey] += 1
            else:
                self.counts[partkey] = 1

    def query(self, nchunks):
        """ Query point count for chunks tuple """
        try:
            return self.counts[nchunks]
        except KeyError:
            return -1

    # FIXME: convert to a staticmethod ?
    def cell_as_ndbox(self, nchunks):
        """ Return ndbox (i.e. space covered) by chunks tuple """
        cur_level = len(nchunks)
        lo = _nchunks_coord(nchunks, self.ndims)
        # side_at_depth 
        # how large side of the nd-cube is at this depth
        side_at_depth = 2**(self.mbits - cur_level)
        lo = tuple(map(lambda x: int(x * side_at_depth), lo))
        hi = tuple(map(lambda x: int(x + side_at_depth), lo))
        return ndbox(lo, hi)


if __name__ == "__main__":
#    import pprint
    from pysfc.ndgeom import visualize_box_2d
    from pysfc.pysfc import nenc, henc

    side = 64
    mbits = _determine_bits(side-1, 2)
    ndims = 2

    use_hilbert = True
#    use_hilbert = False

#    print('mbits:', mbits, 'ndims:', ndims, 'hilbert?:', use_hilbert)
    curve_len = pow(side, ndims)
    #
    hist = Histogram(mbits, ndims, max_depth=3, is_hilbert=use_hilbert)

    if use_hilbert:
        enc = henc
    else:
        enc = nenc

    for x in range(side):
        for y in range(side):
            pt = (x, y)
            key = enc(pt)
#            print(key)
            hist.add_key(key)

    for x in range(16, 32):
        for y in range(16, 32):
            pt = (x, y)
            key = enc(pt)
#            print(key)
            hist.add_key(key)

    for x in range(8, 16):
        for y in range(8, 16):
            pt = (x, y)
            key = enc(pt)
#            print(key)
            hist.add_key(key)

    for i in range(4):
        for x in range(48, 64):
            for y in range(32, 48):
                pt = (x, y)
                key = enc(pt)
    #            print(key)
                hist.add_key(key)

    
#    # add 'points' to the histogram by means of their key
#    for key in range(curve_len):
#        hist.add_key(key)
#    # 2 times as many points in second half of the curve
#    for key in range(curve_len//2, curve_len):
#        hist.add_key(key)
#    # in cell 7 add extra points -> we add in total 4 times points in this cell
#    for _ in range(3):
#        for key in range(curve_len//16 * 7, curve_len//16 * 8):
#            hist.add_key(key)

##    print('query result')
##    for i in range(2**ndims):
##        for j in range(2**ndims):
##                result = hist.query( (i, j) )
##                print((i, j), ' -> ', result)

##            for k in range(2**ndims):
##                result = hist.query( (i, j, k) )
##                print((i, j, k), ' -> ', result)
##                for l in range(2**ndims):   
##                    result = hist.query(tuple([i, j, k, l]))
##                    print((i, j, k, l), ' -> ', result)

#    # the hilbert and morton curve behave differently with respect to
#    # the second 'half' (0+1 vs. 2+3) of the cells in 2D
#    # -- see if the counts are mapped properly
##    if ndims == 2 and side == 64:
#    if not use_hilbert:
#        assert hist.query((0,)) == 4096 / 4 
#        assert hist.query((2,)) == 4096 / 4 * 2
#        assert hist.query((3,)) == 4096 / 4 * 2
#        assert hist.query((1, 2)) == 256
#        assert hist.query((1, 3)) == 1024
#    elif use_hilbert:
#        assert hist.query((0,)) == 4096 / 4
#        assert hist.query((1, 2)) == 1024
#        assert hist.query((1, 3)) == 256
#        assert hist.query((2,)) == 4096 / 4 * 2
#        assert hist.query((3,)) == 4096 / 4 * 2

    for i in range(2**ndims):
        for j in range(2**ndims):
#            for k in range(2**ndims):
                chunks = (i, j)
                print(i, '\t', j, '\t', visualize_box_2d(hist.cell_as_ndbox(chunks)), '\t', hist.query( chunks ))

#    for chunk in range(pow(2, ndims)):
#        result = hist.query(tuple([chunk]))
#        print(chunk, result)

