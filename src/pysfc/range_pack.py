"""Reduce the number of ranges, by looking at the gaps between them.
"""


def gaps(ranges):
    """Get a list with the size of the gaps between the ranges

    """
    gaps = []
    for cur, nxt in zip(ranges[:-1], ranges[1:]):
        gaps.append(nxt[0] - cur[1])    # 1: end, 0: start
    return gaps


def sizes(ranges):
    """Get a list with the size of the ranges"""
    sizes = []
    for r in ranges:
        sizes.append(r[1] - r[0])       # 1: end, 0: start
    return sizes


#def gap_size_index_sorted(ranges):
#    """Decorate each range with gap size to next range and its own size
#    """
#    out = []
#    for i, r in enumerate(ranges):
#        if i+1 == len(ranges):
#            next_start = None
#            gap = None
#        else:
#            next_start = ranges[i + 1][0]
#            gap = next_start - r[1]
#        # gap, size, current index, range
#        tup = (r, gap, r[1] - r[0], i)
#        out.append(tup)
#    return out


def filter_ranges_by_gapsize(ranges, max_gap_size=0):
    """ 
    Go over all ranges and those that should form 1 longer range, glue them

    E.g.
        [(0, 16), (16, 20), (24, 28), (32, 36), (36, 40), (48, 52)]
    becomes:
        [(0, 20), (24, 28), (32, 40), (48, 52)]

    Note, this code does the same as Oscar's merge_consecutive_ranges function
    but is less specialized (works also with gaps <> 0)

    """
    if len(ranges) == 0:
        return []
    else:
        filtered = []
        (range_start, range_end) = ranges[0]
        for idx in range(1, len(ranges)):
            cur_range = ranges[idx]
            gap_size = cur_range[0] - range_end
            if gap_size <= max_gap_size:
                range_end = cur_range[1]
            else:
                filtered.append((range_start, range_end))
                (range_start, range_end) = cur_range
        filtered.append((range_start, range_end))
        return filtered



if __name__ == "__main__":

    from query import query
    from hquery import hquery

    # ranges = [(0, 16), (16, 20), (24, 28), (32, 36), (36, 40), (48, 52)]

    #ranges = [(112, 116), (116, 120), (142, 143), (143, 144)]

    ranges = query()
    print ranges

#    ranges = filter_ranges_by_gapsize(ranges, 0)
#    print ranges


    # sorted gaps (differences between the ranges)
    gapsizes = sorted(gaps(ranges))
    print gapsizes # << will be empty if not enough ranges (e.g. if we only have 1 range to get)

    # nr of gaps = nr of ranges - 1
    
    # when we have the gaps sorted:
    # we can get an estimate on how many ranges will
    # be left after merging

    # [0, 0, 0, 4, 4, 8, 16] # <-- sorted gap size list (for 7 gaps)
    #    7, 6, 5, 4, 3, 2, 1 # nr of ranges that will remain, 
    #                          when earlier gaps are closed
    #
    # Note, we can close a bit more gaps due to equal sized gaps!
    # FIXME: Peter suggested to *not* close *all* equal sized gaps
    # It might be good to check in the sorted gaps list how many equals you have
    # And keep a counter before the index and only close these
    # An example, with gaps like this:
    # [0, 0, 0, 3, 3, 4, 4, 4, 4, 5]
    # if we want to have 4 ranges remaining we have to close the gaps:
    # [0, 0, 0, 3, 3, 4, 4], so we have to keep track that we want to close
    # all gaps < 4 anyway and we have to close gaps of size 4 only 2 times
    # this then would keep the gaps [4, 4, 5] and you get exactly 4 ranges

    # The max gap that we use to filter the ranges is then:
    # >>> gap_sizes[-(remain-1)]
    # where remain is 1-based and used as index into the gaps list,
    # this then gives the gap size that needs to be used for filtering
    # we need to be careful that it does not go out-of-bounds, though!

    # to prevent out-of-bounds we make sure that the maxranges
    # value the user will put will be according to:
    # 1 <= maxrange < len(gaps)+1 ?? 
    # *and* we have to make sure that gaps does have some gaps :-)
    # otherwise there is no need to filter...

    maxranges = 10
    if gapsizes and maxranges <= len(gapsizes):
        # make sure that index stays within bounds
        idx = maxranges
        if idx < 1:
            idx = 1
        if idx > len(gapsizes):
            idx = len(gapsizes)
        # take the gap size that we should use to filter
        # (note, we may get less ranges than specified)
        print -idx
        filter_gap_size = gapsizes[-idx]
        print "filter gaps that are <=", filter_gap_size
        print filter_ranges_by_gapsize(ranges, filter_gap_size)
    else:
        print "only filter 0-wide gaps -- number of ranges is already below", maxranges
        ranges = filter_ranges_by_gapsize(ranges, 0)
        print ranges

    # post condition that should hold: length of ranges <= maxranges !

