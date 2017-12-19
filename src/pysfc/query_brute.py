"""Brute force querying
"""

from hilbert import encode as hencode
from morton_norder import encode as nencode
from morton_zorder import encode as zencode

from relate import ndbox

from itertools import product  # for cartesian product


def compute_deltas(query):
    """Returns how large every dimension is
    """
    dims = query.dims
    diff = [None] * dims
    for d in range(dims):
        diff[d] = query.hi[d] - query.lo[d]
    return diff


def brute(query, tp = 'h'):
    """A way of specifying ranges to search with in a brute force way

    Note, this function is only intended to test whether a more elegant
    way of querying performs in a correct way.
    """
    # the type of the SFC (h: hilbert, n: n-order, z: z-order)
    if tp == 'h':
        enc = hencode
    elif tp == 'n':
        enc = nencode
    elif tp == 'z':
        enc = zencode
    else:
        raise ValueError("Unknown encoding given")

    dims = query.dims
    diff = compute_deltas(query)
    # materialize the values that are there along each dimension
    along_dims = []
    for d in range(dims):
        start_of_dim = query.lo[d]
        along_d = []
        for val in range(diff[d]):
            along_d.append(start_of_dim + val)
        along_dims.append(along_d)

    # we obtain the Cartesian product of all the dimensions
    # and for each value we map the coordinate to the hilbert value
    ranges =[]
    for c in product(*along_dims):
        start = enc(c)
        ranges.append((start, start+1))

    # now sort the ranges we obtained (necessary for hilbert)
    ranges.sort()

    return ranges

    # we could go also with a linear mapping of index to subscripts as well
    # see numpy's *unravel_index*, with how many values there are along each
    # dimension

    #total = diff[0]
    #for ct in diff[1:]:
    #    total *= ct
    #print total

    # could use *reduce*
    # this is equivalent with:
    #
    # from operator import mul
    # reduce(mul, diff)
    # or:
    # reduce(lambda x,y: x*y, diff)


def howmany(diff):
    """How many cells are there for a given set of deltas along dimensions
    """
    total = diff[0]
    for ct in diff[1:]:
        total *= ct
    return total


def _test_small():
    query = ndbox([7, 0], [9, 16])
    print brute(query, 'h')
    print brute(query, 'n')
    print brute(query, 'z')


if __name__ == "__main__":

    query = ndbox([1066052,1642769,1899], [1083529,1677722,2057])
    diffs = compute_deltas(query)
    print diffs
    print howmany(diffs)

