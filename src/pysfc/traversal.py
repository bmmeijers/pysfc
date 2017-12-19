from operator import itemgetter
from hilbert import pack_index


#def nth_child_coordinates_zorder(n, parent_box):
#    """
#    calculate coordinates of the nth-child given coordinates of the parent nd-box
#    where n is given as zorder integer 
#    (to be able to derive which of the nary childs you get)

#    An example for 3D

#    Note that 0 = at low side, and 1 = at high side of dimension
#    for n=0 -> z=0, y=0, x=0 means the child is located with respect to center:
#    (at low side, at low side, at low side)

#        n    zyx
#        -  -----
#        0  0b000
#        1  0b001
#        2  0b010
#        3  0b011

#        4  0b100
#        5  0b101
#        6  0b110
#        7  0b111

#    """
#    c = [l + ((h-l) / 2) for h, l in zip(parent_box.hi, parent_box.lo)]
#    l = parent_box.lo[:]
#    h = parent_box.hi[:]
#    dim = parent_box.dims
#    #    print
#    #    print n, binstr(n, dim)
#    for d in range(dim):
#        rshifted = (n >> d)  # bit shift the child index number to the right by d bits
#        at_high_side_for_dim_d = bool((rshifted) & 1)
#        #        print "", names[d],  ":", binstr(rshifted, dim), "&", binstr(1, dim), "=", at_high_side_for_dim_d
#        if at_high_side_for_dim_d:
#            # keep the high, shift the low ordinate to the center
#            l[d] = c[d]
#        else:
#            # keep the low, shift the high ordinate to the center
#            h[d] = c[d]
#    #    print " >>>", l, h
#    # return the coordinates of this child box
#    return ndbox(l, h)


#def nth_child_coordinates_norder(n, parent_box):
#    """
#    calculate coordinates of the nth-child given coordinates of the parent nd-box
#    where n is given as norder integer 
#    (to be able to derive which of the nary childs you get)

#    An example for 3D

#    Note that 0 = at low side, and 1 = at high side of dimension
#    for n=0 -> z=0, y=0, x=0 means the child is located with respect to center:
#    (at low side, at low side, at low side)

#    n-order:

#        n     xy
#        -  -----
#        0     00
#        1     01
#        2     10
#        3     11

#    z-order (note: flipped xy on top):

#        n     yx
#        -  -----
#        0     00
#        1     01
#        2     10
#        3     11

#    """
#    c = [l + ((h-l) / 2) for h, l in zip(parent_box.hi, parent_box.lo)]
#    l = list(parent_box.lo[:])
#    h = list(parent_box.hi[:])
#    dim = parent_box.dims
#    for d in range(dim): # 2 = x, y = 1, z = 0
##        print n, ">>", d
#        rshifted = (bitreverse(n, dim) >> d)  # bit shift the child index number to the right by d bits AFTER THE BITS OF N ARE REVERSED
#        at_high_side_for_dim_d = bool((rshifted) & 1)
#        #        print "", names[d],  ":", binstr(rshifted, dim), "&", binstr(1, dim), "=", at_high_side_for_dim_d
#        if at_high_side_for_dim_d:
#            # keep the high, shift the low ordinate to the center
#            l[d] = c[d]
#        else:
#            # keep the low, shift the high ordinate to the center
#            h[d] = c[d]
#    #    print " >>>", l, h
#    # return the coordinates of this child box
#    return ndbox(l, h)


def nary(dim):
    # how many childs does a nd-box have given the number of dimensions
    return 2 ** dim


def morton_zcode_ncode(dim, sort_by_ncode=False):
    ary = nary(dim)
    codes = [(zcode, bitreverse(zcode, dim)) for zcode in range(ary)]
    if sort_by_ncode:
        codes.sort(key=itemgetter(1))
    return codes


def path_to_morton(path, dim, order):
    """Given a traversed path in the tree (with integers in z-order),
    determine the morton code that starts in the last box of the path at this level
    """
    # FIXME: check whether geometry (coordinates) correspond with the Morton code needed...
    # diff: how far are we from the bottom of the tree?
    diff = order - 1 # order gives depth of the tree, starting from 1 (with the level that has nary boxes)
    start = 0
    for zorder in path:
        width = (diff * dim)
        start += zorder << width
        diff -= 1    # we move 1 closer to the bottom of the tree
    # width: the size of the range belonging to this node at this level
    level = len(path)
    width = 2**(dim*(order-level)) 
    # the width of the range belonging to this node at this level
    return start, start + width


def path_to_hilbert(path, nD, order):
    """Given a traversed path in the tree (with integers that represent h-codes),
    determine the hilbert code that starts in this box of the path
    """
    level = len(path)
    # width: the size of the range belonging to this node at this level
    width = 2**(nD*(order-level))
    hcode = pack_index(path, nD) * width
    return hcode, hcode + width 
