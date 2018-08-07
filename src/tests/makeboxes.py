from pysfc.relate import ndbox

def make_boxes(mbits=4, ndims=3):
    # 2**M = max side of the cube
    assert mbits > 2
    assert ndims > 2
    boxes = []

    # generate a set of nd-boxes for querying
    # the boxes are positioned around the center of the cube
    # taking a full slice of 1 high in 1 dimension,
    # while using all other dimensions to their full extents

    for n in range(2, ndims):
        for m in range(2, mbits):
            val = 2**m
            half = (val - 1) // 2
            full_range = (0, val)
            below_half_range = (half, half + 1)
            over_half_range = (half + 1, half + 2)
            for where in below_half_range, over_half_range:
                for m in range(0, n):
                    lo, hi = [], []
                    for _ in range(n):
                        if _ != m:
                            r = full_range
                        else:
                            r = where
                        lo.append(r[0])
                        hi.append(r[1])
                    box = ndbox(lo, hi)
                    boxes.append(box)
    return boxes
