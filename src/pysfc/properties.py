from pysfc.encode_decode import nenc, ndec, henc, hdec

n_dims = 2
max_bits = 3
branching_factor = 2 ** n_dims
total_number_of_sfc_keys = branching_factor ** max_bits
last_valid_sfc_key_value = total_number_of_sfc_keys - 1
print("""

""")
print("# curve properties")
print("n_dims", n_dims)
print("m_bits", max_bits)
print("total_number_of_sfc_keys", total_number_of_sfc_keys)
print("last_valid_sfc_key_value", last_valid_sfc_key_value)
print()

levels = {}
resolutions = {}
side_lengths = {}
for level in range(max_bits + 1):
    childs_at_level = pow(branching_factor, level)
    resolution_at_level = total_number_of_sfc_keys // childs_at_level
    side_length_at_level = int(pow(resolution_at_level, 1.0 / n_dims))
    levels[level] = childs_at_level
    resolutions[level] = resolution_at_level
    side_lengths[level] = side_length_at_level

print("# number of child 'nodes' in n-ary tree (per level)")
for level in levels:
    print(level, levels[level])
print()

print("# sfc keys contained in a node (resolution of a n-ary node)")
for level in resolutions:
    print(level, resolutions[level])
print()

print("# side_length of an n-ary node (cube) -> node contains number of keys")
for level in side_lengths:
    print(level,"", side_lengths[level], "->", pow(side_lengths[level], n_dims))
print()

def words(s, length):
    """Split a bit string _s_ into words (groups of characters) of length _length_"""
    # s should be output of bin(val)[2:], where val is positive int (incl zero)
    from collections import deque
    words = deque([])
    while s:
        words.appendleft(s[-length:])
        s = s[:-length]
    # return "-".join(words)
    return tuple(int(n, 2) for n in words)


print()
print(f"# level, key, bin, {n_dims}-words, {n_dims}-words upto level, node range")

def details_for_node(level, sfc_key, end, enc, dec):
    bits = bin(sfc_key)[2:].zfill(max_bits * n_dims)
    # sfc_key.bit_length()
    wrds = words(bits, n_dims)
    f = str(len(str(last_valid_sfc_key_value))) + 'd'
    print(f"{level}, {sfc_key:{f}}, {bits}, {wrds}, {wrds[:level]}, {sfc_key:{f}}—{end-1:{f}}, ")
    coord = dec(sfc_key, ndims=n_dims)
    enc(coord)




from collections import deque

for enc, dec in [(nenc, ndec), (henc, hdec)]:
    print()
    print()
    print("## ", enc, dec)
    max_level = min(2, max_bits)
    level = 0
    start = 0
    size = resolutions[level]
    stack = deque([(level, start, size)])
    prev_level = None
    while stack:
        level, start, size = stack.popleft()
        # print(level, start, size)
        end = start + size
        if prev_level != level:
            print()
        # else:
        #     print("    ⋮")
        details_for_node(level, start, end, enc, dec)
        prev_level = level
        if level < max_level:
            next_size = resolutions[level+1]
            # FIXME: 
            # we could also stack the node-id's: [0..branching_factor)
            # iterate over *number* of n-ary nodes
            # and multiply with size?
            for new_start in range(start, end, next_size):
                stack.append((level + 1, new_start, next_size))
