import relate
import query_hilbert
#import pprint

bx = relate.ndbox([1066052,1642769,1899],[1083529,1677722,13291])

r = query_hilbert.hquery(bx, maxdepth=15)
