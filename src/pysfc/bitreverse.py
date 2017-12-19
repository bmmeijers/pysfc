"""Reverse an integer bit-by-bit
"""

#From: http://www.katjaas.nl/bitreversal/bitreversal.html
# (with nice explanation)
#
#unsigned int bitrev(unsigned int n, unsigned int bits)
#{
#    unsigned int i, nrev;   // nrev will store the bit-reversed pattern
#    N = 1<<bits;            // find N: shift left 1 by the number of bits
#    nrev = n;  

#    for(i=1; i<bits; i++)
#    {
#        n >>= 1;
#        nrev <<= 1;
#        nrev |= n & 1;   // give LSB of n to nrev
#    }
#    nrev &= N-1;         // clear all bits more significant than N-1
#   
#    return nrev;
#}

def memoize(function):
    """Decorator for caching function calls"""
    memo = {}
    def wrapper(*args):
        if args in memo:            # return from cache
            return memo[args]
        else:                       # compute
            rv = function(*args)
            memo[args] = rv
            return rv
    return wrapper


@memoize
def bitreverse(n, bits):
    N = 1 << bits           # find N: shift left 1 by the number of bits
    nrev = n                # nrev will store the bit-reversed pattern
    for i in range(1, bits):
        n >>= 1
        nrev <<= 1
        nrev |= n & 1       # give LSB of n to nrev
    nrev &= N - 1           # clear all bits more significant than N-1
    return nrev


if __name__ == "__main__":
    for j in range(100):
        for n in range(16):
            n, bitreverse(n, 4)

