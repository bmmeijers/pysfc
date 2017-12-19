from query_brute import brute

from query_hilbert import hquery
from query_norder import nquery

from relate import ndbox

def test_query():
    """Check if querying bruteforce gives same result
    as performing the search in a more elegant way.
    """
    query = ndbox([7, 0, 1], [9, 16, 1])

    assert brute(query, 'h') == hquery(query)
    assert brute(query, 'n') == nquery(query)

if __name__ == "__main__":
    test_query()
