import numpy as np

def numIntg(f, x0, x1, eps=1e-5):
    """Do numeric intergral via automatic precision control"""

def build_knots(n, start, end):
    pass

def solve_h1(n, f, x):
    """Solve using V_{h1} finite element space
    
    Uses n+1 knots in total, namely 0, ..., n
    """

    # Check x
    assert(len(x) == n + 1)


    # Build K
    K = np.zeros((n+1, n+1), dtype=np.double)
    for i in range(0, n+1):
        for j in range(i, n+1):
            if j >= i + 2:
                K[i, j] = K[j, i] = 0
            elif j == i + 1:
                K[i, j] = K[j, i] = -( 1.0 / ())