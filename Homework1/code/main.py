import numpy as np

def numIntg(f, x0, x1, eps=1e-5):
    """Do numeric intergral via automatic precision control"""

class PiecewiseLinearFEM:
    def __init__(self, intg = numIntg):
        """
        @param intg: Numerical integrater to use
        """
        self.numIntg = numIntg

    def build_knots(self, n, start, end):
        """Build uniform knots in generation"""
        return np.

    def solve(self, n, f, x):
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
                    K[i, j] = K[j, i] = -( 1.0 / (x[i+1] - x[i])) + (1.0 / (x[j] - x[j-1]))
                elif j == i:
                    K[i, j] = K[j, i] = (1.0 / ((x[i+1] - x[i])**2)) + (1.0 / ((x[i] - x[i-1]) ** 2))
                else:
                    assert(False)
        
        # Build F
        F = np.zeros((n+1, ), dtype=np.double)
        for i in range(0, n+1):
            F[i] = (1.0 / (x[i] - x[i-1])) * numIntg(f, x[i-1], x[i]) \
                - (1.0 / (x[i+1] - x[i])) * numIntg(f, x[i], x[i+1])
        
        # np.linalg.solve() uses _gesv LAPACK routines
        # which uses LU factorization indeed
        U = np.linalg.solve(K, F)

        # Check result
        assert(np.allclose(K @ U, F))

        return U


