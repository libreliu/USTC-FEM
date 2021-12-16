import numpy as np
import matplotlib.pyplot as plt
import math
#import scipy.integrate as integrate

def simpsonIntg(f, x0, x1):
    """Do numeric intergral via Simpson's formula"""
    return (x1 - x0) * ((1.0/6) * f(x0) + (4.0/6) * f((x0 + x1) / 2) + (1.0/6) * f(x1))

def linearIntg(f, x0, x1):
    return (x1 - x0) * (0.5 * f(x0) + 0.5 * f(x1))

# Two point gauss integration
def gaussIntg(f, x0, x1):
    c1 = (x1 - x0) / 2
    c2 = (x1 + x0) / 2
    c3 = 1 / math.sqrt(3)
    c4 = -c3
    return c1 * f(c2 + c1 * c3) + c1 * f(c2 + c1 * c4)

def preciseIntg(f, x0, x1):
    res = integrate.quad(f, x0, x1)
    #print(f"Integrate from {x0} to {x1} yield {res[0]}, error={res[1]}")
    return res[0]

def intgWrapper(useIntg, refIntg):
    def intg(f, x0, x1):
        use = useIntg(f, x0, x1)
        ref = refIntg(f, x0, x1)
        print(f"Integrate error: {use - ref}")
        return use
    return intg

class linearFEM:
    def __init__(self, intg: 'function'):
        """
        @param intg: Numerical integrater to use
        """
        self.numIntg = intg
        self.U = None

    def build_knots(self, n: int):
        """Build uniform knots in generation, total n+1 knots"""
        self.X = [np.linspace(0, 1, num=n+1)]
        self.n = [n]
        self.U = [None]
    
    def build_knots_multigrid(self, numLayers, nStart):
        """
        numLayers: number of multigrid layers
        nStart: starting numInterval (total nStart+1 knots)
        """
        self.X = [np.linspace(0, 1, num=nStart + 1)]
        self.n = [nStart]
        for i in range(1, numLayers):
            self.n.append(self.n[i-1] * 2)
            self.X.append(
                np.linspace(0, 1, num=self.n[i] + 1)
            )

        self.U = [None for i in range(0, numLayers)]
        
    def phiDeriv(self, k, x, i):
        if self.X[k][i-1] <= x <= self.X[k][i]:
            h_i = self.X[k][i] - self.X[k][i-1]
            return 1.0 / h_i
        elif self.X[k][i] <= x <= self.X[k][i+1]:
            h_iP1 = self.X[k][i+1] - self.X[k][i]
            return -1.0 / h_iP1
        else:
            return 0

    def phi(self, k, x, i):
        assert(1 <= i <= self.n[k] - 1)
        if self.X[k][i-1] <= x <= self.X[k][i]:
            h_i = self.X[k][i] - self.X[k][i-1]
            return (x - self.X[k][i-1]) / h_i
        elif self.X[k][i] <= x <= self.X[k][i+1]:
            h_iP1 = self.X[k][i+1] - self.X[k][i]
            return (self.X[k][i+1] - x) / h_iP1
        else:
            return 0

    def phiDerivPhiDerivIntg(self, k, i, j):
        """Calculate \int_0^1 \phi'_i(x) \phi'_j(x) dx
        Here we use i = 1, ..., n-1
        """
        if i > j:
            # This is symmetry so no problem
            return self.phiDerivPhiDerivIntg(k, j, i)

        assert(i <= j)
        if j - i >= 2:
            return 0
        elif j - i == 1:
            return self.numIntg(
                lambda x: self.phiDeriv(k, x, i) * self.phiDeriv(k, x, j),
                self.X[k][i],
                self.X[k][i+1]
            )
        else:
            assert(j == i)
            return self.numIntg(
                lambda x: self.phiDeriv(k, x, i) * self.phiDeriv(k, x, j),
                self.X[k][i-1],
                self.X[k][i]
            ) + self.numIntg(
                lambda x: self.phiDeriv(k, x, i) * self.phiDeriv(k, x, j),
                self.X[k][i],
                self.X[k][i+1]
            )
        

    def phiDerivPhiIntg(self, k, i, j):
        """
        Calculate \int_0^1 \phi'_i(x) \phi_j(x) dx
        Here we use i = 1, ..., n-1
        """
        if i <= j:
            if j - i >= 2:
                return 0
            elif j - i == 1:
                return self.numIntg(
                    lambda x: self.phiDeriv(k, x, i) * self.phi(k, x, j),
                    self.X[k][i],
                    self.X[k][i+1]
                )
            else:
                assert(j == i)
                return self.numIntg(
                    lambda x: self.phiDeriv(k, x, i) * self.phi(k, x, j),
                    self.X[k][i-1],
                    self.X[k][i]
                ) + self.numIntg(
                    lambda x: self.phiDeriv(k, x, i) * self.phi(k, x, j),
                    self.X[k][i],
                    self.X[k][i+1]
                )
        else:
            assert(i > j)
            if i - j >= 2:
                return 0
            elif i - j == 1:
                return self.numIntg(
                    lambda x: self.phiDeriv(k, x, i) * self.phi(k, x, j),
                    self.X[k][j],
                    self.X[k][j+1]
                )

    def phiFIntg(self, k, i, f):
        """
        Calculate \int_0^1 f(x) \phi_i(x) dx
        Here we use i = 1, ..., n-1
        """
        return self.numIntg(
            lambda x: f(x) * self.phi(k, x, i),
            self.X[k][i-1],
            self.X[k][i]
        ) + self.numIntg(
            lambda x: f(x) * self.phi(k, x, i),
            self.X[k][i],
            self.X[k][i+1]
        )

    def solve(self, k, f: 'function'):
        """Solve using V_{h1} finite element space
        
        Uses n+1 knots in total, namely 0, ..., n
        and we have n-1 basis functions in total
        """
        X = self.X[k]
        n = self.n[k]

        # Check x
        assert(len(X) == n + 1)

        # Build K
        K = np.ndarray((n-1, n-1), dtype=np.double)
        for i in range(0, n-1):
            for j in range(0, n-1):
                K[i, j] = self.phiDerivPhiDerivIntg(k, i + 1, j + 1)
        
        #np.set_printoptions(precision=3, linewidth=150)
        #print(K)

        # Build F
        F = np.zeros((n-1, ), dtype=np.double)
        for i in range(0, n-1):
            F[i] = self.phiFIntg(k, i+1, f)
        
        #print(F)

        # np.linalg.solve() uses _gesv LAPACK routines
        # which uses LU factorization indeed
        self.U[k] = np.linalg.solve(K, F)

        # Check result
        chk = K @ self.U[k]
        assert(np.allclose(chk, F))

        return self.U[k]
    
    # TODO: test if it works
    def nonZeroEntries(self, k):
        """Returns an iterator that gives index for non-zero entries.
        
        O(num_of_nonzero_entries) complexity"""
        n = self.n[k]
        for i in range(0, n-1):
            # candidate: i-1, i, i+1
            for j in range(i-1, i+2):
                if j >= 0 and j <= n-2:
                    yield (i, j)

    def computeInnerProduct(self, k, left, right):
        """Compute left.T * A * right"""
        n = self.n[k]
        assert(len(left) == n-1 and len(right) == n-1)

        # only do on non-sparse items
        # TODO: use symmetry
        res = np.double(0.0)
        for i, j in self.nonZeroEntries(k):
            a_ij = self.phiDerivPhiDerivIntg(k, i + 1, j + 1)
            res += a_ij * left[i] * right[j]

        return res
    
    def computeMatrixProduct(self, k, right):
        """Compute A @ right"""
        n = self.n[k]
        assert(len(right) == n-1)

        res = np.zeros((n-1,), dtype=np.double)
        for i, j in self.nonZeroEntries(k):
            a_ij = self.phiDerivPhiDerivIntg(k, i + 1, j + 1)
            res[i] += a_ij * right[j]
        
        return res

    def solveCG(self, k, f: 'function', threshold=1e-5):
        n = self.n[k]
        X = self.X[k]

        assert(len(X) == n + 1)
        U = np.zeros((n-1,), dtype=np.double)
        F = np.zeros((n-1,), dtype=np.double)

        for i in range(0, n-1):
            F[i] = self.phiFIntg(k, i+1, f)

        r = -F
        d = -r

        while np.linalg.norm(r, ord=2) > threshold:
            dT_A_d = self.computeInnerProduct(k, d, d)
            U_correction = np.dot(r, r) / dT_A_d * d
            U += U_correction

            r = self.computeMatrixProduct(k, U) - F
            d = -r + (self.computeInnerProduct(k, r, d) / dT_A_d * d)
        
        self.U[k] = U

    def solveMG(self, f: 'function'):
        pass

    def MG(self, k, z_0, g):
        """k'th MultiGrid; z_0: """
        if k == 1:
            pass

    def genSample(self, k, n):
        T = np.linspace(0, 1, n)
        res = np.zeros((n,), dtype=np.double)
        for i in range(0, n):
            t = T[i]
            for j in range(1, self.n[k]):
                res[i] += self.U[k][j - 1] * self.phi(k, t, j)
        
        return (T, res)

    def plot(self, k, n, ax):
        T = self.X[k]
        res = np.zeros((self.n[k]+1,), dtype=np.double)
        for i in range(0, self.n[k] + 1):
            if i == 0 or i == self.n[k]:
                res[i] = 0
            else:
                res[i] = self.U[k][i - 1]

        ax.scatter(T, res)

        # only works for uniform spacing
        T, res = self.genSample(k, n)
        ax.plot(T, res)


def plot_ref(ax, ref_func):
    x = np.linspace(0, 1, 100)
    ref = [val for val in map(ref_func, x)]

    ax.plot(x, ref)

if __name__ == '__main__':
    f = lambda x: (x-1) * math.sin(x)
    u_ref = lambda x: (x-1) * math.sin(x) + 2 * math.cos(x) + (2 - 2 * math.cos(1))* x - 2

    fig, ax = plt.subplots()

    CGFEM = linearFEM(gaussIntg)
    CGFEM.build_knots(5)
    CGFEM.solveCG(0, f)

    CGFEM.plot(0, 50, ax)

    # directFEM = linearFEM(gaussIntg)
    # directFEM.build_knots(5)
    # directFEM.solve(0, f)

    # directFEM.plot(0, 50, ax)

    plot_ref(ax, u_ref)

    plt.show()