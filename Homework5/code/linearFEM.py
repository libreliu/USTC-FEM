import numpy as np
import matplotlib.pyplot as plt
import math, time
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

        numIter = 0
        while np.linalg.norm(r, ord=2) > threshold:
            dT_A_d = self.computeInnerProduct(k, d, d)
            U_correction = np.dot(r, r) / dT_A_d * d
            U += U_correction

            r = self.computeMatrixProduct(k, U) - F
            d = -r + (self.computeInnerProduct(k, r, d) / dT_A_d * d)
            numIter += 1
        
        print(f"solveCG: done with {numIter} iterations")
        
        self.U[k] = U

    def iterateCG(self, k, z_0, g, numIter):
        """compute A_k z = g using CG; use z_0 as initial value"""
        r = -g + self.computeMatrixProduct(k, z_0)
        d = -r
        Z = z_0

        for i in range(0, numIter):
            dT_A_d = self.computeInnerProduct(k, d, d)
            Z_correction = np.dot(r, r) / dT_A_d * d
            Z += Z_correction

            r = self.computeMatrixProduct(k, Z) - g

            # Exit in time, since r == 0 will lead to d = 0,
            # and Z_correction will blow up
            if np.allclose(r, np.zeros_like(r)):
                return Z

            d = -r + (self.computeInnerProduct(k, r, d) / dT_A_d * d)
        
        return Z

    def coarseToFine(self, k_coarse, z_coarse):
        """Bring k_coarse -> k_coarse + 1"""
        assert(len(self.n) > k_coarse + 1)
        nCoarse = self.n[k_coarse]
        nFine = self.n[k_coarse + 1]

        assert(len(z_coarse) == nCoarse - 1)
        z_fine = np.zeros((nFine-1,), dtype=np.double)
        
        for i in range(0, nFine-1):
            if i % 2 == 0:
                rightNode = i // 2
                leftNode = rightNode - 1
                rightVal = z_coarse[rightNode] \
                    if rightNode < nCoarse - 1 \
                    else 0
                leftVal = z_coarse[leftNode] \
                    if leftNode >= 0         \
                    else 0
                
                z_fine[i] = 0.5 * rightVal + 0.5 * leftVal
            else:
                z_fine[i] = z_coarse[(i - 1) // 2]
        
        return z_fine

    def fineToCoarse(self, k_fine, z_fine):
        """Bring k_fine -> k_fine - 1"""
        assert(k_fine >= 1)
        nFine = self.n[k_fine]
        nCoarse = self.n[k_fine - 1]

        assert(len(z_fine) == nFine - 1)
        z_coarse = np.zeros((nCoarse-1), dtype=np.double)

        for i in range(0, nCoarse-1):
            centerNode = i * 2 + 1
            leftNode = centerNode - 1
            rightNode = centerNode + 1

            # TODO: figure out correct interpolation operator
            z_coarse[i] = z_fine[centerNode]
            z_coarse[i] += 0.5 * z_fine[leftNode]
            z_coarse[i] += 0.5 * z_fine[rightNode]
        
        return z_coarse

    def solveMG(self, f: 'function', r=1, p=1):
        kMax = len(self.n)
        # or we'd not having any basis functions
        assert(self.n[0] >= 2)

        self.U[0] = np.zeros((self.n[0]-1), dtype=np.double)
        for k in range(0, kMax):
            if k != 0:
                self.U[k] = self.coarseToFine(k-1, self.U[k-1])
            
            # construct F_k
            F_k = np.zeros((self.n[k]-1, ), dtype=np.double)
            for i in range(0, self.n[k]-1):
                F_k[i] = self.phiFIntg(k, i+1, f)

            self.U[k] = self.MG(k, self.U[k], F_k, p)
            # self.U[k] = self.MGDirect(k, self.U[k], F_k, p)
            # self.U[k] = self.MGCGDirect(k, self.U[k], F_k, p)

        return self.U[kMax - 1]

    def MG(self, k, z_0, g, p=1):
        """k'th MultiGrid; z_0: """
        n = self.n[k]
        assert(len(z_0) == n-1)
        assert(len(g) == n-1)
        z = z_0

        if k == 0:
            # A_0 ans = g
            # solve with CG
            A = np.zeros((n-1, n-1), dtype=np.double)
            for i, j in self.nonZeroEntries(k):
                A[i, j] = self.phiDerivPhiDerivIntg(k, i + 1, j + 1)
            res = np.linalg.solve(A, g)
            assert(np.allclose(A @ res, g))
            return res

        # Presmoothing
        z = self.iterateCG(k, z, g, 1)

        # Error correction
        g_bar = g - self.computeMatrixProduct(k, z)
        g_bar = self.fineToCoarse(k, g_bar)

        q = np.zeros((self.n[k-1] - 1,), dtype=np.double)
        for i in range(0, p):
            q = self.MG(k-1, q, g_bar, p)
        
        q_fine = self.coarseToFine(k-1, q)
        z += q_fine
        
        # Postsmoothing
        z = self.iterateCG(k, z, g, 1)

        return z
    
    def MGDirect(self, k, z_0, g, p=1):
        """Solve using direct method, as reference"""
        n = self.n[k]
        A = np.zeros((n-1, n-1), dtype=np.double)
        for i, j in self.nonZeroEntries(k):
            A[i, j] = self.phiDerivPhiDerivIntg(k, i + 1, j + 1)
        
        z = np.linalg.solve(A, g)
        assert(np.allclose(A @ z, g))
        return z

    def MGCGDirect(self, k, z_0, g, p=1):
        """Solve using iterateCG, as reference"""
        z = z_0
        z = self.iterateCG(k, z, g, 10)

        return z

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
    
    def estimateError(self, k, u_ref, err_samples=1000):
        T, res = self.genSample(k, err_samples)
        resReal = np.zeros((err_samples,), dtype=np.double)

        for i in range(0, err_samples):
            resReal[i] = u_ref(T[i])
        
        err = np.abs(resReal - res)
        
        err_L1 = 0
        err_Linf = 0
        for i in range(1, err_samples):
            h = T[i] - T[i-1]
            # integrate using trapz scheme
            err_L1 += h * abs(0.5 * err[i] + 0.5 * err[i-1])
            if err_Linf < abs(err[i]):
                err_Linf = abs(err[i])
        
        return (err_L1, err_Linf)



def plot_ref(ax, ref_func):
    x = np.linspace(0, 1, 100)
    ref = [val for val in map(ref_func, x)]

    ax.plot(x, ref)

if __name__ == '__main__':
    f = lambda x: (x-1) * math.sin(x)
    u_ref = lambda x: (x-1) * math.sin(x) + 2 * math.cos(x) + (2 - 2 * math.cos(1))* x - 2

    fig, ax = plt.subplots()

    # CGFEM = linearFEM(gaussIntg)
    # CGFEM.build_knots(5)
    # CGFEM.solveCG(0, f)
    # CGFEM.plot(0, 50, ax)

    # directFEM = linearFEM(gaussIntg)
    # directFEM.build_knots(5)
    # directFEM.solve(0, f)
    # directFEM.plot(0, 50, ax)

    previous_n = None
    previous_L1 = None
    previous_Linf = None
    previous_mg_solve_time = None
    MG_n_used = []
    for numLayers in range(2, 11):
        MGFEM = linearFEM(gaussIntg)
        MGFEM.build_knots_multigrid(numLayers, 2)
        MGFEM_nMax = MGFEM.n[-1]
        MGFEM_kMax = len(MGFEM.n) - 1
        assert(numLayers == MGFEM_kMax + 1)

        mg_solve_start = time.perf_counter()
        MGFEM.solveMG(f, 1, 1)
        mg_solve_end = time.perf_counter()
        mg_solve_time = mg_solve_end - mg_solve_start
        
        err_L1, err_Linf = MGFEM.estimateError(MGFEM_kMax, u_ref, MGFEM_nMax * 3 * 7)
        #MGFEM.plot(MGFEM_kMax, 1000, ax)

        err_L1_order = 0 \
            if previous_L1 is None \
            else math.log(previous_L1 / err_L1, MGFEM_nMax / previous_n)
        err_Linf_order = 0 \
            if previous_Linf is None \
            else math.log(previous_Linf / err_Linf, MGFEM_nMax / previous_n)
        mg_solve_order = 0 \
            if previous_mg_solve_time is None \
            else math.log(mg_solve_time / previous_mg_solve_time, MGFEM_nMax / previous_n)

        print(f"| {numLayers} | {MGFEM_nMax} | {err_L1:.3e} | {err_Linf:.3e} | {err_L1_order:.3f} | {err_Linf_order:.3f} | {mg_solve_time:.3e} | {mg_solve_order:.3f} |")
        previous_n = MGFEM_nMax
        previous_L1 = err_L1
        previous_Linf = err_Linf
        previous_mg_solve_time = mg_solve_time
        MG_n_used.append(MGFEM_nMax)

    previous_n = None
    previous_L1 = None
    previous_Linf = None
    previous_cg_solve_time = None
    for n in MG_n_used:
        CGFEM = linearFEM(gaussIntg)
        CGFEM.build_knots(n)

        cg_solve_start = time.perf_counter()
        CGFEM.solveCG(0, f, 1e-7)
        cg_solve_end = time.perf_counter()
        cg_solve_time = cg_solve_end - cg_solve_start

        err_L1, err_Linf = CGFEM.estimateError(0, u_ref, n * 3 * 7)
        #MGFEM.plot(MGFEM_kMax, 1000, ax)

        err_L1_order = 0 \
            if previous_L1 is None \
            else math.log(previous_L1 / err_L1, n / previous_n)
        err_Linf_order = 0 \
            if previous_Linf is None \
            else math.log(previous_Linf / err_Linf, n / previous_n)
        cg_solve_order = 0 \
            if previous_cg_solve_time is None \
            else math.log(cg_solve_time / previous_cg_solve_time, n / previous_n)

        print(f"| {n} | {err_L1:.3e} | {err_Linf:.3e} | {err_L1_order:.3f} | {err_Linf_order:.3f} | {cg_solve_time:.3e} | {cg_solve_order:.3f} |")
        previous_n = n
        previous_L1 = err_L1
        previous_Linf = err_Linf
        previous_cg_solve_time = cg_solve_time

    #plot_ref(ax, u_ref)

    #plt.show()