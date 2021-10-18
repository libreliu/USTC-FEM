import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

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

class PiecewiseLinearFEM:
    def __init__(self, intg: 'function'):
        """
        @param intg: Numerical integrater to use
        """
        self.numIntg = intg
        self.U = None

    def build_knots(self, n: int):
        """Build uniform knots in generation, total n+1 knots"""
        self.X = np.linspace(0, 1, num=n+1)
        self.n = n
        return self.X

    def phiDerivIntg(self, i, j, d):
        """Calculate \int_0^1 d(x) \phi'_i(x) \phi'_j(x) dx
        Here we use i = 0, ..., n-1, as opposed to what we use in the document (starting from 1)
        """
        assert(i <= j)
        if j - i >= 2:
            return 0
        elif j - i == 1:
            return self.numIntg(
                lambda x: d(x) * (1.0 / (self.X[i+1] - self.X[i+2])) * (1.0 / (self.X[j+1] - self.X[j])),
                self.X[i+1],
                self.X[j+1]
            )
        else:
            assert(j == i)
            return self.numIntg(
                lambda x: d(x) * (1.0 / (self.X[i+1] - self.X[i])) * (1.0 / (self.X[i+1] - self.X[i])),
                self.X[i],
                self.X[i+1]
            ) + self.numIntg(
                lambda x: d(x) * (1.0 / (self.X[i+1] - self.X[i+2])) * (1.0 / (self.X[i+1] - self.X[i+2])),
                self.X[i+1],
                self.X[i+2]
            )
        

    def phiIntg(self, i, j, c):
        """
        Calculate \int_0^1 c(x) \phi_i(x) \phi_j(x) dx
        Here we use i = 0, ..., n-1, as opposed to what we use in the document (starting from 1)
        """
        assert(i <= j)
        if j - i >= 2:
            return 0
        elif j - i == 1:
            return self.numIntg(
                lambda x: c(x) * ((self.X[i+2] - x) / (self.X[i+2] - self.X[i+1])) * ((x - self.X[j]) / (self.X[j+1] - self.X[j])),
                self.X[i+1],
                self.X[j+1]
            )
        else:
            assert(j == i)
            return self.numIntg(
                lambda x: c(x) * ((x - self.X[i]) / (self.X[i+1] - self.X[i])) * ((x - self.X[i]) / (self.X[i+1] - self.X[i])),
                self.X[i],
                self.X[i+1]
            ) + self.numIntg(
                lambda x: c(x) * ((self.X[i+2] - x) / (self.X[i+2] - self.X[i+1])) * ((self.X[i+2] - x) / (self.X[i+2] - self.X[i+1])),
                self.X[i+1],
                self.X[i+2]
            )

    def phiFIntg(self, i, f):
        """
        Calculate \int_0^1 f(x) \phi_i(x) dx
        Here we use i = 0, ..., n-1, as opposed to what we use in the document (starting from 1)
        """
        return self.numIntg(
            lambda x: f(x) * (x - self.X[i]) / (self.X[i+1] - self.X[i]),
            self.X[i],
            self.X[i+1]
        ) + self.numIntg(
            lambda x: f(x) * (self.X[i+2] - x) / (self.X[i+2] - self.X[i+1]),
            self.X[i+1],
            self.X[i+2]
        )

    def solve(self, f: 'function', c: 'function', d: 'function'):
        """Solve using V_{h1} finite element space
        
        Uses n+1 knots in total, namely 0, ..., n
        and we have n-1 basis functions in total
        """
        X = self.X
        n = self.n

        # Check x
        assert(len(X) == n + 1)

        # Build K
        K = np.ndarray((n-1, n-1), dtype=np.double)
        for i in range(0, n-1):
            for j in range(i, n-1):
                K[j, i] = K[i, j] = self.phiIntg(i, j, c) + self.phiDerivIntg(i, j, d)
        
        np.set_printoptions(precision=3, linewidth=150)
        print(K)

        # Build F
        F = np.zeros((n-1, ), dtype=np.double)
        for i in range(0, n-1):
            F[i] = self.phiFIntg(i, f)
        
        print(F)

        # np.linalg.solve() uses _gesv LAPACK routines
        # which uses LU factorization indeed
        self.U = np.linalg.solve(K, F)

        # Check result
        chk = K @ self.U
        assert(np.allclose(chk, F))

        return self.U
    
    def genSample(self, n):
        T = np.linspace(0, 1, n)
        res = np.zeros((n,), dtype=np.double)
        h = self.X[1] - self.X[0]
        for i in range(0, n):
            t = T[i]
            gridCoord = t / h
            leftEndpoint = math.floor(gridCoord)
            rightEndpoint = math.ceil(gridCoord)

            if leftEndpoint > 0 and leftEndpoint < self.n:
                # right arm of left nodal basis
                res[i] += self.U[leftEndpoint - 1] * ((self.X[leftEndpoint + 1] - t) / h)

            if rightEndpoint > 0 and rightEndpoint < self.n and leftEndpoint != rightEndpoint:
                # left arm of right nodal basis
                res[i] += self.U[rightEndpoint - 1] * ((t - self.X[rightEndpoint - 1]) / h)
        
        return (T, res)

    def plot(self, n, ax):
        T = self.X
        res = np.zeros((self.n+1,), dtype=np.double)
        for i in range(0, self.n + 1):
            if i == 0 or i == self.n:
                res[i] = 0
            else:
                res[i] = self.U[i - 1]

        ax.scatter(T, res)

        # only works for uniform spacing
        T, res = self.genSample(n)
        ax.plot(T, res)

def plot_ref(ax, ref_func):

    x = np.linspace(0, 1, 100)
    ref = [val for val in map(ref_func, x)]

    ax.plot(x, ref)

def evaluate():
    f = lambda x: x*(x-1)*(x**2+1) - 2*(math.sin(x) + 2) - (2*x-1) * math.cos(x)
    #f = lambda x: - math.sin(2*x + 1) - 2 * (2 + math.sin(x)) + (x-1) *x* (x**2+1)
    d = lambda x: math.sin(x) + 2
    c = lambda x: x**2 + 1
    ref_func = lambda x: x * (x-1)

    errs = []
    #n_list = [3, 5, 10, 20, 40, 80]
    n_list = [3, 5, 10]
    fig, axs = plt.subplots(1, len(n_list))
    for idx, n in enumerate(n_list):
        linearFEM = PiecewiseLinearFEM(gaussIntg)
        linearFEM.build_knots(n)
        linearFEM.solve(f, c, d)

        axs[idx].set(title=f'n = {n}')
        linearFEM.plot(50, axs[idx])
        plot_ref(axs[idx], ref_func)

        # calculate error w.r.t true value
        err_samples = 1000
        T, res = linearFEM.genSample(err_samples)
        resReal = np.zeros((err_samples,),dtype=np.double)
        for i in range(0, err_samples):
            resReal[i] = ref_func(T[i])

        err = np.abs(resReal - res)
        errs.append(err)

        err_L1 = 0
        err_Linf = 0
        for i in range(1, err_samples):
            h = T[i] - T[i-1]
            # integrate using trapz scheme
            err_L1 += h * abs(0.5 * err[i] + 0.5 * err[i-1])
            if err_Linf < abs(err[i]):
                err_Linf = abs(err[i])
        
        print(f"linear n={n}, err_L1: {err_L1}, err_Linf: {err_Linf}")
    
    for ax in axs.flat:
        ax.label_outer()
    plt.show()

    fig, ax = plt.subplots()
    for err in errs:
        ax.plot(T, err)
    plt.show()


if __name__ == '__main__':
    #plot_ref()
    evaluate()