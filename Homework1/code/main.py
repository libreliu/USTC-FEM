import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

def simpsonIntg(f, x0, x1):
    """Do numeric intergral via Simpson's formula"""
    return (x1 - x0) * ((1.0/6) * f(x0) + (4.0/6) * f((x0 + x1) / 2) + (1.0/6) * f(x1))

def linearIntg(f, x0, x1):
    return (x1 - x0) * (0.5 * f(x0) + 0.5 * f(x1))

def preciseIntg(f, x0, x1):
    res = integrate.quad(f, x0, x1)
    #print(f"Integrate from {x0} to {x1} yield {res[0]}, error={res[1]}")
    return res[0]

class PiecewiseQuaderaticFEM:
    def __init__(self, intg: 'function'):
        self.numIntg = intg
        self.U = None
    
    def build_knots(self, n: int):
        self.X = np.linspace(0, 1, num=n+1)
        self.n = n
        self.h = self.X[1] - self.X[0]
        return self.X
    
    def solve(self, f: 'function', verbose=False):
        X = self.X
        n = self.n

        assert(len(X) == n+1)
        h = X[1] - X[0]
        K = np.zeros((2*n-1, 2*n-1), dtype=np.double)
        for i in range(0, 2*n-1):
            for j in range(i, 2*n-1):
                if i == j and i % 2 == 0:
                    K[i, j] = 16 / (3 * h)
                elif i == j and i % 2 != 0:
                    K[i, j] = 14 / (3 * h)
                elif j == i + 1:
                    K[i, j] = K[j, i] = - 8.0 / (3 * h)
                elif j == i + 2 and i % 2 == 0:
                    K[i, j] = K[j, i] = 0
                elif j == i + 2 and i % 2 != 0:
                    K[i, j] = K[j, i] = 1.0 / (3 * h)
                else:
                    # matrix already initialized
                    pass
        
        F = np.zeros((2*n-1,), dtype=np.double)
        for i in range(0, 2*n-1):
            if i % 2 == 0:
                iPsi = i // 2
                print(iPsi)
                F[i] = self.numIntg(
                    lambda t, i=iPsi: f(t) * (-4*(t-X[i])*(t-X[i+1]) / (h**2)),
                    X[iPsi],
                    X[iPsi + 1]
                )
            else:
                iPhi = (i + 1) // 2
                print(iPhi)
                left = self.numIntg(
                    lambda t: f(t) * ((2*t-X[iPhi]-X[iPhi-1])*(t-X[iPhi-1])/(h**2)),
                    X[iPhi - 1],
                    X[iPhi]
                )
                right = self.numIntg(
                    lambda t: f(t) * ((2*t-X[iPhi]-X[iPhi+1])*(t-X[iPhi+1])/(h**2)),
                    X[iPhi],
                    X[iPhi + 1]
                )
                print(f"left={left} right={right}")
                F[i] = left + right

        print("K eigs:")
        print(np.linalg.eigvals(K))
        self.U = np.linalg.solve(K, F)

        if verbose:
            print("K:")
            print(K)
            
            print("F:")
            print(F)

            print("U:")
            print(self.U)

        chk = K @ self.U
        assert(np.allclose(chk, F))

        return self.U

    def genSample(self, n):
        T = np.linspace(0, 1, n)
        res = np.zeros((n,), dtype=np.double)
        for i in range(0, n):
            t = T[i]
            gridCoord = t / self.h
            leftEndpoint = math.floor(gridCoord)
            rightEndpoint = math.ceil(gridCoord)

            # right phi arm
            if leftEndpoint < self.n and leftEndpoint >= 1:
                res[i] += self.U[leftEndpoint * 2 - 1] * ((2*t-self.X[leftEndpoint]-self.X[leftEndpoint+1])*(t-self.X[leftEndpoint+1])/(self.h**2))

            # left phi arm
            if leftEndpoint != rightEndpoint and rightEndpoint > 0 and rightEndpoint < self.n:
                # right phi
                res[i] += self.U[rightEndpoint * 2 - 1] * ((2*t-self.X[rightEndpoint]-self.X[rightEndpoint-1])*(t-self.X[rightEndpoint-1])/(self.h**2))

            # psi in [left, right]
            if leftEndpoint < self.n:
                res[i] += self.U[leftEndpoint * 2] * (-4*(t-self.X[leftEndpoint])*(t-self.X[leftEndpoint+1]) / (self.h**2))
        
        return (T, res)

    def plot(self, n, ax, verbose=False):
        T = np.linspace(0, 1, 2*self.n+1)
        res = np.zeros((self.n*2+1,), dtype=np.double)
        for i in range(0, self.n*2+1):
            if i == 0 or i == self.n*2:
                res[i] = 0
            else:
                res[i] = self.U[i - 1]

        ax.scatter(T, res)

        # plot the curve with n points
        T, res = self.genSample(n)
        
        if verbose:
            print("T:") 
            print(T)
            print("res:")
            print(res)

        ax.plot(T, res)


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

    def solve(self, f: 'function'):
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
                if j >= i + 2:
                    K[i, j] = K[j, i] = 0
                elif j == i + 1:
                    # No outbound access since we only calculates upper triangular
                    K[i, j] = K[j, i] = -( 1.0 / ((X[i+1 + 1] - X[i + 1])))
                elif j == i:
                    # handle cases for outbound access
                    K[i, j] = (1.0 / ((X[i+1 + 1] - X[i + 1]))) + (1.0 / ((X[i + 1] - X[i-1 + 1])))
                else:
                    assert(False)
        
        # Build F
        F = np.zeros((n-1, ), dtype=np.double)
        for i in range(0, n-1):
            F[i] = self.numIntg(
                lambda t, i=i: f(t) * (t - X[i-1 + 1]) / (X[i + 1] - X[i-1 + 1]),
                X[i-1 + 1],
                X[i + 1]
            ) + self.numIntg(
                lambda t, i=i: f(t) * (X[i+1 + 1] - t) / (X[i+1 + 1] - X[i + 1]),
                X[i + 1],
                X[i+1 + 1]
            )

        
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

def plot_ref(ax):
    ref_func = lambda x: (x-1) * math.sin(x) + 2 * math.cos(x) + (2 - 2 * math.cos(1)) * x - 2

    x = np.linspace(0, 1, 100)
    ref = [val for val in map(ref_func, x)]

    ax.plot(x, ref)

def evaluate():
    f = lambda x: (x - 1) * math.sin(x)
    ref_func = lambda x: (x-1) * math.sin(x) + 2 * math.cos(x) + (2 - 2 * math.cos(1)) * x - 2

    errs = []
    n_list = [3, 5, 10, 20, 40, 80]
    fig, axs = plt.subplots(1, len(n_list))
    for idx, n in enumerate(n_list):
        linearFEM = PiecewiseLinearFEM(simpsonIntg)
        linearFEM.build_knots(n)
        linearFEM.solve(f)

        axs[idx].set(title=f'n = {n}')
        linearFEM.plot(50, axs[idx])
        plot_ref(axs[idx])

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

    errs = []

    # for n in [3, 5, 8, 10, 20, 40, 80]:
    n_list = [1, 2, 3, 5, 8]
    fig, axs = plt.subplots(1, len(n_list))
    for idx, n in enumerate(n_list):
        quadFEM = PiecewiseQuaderaticFEM(simpsonIntg)
        quadFEM.build_knots(n)
        quadFEM.solve(f, verbose=True)

        axs[idx].set(title=f'n = {n}')
        quadFEM.plot(50, axs[idx])
        plot_ref(axs[idx])

        # calculate error w.r.t true value
        # err_samples = 1000
        # T, res = quadFEM.genSample(err_samples)
        # resReal = np.zeros((err_samples,),dtype=np.double)
        # for i in range(0, err_samples):
        #     resReal[i] = ref_func(T[i])

        # err = np.abs(resReal - res)
        # errs.append(err)
    for ax in axs.flat:
        ax.label_outer()
    plt.show()

    # fig, ax = plt.subplots()
    # for err in errs:
    #     ax.plot(T, err)
    # plt.show()


if __name__ == '__main__':
    #plot_ref()
    evaluate()