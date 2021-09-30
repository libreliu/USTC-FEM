import numpy as np
import matplotlib.pyplot as plt
import math

def simpsonIntg(f, x0, x1):
    """Do numeric intergral via Simpson's formula"""
    return (x1 - x0) * ((1.0/6) * f(x0) + (4.0/6) * f((x0 + x1) / 2) + (1.0/6) * f(x1))

def linearIntg(f, x0, x1):
    return (x1 - x0) * (0.5 * f(x0) + 0.5 * f(x1))

class PiccewiseQuaderaticFEM:
    def __init__(self, intg: 'function'):
        self.numIntg = intg
        self.U = None
    
    def build_knots(self, n: int):
        self.X = np.linspace(0, 1, num=n+1)
        self.n = n
        self.h = self.X[1] - self.X[0]
        return self.X
    
    def solve(self, f: 'function'):
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
                    K[i, j] = K[j, i] = -1.0 / (3 * h)
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
                F[i] = self.numIntg(
                    lambda t, i=iPhi: f(t) * ((2*t-X[i]-X[i-1])*(t*X[i-1])/(h**2)),
                    X[iPhi - 1],
                    X[iPhi]
                ) + self.numIntg(
                    lambda t, i=iPhi: f(t) * ((2*t-X[i]-X[i+1])*(t*X[i+1])/(h**2)),
                    X[iPhi],
                    X[iPhi + 1]
                )

        self.U = np.linalg.solve(K, F)

        chk = K @ self.U
        assert(np.allclose(chk, F))

        return self.U

    def plot(self, n, ax):
        T = np.linspace(0, 1, 2*self.n+1)
        res = np.zeros((self.n*2+1,), dtype=np.double)
        for i in range(0, self.n*2+1):
            if i == 0 or i == self.n*2:
                res[i] = 0
            else:
                res[i] = self.U[i - 1]

        ax.scatter(T, res)

        # plot the curve with n points
        T = np.linspace(0, 1, n)
        res = np.ndarray((n,), dtype=np.double)
        for i in range(0, n):
            gridCoord = i / self.h
            leftEndpoint = math.floor(gridCoord)
            rightEndpoint = math.ceil(gridCoord)

            


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
    
    def plot(self, n, ax):
        T = self.X
        res = np.zeros((self.n+1,), dtype=np.double)
        for i in range(0, self.n + 1):
            if i == 0 or i == self.n:
                res[i] = 0
            else:
                res[i] = self.U[i - 1]

        ax.plot(T, res)

def plot_ref(ax):
    ref_func = lambda x: (x-1) * math.sin(x) + 2 * math.cos(x) + (2 - 2 * math.cos(1)) * x - 2

    x = np.linspace(0, 1, 100)
    ref = [val for val in map(ref_func, x)]

    ax.plot(x, ref)
    plt.show()

def evaluate():
    f = lambda x: (x - 1) * math.sin(x)

    for n in [3]:
        linearFEM = PiecewiseLinearFEM(linearIntg)
        linearFEM.build_knots(n)
        linearFEM.solve(f)

        fig, ax = plt.subplots()
        ax.set(title=f'Piecewise linear, n = {n}')
        linearFEM.plot(n, ax)
        plot_ref(ax)
        
        plt.show()

    #for n in [3, 5, 8, 10, 20, 40, 80]:
    for n in [3]:
        quadFEM = PiccewiseQuaderaticFEM(simpsonIntg)
        quadFEM.build_knots(n)
        quadFEM.solve(f)

        fig, ax = plt.subplots()
        ax.set(title=f'Piecewise Quadratic, n = {n}')
        quadFEM.plot(n, ax)
        plot_ref(ax)
        
        plt.show()


if __name__ == '__main__':
    #plot_ref()
    evaluate()