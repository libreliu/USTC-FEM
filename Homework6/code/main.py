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

    def build_shishkin_knots(self, N: int, epsilon: float, n_respective: int):
        """Build shinshkin knots"""
        tau = 1 - 2 * epsilon * math.log(N)
        assert(tau < 1 and tau > 0)
        leftGrid = np.linspace(0, tau, num=n_respective+1)
        rightGrid = np.linspace(tau, 1, num=n_respective+1)
        self.X = np.hstack([leftGrid, rightGrid[1:]])
        self.n = 2 * n_respective
        assert(len(self.X) == self.n + 1)
        
    def phiDeriv(self, x, i):
        if self.X[i-1] <= x <= self.X[i]:
            h_i = self.X[i] - self.X[i-1]
            return 1.0 / h_i
        elif self.X[i] <= x <= self.X[i+1]:
            h_iP1 = self.X[i+1] - self.X[i]
            return -1.0 / h_iP1
        else:
            return 0

    def phi(self, x, i):
        assert(1 <= i <= self.n - 1)
        if self.X[i-1] <= x <= self.X[i]:
            h_i = self.X[i] - self.X[i-1]
            return (x - self.X[i-1]) / h_i
        elif self.X[i] <= x <= self.X[i+1]:
            h_iP1 = self.X[i+1] - self.X[i]
            return (self.X[i+1] - x) / h_iP1
        else:
            return 0

    def phiDerivPhiDerivIntg(self, i, j):
        """Calculate \int_0^1 \phi'_i(x) \phi'_j(x) dx
        Here we use i = 1, ..., n-1
        """
        if i > j:
            # This is symmetry so no problem
            return self.phiDerivPhiDerivIntg(j, i)

        assert(i <= j)
        if j - i >= 2:
            return 0
        elif j - i == 1:
            return self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i],
                self.X[i+1]
            )
        else:
            assert(j == i)
            return self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i-1],
                self.X[i]
            ) + self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i],
                self.X[i+1]
            )
        
    def phiDerivDeltaPhiDerivIntg(self, i, j):
        """Calculate \int_0^1 \delta \phi'_i(x) \phi'_j(x) dx
        Here we use i = 1, ..., n-1
        """
        if i > j:
            # This is symmetry so no problem
            return self.phiDerivDeltaPhiDerivIntg(j, i)

        assert(i <= j)
        if j - i >= 2:
            return 0
        elif j - i == 1:
            return (self.X[i+1] - self.X[i]) * self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i],
                self.X[i+1]
            )
        else:
            assert(j == i)
            return (self.X[i] - self.X[i-1]) * self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i-1],
                self.X[i]
            ) + (self.X[i+1] - self.X[i]) * self.numIntg(
                lambda x: self.phiDeriv(x, i) * self.phiDeriv(x, j),
                self.X[i],
                self.X[i+1]
            )

    def phiDerivPhiIntg(self, i, j):
        """
        Calculate \int_0^1 \phi'_i(x) \phi_j(x) dx
        Here we use i = 1, ..., n-1
        """
        if i <= j:
            if j - i >= 2:
                return 0
            elif j - i == 1:
                return self.numIntg(
                    lambda x: self.phiDeriv(x, i) * self.phi(x, j),
                    self.X[i],
                    self.X[i+1]
                )
            else:
                assert(j == i)
                return self.numIntg(
                    lambda x: self.phiDeriv(x, i) * self.phi(x, j),
                    self.X[i-1],
                    self.X[i]
                ) + self.numIntg(
                    lambda x: self.phiDeriv(x, i) * self.phi(x, j),
                    self.X[i],
                    self.X[i+1]
                )
        else:
            assert(i > j)
            if i - j >= 2:
                return 0
            elif i - j == 1:
                return self.numIntg(
                    lambda x: self.phiDeriv(x, i) * self.phi(x, j),
                    self.X[j],
                    self.X[j+1]
                )

    def RHSIntg(self, i, f):
        """
        Calculate \int_0^1 f(x) (\phi_i(x) + \delta * \phi_i'(x)) dx
        Here we use i = 1, ..., n-1
        """
        return self.numIntg(
            lambda x: f(x) * (self.phi(x, i) + (self.X[i] - self.X[i-1]) * self.phiDeriv(x, i)),
            self.X[i-1],
            self.X[i]
        ) + self.numIntg(
            lambda x: f(x) * (self.phi(x, i) + (self.X[i+1] - self.X[i]) * self.phiDeriv(x, i)),
            self.X[i],
            self.X[i+1]
        )

    def solve(self, f: 'function', epsilon: float):
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
            for j in range(0, n-1):
                K[i, j] = epsilon * self.phiDerivPhiDerivIntg(i + 1, j + 1) + self.phiDerivPhiIntg(j + 1, i + 1) + self.phiDerivDeltaPhiDerivIntg(j + 1, i + 1)
        
        #np.set_printoptions(precision=3, linewidth=150)
        #print(K)

        # Build F
        F = np.zeros((n-1, ), dtype=np.double)
        for i in range(0, n-1):
            F[i] = self.RHSIntg(i+1, f)
        
        #print(F)

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
        for i in range(0, n):
            t = T[i]
            for j in range(1, self.n):
                res[i] += self.U[j - 1] * self.phi(t, j)
        
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

def evaluate(use_shishkin: bool, n_list, epsilon):
    f = lambda x: x
    if epsilon >= 1e-2:
        ref_func = lambda x: (2*epsilon+1)/(2-2*math.exp(1/epsilon))*(math.exp(x/epsilon)-1) + 0.5 * x**2 + epsilon*x
    else:
        # TODO: make small epsilon approximation
        ref_func = lambda x: -math.exp((x-1)/epsilon)*(epsilon+1/2) + 0.5 * x**2 + epsilon*x

    fig, axs = plt.subplots(1, len(n_list), figsize=(17,3))
    errs = []
    err_L1s = []
    err_Linfs = []
    for idx, n in enumerate(n_list):
        linearFEM = PiecewiseLinearFEM(gaussIntg)
        if use_shishkin:
            linearFEM.build_shishkin_knots(n, epsilon, n // 2)
        else:
            linearFEM.build_knots(n)
        linearFEM.solve(f, epsilon)

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
        
        #print(f"{PiecewiseLinearFEM.__name__}: n={n}, err_L1: {err_L1:.3e}, err_Linf: {err_Linf:.3e}")
        err_L1s.append(err_L1)
        err_Linfs.append(err_Linf)

    print("| N | $ L^1 $ error | order | $ L^\infty $ error | order |")    
    print("| - | ------------- | ----- | ------------------ | ----- |")
    for idx, n in enumerate(n_list):
        err_Linf_order = 0 if idx == 0 else math.log(err_Linfs[idx - 1] / err_Linfs[idx], n_list[idx] / n_list[idx-1])
        err_L1_order = 0 if idx == 0 else math.log(err_L1s[idx - 1] / err_L1s[idx], n_list[idx] / n_list[idx-1])

        print(f"|{n}|{err_L1s[idx]:.3e}|{err_L1_order:.5e}|{err_Linfs[idx]:.3e}|{err_Linf_order:.5e}|")
    
    print("")
    print("下面为此时解的情况：")
    print(f"![](Homework6.assets/eps_{eps}_{'shishkin' if use_shishkin else 'uniform'}.png)")
    print("下面为此时误差的情况：")
    print(f"![](Homework6.assets/eps_{eps}_{'shishkin' if use_shishkin else 'uniform'}_errors.png)")


    for ax in axs.flat:
        ax.label_outer()
    #plt.show()
    plt.savefig(f"eps_{eps}_{'shishkin' if use_shishkin else 'uniform'}.png")

    fig, ax = plt.subplots()
    for err in errs:
        ax.plot(T, err)
    #plt.show()
    plt.savefig(f"eps_{eps}_{'shishkin' if use_shishkin else 'uniform'}_errors.png")

if __name__ == '__main__':
    #plot_ref()
    
    for eps in [1e-1, 1e-3, 1e-5, 1e-7]:
        print(f"### $ \epsilon = {eps} $, Uniform")
        evaluate(False, [2,4,8,16,32, 64, 128], eps)
        print("")
        print(f"### $ \epsilon = {eps} $, shishkin")
        evaluate(True, [2,4,8,16,32, 64, 128], eps)
        print("")
