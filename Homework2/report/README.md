# 有限元分析 第二次程序作业

刘紫檀 SA21229063

## 问题描述

给定耗散-反应方程
$$
\left\{
\begin{aligned}
&-(d(x)u'(x))' + c(x)u(x) = f(x) & x \in [0, 1] \\
& u(0) = u(1) = 0
\end{aligned}
\right. \\
\text{where} \quad
\left\{
\begin{aligned}
& 0 < \beta < d(x) \\
& c(x) < \alpha < \infty
\end{aligned}
\right.
\quad
\text{for} \quad x \in (0, 1)
$$
分别下面两种有限维空间中进行方程的数值求解
$$
V_{h1} = \{v \in C^0([0, 1]) \quad v \vert _{I_j} = P^1(I_j), v(0) = v(1) = 0\} \\
V_{h2} = \{v \in C^0([0, 1]) \quad v \vert _{I_j} = P^2(I_j), v(0) = v(1) = 0\}
$$

测试时带入
$$
\left\{
\begin{aligned}
&d(x) = \sin x + 2 \\
&c(x) = x^2 + 1 \\
&u(x) = x(x-1) \\
&f(x) = x (x - 1) (x^2 + 1) - 2 (\sin x + 2) - (2 x - 1) \cos x
\end{aligned}
\right.
$$
验证即可

## 问题求解

### 变分·约定

首先选取原问题的合适的解函数空间。

因为涉及解函数的二阶导数，所以选择试函数 $v$（和解函数 $u $ 的）空间 $ V= W_2^1([0, 1]) $，否则不一定可以做分部积分。

> 其实没那多事，毕竟 $ P^1 $ 和 $ P^2 $ 都好的很。

$$
\begin{aligned}
 &-(d(x)u'(x))' + c(x)u(x) = f(x) \\
\Rightarrow & \int_0^1 -(d(x)u'(x))'v(x) dx + \int_0^1 c(x)u(x) v(x) dx = \int_0^1 f(x) v(x) dx \\
\Rightarrow & -d(x)u'(x)v(x) \bigg | _0^1  + \int_0^1 d(x) u'(x) v'(x) dx + \int_0^1 c(x)u(x) v(x) dx  = \int_0^1 f(x) v(x) dx \\
\Rightarrow & \int_0^1 d(x) u'(x) v'(x) dx + \int_0^1 c(x)u(x) v(x) dx  = \int_0^1 f(x) v(x) dx

\end{aligned}
$$

然后我们有
$$
a(u, v) = \int_0^1 d(x) u'(x) v'(x) dx + \int_0^1 c(x) u(x) v(x) dx \\
F(v) = \int_0^1 f(x) v(x) dx
$$
观察到
$$
\left\{
\begin{aligned}
a(u, v) &= a(v, u) \\
a(u, v) &= \int_0^1 d(x) u'(x) v'(x) dx + \int_0^1 c(x) u(x) v(x) dx \\
&= \int_0^1 -(d(x)u'(x))'v(x) dx + \int_0^1 c(x) u(x) v(x) dx \\
&= \int_0^1 |(d(x)u'(x))'v(x)| dx + \int_0^1 |c(x) u(x) v(x)| dx \\
&\le ||(du')'||_{L^2} ||v||_{L^2} + ||c||_{L^2} ||u||_{L^2} ||v||_{L^2}

\end{aligned}
\right.
$$





### $ V_{h1} $ 空间下的矩阵形式



### $ V_{h2} $ 空间下的矩阵形式



### 代码实现

参考 `main.py` ，`PiecewiseLinearFEM` 为 $ V_{h1} $ 空间下的有限元求解，`PiecewiseQuaderaticFEM` 为 $ V_{h2} $ 空间下的有限元求解。

积分采用 Simpson 积分公式，解线性方程组采用 `np.linalg.solve`，其内部采用 LAPACK 的 `?gesv` 函数，通用实现为 LU 分解后直接法计算得到。

#### 运行说明

1. 确保 Python 不要太旧（我用的是 Python 3.8.3），如果没有 Python 的话，Windows 用户可以去 Microsoft Store 安装一个。

2. 安装 `pip` （必须）并安装 `virtualenv`（可选）。

3. 创建 `virtualenv` 并进入创建好的虚拟环境（可选）。

4. 安装 `scipy`, `numpy` 和 `matplotlib` 三个包

   ```bash
   # 如果没有安装 virtualenv 或没有启用，推荐装在用户目录下
   python -m pip install --user scipy numpy matplotlib
   # 如果安装并创建了虚拟环境，激活之，并在虚拟环境下安装即可
   pip install scipy numpy matplotlib
   ```

5. 运行 `main.py`: `python ./code/main.py`

#### $ V_{h1} $ 实现


#### $ V_{h2} $ 实现



## 一些有趣的问题

