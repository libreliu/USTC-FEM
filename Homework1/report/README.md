# 有限元分析 第一次程序作业

刘紫檀 SA21229063

## 问题描述

给定方程
$$
\left\{
\begin{aligned}
&-u'' = f \\
& u(0) = u(1) = 0
\end{aligned}
\right.
$$
分别下面两种有限维空间中进行方程的数值求解
$$
V_{h1} = \{v \in C^0([0, 1]) \quad v \vert _{I_j} = P^1(I_j), v(0) = v(1) = 0\} \\
V_{h2} = \{v \in C^0([0, 1]) \quad v \vert _{I_j} = P^2(I_j), v(0) = v(1) = 0\}
$$

## 问题求解

### 真实解

$$
-u'' = f \quad \Rightarrow \quad -u(x) = \int_{x_0}^\eta\int_{x_0}^\mu f(t) dt d\eta + Cx + D
$$

对 $ f(x) = (x-1) \sin x $ 我们有
$$
\begin{aligned}
-u(x) &= \int_1^\eta\int_1^\mu f(t) dt d\eta + Cx + D \\
     &= -x \sin x + \sin x - 2 \cos x + Cx + D
\end{aligned}
$$

带入定解条件，我们得到
$$
\left\{
\begin{aligned}
&u(0) = -2 + D = 0 \Rightarrow D = 2 \\
&u(1) = - 2 \cos 1 + C + D = 0 \Rightarrow C = -2 + 2 \cos 1\\
\end{aligned}
\right.
$$
所以解为
$$
u(x) = (x-1) \sin x + 2 \cos x + (2 - 2 \cos 1) x - 2
$$



### 变分·约定

我们首先构造变分问题
$$
\begin{aligned}
\int_0^1 -u''(x) v(x) dx &=  -u'(x)v(x)  \bigg \vert^1_0  + \int_0^1 u'(x)v'(x)dx \\
&=  \int_0^1 u'(x)v'(x)dx \\
&= \int_0^1 f(x)v(x) dx \\
\end{aligned}
$$


我们约定对解空间离散后的节点分别为 $ 0 = x_0 < x_1 \dots < x_{n-1} < x_n = n $，共 $ n + 1 $ 个节点。记 $ h_i =  x_{i} - x_{i-1} $

对于有限维空间 $ V_{h} $，可以选择一组基函数 $ \{ \phi_i(x) \}_{i=0}^n $ 。此时，变分问题可以转换为有限维变分问题。设 $ u_H, \in V_{h1} $ ，我们有
$$
u_H = \sum_j u_j \phi_j \\
$$
同时，由于 $ v_H \in V_{h1}$  可以表示成 $ \{ \phi_i(x) \}_{i=0}^n $  的线性组合，所以不妨直接将 $ v_H $ 带入 $ \phi_0, \phi_1, ..., \phi_n $ ，即可得到下面的形式
$$
\begin{aligned}
\mathrm{LHS} 
&= \int_0^1 u_H'(x) \phi_i'(x) dx =
\int_0^1 (\sum_j u_j \phi_j'(x)) \phi_i'(x) dx \\
&= \sum_j u_j \int_0^1 \phi_j'(x) \phi_i'(x) dx \\

\\
\mathrm{RHS}
&= \int_0^1 f(x) \phi_i(x) dx
\end{aligned}
$$
记 $ K_{ij} = \int_0^1 \phi_j'(x) \phi_i'(x) dx $，$ U_j = u_j $ ，$ F_i =  \int_0^1 f(x) \phi_i(x) dx $，则可以写成矩阵形式
$$
KU = F
$$

### $ V_{h1} $ 空间下的矩阵形式

对于有限维空间 $ V_{h1} $，可以选择基函数 $ \{ \phi_i(x) \}_{i=0}^n $ 如下
$$
\phi_i(x) = \left\{
\begin{aligned}
&\frac{x-x_{i-1}}{x_{i}-x_{i-1}}, &x \in [x_{i-1}, x_{i}] \\
&\frac{x_{i+1} - x}{x_{i+1} - x_{i}}, &x \in [x_i, x_{i+1}] \\
&0, & x \notin [x_{i-1}, x_{i+1}]
\end{aligned}

\right.
$$
那么，我们可以计算 $\phi_i(x)$ 和 $ K_{ij} $ 如下
$$
\phi'_i(x) = \left\{
\begin{aligned}
&\frac{1}{x_i - x_{i-1}} = \frac{1}{h_i}  , & x \in (x_{i-1}, x_i)  \\
&\frac{1}{x_{i} - x_{i+1}} = -\frac{1}{h_{i+1}}, & x \in (x_i, x_{i+1}) \\
\end{aligned}
\right.

\\

K_{ij} = \int_0^1 \phi'_i(x) \phi'_j(x) dx = 
\left\{

\begin{aligned}
& 0 & | j - i | \ge 2 \\
& -\frac{1}{h_{i+1} h_j} & j = i + 1 \\
&\frac{1}{h_{i+1}^2} + \frac{1}{h_{i}^2} & j = i \cup j \ne 0 \cup j \ne n \\
&\frac{1}{h_{i+1}^2} & j = i = 0 \\
&\frac{1}{h_{i}^2} & j = i = n

\end{aligned}

\right.
$$

> 其实 $ \phi_i'(x) $ 是选定的 Sobolev 空间下的弱导数 $ D_w^1 \phi $

$ F_i $ 的计算通过数值积分进行。

### $ V_{h2} $ 空间下的矩阵形式



### 误差分析



### 代码实现



