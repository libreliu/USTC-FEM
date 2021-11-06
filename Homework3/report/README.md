# 有限元分析 第三次程序作业

刘紫檀 SA21229063

## 问题描述

使用有限元方法求解如下的 PDE

$$
\begin{aligned}
&\left\{
\begin{aligned}
&-\epsilon u''(x) + u'(x) = x \\
&u(0) = u(1) = 0
\end{aligned}
\right.
\\

&\text{Using}\ V_h = \{ v \in H^1_0([0, 1]) \ v |_{i_j} \in P^1(I_j) \}
\end{aligned}
$$

分别在 $ \epsilon = 10^{-1} $ 和 $ \epsilon = 10^{-7} $ 上测试。

其中，网格分别选择如下两种：
1. 均匀网格
2. Shishkin 网格： 首先将 $[0, 1]$ 分为 $ [0, \tau] $ 和 $ [\tau, 1] $ 两部分（$ \tau = 1 - 2 \epsilon \ln N $），然后在两个区间分别使用均匀网格

## 原理·计算

### 真实解

用常数变异法可以得到如下的解

$$
u(x) = C_1 \epsilon e^{x/ \epsilon} + \frac{1}{2} x^2 + \epsilon x + C_2
$$

带入定解条件得到

$$
\begin{aligned}
u(x) &= - \frac{1/2 + \epsilon}{\epsilon(e^{1/\epsilon} - 1)} \epsilon e^{x/ \epsilon} + \frac{1}{2} x^2 + \epsilon x + \frac{1/2 + \epsilon}{e^{1/\epsilon} - 1} \\
&= \frac{1}{2} x^2 + \epsilon x - \frac{(1 + 2 \epsilon)(1-e^{x/\epsilon})}{2(1-e^{1/\epsilon})}
\end{aligned}
$$

### 当 $ \epsilon $ 很小时候的处理

$\epsilon$ 很小的时候， $e^{1/\epsilon}$ 很容易就会变成 NaN..

不过，

$$
\lim_{\epsilon \to 0} \frac{(1-e^{x/\epsilon})}{(1-e^{1/\epsilon})} = \lim_{\epsilon \to 0} x e^{\frac{x-1}{\epsilon}} = 0 \quad ( x \le 1 )
$$

所以不妨进行一些近似：

$$
\begin{aligned}
\frac{(1 + 2 \epsilon)(1-e^{x/\epsilon})}{2(1-e^{1/\epsilon})} &= \frac{(1 + 2 \epsilon)}{2} \frac{(1-e^{x/\epsilon})}{(1-e^{1/\epsilon})}\\
&\approx \frac{(1 + 2 \epsilon)}{2} \frac{e^{x/\epsilon}}{e^{1/\epsilon}} \\
&= \frac{(1 + 2 \epsilon)}{2} e^{(x-1)/\epsilon}

\end{aligned}
$$

### 变分

$$
\begin{aligned}
a(u, v) &= - \int_0^1 \epsilon u'' v dx + \int_0^1 u' v  dx \\
&= - \epsilon [u'v \bigg |_0^1 - \int_0^1 u'v' dx  ] + \int_0^1 u'v dx \\
&= \epsilon \int_0^1 u'v' dx + \int_0^1 u'v dx \\
&= \int_0^1 f v dx
\end{aligned}
$$

这个内积不是对称的，需要验证 Lax-Milgram 定理的条件以确保有限元唯一解。


### 有限元离散化

我们约定对解空间离散后的节点分别为 $ 0 = x_0 < x_1 \dots < x_{n-1} < x_n = n $，共 $ n + 1 $ 个节点。记 $ h_i =  x_{i} - x_{i-1} $

对于有限维空间 $ V_{h} $，可以选择一组基函数 $ \{ \phi_i(x) \}_{i=1}^{n-1} $ 。此时，变分问题可以转换为有限维变分问题。设 $ u_H \in V_{h} $ ，我们有
$$
u_H = \sum_j u_j \phi_j \\
$$
同时，由于 $ v_H \in  V_{h} $  可以表示成 $ \{ \phi_i(x) \}_{i=0}^n $  的线性组合，所以不妨直接将 $ v_H $ 带入 $ \phi_0, \phi_1, ..., \phi_n $ ，即可得到下面的形式
$$
\begin{aligned}
\mathrm{LHS} 
&= \epsilon \int_0^1 u_H'(x) \phi_i'(x) dx + \int_0^1 u_H'(x) \phi_i(x) dx \\
&= \epsilon \int_0^1 (\sum_j u_j \phi_j'(x)) \phi_i'(x) dx + \int_0^1 (\sum_j u_j \phi_j'(x)) \phi_i(x) dx  \\
&= \sum_j u_j \left( \epsilon\int_0^1 \phi_j'(x) \phi_i'(x) dx + \int_0^1 \phi_j'(x) \phi_i(x) dx\right) \\

\\
\mathrm{RHS}
&= \int_0^1 f(x) \phi_i(x) dx
\end{aligned}
$$
记 $ K_{ij} = \int_0^1 \phi_j'(x) \phi_i'(x) dx $，$ U_j = u_j $ ，$ F_i =  \int_0^1 f(x) \phi_i(x) dx $，则可以写成矩阵形式
$$
KU = F
$$

对于有限维空间 $ V_{h}  $，可以选择基函数 $ \{ \phi_i(x) \}_{i=1}^{n-1} $ 如下
$$
\phi_i(x) = \left\{
\begin{aligned}
&\frac{x-x_{i-1}}{x_{i}-x_{i-1}}, &x \in [x_{i-1}, x_{i}] \\
&\frac{x_{i+1} - x}{x_{i+1} - x_{i}}, &x \in [x_i, x_{i+1}] \\
&0, & x \notin [x_{i-1}, x_{i+1}]
\end{aligned}

\right.
$$
那么，我们可以计算 $\phi'_i(x)$ 和 $ K_{ij} $ 如下
$$
\phi'_i(x) = \left\{
\begin{aligned}
&\frac{1}{x_i - x_{i-1}} = \frac{1}{h_i}  , & x \in (x_{i-1}, x_i)  \\
&\frac{1}{x_{i} - x_{i+1}} = -\frac{1}{h_{i+1}}, & x \in (x_i, x_{i+1}) \\
\end{aligned}
\right.

\\

K_{ij} = \epsilon\int_0^1 \phi_j'(x) \phi_i'(x) dx + \int_0^1 \phi_j'(x) \phi_i(x) dx
$$

## （面向初学者的超详细）运行指南

