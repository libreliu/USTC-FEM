# 有限元分析 第四次程序作业

刘紫檀 SA21229063

## 问题描述

给定 PDE

$$
\begin{aligned}
&\left\{
\begin{aligned}
&-\Delta u(x, y) = f & \text{in}\ \Omega=[0, 1] \times [0, 1] \\
&u = 0 & \text{on} \ \partial \Omega
\end{aligned}
\right.
\\

&\text{Using}\ V_h = \{ v \in H^1_0(\Omega) \quad v |_{K} \in P^1(K), \forall k \in \mathcal{T}_k \}
\end{aligned}
$$
并且使用非均匀的三角网格（i.e. 不是 criss cross）的协调有限元方法求解。



## 问题求解

$$
\begin{aligned}
a(u, v) &= \int_\Omega (-\Delta u) v  dx = \int_\Omega \nabla u \cdot \nabla v dx - \int_{\partial \Omega} \frac{\partial u}{\partial \vec n} v ds = \int_\Omega \nabla u \cdot \nabla v dx \\
F(v) &= \int_\Omega f v dx

\end{aligned}
$$

$ a(u, v) $ 为对称双线性形式，同时可以证明其强制性，这样知其解唯一。

### 离散化

令 $ v(x, y) = \sum v_i \phi_i(x, y)$ ，$ u(x, y) = \sum u_j \phi_j(x, y)$ ，则得
$$
\begin{aligned}
a(u, v) & = \int_\Omega \nabla (\sum u_j \phi_j(x, y)) \cdot \nabla(\sum v_i \phi_i(x, y)) dxdy \\
&= \int_\Omega  (\sum u_j \nabla\phi_j(x, y)) \cdot (\sum v_i \nabla \phi_i(x, y)) dxdy\\
&= \int_\Omega  (\sum u_j \nabla\phi_j(x, y)) \cdot (\sum v_i \nabla \phi_i(x, y)) dxdy\\
&= \sum \sum u_j v_i  \int_\Omega \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy\\
F(v) &= \int_\Omega f(x, y) (\sum v_i \phi_i(x, y)) dxdy = \sum v_i \int_\Omega f(x, y) \phi_i(x, y) dxdy

\end{aligned}
$$
分别取自然基函数得到如下方程共 $ i $ 个
$$
\sum u_j  \int_\Omega \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy =  \int_\Omega f(x, y) \phi_i(x, y) dxdy
$$
则只要计算 $ K_{ij} = \int_\Omega \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy $ 和 $ F_i = \int_\Omega f(x, y) \phi_i(x, y) dx dy$ 再拼装解方程组即可。

### 二维分片一次多项式

我们取线性 Lagrange 元作为我们的有限元。对于每个三角形 $ K \in \mathcal T_h $ ，其三个顶点分别记为 $ (x_i, y_i) $ （i=1,2,3 )，我们定义其内部的**重心坐标** $ (\xi_1, \xi_2, \xi_3 ) $ 满足如下关系：
$$
\left\{
\begin{aligned}
x &= x_1 \xi_1 + x_2 \xi_2 + x_3 \xi_3 \\
y &= y_1 \xi_1 + y_2 \xi_2 + y_3 \xi_3\\
1 &= \xi_1 + \xi_2 + \xi_3\\
\end{aligned}
\right.
$$

用矩阵逆可以求出给定 $ (x, y) $ 时 $ (\xi_1, \xi_2, \xi_3) $ 的值
$$
\begin{bmatrix}
\xi_1 \\
\xi_2 \\
\xi_3 \\
\end{bmatrix} =
\frac{1}{ \det A}
\left[
\begin{array}{ccc}
 {y_2}-{y_3} & {x_3}-{x_2} & {x_2} {y_3}-{x_3} {y_2} \\
 {y_3}-{y_1} & {x_1}-{x_3} & {x_3} {y_1}-{x_1} {y_3} \\
 {y_1}-{y_2} & {x_2}-{x_1} & {x_1} {y_2}-{x_2} {y_1} \\
\end{array}
\right]
\begin{bmatrix}
x \\
y \\
1 \\
\end{bmatrix} \\
\text{where} \ 
\det A = 
\left |
\begin{array}{ccc}
 x_1 & x_2 & x_3 \\
 y_1 & y_2 & y_3 \\
 1   & 1   & 1 \\
\end{array}
\right|
 =x_1 y_2 - x_1 y_3 - x_2 y_1 + x_2 y_3 + x_3 y_1 - x_3 y_2
$$
同时取 $ \xi_1, \xi_2 $ 为独立变量时容易得到
$$
\begin{bmatrix}
\frac{\partial}{\partial \xi_1} \\ 
\frac{\partial}{\partial \xi_2}
\end{bmatrix}

= 
\begin{bmatrix}
x_1 - x_3 & y_1 - y_3 \\ 
x_2 - x_3 & y_2 - y_3
\end{bmatrix}

\begin{bmatrix}
\frac{\partial}{\partial x} \\ 
\frac{\partial}{\partial y}
\end{bmatrix}
$$
我们观察到 $ (x, y) $ 到 $ (\xi_1, \xi_2) $ 的坐标变换为非退化的，只要三角形 $ K $ 非退化，则我们可以在重心坐标 $ (\xi_1, \xi_2, 1 - \xi_1 - \xi_2 ) $ 上定义线性 Lagrange 元的基函数如下

$$
\phi_i(\xi_1, \xi_2) = \xi_i \quad (\text{where} \ \xi_3 = 1 - \xi_1 - \xi_2)
$$

注意此定义为一局部定义，全局的情形是由各个局部情形拼接而成的，所以严格来说需要把全局和局部的关系在详细写出，此处略。

#### $ K_{ij} $ 计算

变分后的问题需要我们计算 $  \int_\Omega \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dx $ ，而由 $ \Omega = \cup K $ 我们只需要处理在每个 $ K $ 上的积分在求和即可。

对于某个 $ K $ ，我们设其三个顶点坐标 $ v_1, v_2, v_3 $ 分别为 $  (x_1, y_1), (x_2, y_2), (x_3, y_3) $。

我们设对于某个全局编号为 $ G_i $ 的顶点在此三角形中顶点的编号为 $ i $ ，则对于使用三角形内部标号的节点基函数 $ \phi_i(x, y) $ 和 $ \phi_j (x, y) $ ，我们有
$$
\begin{aligned}
(\nabla \phi_i (x, y))_s &= \sum_m\partial_{\xi_m} \phi_i \partial_{x_s} \xi_m \\
&= \partial_{x_s} \xi_i \qquad {(\because \phi_i(\xi_m) = \delta_{im})} \\
&= \frac{1}{\det A}\left\{

\begin{aligned}
& y_2 - y_3 & i=1,s=1 \\
& x_3 - x_2 & i=1,s=2 \\
& y_3 - y_1 & i=2,s=1 \\
& x_1 - x_3 & i=2,s=2 \\
& y_1 - y_2 & i=3,s=1 \\
& x_2 - x_1 & i=3,s=2 \\

\end{aligned}

\right.
\end{aligned}
\\
\text{where} \ |A| = x_1 y_2 - x_1 y_3 - x_2 y_1 + x_2 y_3 + x_3 y_1 - x_3 y_2\ \text{and}\ x_1 := x,\ x_2 := y
\\

\begin{aligned}
\int_K \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy &= \int_K \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy \\
&= \int_K (\nabla\phi_j(x, y))_1 (\nabla \phi_i(x, y))_1 + (\nabla\phi_j(x, y))_2 (\nabla \phi_i(x, y))_2 dxdy \\ 
&= \frac{S_K}{(\det A)^2} \left\{
\begin{aligned}
& (y_2 - y_3)^2 + (x_3 -x_2)^2 & j=1,i=1 \\
& (y_2 - y_3)(y_3 - y_1) + (x_3 - x_2)(x_1 - x_3) & j=1,i=2 \\ 
& (y_2 - y_3)(y_1 - y_2) + (x_3 - x_2)(x_2 - x_1) & j=1,i=3 \\
& (y_3 - y_1)^2 + (x_1 -x_3)^2 & j=2,i=2 \\
& (y_3 - y_1)(y_1 - y_2) + (x_1 - x_3)(x_2 - x_1) & j=2,i=3 \\
& (y_1 - y_2)^2 + (x_2 -x_1)^2 & j=3,i=3 \\
\end{aligned}
\right.
\end{aligned}
\\

\text{where}\ S_k = \frac{1}{2}|\det A| = \frac{1}{2} |x_1 y_2 - x_1 y_3 - x_2 y_1 + x_2 y_3 + x_3 y_1 - x_3 y_2| \qquad \text{(Shoelace Formula)}
$$

#### $ F_i $ 计算

记号同上。
$$
\begin{aligned}
\int_K f(x, y) \phi_i(x, y) dxdy 
&= \int_K f(x(\xi_1, \xi_2), y(\xi_1, \xi_2)) \phi_i(\xi_1, \xi_2) | \frac{\partial(x, y)}{\partial(\xi_1, \xi_2)} | d\xi_1 d\xi_2 \\
&= \int_K f(x(\xi_1, \xi_2), y(\xi_1, \xi_2)) \phi_i(\xi_1, \xi_2) \left| 
\begin{matrix}
x_1 - x_3 & x_2 - x_3 \\
y_1 - y_3 & y_2 - y_3
\end{matrix}

\right| d\xi_1 d\xi_2 \\
&= \left| 
\begin{matrix}
x_1 - x_3 & x_2 - x_3 \\
y_1 - y_3 & y_2 - y_3
\end{matrix}

\right|
 \left\{ \begin{aligned}
&\int_K f(\xi_1, \xi_2) \xi_i d \xi_1 d\xi_2 & i = 1, 2 \\
&\int_K f(\xi_1, \xi_2) (1-\xi_1-\xi_2) d \xi_1 d\xi_2 & i = 3
\end{aligned}
\right.
\end{aligned}
$$

容易观察到积分区域 K 在 $(\xi_1, \xi_2)$ 坐标下是 $ \xi_2 = -\xi_1 + 1 $ 与坐标轴围成的三角形，则在此区域上进行数值积分即可。

### 算法步骤

二维问题思考稍显复杂，于是略记如下以飨读者。

1. 计算总自由度 = 所有顶点 - 在边界上的顶点；对于每一个自由度（i.e. 每个自由顶点），都对应着一个基函数 $ \phi_i(x, y) $

2. 维护一个双射 $ M: \text{vertex}_i \to i_{M} $  ，$i_M \in \{1, ..., N_f\}$ ，其中 $ N_f $ 为总自由度

   之后更新矩阵的时候，下标 $ K_{ij} $  中的 $ i $ 和 $ j $ 分别为  $ M(\text{vertex}_i) $ 和   $ M(\text{vertex}_j) $ 

3. 初始化 $ K_{ij} = 0, \forall i, j $

4. 对每个自由顶点 $ \text{vertex}_i $ 和其对应的基函数 $ \phi_i(x,y) $ 
   1. $ F_i \leftarrow 0$
   2. 确定影响范围，即找到所有与该顶点临接的三角面，赋值到列表 $ F $
   3. 对于每个 $ F $ 中的面 $ A $
      1.  $ F_i \leftarrow F_i + \int_A f(x, y) \phi_i(x, y) dx $
      2. 寻找面 $A$ 中不为 $ \text{vertex}_i $ 的两个点 $ \text{vertex}_j$ 和  $ \text{vertex}_k$ 
         1. 若 $ i < j $，或者 $ j $ 非自由顶点，则跳过计算；否则 $ K_{ij} \leftarrow K_{ij} + \int_A \nabla\phi_j(x, y) \cdot \nabla \phi_i(x, y) dxdy $ 
         2. 若 $ i < k $，或者 $ k $ 非自由顶点，则跳过计算；否则 $ K_{ik} \leftarrow K_{ik} + \int_A \nabla\phi_k(x, y) \cdot \nabla \phi_i(x, y) dxdy $ 
      3. $ K_{ii} \leftarrow K_{ii} + \int_A \nabla\phi_i(x, y) \cdot \nabla \phi_i(x, y) dxdy$ 

5. 利用对称性把没有填好的 $ K_{ij} $ 中的元素填好

6. 解线性方程组 $ \mathrm{KU}=\mathrm{F} $ 得到 $ \mathrm U $ 