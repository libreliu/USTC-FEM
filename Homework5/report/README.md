# 有限元分析 第五次程序作业

## 方法整理

### 共轭梯度法

对于正定矩阵 $ A $， $ Ax = b $ 可以用梯度方法求解。

考虑 $ f(z) = \frac{1}{2} z^{T} A z - z^T b $ ，则正定时 $ f(x) \le f(z) ,\forall z $。

迭代公式如下

$$
\begin{aligned}
&x_{k+1} = x_{k} + \frac{r_k^T r_k}{d_k^T A d_k} d_k \\
&r_{k+1} = A x_{k+1} -b \\
&d_{k+1} = -r_{k+1} + \frac{r_{k+1}^T A d_k}{d_k^T A d_k} d_k \\

\end{aligned}
$$

其中
$$
\begin{aligned}
&r_0 = Ax_0 - b \\
&d_0 = -r_0
\end{aligned}
$$

对于共轭梯度方法，如果可以直接计算 $ x^T A y $，那么就可以略去显式构造 $ A $ 这一步。

事实上，

$$
(x^T A y) = \sum_i x_i \cdot (\sum_k a_{ik} y_{k}) = \sum_{ij} a_{ij} x_i y_j
$$

所以只要可以以 $ O(1) $ 的时间和空间复杂度计算 $ a_{ij} $ 并求和，就可以以合适的复杂度得到 $ x^T A y $ 的结果。

### 多重网格方法

存在 $ k \ge 1 $ 层网格，其中第 $ m+1 $ 层是第 $ m $ 层采用 Tessellation 得到的。

对于第 $ k $ 层的近似解 $ \hat u_k $（$ A \hat u_k \approx f $），其计算方法如下：
$$
\begin{aligned}
u_0^k &= I_{k-1}^k \hat u_{k-1} \\
u_l^k &= MG(k, u_{l-1}^k, f_k), 1 \le l \le r \\
\hat u^k &= u_r^k

\end{aligned}
$$

其中 $ ans = MG(k, z_0, g) $ 的计算如下：
> $ z_0 $ 为初始近似解，$ g $ 为右端项，$ k $ 为使用的网格的精度
- 若 $ k = 1 $，$ ans \leftarrow A_1 ^{-1} g $
- 若 $ k \ge 2 $
  1. (Presmoothing Step) 进行 $ m_1 $ 次如下操作：

     $$ z_l = z_{l-1} + \frac{1}{\Lambda_k} (g - A_k z_{l-1}) $$

     其中 $ \Lambda_k $ 是满足 $ \Lambda_k \le C h_k^{-2} $ 的 $ A_k $ 的谱半径上界估计。
  
  2. (Error Correction Step) 首先将解 $ z_{m_1} $ 的残差放到 $ V_{k-1} $ 中去。记
    
     $$ \bar g := I_k^{k-1} (g-A_k z_{m_1} ) $$

     然后取 $ q_0 = 0 $ 并迭代 $ p $ 次

     $$ q_i = MG(k-1, q_{i-1}, \bar g) $$

     再把修正项放回来

     $$ z_{m_1 + 1} := z_{m_1} + I_{k-1}^{k} q_p $$

  3. (Postsmoothing Step) 进行 $ m_2 $ 次如下操作：

     $$ z_l = z_{l-1} + \frac{1}{\Lambda_k} (g - A_k z_{l-1}) $$

  4. $ ans \leftarrow z_{m_1 + m_2 + 1} $  

其中有几个问题需要注意：
1. $ \Lambda_k $ 的确定
   
2. 

## 问题描述

$$
\left\{
\begin{aligned}
&-u'' = f \\
&u(0) = u(1) = 0
\end{aligned}
\right.
$$

在分片 $ P^1 $ 多项式元下求解，给出共轭梯度方法和多重网格方法的时间。

