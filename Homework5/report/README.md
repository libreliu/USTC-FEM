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
   
   > TODO: 严格证明
   
   此处使用 CG 方法迭代一轮。
   
2. $ m_1 , m_2 $ 的选择

   此处均使用 1。

3. Coarse to fine 和 fine to coarse 算子的选择

   ```python
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
   ```

   coarseToFine 采用自然嵌入，fineToCoarse 采用中心点权重 1，周围权重 1/2 的方式进行。

   > 周围权重是 0 的话误差会到 1e-4 下不去，为什么？

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

## 实验结果

使用 

$$
\begin{aligned}
&u(x) = (x-1) \sin x + 2 \cos x + (2 - 2 \cos 1) x - 2 \\
&f(x) = (x-1) \sin x
\end{aligned}
$$

进行测试，结果如下

> MG 方法的 r = 1, p = 1，详情参考代码

### MG 方法

| kMax | n | $ L^1 $ error | $ L^\infty $ error | $ L^1 $ order | $ L^\infty $ order | MG 求解时间 (sec) | MG 求解时间阶估计 |
| - | -- | ----- | ------------------ | ----- | ----- | ----- | ----- |
| 2 | 4 | 1.017e-03 | 2.041e-03 | N/A | N/A | 9.979e-04 | 0.000 |
| 3 | 8 | 2.286e-04 | 4.867e-04 | 2.154 | 2.069 | 3.529e-03 | 1.822 |
| 4 | 16 | 7.477e-05 | 1.626e-04 | 1.612 | 1.581 | 9.054e-03 | 1.359 |
| 5 | 32 | 1.597e-05 | 3.384e-05 | 2.227 | 2.265 | 2.174e-02 | 1.264 |
| 6 | 64 | 3.786e-06 | 8.616e-06 | 2.076 | 1.973 | 4.804e-02 | 1.144 |
| 7 | 128 | 7.987e-07 | 1.848e-06 | 2.245 | 2.221 | 1.019e-01 | 1.085 |
| 8 | 256 | 1.586e-07 | 3.992e-07 | 2.333 | 2.211 | 2.136e-01 | 1.067 |
| 9 | 512 | 4.639e-08 | 1.102e-07 | 1.773 | 1.857 | 4.303e-01 | 1.011 |
| 10 | 1024 | 1.070e-08 | 2.660e-08 | 2.116 | 2.050 | 8.042e-01 | 0.902 |

### CG 方法

| n | $ L^1 $ error | $ L^\infty $ error | $ L^1 $ order | $ L^\infty $ order | CG 求解时间 (sec) | CG 求解时间阶估计 |
| -- | ----- | ------------------ | ----- | ----- | ----- | ----- |
| 4 | 8.352e-04 | 1.765e-03 | 0.000 | 0.000 | 9.782e-04 | 0.000 |
| 8 | 2.071e-04 | 4.640e-04 | 2.012 | 1.927 | 4.031e-03 | 2.043 |
| 16 | 5.166e-05 | 1.171e-04 | 2.003 | 1.987 | 1.667e-02 | 2.048 |
| 32 | 1.291e-05 | 2.930e-05 | 2.001 | 1.999 | 6.778e-02 | 2.023 |
| 64 | 3.227e-06 | 7.327e-06 | 2.000 | 2.000 | 2.777e-01 | 2.034 |
| 128 | 8.067e-07 | 1.832e-06 | 2.000 | 2.000 | 1.095e+00 | 1.980 |
| 256 | 2.017e-07 | 4.580e-07 | 2.000 | 2.000 | 4.325e+00 | 1.982 |
| 512 | 5.042e-08 | 1.145e-07 | 2.000 | 2.000 | 1.755e+01 | 2.021 |
| 1024 | 1.261e-08 | 2.863e-08 | 2.000 | 2.000 | 6.972e+01 | 1.990 |

### 总结

可以看到，共轭梯度为 $ O(n^2) $ 的时间复杂度，MG 方法为 $ O(n) $ 的时间复杂度；收敛阶上两者均近似为二阶方法。