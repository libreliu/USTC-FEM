# 《有限元分析》 总结与回顾

## 章节一览

- Ch0. 基本概念
  - Ritz-Galerkin 近似
  - 与有限差分法的比较
    > 收敛性证明容易; 编写容易，只需要替换无穷维到有限维;
  - 泛函和实分析超快速回顾 
- Ch1. Sobolev 空间
  - 弱导数
  - Sobolev 范数和 Sobolev 空间
  - Sobolev 嵌入定理
  - 迹定理
  - 负范数和对偶性
- Ch2. 椭圆边值问题的变分形式
  - 对称变分问题
  - 非对称变分问题
  - Lax-Milgram 定理
  - 高维情形
- Ch3. 有限元空间构造
  - 有限元的要求
  - 三角元
  - 矩形元
- Ch4. Sobolev 空间上的多项式近似理论
- Ch5. n 维变分问题
- Ch6. 多重网格方法

## 编写计划

ddl 好多..争取每周一章？要不然都堆到期末肯定要 gg，我们还有两个月就期末了

## 第零章

## 第一章


### Minkowski's Inequality
对于 $ 1 \le p \le \infty $ 和 $ f, g \in L^p(\Omega) $ 有  

$$
|| f+g || _{L^p(\Omega)} \le || f || _{L^p(\Omega)} + || g || _{L_p(\Omega)}
$$

证明：

### Holder's Inequality
对于 $ 1 \le p,q \le \infty $ 使得 $ 1 = 1/p + 1/q $，如果 $ f \in L^p(\Omega) $ 且 $ g \in L^q(\Omega) $，则有

$$
|| fg || _L^1(\Omega) \le || f || _{L^p(\Omega)} || g || _{L^q(\Omega)}
$$

证明：

### Schwarz's Inequality

就是 $ p = q = 2 $ 的 Holder 不等式

$$
\int_{\Omega} | f(x) g(x) | dx \le || f || _{L^2(\Omega)} || g || _{L^2(\Omega)}
$$

> 其实是内积空间的常用结论 $ |(x, y)|^2 \le (x, x) (y, y) $，如果利用 $ L^2 $ 范数定义 $ L^2 $ 内积的话

证明：

### $ L^2 $ 范数诱导的 $ L^2 $ 内积

首先，我们知道

### 内积诱导的范数

## 第二章



## 第三章

## 第四章

## 第五章

## 第六章