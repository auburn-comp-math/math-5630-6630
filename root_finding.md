---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Root Finding

Root finding is a fundamental problem in numerical analysis and has many applications in science and engineering such as solving nonlinear equations, optimization problems, and differential equations. Usually, a closed form of the root is not available, and we need to compute the root numerically. In this chapter, we will discuss some of the most common methods for root finding.

## Bracket Methods

If $f$ is a continuous function, and $f(a)$ and $f(b)$ have opposite signs, then by the Intermediate Value Theorem, there exists a root of $f$ on the interval $[a, b]$. The bracket method is based on this fact and iteratively locates the pair of points $a$ and $b$ such that $f(a)$ and $f(b)$ have opposite signs. The most common bracket methods are the bisection method and the false position method.

### Bisection Method

The simplest bracket method is the **bisection method**. Once $f(a)f(b) < 0$, one can select the midpoint $c = \frac{a + b}{2}$ and check the sign of $f(a) f(c)$.

- If $f(a)f(c) < 0$, then the root is in the interval $[a, c]$.
- If $f(a)f(c) > 0$, then the root is in the interval $[c, b]$.
- If $f(a)f(c) = 0$, then $c$ is the root.

```{margin} Iterative vs Recursive
The bisection method can be implemented either iteratively or recursively. The recursive program usually is more compact, but may suffer from inefficiency (more instructions) and the risk of stack overflow.

For a variety of programming languages, recursive implementation can be made more efficient by using the so called **tail recusion optimization**. This is a compiler feature that allows the recursive program to be executed with the same efficiency as the iterative one. However, neither ``Python`` nor ``MATLAB`` natively supports tail recursion optimization. : )
```

For the first two cases, we can repeat the process with the new interval until certain stop criteria are met. Each iteration reduces the size of the interval by half (gaining one bit each iteration), the total number of iterations required to reduce the interval to a certain size is $\lceil\log_2\left(\frac{b - a}{\epsilon}\right)\rceil$, where $\epsilon$ is the desired tolerance.

Once the function $f$ has a sign change over the bracket $[a, b]$, the bisection method is guaranteed to converge to a root, but it is not very efficient. It is usually used to obtain a rough estimate of the root, which is then taken as an initial guess for a more efficient method.

### False Position Method

The bisection method only uses $\text{sgn}(f(a))$ and $\text{sgn}(f(b))$ instead of the function values. The **false position method** (*Regula falsi* in Latin) improves the bisection method by taking the function values into account. Instead of selecting the midpoint $c = \frac{a + b}{2}$, the false position method selects the point $c\in[a, b]$ that lies on the line connecting $(a, f(a))$ and $(b, f(b))$, that is

$$c = \frac{f(b)}{f(b) - f(a)} a + \frac{-f(a)}{f(b) - f(a)} b = \frac{f(b) a -  f(a) b}{f(b) - f(a)}.$$

The false position method is also guaranteed to converge to a root if $f(a)f(b) < 0$ and the implementation is quite similar to the bisection method. It usually converges faster than the bisection method, but sometimes exceptions occur.

````{note}
When $f''$ keeps the same sign over $[a, b]$, it is not hard to show that only one side of the bracket is updating. The bracket size will never decrease to zero, which is different from the bisection method. In the following, we use an example to illustrate this. The function $f(x) = x^2 - 1$ over the initial bracket $[0, 2]$. The left endpoint is updating to the root while the right endpoint is fixed at $2$.

```{image} images/doc/root_finding_img_0.png
:alt: false position method
:class: bg-primary 
:width: 600px
:align: center
```
````

It is actually easy to improve the false position method by forcing more weight towards the other endpoint. This is called the **Illinois method**. Once the same side is to update in two consecutive iterations, the Illinois method will adjust $c$ using a slightly different formula.

$$
{c} = \frac{\lambda_b f(b) a - \lambda_a f(a) b}{\lambda_b f(b) - \lambda_a f(a)}.
$$

where the weights $\lambda_a$ and $\lambda_b$ are the weights. The weights are initially set to $1$. If the same side is to update, the weight of the other side will be halved. If the side changes, the weight on the new bracket point will be reset to $1$. See {prf:ref}``AL-ILLINOIS``.

```{margin}
The choice of decay factor $\frac{1}{2}$ is optimal if it has to be a constant. The factor can be replaced with other variable values. Two usual replacements are **Pegasus** method and **Anderson-Bjorck** method.

The essential changes in the algorithm are:

- **Pegasus**: replace $\lambda_b \gets \lambda_b/2$ with $\lambda_b\gets \lambda_b\frac{f(a)}{f(a) + f(c)}$ and replace $\lambda_a \gets \lambda_a/2$ with $\lambda_a\gets\lambda_a\frac{f(b)}{f(b) + f(c)}$.
- **Anderson-Bjorck**: replace $\lambda_b \gets \lambda_b/2$ with $\lambda_b\gets \lambda_b m_b$ and replace $\lambda_a \gets \lambda_a/2$ with $\lambda_a\gets\lambda_a m_a$, where 

$$m_a = \begin{cases}1 - \frac{f(c)}{f(b)} &\text{if positive}\\ \frac{1}{2} &\text{otherwise}\end{cases}$$

and 

$$m_b = \begin{cases}1 - \frac{f(c)}{f(a)} &\text{if positive}\\ \frac{1}{2} &\text{otherwise}\end{cases}$$
```

```{prf:algorithm} Illinois Method
:label: AL-ILLINOIS

**Inputs** $f$, $a, b$, $\epsilon$ 

**Outputs** $c$ 

1. $\lambda_a \gets 1$, $\lambda_b \gets 1$ and $c \gets \frac{\lambda_b f(b) a - \lambda_a f(a) b}{\lambda_b f(b) - \lambda_a f(a)}$.  $s \gets 0$. //initialization
2. While True do:
    1. If $f(a)f(c) < 0$, then 
        - If $s \le 0$, $\lambda_a \gets 
    \lambda_a/2$; Else $\lambda_b \gets 1$.
        -  $b \gets c$, $s\gets -1$.
    2. Else If $f(a)f(c) > 0$, then 
        - If $s \ge 0$, $\lambda_b \gets \lambda_b/2$; Else $\lambda_a\gets 1$. 
        - $a \gets c$, $s\gets 1$.

    3. $c \gets \frac{\lambda_b f(b) a - \lambda_a f(a) b}{\lambda_b f(b) - \lambda_a f(a)}$; //compute $c$
    4. If $ |f(c)| < \epsilon$, then return $c$; //check stopping criteria

```

````{note}
We use the previous example to illustrate the difference between the false position method and the Illinois method. At the second iteration, the Illinois method finds the updating is still on left side, so it modifies right endpoint $f(b)$ into $\frac{1}{2} f(b)$ to compute the new $c$, which makes the selected point ${c}$ closer to the right endpoint than the false position method.

```{image} images/doc/root_finding_img_1.png
:alt: Illinois method
:class: bg-primary 
:width: 600px
:align: center
```
````

```{margin} Stopping Criteria
There are several types of stopping criteria to terminate the iteration. Common ones include:

- $|f(c)| < \texttt{ftol}$, tolerance on the function value.
- $|c_n - c_{n-1}| < \texttt{atol}$, absolute tolerance.
- $|c_n - c_{n-1}| < \texttt{rtol} |c_n|$, relative tolerance.
```

``````{prf:example}
Let us try the aforementioned methods to find the root of $f(x) = x^3 - 2x^2 - 4$ on the interval $[1, 3]$. The root $x^{\ast}$ can be computed analytically through cubic root formula, which is roughly ``2.5943130163548496``.

Using the previous methods, we obtain the sequence of selection $c_n$ and the error $|c_n - x^{\ast}|$. The tolerance is set to $|f(c)| < 10^{-6}$. The results are shown in the following table.

`````{tab-set}
````{tab-item} Bisection Method
```{code-block}
iter  1 | 2.0000000000000000000000000 | 5.94e-01 
iter  2 | 2.5000000000000000000000000 | 9.43e-02 
iter  3 | 2.7500000000000000000000000 | 1.56e-01 
iter  4 | 2.6250000000000000000000000 | 3.07e-02 
iter  5 | 2.5625000000000000000000000 | 3.18e-02 
iter  6 | 2.5937500000000000000000000 | 5.63e-04 
iter  7 | 2.6093750000000000000000000 | 1.51e-02 
iter  8 | 2.6015625000000000000000000 | 7.25e-03 
iter  9 | 2.5976562500000000000000000 | 3.34e-03 
iter 10 | 2.5957031250000000000000000 | 1.39e-03 
iter 11 | 2.5947265625000000000000000 | 4.14e-04 
iter 12 | 2.5942382812500000000000000 | 7.47e-05 
iter 13 | 2.5944824218750000000000000 | 1.69e-04 
iter 14 | 2.5943603515625000000000000 | 4.73e-05 
iter 15 | 2.5942993164062500000000000 | 1.37e-05 
iter 16 | 2.5943298339843750000000000 | 1.68e-05 
iter 17 | 2.5943145751953125000000000 | 1.56e-06 
iter 18 | 2.5943069458007812500000000 | 6.07e-06 
iter 19 | 2.5943107604980468750000000 | 2.26e-06 
iter 20 | 2.5943126678466796875000000 | 3.49e-07 
iter 21 | 2.5943136215209960937500000 | 6.05e-07 
iter 22 | 2.5943131446838378906250000 | 1.28e-07 
iter 23 | 2.5943129062652587890625000 | 1.10e-07 
iter 24 | 2.5943130254745483398437500 | 9.12e-09 
```
````

````{tab-item} False Position Method
```{code-block}
iter  1 | 2.0000000000000000000000000 | 5.94e-01 
iter  2 | 2.4444444444444446418174266 | 1.50e-01 
iter  3 | 2.5621621621621617492792211 | 3.22e-02 
iter  4 | 2.5876913365185605364615640 | 6.62e-03 
iter  5 | 2.5929610854818996301673906 | 1.35e-03 
iter  6 | 2.5940374914642010395482430 | 2.76e-04 
iter  7 | 2.5942568846837747997824408 | 5.61e-05 
iter  8 | 2.5943015817106331866170876 | 1.14e-05 
iter  9 | 2.5943106870264029950590157 | 2.33e-06 
iter 10 | 2.5943125418534931370118102 | 4.75e-07 
iter 11 | 2.5943129196954899384763849 | 9.67e-08 
```
````

````{tab-item} Illinois Method
```{code-block}
iter  1 | 2.0000000000000000000000000 | 5.94e-01 
iter  2 | 2.6153846153846154187760931 | 2.11e-02 
iter  3 | 2.5847750865051901669744439 | 9.54e-03 
iter  4 | 2.5941951587569969106539247 | 1.18e-04 
iter  5 | 2.5944267005726100450146987 | 1.14e-04 
iter  6 | 2.5943130084597889606357057 | 7.90e-09 
```
````

`````
``````

```{prf:remark}
:label: rmk:bracket_methods
The bracket methods need to first locate an interval $[a, b]$ such that $f(a)f(b) < 0$. A common approach is to sample a few equidistant points in a large interval and then use the sign of the function values to identify the bracket. This is a simple and robust approach, but it may require a large number of function evaluations.
```

### Order of Convergence

The order of convergence quantifies how fast the sequence approximates the limiting value.

```{prf:definition}
:label: dfn-order_of_convergence
The order of convergence of a sequence $\{x_n\}$ is $p > 0$ if 

$$
\lim_{n\to\infty}\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p} = \rho,
$$

The constant $\rho$ is the rate of convergence. If $p = 1$, the sequence is said to have linear convergence. If $p = 2$, the sequence is said to have quadratic convergence. 

If the limit does not exist while the upper bound exists for sufficiently large $n$:

$$
\limsup_{n\to\infty}\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p} \le \rho,
$$

then the order of convergence is **at least** $p$ and the rate of convergence is **at most** $\rho$.
```

```{margin} More About Convergence
```

```{margin} More About Convergence
In practice, the limit or even upper bound may not exist for $\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p}$, it is possible to consider the convergence rate in weaker sense. For instance, for sufficiently large $n$ that the inequality 

$$
\lim_{k\to\infty}\sqrt[k]{\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|^{p^k}}}= \rho
$$

holds for certain $p > 0$ and $\rho > 0$, then the sequence has a mean order of convergence is $p$ and a mean convergence rate $\rho$. 

On averaged, each iteration contributes an order of $p$ and a rate of $\rho$.
```

````{prf:theorem}
:label: thm-bisection_convergence

The bisection method has a linear convergence rate. 

````

````{prf:proof}
:label: prf-bisection_convergence

Without loss of generality, we may assume the initial bracket is on $[0, 1]$. Let the root $x^{\ast} = 0.b_1 b_2\cdots$ be the binary representation, then the sequence of bisection method can be written as 

$$
x_n = 0.b_1 b_2\cdots b_{n-1} 1,
$$

where $b_i$ is the $i$-th bit of the binary representation. The error at the $n$-th iteration is

$$|x_n - x^{\ast}| = |2^{-n} - \sum_{j\ge n} 2^{-j} b_j|=\begin{cases}\sum_{j>n} 2^{-j}b_j & \text{if } b_n=1\\\sum_{j>n} 2^{-j}(1-b_j) &\text{if }b_n=0\end{cases}$$

For each $n$ such that $b_n = 1$, we can find $s\in\N$ such that $b_{s+n} = 1$, otherwise we arrive at the exact solution. 

$$
\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|} = \frac{\sum_{j>n+k} 2^{-j} b_j }{\sum_{j>n} 2^{-j}b_j} \le \frac{2^{-(n+k)}}{2^{-(n+s)}} =2^{s - k} $$

Therefore, the geometric mean of the convergence rate is bounded by $\frac{1}{2}$.

$$
\rho =\lim_{k\to\infty} \sqrt[k]{\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|}} \le \lim_{k\to\infty}2^{s/k}\frac{1}{2} =\frac{1}{2}.
$$
````

## Iterative Methods

### Newton-Raphson Method

### Secant Method

## Applications in Optimization

## Exercises

### Theoretical Part

```{admonition} Problem 1
Derive the convergence rate for the Illinois method.
```

### Computational Part
