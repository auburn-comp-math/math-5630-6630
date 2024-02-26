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
When $f''$ keeps the same sign over $[a, b]$, it is not hard to show that only one side of the bracket is updating. The bracket size will never decrease to zero, which is different from the bisection method. In the following, we use an example to illustrate this, see {numref}`root-finding-false-img`. The function $f(x) = x^2 - 1$ over the initial bracket $[0, 2]$. The left endpoint is updating to the root while the right endpoint is fixed at $2$.

```{figure} images/doc/root_finding_img_0.png
---
name: root-finding-false-img
scale: 80%
align: center
---
False position method
```

````

It is actually easy to improve the false position method by forcing more weight towards the other endpoint. This is called the **Illinois method**. Once the same side has updated in two consecutive iterations, the Illinois method will adjust $c$ using a slightly different formula.

$$
{c} = \frac{\lambda f(b) a -  f(a) b}{\lambda f(b) - f(a)}
$$

or

$$
{c} = \frac{ f(b) a - \lambda f(a) b}{f(b) - \lambda f(a)},
$$

where the weight $\lambda$ controls the position of $c$. The weight is initially set to $1$ which corresponds to the standard false position method. If the same side is about to update twice, the weight of the other side will be halved. If not, the weight on the new bracket point will be reset to $1$. See {prf:ref}``AL-ILLINOIS``.

```{margin}
The choice of decay factor $\frac{1}{2}$ is optimal if it has to be a constant (explain later). The factor can be replaced with other variable values. A usual replacement is the **Pegasus** method.

Essentially, the **Pegasus** method replaces $\lambda = 1/2$ with $\lambda=\frac{f_1}{f_1 + f_2}$.

<!-- - **Anderson-Bjorck**: it is slightly different, instead of detecting two consecutive false position iterations on the same side, it will try to prevent two consecutive false position iterations on the same side by interchanging endpoints. Otherwise, it replace $\lambda_b \gets \lambda_b/2$ with $\lambda_b\gets \lambda_b m_b$ and replace $\lambda_a \gets \lambda_a/2$ with $\lambda_a\gets\lambda_a m_a$, where 

    $$m_a = \begin{cases}1 - \frac{f(c)}{f(b)} &\text{if positive}\\ \frac{1}{2} &\text{otherwise}\end{cases}$$

    $$m_b = \begin{cases}1 - \frac{f(c)}{f(a)} &\text{if positive}\\ \frac{1}{2} &\text{otherwise}\end{cases}$$
-->
```

```{prf:algorithm} Illinois Method
:label: AL-ILLINOIS

**Inputs** $f$, $a, b$, $\epsilon$ 

**Outputs** $c$ 

1. $x_0\gets a$, $x_1\gets b$, $f_0\gets f(x_0)$, $f_1\gets f(x_1)$.  //initialization.
2. While True do:
    1. $x_2 \gets \frac{f_1 x_0 - f_0 x_1}{f_1 - f_0}$, $f_2 \gets f(x_2)$. //standard false position step, $x_1$ and $x_2$ are two latest iterations.
    2. If $ |f_2| < \epsilon$, then return $x_2$; //check stopping criteria
    3. While $f_1 f_2 > 0$, then // adjust until the sign changes
        - $(x_0, f_0)\gets (x_0, \lambda f_0)$, where $\lambda = \frac{1}{2}$
        - $(x_1, f_1)\gets (x_2, f_2)$ and  $x_2 \gets \frac{f_1 x_0 - f_0 x_1}{f_1 - f_0}$, $f_2 \gets f(x_2)$
    4. If $f_1 f_2 < 0$, then // perform a false position step
        - $(x_0, f_0)\gets (x_1, f_1)$ and $(x_1, f_1)\gets (x_2, f_2)$. 

```

````{note}
We use the previous example to illustrate the difference between the false position method and the Illinois method. At the second iteration, the Illinois method finds the updating is still on left side, so it modifies right endpoint $f(b)$ into $\frac{1}{2} f(b)$ to compute the new $c$, which makes the selected point ${c}$ closer to the right endpoint than the false position method, see {numref}`root-finding-illinois-img`.

```{figure} images/doc/root_finding_img_1.png
---
name: root-finding-illinois-img
scale: 80%
align: center
---
Illinois method
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

````{tab-item} More Digits
Using ``bigfloat`` package, it is possible to extract more digits to observe the dynamics of the convergence. For instance, by setting the precision to 900 bits, the Illinois method's convergence is shown in the following table.
  
```{code-block}
iter  1 | 2.0000000000000000000000000 | 5.94e-01 
iter  2 | 2.6153846153846154187760931 | 2.11e-02 
iter  3 | 2.5847750865051901669744439 | 9.54e-03 
iter  4 | 2.5941951587569973547431346 | 1.18e-04 
iter  5 | 2.5944267005726091568362790 | 1.14e-04 
iter  6 | 2.5943130084597889606357057 | 7.90e-09 
iter  7 | 2.5943130163543197674869134 | 5.29e-13 
iter  8 | 2.5943130163553775879847763 | 5.29e-13 
iter  9 | 2.5943130163548486777358448 | 1.65e-25 
iter 10 | 2.5943130163548486777358448 | 5.13e-38 
iter 11 | 2.5943130163548486777358448 | 5.13e-38 
iter 12 | 2.5943130163548486777358448 | 1.55e-75 
iter 13 | 2.5943130163548486777358448 | 4.70e-113 
iter 14 | 2.5943130163548486777358448 | 4.70e-113 
iter 15 | 2.5943130163548486777358448 | 1.30e-225 
```
````

`````
``````

```{prf:remark}
:label: rmk:bracket-methods
The bracket methods need to first locate an interval $[a, b]$ such that $f(a)f(b) < 0$. A common approach is to sample a few equidistant points in a large interval and then use the sign of the function values to identify the bracket. This is a simple and robust approach, but it may require a large number of function evaluations.
```

### Order of Convergence

The order of convergence quantifies how fast the sequence approximates the limiting value.

```{margin} More About Convergence
In practice, the limit or even upper bound may not exist for $\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p}$, it is possible to consider the convergence rate in weaker sense. For instance, for sufficiently large $n$ that the inequality 

$$
\lim_{k\to\infty}\sqrt[k]{\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|^{p^k}}}= \rho
$$

holds for certain $p > 0$ and $\rho > 0$, then the sequence has a mean order of convergence is $p$ and a mean convergence rate $\rho$. 

On averaged, each iteration contributes an order of $p$ and a rate of $\rho$.
```

````{prf:definition}
:label: dfn-order-of-convergence
The order of convergence of a sequence $\{x_n\}$ is $p > 0$ if 

$$
\lim_{n\to\infty}\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p} = \rho,
$$

The constant $\rho$ is the rate of convergence. If $p = 1$, the sequence is said to have linear convergence. If $p = 2$, the sequence is said to have quadratic convergence. 

If the limit does not exist while the upper bound exists for sufficiently large $n$:

$$
\frac{|x_{n+1} - x^{\ast}|}{|x_n - x^{\ast}|^p} \le \rho,
$$

then the order of convergence is **at least** $p$ and the rate of convergence is **at most** $\rho$.
````

````{prf:theorem}
:label: thm-bisection-convergence

The bisection method has a linear convergence rate. 

````

````{prf:proof}
Without loss of generality, we may assume the initial bracket is on $[0, 1]$. Let the root $x^{\ast} = 0.b_1 b_2\cdots$ be the binary representation, then the sequence of bisection method can be written as 

$$
x_n = 0.b_1 b_2\cdots b_{n-1} 1,
$$

where $b_i$ is the $i$-th bit of the binary representation. The error at the $n$-th iteration is

$$|x_n - x^{\ast}| = |2^{-n} - \sum_{j\ge n} 2^{-j} b_j|=\begin{cases}\sum_{j>n} 2^{-j}b_j & \text{if } b_n=1\\\sum_{j>n} 2^{-j}(1-b_j) &\text{if }b_n=0\end{cases}$$

For each $n$ such that $b_n = 1$ (similar argument holds for $b_n=0$), we can find $s\in\N$ such that $b_{s+n} = 1$, otherwise we arrive at the exact solution. 

$$
\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|} = \frac{\sum_{j>n+k} 2^{-j} b_j }{\sum_{j>n} 2^{-j}b_j} \le \frac{2^{-(n+k)}}{2^{-(n+s)}} =2^{s - k} $$

Therefore, the geometric mean of the convergence rate is bounded by $\frac{1}{2}$.

$$
\rho =\lim_{k\to\infty} \sqrt[k]{\frac{|x_{n+k} - x^{\ast}|}{|x_n - x^{\ast}|}} \le \lim_{k\to\infty}2^{s/k}\frac{1}{2} =\frac{1}{2}.
$$
````

A common technique to study the order of convergence is to use the Taylor expansion. Let us use the Illinois method as an example.

````{prf:theorem}
:label: thm-illinois-convergence
The order of convergence of the Illinois method is $\sqrt[3]{3}$.
````

```{margin} Big $O$ Notation

The big $O$ notation is used to describe the upper bound of a function. For instance, $f(x) = O(g(x))$ means that there exists a constant $C$ such that $|f(x)| \le C|g(x)|$ for sufficiently large $x$. 

Another notation is the big $\Theta$ notation, which is used to describe the upper and lower bounds of a function. For instance, $f(x) = \Theta(g(x))$ means that there exists constants $C_1$ and $C_2$ such that $C_1|g(x)| \le |f(x)| \le C_2|g(x)|$ for sufficiently large $x$.
```

```{margin} Selection of Decay Factor
In theory, any decay factor smaller than one would produce the same order of convergence. The decay factor $\frac{1}{2}$ is optimal in the sense that any smaller factor will make the adjustment iteration move away from the root. While having a larger factor will not be efficient to reduce the bracket size.

```

```{prf:proof}
Let $f(x)\in C^2[a, b]$ and $x^{\ast}$ be the root. For simplicity, we assume the root is simple. Let $\ell(x)$ be the linear interpolation of $f$ at the bracket $[a, b]$, the slope is $f'(\zeta)=\frac{f(b) - f(a)}{b -a}$, then 

$$
f(x) - \ell(x) = \frac{f''(\xi)}{2}(x - a)(x - b),
$$
therefore, the root satisfies $-\ell(x^{\ast}) = \frac{f''(\xi)}{2}(x^{\ast} - a)(x^{\ast} - b)$. Let $c$ be the next iteration, then 

$$
-\ell(x^{\ast}) + \ell(c) = -(x^{\ast} - c) f'(\zeta) =  \frac{f''(\xi)}{2}(x^{\ast} - a)(x^{\ast} - b),
$$

which imples $(x^{\ast} - c) = -\frac{f''(\xi)}{2 f'(\zeta)} (x^{\ast} - a)(x^{\ast} - b)$. When $a$ and $b$ are already near the root, the right-hand side has a fixed sign, which implies that using the unadjusted false position method, the new point $c$ will always fall into a fixed side of the root.

Now, we consider the adjustment of the Illinois method. Without loss of generality, we assume that $c$ falls into the left side of the root. Then the next iteration will be using the bracket $[c, b]$ with adjusted weight for $b$. The new point $c'$ satisfies 

$$
\begin{aligned}
c' - x^{\ast} &= \frac{\frac{1}{2} f(b) (c - x^{\ast}) - f(c) (b - x^{\ast})}{\frac{1}{2} f(b) - f(c)}  = (c - x^{\ast}) - \frac{f(c)(b-x^{\ast})}{\frac{1}{2}f(b) - f(c)}
+ \frac{f(c)(c - x^{\ast})}{\frac{1}{2}f(b) - f(c)}\\
&= (c - x^{\ast}) - \frac{2f(c)}{f(b)}(b - x^{\ast}) \frac{1}{1 - \frac{2f(c)}{f(b)}} + \frac{2f(c)}{f(b)}(c - x^{\ast}) \frac{1}{1 - \frac{2f(c)}{f(b)}}
\end{aligned}
$$

Using Taylor expansion at $x^{\ast}$, since $f'$ is bounded from below near the root, 

 $$
 \begin{aligned}
 \frac{f(c)}{f(b)} &= \frac{(c - x^{\ast}) + O(|c - x^{\ast}|^2)}{(b - x^{\ast}) + O(|b - x^{\ast}|^2)} = \frac{(c - x^{\ast})}{(b - x^{\ast})}\left(1 + \frac{O(|b-x^{\ast}|) + O(|c - x^{\ast}|)}{1 + O(|b-x^{\ast}|)}\right) \\&= \frac{(c - x^{\ast})}{(b - x^{\ast})}\left(1 + O(|b-x^{\ast}|)\right) = \frac{f''(\xi)}{2 f'(\zeta)} (a - x^{\ast}) \left(1 + O(|b-x^{\ast}|)\right)  .
 \end{aligned}
 $$ 

Therefore, we obtain

$$
\begin{aligned}
c' - x^{\ast} &= (c- x^{\ast}) + 2(x^{\ast} - c)\left(1 + O(|a - x^{\ast}|)\right) + O(|a - x^{\ast}||c - x^{\ast}|) \\
&= (x^{\ast} - c)(1 + O(|a - x^{\ast}|)),
\end{aligned}
$$
which moves the new point to the other side of the root by almost reflection. The next bracket becomes $[c, c']$ such that $|c - x^{\ast}|\approx |c' - x^{\ast}|$. Therefore, if denote the last two points as $c = x_{n}$ and $c'=x_{n+1}$, then $|x_{n} - x^{\ast}|\approx |x_{n+1} - x^{\ast}|$, currently $c'$ is on the right side of the root. The next three iterations satisfy

-  $|x_{n+2} - x^{\ast}| = \Theta(|x_n - x^{\ast}| |x_{n+1} - x^{\ast}|) = \Theta(|x_{n+1} - x^{\ast}|^2) $, the iteration $x_{n+2}$ is on left side. 
-  $|x_{n+3} - x^{\ast}| = \Theta(|x_{n+2} - x^{\ast}||x_{n+1} - x^{\ast}|) = \Theta(|x_{n+1} - x^{\ast}|^3)$, the iteration $x_{n+3}$ is on left side, now two consecutive iterations on left side
- $|x_{n+4} - x^{\ast}| = \Theta(|x_{n+3} - x^{\ast}|) = \Theta(|x_{n+1} - x^{\ast}|^3)$, adjusted by Illinois method, the iteration $x_{n+4}$ is on right side, which completes a cycle.

By the previous definition of order of convergence, the Illinois method has an order of convergence $\sqrt[3]{3}$.
```

```{prf:remark}
:label: rmk-illinois-convergence
A quick way of proof uses the asymptotic analysis. 

If we denote $\epsilon_{n} =x_n - x^{\ast}$, then if $\epsilon_{i} > 0$ and $\epsilon_{i-1} < 0$, then $\epsilon_{i+1} \simeq C \epsilon_{i}\epsilon_{i-1}$, where $C\approx \frac{f''(x^{\ast})}{2f'(x^{\ast})}$, here we may assume $C > 0$, then $\epsilon_{i+1} < 0$. 

The next iteration applies the same rule, which is $\epsilon_{i+2} \simeq C \epsilon_{i+1}\epsilon_{i} = C^2 \epsilon_i^2 \epsilon_{i-1} < 0$. 

Now, we need an adjustment step, which gives $\epsilon_{i+3} = -C^2 \epsilon_i^2 \epsilon_{i-1} > 0$, which completes a cycle. Every three iterations make a cycle, to derive the order of convergence, we continue to derive $\epsilon_{i+6}$, which equals to $-C^8\epsilon_{i}^6\epsilon_{i-1}^3 = C^2\epsilon_{i+3}^3$. 

Therefore, the order of convergence is $\sqrt[3]{3}$.
```

By taking the decay factor as a constant, the order of convergence is at most $\sqrt[3]{3}$. However, if the decay factor can be chosen as a variable one, the order of convergence can be improved. Next, we apply the same technique to analyze the Pegasus method, which is a variant of the Illinois method. It is interesting that the Pegasus method was initially discovered in a subroutine for the **Ferranti Pegasus** computer, but no author information is included. The convergence analysis was given by {cite}`dowell1972pegasus` after its discovery.

```{prf:theorem} Pegasus Method
:label: thm-pegasus-convergence
The Pegasus method is a variant of Illinois method, whose updating scheme for the decay factors is replaced by 

$$\lambda=\frac{f_1}{f_1+f_2}$$ 

in {prf:ref}`AL-ILLINOIS`. Then the order of convergence is $\sqrt[4]{\frac{7+\sqrt{57}}{2}}$.
```

```{prf:proof}
Similar to the previous {prf:ref}`rmk-illinois-convergence`, we assume $\epsilon_{i-1} < 0$ and $\epsilon_{i} > 0$ such that $|\epsilon_i|\ll |\epsilon_{i-1}|$, the constant $C = \frac{f''}{2f'}|_{x=x^{\ast}} > 0$.
Similar to the previous analysis, we can derive more terms in the asymptotic form

$$
\epsilon_{i+1}\simeq C \epsilon_{i}\epsilon_{i-1} + D \epsilon_{i} \epsilon_{i-1}(\epsilon_{i} + \epsilon_{i-1})< 0
$$

where $D = -C^2 + \frac{f'''}{6f'}|_{x=x^{\ast}}$. Then due to a different sign for $\epsilon_{i}$ and $\epsilon_{i+1}$, we can derive 

$$\epsilon_{i+2}\simeq C \epsilon_{i+1}\epsilon_{i} + D \epsilon_{i+1} \epsilon_{i}(\epsilon_{i+1} + \epsilon_{i}) < 0.$$

Now use the adjustment step, we obtain (needs some calculation)

$$
\begin{aligned}
\epsilon_{i+3} &= \frac{\epsilon_{i+2}\lambda f(x_i) - \epsilon_i f(x_{i+2})}{\lambda f(x_{i}) - f(x_{i+2})},\quad \lambda = \frac{f(x_{i+1})}{f(x_{i+2}) + f(x_{i+1})}\\
&\approx C^2 \epsilon_{i}^2\epsilon_{i+2} - D \epsilon_{i} \epsilon_{i+1} \epsilon_{i+2}, 
\end{aligned}
$$
it implies $\epsilon_{i+3} < 0$ as well. Therefore, another adjustment step is needed and (after some more calculations, see {prf:ref}`rmk-pegasus-asymptotic`)

$$
\begin{aligned}
\epsilon_{i+4} &=  \frac{\epsilon_{i+3}\lambda f(x_i) - \epsilon_i f(x_{i+3})}{\lambda f(x_{i}) - f(x_{i+3})},\quad \lambda = \frac{f(x_{i+2})}{f(x_{i+3}) + f(x_{i+2})}\frac{f(x_{i+1})}{f(x_{i+2}) + f(x_{i+1})}\\
&\approx C^5 \epsilon_{i}^4\epsilon_{i+1}^2 > 0.
\end{aligned}
$$

Therefore, $\epsilon_{i+4}\approx C \epsilon_{i+3} \epsilon_{i+2}$. A full cycle consists of 4 iterations and 

$$
\epsilon_{i+8} \simeq C^7 \epsilon_{i+4}^6 \epsilon_{i+3}^2 \simeq C^{8}\epsilon_{i+4}^7 \epsilon_i^2.
$$

The order of convergence $p$ solves $p^2 - 7p - 2 = 0$, which gives $p = \sqrt[4]{\frac{7+\sqrt{57}}{2}}$.
```

```{prf:remark}
:label: rmk-pegasus-asymptotic
Actually, the asymptotic expansion of $\epsilon_{i+4}$ is (expanded in $\epsilon_{i+1}$ first )

$$
\epsilon_{i+4} = C^5 \epsilon_{i}^4\epsilon_{i+1}^2 + C^7 \epsilon_i^7\epsilon_{i+1}
+ O(\epsilon_{i}^8 \epsilon_{i+1} + \epsilon_{i}^5 \epsilon_{i+1}^2 + \epsilon_{i}^3 \epsilon_{i+1}^3)$$

It is at first not clear why we can retain the first term and drop the rest, because it requires the following inequalities to hold

$$\epsilon_{i}^3 \ll \epsilon_{i+1}\ll \epsilon_i.$$

The latter one is correct because of the relation $\epsilon_{i+1} \simeq C \epsilon_{i}\epsilon_{i-1} \ll \epsilon_i$ once the iterations are close to the root. The former one is equivalent to $\epsilon_{i}^2 \ll \epsilon_{i-1}$. This can be made into an assumption because $\epsilon_{i-1}$ and $\epsilon_{i}$ are the last two iterations in the previous cycle, which correspond to $\epsilon_{i+3}$ and $\epsilon_{i+4}$ in the current cycle, we can use the above asymptotic estimate to find $\epsilon_{i+3}^3 \approx C^6 \epsilon_{i}^9 \epsilon_{i+1}^3 \ll \epsilon_{i+4}$. 

Therefore, if the inequality does not hold, one can use the current cycle as the starting point to perform the same analysis.
```

In the Pegasus method, two consecutive standard false position steps and two adjustment steps are performed in a full cycle. Although the decay over a full cycle is significant, such advantage will be gone if the cycle is long.

In order to make the cycle shorter, we need to drop at least one step. According to the previous analysis (see Illinois method), the standard false position step does not change the sign of the error asymptotically, thus it is preferred to drop one of the standard false position steps.

Let us finish this section with a brief discussion on the **Improved Pegasus** method {cite}`king1973improved`, which takes the advantage of the symmetry in the false position method to avoid consecutive standard false position steps. The similar technique can be also applied to other methods such as **Anderson-Bjorck** method {cite}`anderson1973new`.

```{margin}
The **Improved Pegasus** method can also use the decay factor of **Anderson-Bjorck** method as

$$\lambda = \begin{cases} 1 - \frac{f_2}{f_1}&\text{ if positive},\\
        \frac{1}{2}&\text{ otherwise}.\end{cases}
$$
```

```{prf:algorithm} Improved Pegasus Method
:label: AL-IMPROVED-PEGASUS

**Inputs** $f$, $a, b$, $\epsilon$ 

**Outputs** $c$ 

1. $x_0\gets a$, $x_1\gets b$, $f_0\gets f(x_0)$, $f_1\gets f(x_1)$.  //initialization.
2. While True do:
    1. $x_2 \gets \frac{f_1 x_0 - f_0 x_1}{f_1 - f_0}$, $f_2 \gets f(x_2)$. //standard false position step, $x_1$ and $x_2$ are two latest iterations.
    2. If $ |f_2| < \epsilon$, then return $x_2$; //check stopping criteria
    3. If $f_1 f_2 < 0$, then swap $(x_0, f_0)$ and $(x_1, f_1)$. // avoid false position step
    4. While $f_1 f_2 > 0$, then // adjust until the sign changes
        - $(x_0, f_0)\gets (x_0, \lambda f_0)$, where $\lambda =\frac{f_1}{f_1 + f_2}$.
        - $(x_1, f_1)\gets (x_2, f_2)$ and  $x_2 \gets \frac{f_1 x_0 - f_0 x_1}{f_1 - f_0}$, $f_2 \gets f(x_2)$
    5. If $f_1 f_2 < 0$, then // perform a false position step
        - $(x_0, f_0)\gets (x_1, f_1)$ and $(x_1, f_1)\gets (x_2, f_2)$. 

```

```{prf:theorem}
:label: thm-improved-pegasus-convergence
The **Improved Pegasus** method has an order of convergence at least $\sqrt[3]{5}$.
```

```{prf:proof}

In the same setting as the previous step, we will find the first iteration is the same as the Pegasus method (false position) that

$$\epsilon_{i+1}\simeq C \epsilon_i \epsilon_{i-1} + D \epsilon_i \epsilon_{i-1}(\epsilon_{i} + \epsilon_{i-1}) < 0.$$

For the second iteration, although $x_{i+1}$ and $x_{i}$ are on different sides, but a false position step has been performed in the previous step, thus it will perform adjustment step with $\lambda = \frac{f(x_{i}) - f(x_{i+1})}{f(x_{i})}$, then

$$\epsilon_{i+2} \simeq -D \epsilon_{i+1}\epsilon_{i}\epsilon_{i-1},$$

note the leading term is different from the Pegasus method. There are two options.

- If $D < 0$, then $\epsilon_{i+2} > 0$ which completes a cycle with two iterations, which is very compact. In this case, we find

    $$\epsilon_{i+4} \simeq C^{2}  \epsilon_{i+2}^3,$$

    which implies an order of convergence at $\sqrt{3}$.

- If $D > 0$, then it will perform another adjustment step,

    $$\epsilon_{i+3} =  \frac{\epsilon_{i+2}\lambda f(x_i) - \epsilon_i f(x_{i+2})}{\lambda f(x_{i}) - f(x_{i+2})},\quad \lambda = \frac{f(x_{i+1})-f(x_{i+2})}{f(x_{i+1})}\frac{f(x_{i})-f(x_{i+1}) }{ f(x_{i})}$$

    which will make $\epsilon_{i+3} \simeq -C^3 D \epsilon_{i-1}^3 \epsilon_i^3 > 0$, which completes a cycle with three iterations. In this case, we find

    $$\epsilon_{i+6} \simeq D^{2}  \epsilon_{i+3}^5,$$

    which implies an order of convergence at $\sqrt[3]{5}$.

```

## Iterative Methods

### Newton-Raphson Method

### Secant Method

## Applications in Optimization

## Exercises

### Theoretical Part

```{admonition} Problem 1
The only difference between {prf:ref}`AL-IMPROVED-PEGASUS` and the usual Pegasus method is at the 3rd step in the while loop, which eliminates consecutive false position steps. This change is very simple, but it leads to an improvement in the order of convergence.

Explain why the Illnois method $\lambda = \frac{1}{2}$ cannot be faster by the same technique.
```

### Computational Part

## Extended Reading

```{bibliography}
:filter: docname in docnames
```
