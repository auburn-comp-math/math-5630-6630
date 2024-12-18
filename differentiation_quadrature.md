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

# Differentiation and Quadrature

A vast number of applications such as the calculation of tangent vectors or areas lead to the problem of computing

$$
    \mathcal{D}(f) := \frac{d}{dx} f(x),\quad \mathcal{I}(f) := \int_a^b f(x) dx,
$$

for certain function $f(x)\in C^k([a, b])$. Accurate evaluations would sometimes be difficult if an analytic expression is absent. Especially when the function values of $f$ are only accessible at a finite number of nodes. Therefore, it is important to find simple yet effective methods to approximate the derivatives and integrals.

## Extrapolation

From the previous discussion, we already know that interpolation provides an estimate **within** the original observation range. The extrapolation is similar but aims to produce estimates outside the observation range. However, sometimes extrapolation may be subject to a greater uncertainty (see following example), one should use it only when an overestimate is hardly occurring.

```{code-cell} ipython3
:tags: [remove-input]
:mystnb:
:   image:
:       align: "center"
:   figure:
:       caption: Extrapolation behavior of polynomial fitting on Chebyshev nodes.

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.set_loglevel("critical")

%matplotlib inline
def runge(x):
  return 1 / (1 + x ** 2)

def polyfit(x, y):
  return np.poly1d(np.polyfit(x, y, x.size -1))

def overlay_plot_fit(f, f_fit, x, y, x0=-1, x1=1):
    _x = np.linspace(x0, x1, 500)
    plt.plot(_x, f_fit(_x), label='fitted function')
    plt.plot(_x, f(_x), label='actual function')
    plt.plot(x, y, 'k.')
    # plt.title('Polynomial fitting for Runge function')
    plt.legend(loc='lower center', fontsize=9)

# x = np.linspace(-5, 5, 11) # use chebyshev nodes todo
x = 5 * np.cos(-np.linspace(0, np.pi, 15))
y = runge(x)
runge_fit = polyfit(x, y)

plt.figure(figsize=(6, 4))

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

overlay_plot_fit(runge,runge_fit, x, y, -5.2,5.2)
```

### Richardson Extrapolation

Suppose there is a sequence of estimates $A(h)$ depending on the parameter $h$ smoothly, the limit $A^{\ast} = \lim_{h\to 0^{+}} A(h)$ is the quantity to be computed. In practice, we only have access to $A(h)$ for a few values of $h$. Using these values to estimate $A^{\ast}$ is a typical problem in extrapolation.

The basic idea behind **Richardson extrapolation** is to use polynomial interpolation with a sequence of nodes $h_j\to 0$. Suppose that the function $A(h)$ admits the following asymptotic expansion:

$$
    A(h) = a_0 + a_1 h^{\gamma} + a_2 h^{2\gamma} + \dots + a_k h^{k\gamma} + \cO(h^{(k+1)\gamma})
$$

for any $h > 0$ and $k\ge 0$. Then $A^{\ast} = a_0$ and $A(h) = A^{\ast} + \cO(h^{\gamma})$. Suppose we have access to the values $A(h_0),\dots, A(h_n)$, then this uniquely determines a polynomial $f_n\in\Pi_n$ and
$f_n(h_j^{\gamma}) = A(h_j)$. We will approximate $A(0)\approx f_n(0)$. The computation of $f_n$ follows the construction of the Newton form.

```{prf:lemma} Richardson Extrapolation
:label: lem-richardson

Suppose $h_j$ can be represented as 

$$
h_j = \frac{\hbar}{t_j}
$$

for some adjustable parameter $\hbar$ and scaling constants $ 1 < t_0 < t_1<\dots < t_{n-1}$. Then 

$$
f_n(0) = A^{\ast} + (-1)^n \frac{a_{n+1}}{\prod_{j=0}^n t_j^{\gamma}} \hbar^{(n+1)\gamma} + \cO(\hbar^{(n+2)\gamma}),\quad \text{as }\hbar\to 0.$$
```

```{prf:proof}
We view $A(h)$ as a polynomial with respect to $h^{\gamma}$ of degree $(n+1)$ with an addition perturbation $\cO(h^{(n+2)\gamma})$. Then we have the following. 

$$
      A(h) = p_{n+1}(h^{\gamma}) + \cO(h^{(n+2)\gamma}).
$$

Let $\tilde{f}_n$ be the interpolation polynomial of degree $n$ to $p_{n+1}$, 

$$
    p_{n+1}(x) \equiv \tilde{f}_n(x) + p[x, h_0, h_1, \dots, h_n] \prod_{j=0}(x - h_j^{\gamma}), 
$$

where $p[x, h_0, h_1, \dots, h_n]$ is the coeffimathcal{I}ent of the leading power in $p_{n+1}$, $a_{n+1}$. Thus, 

$$
    A^{\ast} = p_{n+1}(0) = \tilde{f}_n(0) +a_{n+1} \prod_{j=0}(0 - h_j^{\gamma}). 
$$

Here we use the result for the stability of polynomial interpolation in earlier chapters. Therefore, 

$$
    |\tilde{f}_n(0) - f_n(0)|\le \lambda_n(0) \cdot \cO(\hbar^{(n+2)\gamma}).
$$
Here the Lebesgue function at zero $\lambda_n(0)$ is 

$$
    \lambda_n(0) = \sum_{j=0}^n \prod_{k=0, k\neq j}^n \left|\frac{h_k^{\gamma}}{h_k^{\gamma} - h_j^{\gamma}}\right| =  \sum_{j=0}^n \prod_{k=0, k\neq j}^n \left|\frac{1}{1 - (\frac{t_k}{t_j})^{\gamma}}\right|,
$$

which is independent of $\hbar$.

```

The Richardson extrapolation considers the special choice of $t_j = t^j$ for some $t > 1$. The error estimate then is

$$
    f_n(0) = A^{\ast} + \left(\frac{(-1)^n}{t^{n(n+1)\gamma/2}} a_{n+1}\right) \hbar^{(n+1)\gamma} + \cO(\hbar^{(n+2)\gamma}).
$$

There are easier ways to calculate the Richardson extrapolation using the following expansion.

$$
    A(h) - t^{\gamma} A\left(\frac{h}{t}\right) = (1 - t^{\gamma}) A^{\ast} + \cancel{a_1 \left(h^{\gamma} - t^{\gamma} \left(\frac{h}{t}\right)^{\gamma}\right)} + a_2\left(h^{2\gamma} - t^{\gamma}\left(\frac{h}{t}\right)^{2\gamma}\right) +\dots.
$$

Let $A_1(h) = \frac{A(h) - t^{\gamma} A(\frac{h}{t}) }{1 - t^{\gamma}} $, we obtain the first iteration result as

$$
    A^{\ast} \approx A_1(h) + \cO(h^{2\gamma}),
$$

then follow the same idea, we cancel the $\cO(h^{2\gamma})$ term by

$$
    A_1(h) - t^{2\gamma} A_1\left(\frac{h}{t^2}\right) = (1 - t^{2\gamma}) A^{\ast} + \cO(h^{3\gamma}).
$$

Therefore by taking $A_2(h) = \frac{A_1(h)- t^{2\gamma} A_1(\frac{h}{t^{2\gamma}})}{1 - t^{2\gamma}}$, the second iteration satisfies

$$
    A^{\ast}\approx A_2(h) + \cO(h^{3\gamma}).
$$

However, such a process may not constantly refine the approximation due to the potentially fast-growing constant in the $\cO$ notation.

### Aitken Extrapolation

Alexander Aitken rediscovered Aitken extrapolation, which has been used to accelerate a sequence's convergence.

Given a sequence $S = \{s_n\}_{n\ge 0}$, the Aitken extrapolation generates a new sequence

$$AS = \left\{\frac{s_n s_{n+2} - s_{n+1}^2}{s_n + s_{n+2} - 2 s_{n+1}}\right\}_{n\ge 0}. $$

A more stable formulation (why?) is written as 

$$(AS)_n  = s_n - \frac{(\Delta s_n)^2}{\Delta^2 s_n},$$

where $\Delta$ is the forward difference operator that $\Delta s_n = s_{n+1} - s_n$. It is not difficult to see that Aitken extrapolation can accelerate linearly convergent sequences. 

```{prf:theorem}
:label: aitken_extrapolation
Assume that the sequence $S = \{s_n\}_{n\ge 0}$ satisfies 

$$
    \lim_{n\to\infty} \frac{|s_{n+1} - \mu|}{|s_n - \mu|} = \rho \in (0, 1).
$$

Then the accelerated sequence $AS$ converges faster than $S$.
```

```{prf:proof}
Let $\rho_n := \frac{s_{n+1} - \mu}{s_n - \mu}$, without loss of generality, we assume $\lim_{n\to\infty} \rho_n = \rho$.  
By the definition of the acceleration formula, we obtain 

$$
    \frac{ (AS)_{n} - \mu}{ s_{n} - \mu } =  1 - \frac{(\Delta s_{n})^2}{\Delta^2 s_{n}} \frac{1}{s_n -\mu}.
$$

For sufficiently large $n$, we have $\rho_n \approx \rho$ and 

$$
\begin{aligned}
    \Delta s_n &= s_{n+1} - \mu - (s_n -\mu) = (\rho_n - 1) (s_n -\mu)    \\
    \Delta^2 s_n &= s_{n+2} - \mu + s_n -\mu - 2(s_{n+1} -\mu) = (\rho_{n+1}\rho_n-2\rho_{n}+1)(s_n -\mu)
\end{aligned}
$$

Therefore

$$
    \frac{ (AS)_{n} - \mu}{ s_{n} - \mu } = 1 - \frac{(\rho_n-1)^2}{(\rho_{n+1}\rho_n - 2\rho_{n} + 1)}=\frac{\rho_n(\rho_{n+1} - \rho_n)}{(\rho_{n+1}\rho_n - 2\rho_{n} + 1)}  \to 0. 
$$
```

```{prf:remark}
Let $ \eps_n := s_n - \mu$. If the error satisfies the following relation (common in most fixed-point iterations) $$\eps_{n+1} = \eps_n (\rho + c_1 \eps_n + o(\eps_n)),$$ then

$$
    \frac{ (AS)_{n} - \mu}{ s_{n} - \mu }  = \frac{\rho_n(\rho_{n+1} - \rho_n)}{(\rho_{n+1}\rho_n - 2\rho_{n} + 1)} \approx \frac{c_1 \rho}{(1-\rho)^2} \left((\rho - 1)\eps_n + o(\eps_n) \right).
$$

That is, $| (AS)_{n} - \mu| = \cO( |s_{n} - \mu |^2 )$.
```

The Aitken extrapolation can be used to accelerate fixed-point iteration solving the root $x^{\ast}$ for $f(x)$. Steffensen's method is a root-finding algorithm based on such an acceleration technique, the iteration is given by 

$$
        x_{n+1} = x_n - \frac{f(x_n)}{g(x_n)} ,\quad 
        g(x)= \frac{f(x + f(x))}{f(x)} - 1.
$$

If $f$ is twice differentiable, then clearly the function $g(x)\approx f'(x)$, which recovers the Newton-Raphson method asymptotically if $f(x)\approx 0$.

```{prf:theorem}
Suppose $f'(x^{\ast})\in (-1, 0)$, then the order of convergence is 2 for Steffensen's method.
```

```{prf:proof}
Let us first point out the connection between Aitken extrapolation and Steffensen's method. Denote $x_0$ as the starting point, the intermediate values $z_1 = x_0 +f (x_0)$ and $z_2 =z_1 +f(z_1)$. The Aitken extrapolation finds

$$
\begin{aligned}
A x_0 &= x_0 - \frac{(z_1 - x_0)^2}{z_2 + x_0 - 2 z_1 } = x_0 - \frac{|f(x_0)|^2}{z_1 + f(z_1) + x_0 - 2 (x_0 + f(x_0)) }     \\
&=x_0 - \frac{|f(x_0)|^2}{f(z_1) - f(x_0)} = x_0 - \frac{f(x_0)}{g(x_0)} = x_1.
\end{aligned}
$$

The convergence is linear when the iterates $x_0, z_1, z_2$ are close to the root $x^{\ast}$. Therefore, using the same derivation, we find 

$$
    x_1 - x^{\ast} = A x_0 - x^{\ast} = \cO( | x_0 - x^{\ast}|^2 ).
$$
```

```{prf:remark}
It should be noticed that Steffensen's method evaluates $f$ twice during each iteration, the same as Newton's method (one for $f$ and one for $f'$, although the evaluation of $f'$ is more expensive in practice). Each evaluation brings an order of $1$ of convergence on average. From the viewpoint of efficiency, the secant method is preferred, since it only evaluates $f$ once each iteration with an order of $\frac{1+\sqrt{5}}{2}$ of convergence. 
```

### Wynn's Epsilon Method

Wynn's $\eps$ method {cite}`wynn1966convergence,wynn1956device` is another kind of extrapolation algorithm that is recommended as the best all-purpose acceleration method. It has a strong connection with Padé approximation and continued fractions. We will not cover the detailed derivation of the theory in this section. However, Wynn's $\eps$ method still has its limitations if the sequence converges to the desired value too slowly.

The algorithm is stated as follows.

```{prf:algorithm} Wynn's $\eps$ Method
:label: AL-WYNN-EPS
**Input**: $s_0, s_1, \dots, s_n,\dots$, which is a sequence converging to the desired quantity.

**Output**: $\eps_{2l}^{(j)}$, $j, l=0,1, \dots$.

1. For $j = 0, 1, 2,\dots$, set 
  
  $$
  \eps_{-1}^{(j)} = 0 \text{ (guarding elements)}, \quad \eps_{0}^{(j)} = s_j.
  $$ 

2. For $j,k = 0,1,2,\dots$, 

  $$
  \eps_{k+1}^{(j)} = \eps_{k-1}^{(j+1)} + [\eps_{k}^{(j+1)} - \eps_{k}^{j}]^{-1}.
  $$

```

```{prf:example}
It is known that $\frac{\pi}{4}$ can be calculated by the asymptotic expansion:

$$
    \arctan z = z - \frac{z^3}{3} + \frac{z^5}{5} - \dots 
$$

at $z = 1$. Define the function $A(h)$ such that 

$$
    A(h)= \sum_{j=0}^{1/h} \frac{(-1)^j}{2 j + 1} = \frac{\pi}{4} + a_1 h + a_3 h^3 + \dots
$$

Then the approximation error is $\cO(h)$, which is very slow. Taking $h=10^{-3}$ has around $2\times 10^{-4}$ error.  We test with two extrapolation algorithms.

- Richardson Extrapolation. We take $1/h = 250, 500, 1000, 2000$ and calculate that the Richardson extrapolation three times would result in almost machine premathcal{I}sion.
- Wynn's $\eps$ method. We take the sequence 

  $$
  s_k = \sum_{j=0}^k \frac{(-1)^j}{2 j + 1}
  $$

  as the truncated series at $z = 1$. With about 20 terms, we already reached machine precision.
```

```{code-cell} ipython3
:tags: [remove-input]
:mystnb:
:   image:
:       align: "center"
:   figure:
:       caption: Absolute error of $\eps_{2l}^{(j)}$ for Wynn's $\eps$ method, $0\le l\le 5$.

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.set_loglevel("critical")
from bigfloat import *

prec = 900

def wynn(s):
    with precision(prec):
        n = (len(s) - 1) // 2
        N = 2 * n + 1
        e = [ [BigFloat('0', precision(prec))] * (N+1) for _ in range(N+1)]# 2 * (n + 2)


        for j in range(1, N + 1):
            e[j][1] = s[j - 1]

        for k in range(3, N + 2):
            for j in range(3, k + 1):
                e[k-1][j-1] = e[k- 2][j - 3] + BigFloat('1', precision(prec))/(e[k - 1][j - 2] - e[k - 2][j - 2])
        return e

with precision(prec):
    s = []
    ret = 0
    for i in range(20):
        ret = ret + (-1)**i * (BigFloat('1', precision(prec))) / (2 * i + 1)
        s.append(ret)

    t = wynn(s)
    v = const_pi()/4 

    for i in range(len(t)):
        for j in range(len(t[0])):
            t[i][j] =  t[i][j] -  v


    fig = plt.figure()
    ax = fig.add_subplot(111)

    for m in range(6):
        plt.loglog(range(2 * m + 1, len(t)), [(abs(t[j][2 * m + 1])) for j in range(2 * m + 1, len(t))], "-o", label="error at iter {}".format(m))

    plt.legend()
    plt.show()
```

## Quadrature

The numerical quadrature finds the value of an integral

$$\mathcal{I}(f) = \int_a^b f(x) dx $$

from the function values at a finite number of points. We are mostly interested in the following quadrature formula

$$
    \mathcal{I}_n(f) = (b - a) \sum_{j=0}^n w_j f(x_j)
$$

where $x_j$ are the **nodes**and $w_j$ are the **weights**. Similar to the numerical derivatives, we also define some terminologies. The formula $\mathcal{I}_n$ is said to have **degree $k$ accuracy** if $\mathcal{I}_n$ is exact for all polynomials $f\in\Pi_k$. Since the integration formula is linear, the exactness can be rephrased as

$$
    \begin{aligned}
        \mathcal{I}(x^p) &= \mathcal{I}_n(x^p),\quad 0\le p \le k, \\
        \mathcal{I}(x^{k+1}) &\neq \mathcal{I}_n(x^{k+1}).
    \end{aligned}
$$

### Interpolation Based Quadrature

The interpolation-based idea is intuitive. Let $q_n(x)$ be the interpolation polynomial on the nodes $x_j$ with values $f(x_j)$, $j=0,1,\dots, n$, respectively. We define the quadrature by interpolation formula as

$$
   \mathcal{I}_n(f):= \int_a^b q_n(x) dx.
$$

The above quadrature formula is exact for all degree $n$ polynomials $f$, therefore it has at least degree $n$ accuracy.

```{prf:remark}
:label: rem-runge-quad
We have seen that the $L^{\infty}$ error between $f$ and $q_n$ could be large (e.g. Runge phenomenon) as $n\to \infty$. So the nodes would be important as well for quadratures.
```

Using the Lagrange polynomials, we can represent

$$
q_n(x) = \sum_{j=0}^n f(x_j) L_j(x),\quad L_j(x) = \prod_{k=0, k\neq j}^n \frac{x - x_k}{x_j - x_k}.
$$

Then it is not difficult to derive

$$
\begin{aligned}
    \mathcal{I}_n(f)& = \sum_{j=0}^n f(x_j) \int_a^b L_j(x) dx  \\ &=(b-a) \sum_{j=0}^n f(x_j) \int_0^1 \prod_{k=0, k\neq j }^n \frac{t - t_k}{t_j - t_k} dt, \quad t_j = \frac{x_j - a}{b-a}.
\end{aligned}
$$

therefore the weights $w_j =  \int_0^1 \prod_{k=0, k\neq j }^n \frac{t - t_k}{t_j - t_k} dt$.

```{prf:example}
:label: quad-rec-rule
The rectangle rule is the simplest, where we choose $x_0 = \frac{a+b}{2}$ as the middle point. Then the quadrature rule writes

$$\mathcal{I}_{0, rectangle}(f)  = (b-a) f(x_0). $$

Such a rule is exact for any linear function, so it has a degree of accuracy of one. We can see that the degree of exactness could exceed $n$. In the next chapter, we will see that the maximum degree of exactness for such a form is $2n+1$.
```

```{prf:example}
:label: quad-trap-rule
The trapezoid rule takes $x_0 = a$ and $x_1= b$. 

$$\mathcal{I}_{1,trapezoid} (f) = (b-a)\left( \frac{1}{2} f(x_0) + \frac{1}{2} f(x_1) \right).$$

One can check that this rule is exact for $f(x) = 1, x$. Its degree of accuracy is one. It has a slightly larger constant in error estimate than the rectangle rule.
```

```{prf:example}
:label: quad-simp-rule
The Simpson's rule takes $x_0 = a$, $x_1 = \frac{a+b}{2}$, and $x_2 = b$. 

$$\mathcal{I}_{2, Simpson}(f) = (b-a) \left(\frac{1}{6}f(x_0) + \frac{2}{3}f(x_1) + \frac{1}{6} f(x_2) \right).$$

This rule is exact for $f(x) = 1, x, x^2, x^3$, thus the degree of accuracy is three. 
```

### Newton-Cotes Quadrature

### Romberg Integration

### Gaussian Quadrature

### Monte Carlo Integration

## Exercises

### Theoretical Part

<!-- ```{admonition} Problem 1
Let $f(x)$ be a real-valued smooth function with a simple root $x^{\ast}\in [a, b]$. Suppose the bisection method produces a sequence of approximation $\{x_k\}_{k=1}^{\infty}$, show that 

- The Wynn's $\eps$ method applied to $\{x_k\}_{k=1}^n$ in general cannot accelerate the convergence of the sequence to $x^{\ast}$. 

- If the sequence is generated by the false position method, then the Wynn's $\eps$ method can accelerate the convergence.
``` -->

### Computational Part

## Extended Reading

```{bibliography}
:filter: docname in docnames
```
