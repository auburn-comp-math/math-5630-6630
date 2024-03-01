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

# Interpolation

The interpolation (1D) solves problems of the following type:

> Given a set of predefined functions $\mathcal{K}$, find an element $f: \mathbb{I}\mapsto \mathbb{R}$ in $\mathcal{K}$ such that $y_j = f(x_j)$ for all $j=0,\dots, n$.

Here $\mathbb{I}$ denotes a finite or infinite interval such that $x_1,\dots x_n\in \mathbb{I}$. One of the important applications for interpolation is Computer-Aided Design (CAD) which is used extensively in the manufacturing industry. Generally speaking, the interpolation provides a closed form of the function to determine the value of $y$ where the parameter $x$ is not accessible.

## Polynomial Interpolation

The polynomial interpolation considers the set $\mathcal{K} = \Pi_m$, where the set $\Pi_m$ represents the polynomials of with degree $ \le m$. We will seek for a polynomial $f(x)$ with the constraints that

$$
    \begin{cases}
        f\in \mathcal{K} = \Pi_n,&\\
        f(x_k) = y_k &\text{ for } k = 0, 1,\dots, n.
    \end{cases}
$$

The points $x_k$ are called **interpolation nodes**, if $m > n$ (resp. $m < n$), the problem is underdetermined (resp. overdetermined). For the case that $ m = n$, we have

```{prf:theorem}
:label: THM-INTER-UNIQ
There exists a unique polynomial function $f\in \Pi_n$ such that $f(x_j) = y_j$ for $j=0,\dots, n$.
```

```{prf:proof}
**Existence**: In order to construct the polynomial $f$, it is straightforward to consider the general form of polynomial $f(x) = \sum_{j=0}^n a_j x^j$, then we can formulate a linear system for the coefficients $a_j$, $j=0,\dots, n$, which is 

$$
\begin{pmatrix}
    1 & x_0 & x_0^2 & \dots & x_0^n \\
    1 & x_1 & x_1^2 & \dots & x_1^n \\
    \vdots & \vdots &\vdots & \ddots & \vdots\\
    1 & x_n & x_n^2 & \dots & x_n^n
\end{pmatrix} \begin{pmatrix}
    a_0\\a_1\\\vdots\\a_n
\end{pmatrix} = \begin{pmatrix}
    y_0\\y_1\\\vdots \\y_n
\end{pmatrix}.
$$

The matrix 

$$
V = \begin{pmatrix}
    1 & x_0 & x_0^2 & \dots & x_0^n \\
    1 & x_1 & x_1^2 & \dots & x_1^n \\
    \vdots & \vdots &\vdots & \ddots & \vdots\\
    1 & x_n & x_n^2 & \dots & x_n^n
\end{pmatrix} 
$$

is called **Vandermonde matrix**. To determine the coefficients $a_j$, one needs the matrix $V$ be invertible. Its determinant can be computed (as an exercise) as

$$
        \det(V) = \prod_{0\le i\le j\le n}(x_j - x_i).
$$
When $x_j$ are distinct, the determinant is nonzero.

**Uniqueness**:  Suppose there are two distinct polynomials $f, g\in \Pi_n$ satisfying the condition that $f(x_j) = g(x_j) = y_j$, then $f - g$ has $(n+1)$ roots $x_j$, $j=0, \dots, n$. If $f\neq g$, it is clear that $f-g\in\Pi_n$ has at most $n$ roots. Contradiction.

```

In the above proof, the interpolation polynomial can be uniquely determined by solving the linear system

$$        \begin{pmatrix}
    1 & x_0 & x_0^2 & \dots & x_0^n \\
    1 & x_1 & x_1^2 & \dots & x_1^n \\
    \vdots & \vdots &\vdots & \ddots & \vdots\\
    1 & x_n & x_n^2 & \dots & x_n^n
\end{pmatrix} \begin{pmatrix}
    a_0\\a_1\\\vdots\\a_n
\end{pmatrix} = \begin{pmatrix}
    y_0\\y_1\\\vdots \\y_n
\end{pmatrix}.
$$

However, it is generally easier to compute the polynomial $f$ with the **Lagrange polynomial interpolation** (which is somewhat equivalent to compute the inverse of $V$).

### Lagrange Polynomial

```{prf:definition}
:label: DEF-LA-PO
For the given distinct $x_j$, $j = 0, 1, \dots, n$, the $(n+1)$ Lagrange polynomials $L_0, L_1,\dots, L_n\in\Pi_n$ are defined by 

$$
    L_j(x) = \prod_{s = 0, s\neq j}^n \frac{x - x_s}{x_j - x_s}, \quad j = 0, 1,\dots , n.
$$
```

It is clear that these polynomials satisfy the conditions that

$$
    L_j(x_k) = \delta_{jk} := \begin{cases}
        1&\text{for } k=j,\\
        0&\text{for } k\neq j.
    \end{cases}
$$

Therefore these polynomials are linearly independent, which form a basis of the $(n+ 1)$-dimensional space $\Pi_n$.

```{prf:theorem}
:label: THM-UNIQ-LAG-INTER
The unique interpolating polynomial $f$ satisfying $f(x_j) = y_j$, $j=0,1,\dots, n$ can be represented by 

$$
    f(x) = \sum_{j=0}^n y_j L_j(x).
$$
```

```{prf:proof}
It is straightforward to check the interpolation conditions are satisfied.
```

```{prf:remark}
We introduce a preliminary procedure to compute value of the interpolating polynomial $f$ at a point $x$. Let constants $k_j$ and $q(x)$ be defined as 

$$
    k_j = \prod_{s= 0, s\neq j}\frac{1}{x_j - x_s},\quad q(x) = \prod_{j=0}^n (x - x_j),
$$

then 

$$
    f(x) = \sum_{j=0}^n y_j L_j(x) = q(x)  \sum_{j=0}^n k_j y_j \frac{1}{x - x_j}.
$$

One can first compute $k_j$ with $\cO(n^2)$ flops, then $f(x)$ can be computed with $\cO(n)$ flops. The advantage of the above scheme is the constants $k_j$ are independent of $y_j$, therefore evaluating another instance of the interpolating polynomial will not need to re-compute them. The disadvantage is that if we add a new node, the constants $k_j$ have to be updated with an additional cost of $\cO(n)$ flops. Later we will see the Newton's form can overcome this issue.
```

### Interpolation Error

When the data pairs $(x_j, y_j)$, $j=0,1,\dots, n$ are generated by a sufficiently smooth function $h(x)$, it is possible to quantify the error between the interpolating polynomial $f(x)$ and $h(x)$.

```{prf:theorem}
:label: THM-INTERP-ERROR
Let $h: [a, b]\mapsto \mathbb{R}$ be a $(n+1)$-times differentiable function. If $f(x)\in\Pi_n$ is the interpolating polynomial that 

$$
f(x_j) = h(x_j),
$$

for $j=0,1,\dots, n$. Then for each $\overline{x}\in [a, b]$, the error can be represented by 

$$
h(\overline{x}) - f(\overline{x}) = \frac{\omega(\overline{x})}{(n+1)!} h^{(n+1)}(\xi),
$$

where $\xi = \xi(\overline{x})\in [a, b]$ and $\omega(x) = \prod_{j=0}^n (x - x_j)$.
```

```{prf:proof}
The proof is based on the Rolle's Theorem. Select any $\overline{x}\in[a, b]$ such that $\omega(\overline{x})\neq 0$, then let 

$$
\psi(x) = h(x) - f(x) - k\omega (x)
$$

the constant $k$ is chosen such that $\psi(\overline{x}) = 0$. Then $\psi(x) = 0$ at $(n+2)$ points 

$$
x_0, x_1, \dots, x_n, \overline{x}\in [a, b]
$$

By Rolle's Theorem, $\psi^{(n+1)}$ has at least one zero $\xi$ in $[a,b]$. Therefore

$$
    \psi^{(n+1)}(\xi) = h^{(n+1)}(\xi) - 0 - k(n+1)! = 0.
$$
```

```{prf:corollary}
If $h(x)\in C^{\infty}([a, b])$ satisfies that $\max_{x\in[a,b]} |h^{(n)}(x)|\le M <\infty$ for all $n\ge 0$, then the interpolating polynomial approximates $h$ uniformly as the number of nodes $n\to \infty$.
```

```{prf:proof}
Since $|x -x_j|\le b-a$, the error is bounded by $\frac{(b-a)^{n+1}}{(n+1)!} M$, which converges to zero.
```

It is interesting to think about the converse: under what kind of condition the interpolation error is not vanishing as the number of nodes tends to infinity? From the {prf:ref}``THM-INTERP-ERROR``, the error depends on the sizes of three terms.

- The bound of the $(n+1)$-th derivative, $\max_{x\in[a,b]}|h^{(n+1)}(x)|$. This could grow rapidly. For instance, $h(x) = 1/\sqrt{x}$ on $[\frac{1}{2}, \frac{3}{2}]$,

$$
h^{(n+1)}(x) = \frac{(-1)^{n+1}}{2^{n+1}} (2n+1)!! x^{-(2n+3)/2}.
$$

- The function $\omega(x) = \prod_{j=0}^n (x - x_j)$, such product could be large if $x$ and the nodes $x_j$ are not so close.

- The term $\frac{1}{(n+1)!}$, which decays fast.

We can see that for the function $h(x) = 1/\sqrt{x}$ on $[\frac{1}{2}, \frac{3}{2}]$, it is not trivial to show the interpolating polynomial could converge to $h$ anymore (it is still true for certain choices of $x_j$).

Next, we try to provide a better estimate of $\omega$ for the special choice: equally spaced nodes.  Let the nodes $x_j = a + j\Delta$, where $\Delta = \frac{b-a}{n}$. It is not difficult (prove it) to see $\omega(x)$ will be the worst if $x$ is located on the end sub-intervals, $[x_0, x_1]$ and $[x_{n-1}, x_n]$. Without loss of generality, we assume $x$ is located on $[x_0, x_1]$, then

$$|x - x_j|\le j \Delta$$

for $j = 2, \dots, n$, which implies

$$
    |\omega(x)|\le \prod_{j=0}^n |x - x_j|\le n! \Delta^{n-1} \sup_{x\in [x_0, x_1]} |(x - x_0)(x-x_1)|= \frac{n!}{4} \frac{(b-a)^{n+1}}{n^{n+1}}.
$$

Thus the interpolation error is bounded by

$$  \|h  - f\|_{\infty} \le  \frac{\sup_{x\in[a,b]}|h^{(n+1)}(x)|}{4(n+1)} \frac{(b-a)^{n+1}}{n^{n+1}}.$$

Such an estimate is useful to derive uniform convergence.

```{prf:example}
Consider $h(x) = 1/x$ on $[\frac{1}{2}, \frac{3}{2}]$. Then $h^{(n+1)}(x) = \frac{(n+1)!(-1)^{n+1}}{x^{n+2}}$, hence 

$$
    \frac{| h^{(n+1)}(x) |}{4(n+1)}  \frac{(b-a)^{n+1}}{n^{n+1}} \le \frac{1}{4n^{n+1}}\max_{x\in[\frac{1}{2}, \frac{3}{2}]} \left| \frac{1}{x^{n+2}} \right| = \frac{1}{2} \left(\frac{2}{n}\right)^{n+1},
$$

therefore, the interpolation error converges to zero exponentially. It is important to notice that the above method only works for intervals away from the origin. 

```

### Runge's Phenomenon

From the above discussion, we can see there is a possibility that $\max_{x\in[a,b]}|h^{n+1}(x)|\omega(x)$ grows faster than $(n+1)!$, which would lead to divergence. Hence increasing the number of interpolation nodes (at least for equally spaced nodes) is not guaranteed for better approximation. The most famous example is the one made by Carl Runge.

```{math}
:label: EQ-RUNGE-EXAMPLE
h(x) = \frac{1}{1+x^2},\quad x\in [-5, 5].
```

```{code-cell} ipython3
:tags: [remove-input]
:mystnb:
:   image:
:       align: "center"
:   figure:
:       caption: Runge phenomenon with 11 equally spaced nodes

import numpy as np
import matplotlib.pyplot as plt

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
    plt.title('Runge phenomenon')
    plt.legend(loc='upper center', fontsize=9)

x = np.linspace(-5, 5, 11)
y = runge(x)
runge_fit = polyfit(x, y)

plt.figure(figsize=(5, 3))

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

overlay_plot_fit(runge,runge_fit, x, y, -5,5)
```

It can be shown that the interpolation will diverge at around $3.6$ as $n\to \infty$ and the maximum error $\max_{x\in[-5, 5]} |f_n(x) - h(x) |$ grows exponentially, where $f_n$ is the interpolating polynomial with $n+1$ equally spaced nodes.

### Interpolation Remainder Theory

Let $f_n$ be the degree-$n$ polynomial interpolates $h$ at nodes $\{x_j\}_{j=0}^n$. If $h$ is analytic in a domain $T$ (possibly contains holes),  then the interpolation (Lagrange interpolant) can be written as

$$
    f_n(z) = \sum_{j=0}^n \frac{\omega(z) h(x_j)}{(z - x_j) \omega'(x_j)}
$$

Let $\psi(\xi; z) = \frac{(\omega(\xi) - \omega(z)) h(\xi)}{(\xi - z) \omega(\xi)}$, then by the Residue theorem for simple poles, if $z \neq x_j$,  

$$
    \frac{1}{2\pi i}\int_{\partial T} \psi(\xi; z) d\xi = \sum_{j=0}^n \mathrm{Res}(\psi, x_j) = \sum_{j=0}^n \frac{(\omega(x_j) - \omega(z)) h(x_j)}{(x_j - z)\omega'(x_j)} = f_n(z),
$$

which implies that

$$
    h(z) - f_n(z) = \frac{1}{2\pi i}\int_{\partial T} \frac{\omega(z) h(\xi)}{(\xi - z)\omega(\xi)} d\xi.
$$

The error analysis focuses on studying the behavior of $|\omega(z)|$ as $n\to \infty$.

````{prf:lemma}
:label: Lem-2-Ome-Lim
If $\{x_j\}_{j=0}^n$ are equispaced nodes over $[a, b]$, then 

$$
\lim_{n\to\infty} |\omega(z)|^{\frac{1}{n+1}} = \exp\left(\frac{1}{b - a}\int_a^b \log|z-\xi| d\xi \right).
$$
````

````{prf:proof}
Taking $\log$ on $|\omega|^{\frac{1}{n+1}}$, then 

$$
\log |\omega|^{\frac{1}{n+1}} = \frac{1}{n+1}\sum_{j=0}^n \log |z - x_j|\to \frac{1}{b-a}\int_a^b \log|z - \xi| d\xi. 
$$
````

Let $\sigma_n(z):= |\omega(z)|^{\frac{1}{n+1}}$ and define the contour $C_{\rho} = \{z\in \mathbb{C}\mid \sigma_n(z) = \rho\}$. These level sets are concentric closed curves about the midpoint of $[a,b]$.

````{prf:lemma}
:label: Lem-2-Ana-Uni-Con
Suppose the interpolation nodes $\{x_j\}_{j=0}^n$ are enclosed by $C_{\rho}$ and $h$ is analytic inside $C_{\rho}$. Let $z\in C_{\rho'}$ be such that $\rho'<\rho$, then $f_n\to h$ uniformly as $n\to\infty$. 
````

````{prf:proof}
Using the maximum modulus principle, the analytic function $h - f_n$ must attain its maximum modulus at the boundary $C_{\rho}$, thus 

$$
|h(z) - f_n(z)| = \frac{1}{2\pi}\sup_{z\in C_{\rho'}} \left|\int_{C_{\rho}} \frac{\omega(z)}{\omega(\xi)} \frac{h(\xi)}{(\xi - z)} d\xi\right| \le C(\rho, \rho') \sup_{\xi\in C_{\rho}} \frac{|\omega(z)|}{|\omega(\xi)|},
$$

where $C(\rho, \rho')$ is a positive constant independent of $n$. For $n$ sufficiently large, we can find $0 < \delta < \frac{1}{3}(\rho - \rho')$ sufficiently small such that 

$$
\sup_{\xi\in C_{\rho}} \frac{|\omega(z)|}{|\omega(\xi)|} \le \left|\frac{\rho' + \delta}{\rho - \delta}\right|^{n+1} \to 0\; \text{ as }n\to \infty. 
$$
````

When $h$ is not analytic inside $C_{\rho}$, let us consider a generic situation in which there exist isolated simple poles $z_k\in C_{\rho_k}$, $k\in [m]$ with $\rho_k < \rho$, then we select a contour $C_{\rho'}$ that $\rho_k <\rho'<\rho$ for all $k$. For $z\in C_{\rho'}$, we have

$$
\begin{aligned}
    h(z) - f_n(z) &= \frac{1}{2\pi i} \int_{C_{\rho} - \bigcup_{k=1}^m \Gamma_k} \frac{\omega(z)h(\xi)}{(\xi - z) \omega(\xi)} d\xi \\
    &= \frac{1}{2\pi i} \int_{C_{\rho}} \frac{\omega(z)h(\xi)}{(\xi - z) \omega(\xi)} d\xi - \frac{1}{2\pi i}\sum_{k=1}^m \int_{\Gamma_k}  \frac{\omega(z)h(\xi)}{(\xi - z) \omega(\xi)} d\xi,
\end{aligned}
$$

where $\Gamma_k$ is a small path surrounding $z_k$. The first term can be estimated using the lemma~\ref{Lem: 2-Ana-Uni-Con} whose limit goes to zero as $n\to \infty$. The second term is the summation

$$
\sum_{k=1}^m \frac{\omega(z)}{\omega(z_k)}\frac{\mathrm{Res}(h, z_k)}{z_k - z}.
$$
Since for sufficiently large $n$, we can find $0<\delta<\min_{k\in[m]}\frac{1}{3}(\rho'-\rho_k)$,

$$
\left|\frac{\omega(z)}{\omega(z_k)}\right|^{\frac{1}{n+1}} = \frac{\sigma_n(z)}{\sigma_n(z_k)} > \min_{k\in [m]}\frac{\rho'-\delta}{\rho_k + \delta} > 1.
$$

Define the unique set $\mathcal{U} = \{u_1, u_2,\cdots, u_l\}$ that $u_1 < u_2<\cdots < u_l$ which consists of all distinct values from $\{ \rho_1, \rho_2, \cdots,  \rho_m\}$, then the summation can be decomposed into $l$ groups:

$$
\sum_{k=1}^m \frac{\omega(z)}{\omega(z_k)}\frac{\mathrm{Res}(h, z_k)}{z_k - z} = \sum_{s=1}^l \sum_{\substack{k\in [m] \\ \rho_k = u_s}} \frac{\omega(z)}{\omega(z_k)}\frac{\mathrm{Res}(h, z_k)}{z_k - z}.
$$

As $n\to\infty$, if any of the groups does not vanish, then the whole summation must blow up as $n\to\infty$ (why?). Otherwise, the limiting summation should vanish, which violates the maximum modulus principle.

In Runge's example {eq}``EQ-RUNGE-EXAMPLE``, the simple poles $\pm i$ are on the contour $C_{\rho}$ which intersects the real line at $x_c\approx 3.6334$. Therefore, for $|x| < x_c$, the interpolation $f_n$ uniformly converges to $h$ and diverges once $|x| > x_c$.

```{code-cell} ipython3
:tags: [remove-input]
:mystnb:
:   image:
:       align: "center"
:   figure:
:       caption: Contour curve passing through the poles at $\pm i$.

from pylab import *
from scipy.special import p_roots
import numpy as np
import matplotlib.pyplot as plt

%matplotlib inline

N = 300

def gauss1(f,n):
    [x,w] = p_roots(n+1)
    G=sum(w*f(x))
    return G

def gauss(f,n,a,b):
    [x,w] = p_roots(n+1)
    G=0.5*(b-a)*sum(w*f(0.5*(b-a)*x+0.5*(b+a)))
    return G

def my_f(x):
    return log(x**2+1)/20

v = gauss(my_f,N,-5,5)

def target_f(x, a, b):
    return log((x - a)**2 + b**2) / 20

def d_target_f(x, a, b):
    return 2 * b / ((x - a)**2 + b**2) / 20

a = np.linspace(0, 3.63334, N)
y = np.zeros(len(a))
for j in range(len(a)):
    b0 = 0.1
    while True:
        def target_f2(x):
            return target_f(x, a[j], b0)
        def d_target_f2(x):
            return d_target_f(x, a[j], b0)

        u = gauss(target_f2, N, -5, 5)
        w = gauss(d_target_f2, N, -5, 5)

        delta = (u - v)/w

        if np.abs(delta) > 1e-12:
            b0 = b0 - (u - v)/w
        else:
            break

    y[j] = b0

plt.figure(figsize=(5, 3))

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

plt.axhline(y = 0., color = 'r', linestyle = '--') 
plt.axvline(x = 3.63334, color = 'b', linestyle = '--')
plt.plot(a, y,  label='contour curve')
plt.title('Contour curve passing through (0, 1)')
plt.legend(loc='center left', fontsize=9)
plt.show()
```

There exist better choices of interpolation nodes to prevent such a phenomenon. We will discuss this topic in the next Section.

### Chebyshev Interpolation

The Chebyshev interpolation aims to minimize the bound of the interpolation error. The bound of $\omega(x)$ only depends on the choice of the nodes, so a natural question is: what kind of interpolation nodes will minimize
$\max_{x\in [a, b]} \prod_{j=0}^n |x-x_j|$. We first restrict our analysis to the interval $[a, b] = [-1,1]$ for simplicity, the general case will be discussed later.

```{prf:example}
When $n = 1$, $\omega(x) = (x - x_0)(x - x_1)$, this function changes sign over the sub-intervals $[-1, x_0)$, $(x_0, x_1)$, $(x_1, 1]$, then we can compute the maximum of $|\omega(x)|$ on these sub-intervals. Therefore we need to solve 

$$
    \min_{x_0, x_1\in [-1,1]}\max((1 + x_0)(1 + x_1), \frac{(x_1-x_0)^2}{4}, (1 - x_0)(1 - x_1) ),
$$

while we can observe that 

$$
    \frac{1}{2} (1 + x_0)(1 + x_1) +  \frac{(x_1-x_0)^2}{4} + \frac{1}{2}(1 - x_0)(1 - x_1) = 1 + \frac{(x_0 + x_1)^2}{4}\ge 1
$$

holds for any choice of $x_0, x_1$, which means the maximum is at least $\frac{1}{2}$, it occurs when all terms are equal and $x_0 + x_1 = 0$. Hence $x_0 = -\frac{\sqrt{2}}{2}, x_1 = \frac{\sqrt{2}}{2}$.
```

```{prf:definition}
The Chebyshev polynomials of the first kind are defined by: 

$$
T_k(x) = \cos (k\arccos x),\quad x\in[-1,1].
$$
```

```{prf:theorem}
The Chebyshev polynomial satisfies the following: 

- $T_k(\cos\theta) = \cos k\theta, \quad \theta\in [0, \pi]$.
- $T_0 \equiv 1$, $T_1(x) = x$ and $T_{k+1}(x) = 2 x T_{k}(x) - T_{k-1}(x), \quad k\ge 1.$
- $\max_{x\in[-1,1]} |T_k(x)| = 1$.
- The leading coefficient of $T_k(x)$ is $2^{k-1}$.
- $T_k$ has a total of $(k+1)$ extrema $s_j = \cos(\frac{j\pi}{k}), j = 0, 1,\dots, n$ in the interval $[-1,1]$ such that $T_k(s_j) = (-1)^j$.
```

```{prf:proof}
The first three statements are straightforward after replacing the variable $x = \cos\theta$. The fourth statement is an immediate result with induction through the recursion formula 
$T_{k+1}(x) = 2 x T_k(x) - T_{k-1}$. The last statement is trivial.
```

More importantly, the Chebyshev polynomial has the following optimality property.

````{prf:theorem}
The optimal choice of interpolation nodes that minimize $\max |\omega(x)|$ are the extrema of Chebyshev polynomial $T_{{n+1}}$.

```{math}
:label: EQ-CHEBY
\min_{x_j\in[-1,1]} \max_{x\in[-1,1]} |\omega(x)| =   \max_{x\in[-1,1]} \frac{1}{2^n}|T_{n+1}(x)|  = \frac{1}{2^n}
```
````

```{prf:proof}
Let the roots of $T_{n+1}(x)$ be $z_0, z_1, \dots, z_{n}\in [-1, 1]$, then we can write 

$$T_{n+1} = 2^{n}(x - z_0)(x-z_1)\dots (x - z_{n})$$

therefore $\frac{1}{2^n} T_{n+1}(x)$ is a polynomial with leading coefficient as $1$. Since $\max_{x\in[-1,1]} |T_{n+1}(x)| = 1$, it is clear that $\max_{x\in [-1,1]} \frac{1}{2^n}|T_{n+1}(x)| = \frac{1}{2^n}$, which is the second equality in {eq}`EQ-CHEBY`. For the first equality, we try to prove by contradiction. Let $x_0, x_1, \dots, x_n\in [-1, 1]$, such that 

$$\max_{x\in[-1,1]}|\omega(x)| < \frac{1}{2^n},$$

then we define the polynomial $\psi(x) = \frac{1}{2^n}T_{n+1}(x)- \omega(x)$, its degree is at most $n$ due to cancellation, therefore at most have $n$ zeros. On the other hand, because $\frac{1}{2^n}T_{n+1}(s_j) = \frac{1}{2^n}(-1)^j$ at the extrema $s_j = \cos(\frac{j\pi}{n+1})$, $j=0,\dots, (n+1)$, the polynomial $\psi(s_j)$ must share the same sign of $\frac{1}{2^n}T_{n+1}(s_j)$. This means $\psi(x)$ changes sign $(n+1)$ times, hence $(n+1)$ zeros. It is a contradiction.
```

```{prf:definition}
The interpolation nodes $z_j = \cos(\frac{(2j+1)\pi}{2(n+1)})$, $j = 0, 1, \dots, n$ are called **Chebyshev nodes**. These nodes are the zeros of Chebyshev polynomial $T_{n+1}$.
```

Now we can generalize the above theorem to interval $[a, b]$. One can defined the affine transformation $\phi$ mapping $[-1,1]$ to $[a, b]$ by $\phi(x) = \frac{1}{2} (a + b + (b-a)x)$. It is not difficult to prove the following.

```{prf:corollary}
:label: COR-CHEBY
The optimal choice of interpolation nodes that minimize $\max |\omega(x)|$ on $[a, b]$ are $\phi(z_j)$ and 

$$\renewcommand{\eps}{\varepsilon}
\min_{x_j\in [a, b]} \max_{x\in [a, b]} |\omega(x)| = \frac{(b-a)^{n+1}}{2\cdot 4^n}.
$$

This bound is much smaller than the bound for equally spaced nodes.
```

### Stability of Polynomial Interpolation

Suppose there is some perturbation of the data $\tilde{y}_j = y_j + \eps_j$ at the interpolation node $x_j$. Let $\tilde{f}_n(x)$ and $f_n(x)$ be the interpolating polynomials on perturbed data and original data. Then, with Lagrange polynomials,

$$
    \begin{aligned}
        |f_n(x) - \tilde{f_n}(x)| &= |\sum_{j=0}^n (y_j - \tilde{y}_j) L_j(x)| \\
        &\le \left(\max_{j} |\eps_j|\right) \sum_{j=0}^n |L_j(x)|.
    \end{aligned}
$$

Here $\lambda_n(x) := \sum_{j=0}^n |L_j(x)|$ is the Lebesgue function. It is a piecewise polynomial. Its maximum $\Lambda_n$ is the **Lebesgue constant** and only depends on the choice of interpolation nodes.

For the equally spaced nodes, this Lebesgue constant grows exponentially. In fact,{cite}`turetskii1940bounding` proved the following sharp result.

```{prf:lemma}
Let $\{x_j\}_{j=0}^n$ be equispaced nodes on $[0, 1]$, then the Lebesgue constant 

$$
\Lambda_n  =  \frac{2^{n+1}}{e n \log n} \left(1 + o(1)\right), \quad n\to\infty. 
$$
```

```{prf:proof}
Assume $n\ge 3$, we prove the lower bound by construction. Let $f(x) = e^{i n\pi x}$ and define $p_n(x) = \sum_{j=0}^{n} f(x_j) L_j(x)$ as the interpolation polynomial at nodes $x_j = j \Delta$, $\Delta = \frac{1}{n}$, then 

$$
|p_n(x)| = \left| \sum_{j=0}^{n} (-1)^j L_j(x) \right| \le \Lambda_n.

$$
Let $\delta_{+}$ be the forward difference operator

$$
    \delta_{+} p_n(x) := p_n(x+\frac{1}{n}) - p_n(x),
$$

then $ p_n(x_j) =  (1+\delta_{+})^j p_n(0)$ for $j=0,1,\cdots, n$, which means the interpolation polynomial is

$$
    p_n(x) = \sum_{k=0}^n \binom{nx}{k} \delta_+^k p_n(0) = \sum_{k=0}^n \binom{nx}{k} \delta_+^k f(0),
$$

which is exactly the Newton form at the equispaced nodes.Because $\delta_+ f(x) = (e^{i\pi} - 1)f(x) = -2 f(x)$ acts as a multiplicative operator, therefore 

$$
    p_n(x) = \sum_{k=0}^n \binom{nx}{k}  (-2)^k.
$$

Let $x_{\ast} = \frac{1}{n \log n}$ and set $\mu = n x_{\ast}\in(0, \frac{1}{\log(3)})$. Then, for $k\ge 1$, we have 

$$
    \frac{ \binom{\mu}{k}  (-2)^k  }{ \binom{\mu}{k + 1}  (-2)^{k + 1} } = \frac{k+1}{2(k - \mu)} > \frac{1}{2}.
$$

Therefore, 

$$
\begin{aligned}
        \Lambda_n &\ge  |p_n(x_{\ast})| \ge \binom{\mu}{n}  (-2)^n \left(2 - \frac{1}{2^n}\right) - 1 \\
        &=\exp\left( \log\mu - \log n + \sum_{k=1}^{n-1} \log\left(1 - \frac{\mu}{k}\right) \right) (2^{n+1} - 1) - 1.
\end{aligned}
$$

Denote the positive constant $C = 2(1 - \frac{1}{\log(3)})$ and notice $\log(1 - x) \ge -x - \frac{x^2}{C}$ for all $x\in (0, \frac{1}{\log(3)})$, then 

$$
\begin{aligned}
    \Lambda_n &\ge \frac{\mu}{n} \exp\left(- \mu\sum_{k=1}^{n-1}\frac{1}{k}\right) \exp\left(-\frac{\mu^2}{C} \sum_{k=1}^{n-1}\frac{1}{k^2}\right)(2^{n+1} - 1) - 1 \\
    & \ge \frac{\mu}{n} \exp(-\mu (\log n + \gamma))\exp\left(-\frac{\pi^2\mu^2}{6C}\right)(2^{n+1} - 1) - 1 \\
    & =\frac{2^{n+1}}{e n \log n} \left(1 + o(1)\right).
\end{aligned}
$$

We notice that if $\|f\|_{\infty}\le 1$,  then $\|\delta_{+} f\|_{\infty}\le 2$, using~\eqref{EQ: NEWTON-FORM-EQUI}, the upper bound of $\Lambda_n$ can be estimated by

$$
    \Lambda_n \le \sup_{x\in (0, 1)}\sum_{k=0}^n \left|\binom{nx}{k}\right| 2^k.
$$

By symmetry, it is sufficient to consider $x\le \frac{1}{2}$, otherwise, a similar bound can be derived with the backward difference operator $\delta_{-}$. Hence, 

$$
\begin{aligned}
    \Lambda_n &\le \sum_{k=0}^{\lceil \frac{n}{2} \rceil} \binom{\lceil \frac{n}{2} \rceil}{k} 2^k + \sup_{z\in (0, \frac{n}{2})}  \sum_{k=\lceil \frac{n}{2} \rceil + 1}^{n} \binom{z}{k} 2^k \\
    &= 3^{ \lceil \frac{n}{2} \rceil } + \sum_{k = \lceil\frac{n}{2} \rceil + 1}^n  \frac{2^k}{e k \log k}\left(1 + o(1)\right),
\end{aligned}
$$

where the following inequalities are applied:  

$$
    \begin{aligned}
        &\left|\binom{z}{k}\right|\le \binom{\lceil \frac{n}{2} \rceil}{k},\quad \forall z,k\le \lceil \frac{n}{2} \rceil
    \end{aligned}
$$

and $\forall k \ge \lceil\frac{n}{2} \rceil + 1$, 

$$
    \begin{aligned}
        \sup_{z\in (0, \frac{n}{2})}\left|\binom{z}{k}\right| &=  \sup_{z\in (0, \frac{n}{2})} \frac{\sin(\pi(z - l))}{\pi k!} \Gamma(z+1)\Gamma(k - z) \\
        &=\sup_{z\in (0, 1)} \frac{\sin(\pi(z - l))}{\pi k!} \Gamma(z+1)\Gamma(k - z) \quad (\text{by comparison})\\&\le \sup_{z\in (0, 1)}\frac{z}{k}\left(\frac{z+1}{k - z}\right)^z \quad (\text{by Gautschi's inequality})\\
        &= \frac{1}{e k \log k}\left(1 + o(1)\right).
    \end{aligned}
$$

Therefore,  

$$
    \Lambda_n \le 3^{\lceil \frac{n}{2}\rceil}+  \frac{1}{e \log n} \left( \sum_{k=\lceil \frac{n}{2}\rceil+1}^n \frac{2^k}{k}\right) (1 + o(1)) =  \frac{2^{n+1}}{e n \log n}  (1 + o(1)).
$$
```

For the general case, it has been proved by Paul Erdös (1964) that for any set of interpolation nodes,

$$
    \Lambda_n > \frac{2}{\pi}\log(n+1) + \frac{1}{2},\quad  n\ge 0.
$$

As the number of nodes $n\to \infty$, $\Lambda_n \to \infty$. This leads to the result of Faber that, for any choice of nodes, there exists a continuous function not able to be approximated by the interpolating polynomial. The Chebyshev nodes are almost optimal, in the sense that

$$
    \Lambda_{n, \textrm{Chebyshev}} < \frac{2}{\pi}\log(n+1) + 1.
$$

The set of nodes that minimizes $\Lambda_n$ is difficult to compute. A slightly better set of nodes than Chebyshev nodes are the **extended Chebyshev nodes**:

$$
    \tilde{x}_j = \frac{\cos\left(\frac{2j+1}{2(n+1)\pi}\right)}{\cos\left(\frac{\pi}{2(n+1)}\right)}.
$$

### Newton Form

The Newton form is useful when we dynamically add interpolation nodes. Consider the following scenario: we already have found an interpolation polynomial $f_k$ through $(x_0, y_0)$, $(x_1, y_1)$,$\dots, (x_k, y_k)$, then we are provided an addition pair $(x_{k+1}, y_{k+1})$, how to effectively transform $f_k$ to $f_{k+1}$? If we write

$$
f_{k+1}(x) = f_k(x) + c_{k+1} (x - x_0)(x - x_1)\dots (x - x_{k}),
$$

then $f_{k+1}(x_j) = f_k(x_j)$, $j = 0, 1,\dots, k$, automatically. Therefore we only need to take care of $f_{k+1}(x_{k+1}) = y_{k+1}$, which means

```{math}
:label: EQ-CK
c_{k+1} = \frac{y_{k+1} - f_k(x_{k+1})}{\prod_{j=0}^k (x_{k+1} - x_j)}.
```

Such an inductive procedure produces the Newton form:

```{math}
:label: EQ-NEWTON
f_n(x) = c_0 + c_1 ( x - x_0) + c_2 (x - x_0)(x - x_1)+\dots+c_{n}(x-x_0)\dots (x - x_{n-1}).
```

where the constant $c_j$ depends on $x_0, x_1, \dots, x_{j}$ only. The polynomials $\prod_{j=0}^k (x - x_j)$ are called **Newton polynomials**. When the coefficients $c_k$ are known, the Newton form {eq}`EQ-NEWTON` can be evaluated by the famous **Horner's scheme**, which is

```{math}
:label: EQ-HORNER
    f_n(x) = c_0 + (x-x_0)(c_1 + (x-x_1)(c_2 + (x-x_2)(c_3 + \dots))),
```

the evaluation order starts from the innermost part $c_n (x -x_{n-1})$. This formulation has a complexity of $3n$ flops.

```{prf:remark}
The computation of $c_k$ is not cheap from {eq}`EQ-CK`. A naive algorithm with Horner's scheme roughly takes $5n^2/2+\cO(n)$ flops to compute all coefficients. The **divided differences** is a better way to compute $c_k$.
```

```{prf:definition}
Let the interpolation nodes be $\{x_0, x_1, \dots, x_n\}$, the **divided differences** are defined recursively as follows (the square bracket is used to distinguish from the usual bracket): 

$$
    \begin{aligned}
        f[x_j] &:= f(x_j),\\
    f[x_{j}, \dots, x_{j+k}] &:= \frac{f[x_{j+1},\dots, x_{j+k}] - f[x_j,\dots, x_{j+k-1}]}{x_{j+k} - x_{j}},
    \end{aligned}
$$

where $0\le j, k\le n$ and $j+k\le n$.
```

The following example graph is helpful to understand the relationships among the divided differences.

```{math}
:label: EQ-ALG-NEWTON
\begin{aligned}
    f[x_0] &               &\\ 
            &\searrow       &\\ 
    f[x_1] &\to f[x_0, x_1]& \\
            &\searrow        &\searrow&\\ 
    f[x_2] & \to f[x_1, x_2]&\to& f[x_0, x_1, x_2]\\
            &\searrow        &\searrow& &\searrow&\\ 
    f[x_3] &\to f[x_2, x_3] &\to& f[x_1, x_2, x_3] &\to & f[x_0, x_1, x_2, x_3]\\
            &\searrow        &\searrow & &\searrow&  &\searrow &\\ 
    f[x_4] & \to f[x_3, x_4]&\to& f[x_2, x_3, x_4] &\to & f[x_1, x_2, x_3, x_4] &\to&  f[x_0, x_1, x_2, x_3, x_4] 
\end{aligned}
```

It is clear that computing all of the divided differences requires $\frac{3n^2}{2} +\cO(n)$ flops. The following theorem is the main statement for the Newton form.

````{prf:theorem}
The interpolation polynomial $f_n$ in Newton form is given by 

```{math}
:label: EQ-NEWTON-FORM
\begin{aligned}
    f_n(x) =\; &f[x_0] + f[x_0, x_1](x-x_0) + \dots +\\
             & f[x_0, \dots, x_n](x - x_0)(x - x_1)\dots (x - x_{n-1}).
\end{aligned}
```

In other words, $c_k = f[x_0, \dots, x_k]$.
````

```{prf:proof}
We prove this by induction. Assume the statement is true for $n$ and interpolation node and corresponding values $(x_0, f[x_0]), (x_1, f[x_1]), \dots, (x_n, f[x_n])$. For a new node and value $(x_{n+1}, f[x_{n+1}])$, it is known from {eq}`EQ-CK` that $c_{n+1}$ is the coefficient of leading power. Let $g_n$ be the interpolation polynomial in Newton form through nodes $(x_1, f[x_1]), \dots, (x_{n+1}, f[x_{n+1}])$, then 

$$
\psi(x)  := g_n(x)(x - x_0) - f_n(x)(x - x_{n+1})
$$

satisfies that $\psi(x_j) = f[x_j](x_{n+1} - x_0)$ for $0\le j\le {n+1}$. Therefore 

$$
f_{n+1}(x) = \frac{g_n(x)(x - x_0) - f_n(x)(x - x_{n+1})}{x_{n+1} - x_0}
$$

The leading power's coefficient is then 

$$
    \frac{f[x_1, \dots, x_{n+1}] - f[x_0, \dots, x_n]}{x_{n+1} - x_0} = f[x_0, x_1,\dots, x_{n+1}].
$$
```

```{prf:remark}
The divided difference $f[x_j, \dots, x_{j+k}]$ is the coefficient of leading power of the interpolating polynomial through $(x_j, f[x_j]), \dots, (x_{j+k}, f[x_{j+k}])$. It can be shown that 

$$
f[x_j, \dots, x_{j+k}] = \frac{1}{k!}f^{(k)}(\xi)
$$
    
for some $\xi\in [a, b]$. See exercise.
```

```{prf:remark}
The error estimate can be derived as 

$$
    f(x) - f_n(x) = f[x_0, x_1, \dots, x_n, x] (x-x_0)\dots (x - x_n).
$$
```

```{prf:remark}
The Newton form {eq}`EQ-NEWTON-FORM` actually does not require distinct nodes. The divided difference can be defined as a limit for repeated nodes:

$$
    f[x_0, x_0] = \lim_{x_1 \to x_0} \frac{f[x_1] - f[x_0]}{x_1 - x_0} = f'(x_0).
$$

Moreover, using Taylor expansion, $f[\underbrace{x_0,\dots, x_0}_{(k+1)\,\text{times}}] = \frac{1}{k!}f^{(k)}(x_0)$. However, in such case the divided differences are not possible to be computed if the derivative values are not provided. We will discuss this scenario later in Hermite interpolation polynomial.
```

```{prf:remark}
The algorithm to compute the divided difference can be made more efficient with a single column to store the diagonal elements. $\leadsto$ is representing the number is not changing.

$$
\begin{aligned}
\color{red}{f[x_0]} &    \leadsto  \color{green}{f[x_0]}          & \leadsto& \color{cyan}{f[x_0]} &\leadsto& \color{blue}{f[x_0]} &\leadsto& \color{black}{f[x_0]}\\ 
&\searrow       &\\ 
\color{red}{f[x_1]} &\to \color{green}{f[x_0, x_1]}&  \leadsto& \color{cyan}{f[x_0, x_1]}  &\leadsto& \color{blue}{f[x_0, x_1]} &\leadsto& \color{black}{f[x_0, x_1]}\\
&\searrow        &\searrow&\\ 
\color{red}{f[x_2]} & \to \color{green}{f[x_1, x_2]}&\to& \color{cyan}{f[x_0, x_1, x_2]} &\leadsto& \color{blue}{f[x_0, x_1, x_2]}&\leadsto& \color{black}{f[x_0, x_1, x_2]}\\
&\searrow        &\searrow& &\searrow&\\ 
\color{red}{f[x_3]} &\to \color{green}{f[x_2, x_3] }&\to& \color{cyan}{f[x_1, x_2, x_3]} &\to & \color{blue}{f[x_0, x_1, x_2, x_3]}&\leadsto& \color{black}{f[x_0, x_1, x_2, x_3]}\\
&\searrow        &\searrow & &\searrow&  &\searrow &\\ 
\color{red}{f[x_4]} & \to \color{green}{f[x_3, x_4]}&\to& \color{cyan}{f[x_2, x_3, x_4]} &\to & \color{blue}{f[x_1, x_2, x_3, x_4]} &\to&  f[x_0, x_1, x_2, x_3, x_4] 
\end{aligned}
$$

```

### Hermite Polynomial Interpolation

The Lagrange polynomial interpolation only requires the values of the data function $h$ at each node. It can be generalized when the derivative values of $h$ are also available.

Let the tuple $(h(x_j), h^{(1)}(x_j), \dots, h^{(m_j)}(x_j))$ be the provided derivative values at the interpolation node $x_j$, $j=0,\dots, n$ and $m_j\ge 0$. $N = \sum_{j=0}^n (m_j + 1)$ is the total number of constraints. It can be shown that there exists a unique polynomial $H_{N-1}\in \Pi_{N-1}$ satisfies

$$
H_{N-1}^{(k)}(x_j) = y_j^k:= h^{(k)}(x_j),\quad j=0,\dots, n,\quad 0\le k\le m_j.
$$

This polynomial is called **Hermite interpolation polynomial**. The idea to construct the Hermite interpolation polynomial borrows from the Lagrange polynomials, which is to find basis $L_{jk}$ such that

```{math}
:label: EQ-HERMITE-CONSTRUCTION
    \frac{d^p}{d x^p}L_{jk}(x_l) = \begin{cases}
        1, & l = j, k = p\\
        0, & \text{otherwise}.
    \end{cases}
```

Once these polynomials are obtained, the Hermite interpolation is straightforward:

$$H_{N-1}(x) = \sum_{j=0}^n \sum_{k=0}^{m_j} y_j^k L_{jk}(x).$$

Its uniqueness can be concluded from the linearly independence of the basis $L_{jk}$. However, the construction method in {eq}`EQ-HERMITE-CONSTRUCTION` is not the simplest. It is known that the Newton form {eq}`EQ-NEWTON-FORM` works for repeated nodes as long as the diagram's diagonal {eq}`EQ-ALG-NEWTON` can be filled. Therefore, we can arrange the nodes

$$
\underbrace{x_0,\dots, x_0}_{(m_0+1)\,\text{times}}, \quad \underbrace{x_1,\dots, x_1}_{(m_1+1)\,\text{times}}, \quad \dots,\quad  \underbrace{x_n,\dots, x_n}_{(m_n+1)\,\text{times}}
$$

In this way, all of the necessary divided differences can be computed. See exercise.

```{prf:remark}
The error estimate for Hermite polynomial interpolation will be the same as the Newton form case.
```

## Trigonometric Interpolation

Periodic functions occur in many applications, that is, $f(x + T) = f(x)$, $x\in \mathbb{R}$ for some $T > 0$. For example, a closed planar curve can be parameterized as a periodic function naturally. The polynomial interpolation does not suit periodic functions, this is because polynomials will eventually go to infinity as $x\to\infty$. The most used interpolation for the periodic function is the \emph{trigonometric polynomial interpolation}. In the following, we assume the period $T = 2\pi$ without loss of generality.

### Fourier Series

```{prf:definition}
For $n\ge 0$, we defined $F_n$ the space of trigonometric polynomials 

$$
    F_n := \{f(x) \mid f(x) = \frac{a_0}{2} + \sum_{k=1}^n a_k \cos kx + \sum_{k=1}^n b_k \sin kx,\; a_k, b_k\in\mathbb{R}\}.
$$

The coefficients $a_0,\dots, a_n$, $b_1,\dots, b_n$ can be also chosen as complex numbers. $f\in F_n$ is said to be of degree $n$ if $|a_n| + |b_n| > 0$.
```

The concept of **degree** here can be validated by the addition theorem of trigonometric functions. For instance, if $f_1\in F_k$, $f_2\in F_l$, then $f_1 f_2 \in F_{k+l}$. In the next, we discuss the uniqueness of the interpolation with the trigonometric polynomial.

```{prf:lemma}
A trigonometric polynomial $f\in F_n$ that has more than $2n$ zeros in $[0, 2\pi)$ must vanish identically.
```

````{prf:proof}
Rewrite the trigonometric function in the form of 
```{math}
:label: EQ-COMPLEX
    f_n(x) = \sum_{k=-n}^{n} \gamma_k e^{ik x}. 
```
where $\gamma_0 = \frac{1}{2}a_0$ and $\gamma_{k} = \frac{1}{2}(a_k - ib_k)$ and $\gamma_{-k} = \frac{1}{2}(a_k + i b_k)$, $k=1,\dots, n$. Then substitute $z = e^{ix}$ and set 

$$
    p(z) = \sum_{k = -n}^n \gamma_k z^{n + k}, 
$$

one can rewrite $f_n(x) = z^{-n} p(z)$. If $f_n(x)$ has more than $2n$ zeros, then $p(z)$ has more than $2n$ zeros, which is a contradiction since $p(z)$ is a polynomial of degree $2n$.
````

```{prf:remark}
Since $\sin nx\in F_n$ has $2n$ zeros $\frac{\pi j}{n}$, $j=0,\dots, 2n-1$, it means to uniquely determine a trigonometric polynomial in $F_n$, exactly $2n+1$ values are needed. This is also known as the Nyquist sampling theorem.
```

A direct corollary is the linear independence of the functions $1$, $\cos k x$ and $\sin k x$, $k = 1, \dots n$, these $(2n+1)$ functions form a natural basis for the trigonometric polynomial space $F_n$.

```{prf:corollary}
The functions $1, \cos kx, \sin kx$, $k=1,\dots, n$ are linearly independent on $C([0, 2\pi])$, hence $F_n$ is a $(2n+1)$ dimensional space.
```

To determine the coefficients $a_k, b_k$ from $(2n+1)$ data paris $(x_j, y_j)$, $j=0, \dots, 2n$. We simply follow the idea of Lagrange polynomials by creating the basis polynomial $l_k(x)$ such that

$$
l_k(x_j) = \begin{cases}
    1, &j = k,\\
    0,&\text{otherwise}.
\end{cases}
$$

```{prf:remark}
A natural idea is replace $x - x_j$ in the Lagrange basis by $\sin(x - x_j)$ and produce something like 

$$
\prod_{j=0, j\neq k}^{2n}\frac{\sin(x - x_j)}{\sin(x_k - x_j)},
$$

but $\sin(x - x_j)$ has two roots on $[0, 2\pi)$, therefore we need to rescale it to $[0, \pi)$.
```

```{prf:theorem}
:label: THM-TRIG
Let the basis trigonometric polynomial 

$$
l_k(x) =\prod_{j=0, j\neq k}^{2n}\frac{\sin(\frac{x - x_j}{2})}{\sin(\frac{x_k - x_j}{2})} ,
$$

then the interpolation trigonometric polynomial is 

$$
f_n(x) = \sum_{k=0}^{2n} y_k l_k(x).
$$
```

```{prf:proof}
It remains to show $l_k\in F_n$. This can be seen by splitting $l_k$ into $n$ pairs, each pair takes the form of 

$$
    \sin(\frac{x-x_0}{2})\sin(\frac{x-x_1}{2}) = \frac{1}{2}\cos\left( \frac{x_0 - x_1}{2}\right) - \frac{1}{2}\cos\left(\frac{2x - x_0 - x_1}{2}\right)\in F_1.
$$
```

Computationally, we can reuse the previously known barycentric form but there exist better methods. For simplicity, we consider the equal space nodes in the following (non-uniform nodes could achieve the same complexity though).

$$
    x_j = \frac{2\pi j}{2n + 1}, \quad j = 0, \dots, 2n.
$$

We will try to locate the coefficients $\gamma_k$ in the complex form (see {eq}`EQ-COMPLEX`) from the interpolation conditions.

$$
    f_n(x_j) = y_j = \sum_{k=-n}^n \gamma_k e^{i k x_j}.
$$

Use the property of the functions $e^{ikx_j}$ that

$$
    \sum_{k=0}^{2n} e^{ik x_j} = \begin{cases}
        2n+1, & k = 0\\
        0, &\text{otherwise}.
    \end{cases}
$$

It is not difficult to derive

$$
    \sum_{j = 0}^{2n} y_j e^{-im x_j} = \sum_{k=-n}^{n} \gamma_k \sum_{j=0}^{2n} e^{i (k - m) x_j} = \gamma_{m} (2n+1).  
$$

Therefore, we can compute

```{math}
:label: EQ-GAMMA
    \gamma_m = \frac{1}{2n+1}  \sum_{j = 0}^{2n} y_j e^{-im x_j}.
```

When the coefficients $\gamma_m$ are known, a Horner's scheme can be employed to evaluate the trigonometric polynomial in $\cO(n)$ time complexity.  However, naive computing of all of the coefficients $\gamma_k$ will cost $\cO(n^2)$ flops. The fast Fourier transform can reduce the time complexity to $\cO(n\log n)$.

### Fast Fourier Transform

The discrete Fourier transform $\texttt{DFT}$ of a vector $\mathbf{a} = (a_0, \dots, a_{n-1})$ is to evaluate the following vector:

$$
\texttt{DFT}(\mathbf{a})_{k} := \frac{1}{n}\sum_{j=0}^{n-1} a_je^{-2\pi i jk/n},\quad k = 0,\dots, n-1.
$$

This is the exact formula to compute the coefficients for the trigonometric interpolation polynomial. Such transform is most efficiently calculated through the fast Fourier transform ($\texttt{fft}$).
The fast Fourier transform exploits the symmetry in $e^{2\pi i j/n}$ when $n$ is the power of two using divide-and-conquer. Let $\omega = e^{-2\pi i/n}$ and $c_k$ be

$$
    c_k = \frac{1}{n}\sum_{j=0}^{n-1} y_j\omega^{jk},\quad k = 0,\dots, n-1.
$$

Let $m = n/2\in \mathbb{N}$, then $\omega^n = 1$, $\omega^m = -1$. We can separate $c_k$ into two parts with even $j$ and odd $j$.

```{math}
:label: EQ-CK2
    c_k = \frac{1}{2} A_k + \frac{1}{2} B_k \omega^k,\quad c_{k+m} = \frac{1}{2} A_k - \frac{1}{2} B_k \omega^k 
```

where

```{math}
:label: EQ-AK-BK
\begin{aligned}
    A_k = \frac{1}{m} \sum_{j=0}^{m-1} y_{2j} (\omega^2)^{jk},\quad B_k = \frac{1}{m} \sum_{j=0}^{m-1} y_{2j+1} (\omega^2)^{jk},
\end{aligned}
```

both $A_k$ and $B_k$ are in the same form and similar to $c_k$, but with only half of the terms in summation. This implies a recursive algorithm. Suppose $A_k$ and $B_k$, $0\le k\le m-1$ can be computed with $f(m)$ operations each, then

$$
   f(n) =  f(2m) = 2 f(m) + 4m
$$

The second term includes $2m$ multiplications and $2m$ additions in {eq}`EQ-CK2`. Therefore

$$
    \begin{aligned}
        f(n)& = 2f(\frac{n}{2}) + 2n \\
&= 4f(\frac{n}{4}) + 2 n + 2n \\
&=\dots \\
&= n f(1) + \underbrace{2n + \dots + 2n}_{\log_2 n \text{ times}} = 2n \log_2 n.
    \end{aligned}
$$

since $f(1) = 0$, no computation is needed in this case.  The $\texttt{fft}$ is usually a standard routine in modern scientific computing software.

### Interpolation Error of Trigonometric Polynomial

The $L^2$ error estimate will be discussed at a later point. In this part, we only focus on the $L^{\infty}$ error estimate.

```{prf:theorem}
If $f\in C^2(\mathbb{R})$ is a $2\pi$-period function, then the trigonometric interpolation polynomial with $2n+1$ equally spaced nodes converges uniformly as the number of nodes tends to infinity.
```

````{prf:proof}
Since $f$ is continuously differentiable, its Fourier series converges to $f$ uniformly. Let $f(x), g(x)$ be

$$
f(x) = \sum_{s = -\infty}^{\infty} \gamma_s e^{is x},\quad g(x) = \sum_{s = -n}^n \gamma_s e^{is x},
$$

respectively and denote $h(x) = f(x)- g(x) = \sum_{|s| > n} \gamma_s e^{isx}$ the reminder. Use integration by parts twice,

$$
\gamma_s = \frac{1}{2\pi}\int_{0}^{2\pi} f(x) e^{-isx} dx = -\frac{1}{2\pi}\frac{1}{s^2}\int_{0}^{2\pi} f''(x) e^{-isx} dx \le \frac{1}{2\pi s^2}\|f''\|_{\infty}.
$$

On the other hand, the interpolation polynomial 

$$
\begin{aligned}
    f_n(x)& = \sum_{m=-n}^n \left(\frac{1}{2n+1}\sum_{j=0}^{2n} y_j e^{-im x_j}\right) e^{im x}  \\
        &=\sum_{m=-n}^n \left(\frac{1}{2n+1}\sum_{j=0}^{2n} ( g(x_j) + h(x_j) ) e^{-im x_j}\right) e^{im x} \\
        &= \sum_{m=-n}^n \left(\frac{1}{2n+1}\sum_{j=0}^{2n} ( \sum_{s = -n}^n \gamma_s  e^{is x_j }+ h(x_j) ) e^{-im x_j}\right) e^{im x} \\
        &=g(x)+ \sum_{m=-n}^n \left(\frac{1}{2n+1}\sum_{j=0}^{2n} h(x_j) e^{-im x_j}\right)e^{im x} \\
        &=g(x)+ \frac{1}{2n+1} \sum_{j=0}^{2n} h(x_j) \frac{\sin((2n+1)(x-x_j)/2)}{\sin((x-x_j)/2)}
\end{aligned}
$$

Therefore 

```{math}
:label: EQ-DIFF
        |f - f_n| \le \|h\|_{\infty} + \|h\|_{\infty}  \frac{1}{2n+1} \sum_{j=0}^{2n} \left|\frac{\sin((2n+1)(x-x_j)/2)}{\sin((x-x_j)/2)}\right|.
```

It is simple to derive $\|h\|_{\infty} \le \sum_{|s|>n}|\gamma_s| \le \frac{1}{n\pi}\|f''\|_{\infty}$. We only need to estimate 

$$ 
\sum_{j=0}^{2n} \left|\frac{\sin((2n+1)(x-x_j)/2)}{\sin((x-x_j)/2)}\right|
$$
    
Separate the nodes into two groups: The first group with $|x - x_j| < \frac{2\pi}{2n+1}$, the absolute value is bounded by $(2n+1)$, there are at most $3$ nodes lying in this region, thus the contribution is at most $\cO(n)$. The second group is $\pi\ge |x - x_j| \ge \frac{2\pi}{2n+1}$, while the rest is symmetric, then one can estimate 

$$
\left|\frac{\sin((2n+1)(x-x_j)/2)}{\sin((x-x_j)/2)}\right| \le \frac{\pi}{|x - x_j|}
$$

where the inequality $\sin x \ge \frac{2}{\pi}x$ for $0\le x\le \frac{\pi}{2}$, then this part will be at most $\cO(n\log n)$, the total contribution is bounded by $\cO(n\log n)$. Then {eq}`EQ-DIFF` can be bounded by 

$$
    |f - f_n|\le \|f''\|_{\infty} \cO\left(\frac{\log n}{n}\right)\to 0,\quad \text{as } n\to \infty.
$$

````

The above result can be extended to the case of Hölder continuous function, see the work of Dunham Jackson (1913).

## Spline Interpolation

It has been seen that increasing the number of interpolation nodes will not always help to improve the approximation. The spline interpolation is to conquer this issue by using the piecewise low-degree polynomials.

```{prf:definition}
Let $x_0, \dots, x_n$ be the distinct nodes on $[a, b]$ such that $a = x_0 <\dots < x_n = b$. The piecewise defined function $s_k(x)$ on the interval $[a, b]$ is a spline of degree $k$ to the nodes if  

$$
    s_k|_{[x_j, x_{j+1}]} \in \Pi_k, \quad s_k\in C^{k-1}([a, b]).
$$

The spline function $s_k$ is $(k-1)$-times continuously differentiable and piecewise polynomial of degree $k$.
```

Then the space of splines $s_k$ will be $(n + k)$ dimension: each interval has $(k+1)$ dimensions, each interface imposes $k$ constraints, therefore $n (k+ 1) - (n-1) k = n + k$ dimensions. This shows that to determine a spline on the nodes uniquely, we will require $n+1$ interpolation values and $k-1$ additional constraints. Usual choices are

- periodic splines. $s_k^{(m)}(a) = s_k^{(m)}(b)$ for $m = 0, 1, \dots, k-1$.
- natural splines. $s_k^{(l+j)}(a) = s_k^{(l+j)}(b) = 0$, $j = 0, 1,\dots, l-2$ and $k = 2l-1$ with $l\ge 2$.

In the following, we discuss some useful examples of spline.

### Linear Spline

The linear splines are a special case of splines. It uses piecewise linear polynomials on each subinterval and does not impose any derivative continuity. Let $y_j$ be the interpolation values at nodes $x_j$, respectively. The interpolation has an explicit form:

$$
    s_1(x) = y_{j-1} + \frac{x - x_{j-1}}{x_j - x_{j-1}} y_j
$$

on the interval $[x_{j-1}, x_j]$. It can be represented as a linear combination of the ``hat'' basis function $\theta_j(x)$, defined by

$$
    \theta_j(x) = \begin{cases}
        \frac{x - x_{{j-1}}}{x_j - x_{j-1}}, & x\in [x_{j-1}, x_{j}],\quad 1\le j\le n\\
         \frac{x - x_{j+1}}{x_j - x_{j+1}}, & x\in [x_j, x_{j+1}],\quad 0\le j\le n-1\\
         0, & \text{otherwise}
    \end{cases}
$$

then $s_1(x)$ can be written as

$$s_1(x) = \sum_{j=0}^n \theta_j(x) y_j.$$

The interpolation error can be derived directly from the previous theory for two interpolation nodes. Let $f\in C^2([a, b])$, then on $[x_{j-1}, x_j]$, the interpolation error is

$$
    \frac{1}{2!}f''(\xi) (x - x_{j-1})(x - x_j) \le \frac{1}{8} \|f''\|_{\infty} |x_j - x_{j-1}|^2.
$$

Therefore, the interpolation error on $[a, b]$ is $\frac{1}{8} \|f''\|_{\infty} h^2$, where $h = \max|x_j - x_{j-1}|$. Once $f''$ is not uniformly bounded or even $f'$ is not well-defined somewhere, e.g., $f\in C^{0,\alpha}[a, b]$, the interpolation error will be replaced by the modulus of continuity. On $x\in [x_j, x_{j+1}]$, we have

$$
\begin{aligned}
    |s_1(x) - f(x)| &= \left| \frac{ (x_{j+1} - x) y_j+ (x - x_j) y_{j+1}}{x_{j+1} - x_j} - f(x) \right| \\
    &=  \left| \frac{ (x_{j+1} - x) (y_j - f(x)) + (x - x_j) (y_{j+1}-f(x))}{x_{j+1} - x_j}  \right| \\ &\le \omega(f; |x_{j+1} - x_j|),
\end{aligned}
$$

where $\omega(f;\tau)$ is the modulus of continuity of $f$.
