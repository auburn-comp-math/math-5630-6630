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

From the previous discussion, we already know that interpolation provides an estimate **within** the original observation range. The extrapolation is similar but aims to produce estimates outside the observation range. However, sometimes extrapolation may be subject to a greater uncertainty, see {numref}`extrapolation`, one should use it only when an overestimate is hardly occurring.

```{figure} images/doc/extrapolate.png
---
name: extrapolation
scale: 120%
align: center
---
Extrapolation behavior for Chebyshev interpolation with $15$ nodes.
```

### Richardson Extrapolation

Suppose there is a sequence of estimates $A(h)$ depending on the parameter $h$ smoothly, the limit $A^{\ast} = \lim_{h\to 0^{+}} A(h)$ is the quantity to be computed. In practice, we only have access to $A(h)$ for a few values of $h$. Using these values to estimate $A^{\ast}$ is a typical problem in extrapolation.

The basic idea behind \emph{Richardson extrapolation} is to use polynomial interpolation with a sequence of nodes $h_j\to 0$. Suppose that the function $A(h)$ admits the following asymptotic expansion:

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
    p_{n+1}(x) \equiv \tilde{f}_n(x) + p[x, x_0, x_1, \dots, x_n] \prod_{j=0}(x - h_j^{\gamma}), 
$$

where $p[x, x_0, x_1, \dots, x_n]$ is the coefficient of the leading power in $p_{n+1}$, $a_{n+1}$. Thus, 

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

However, such a process can constantly refine the approximation due to the potentially fast-growing constant in the $\cO$ notation.

### Wynn's Epsilon Method

Wynn's $\eps$ method {cite}`wynn1966convergence` is another kind of extrapolation algorithm that is recommended as the best all-purpose acceleration method. It has a strong connection with Padé approximation and continued fractions. We will not cover the detailed derivation of the theory in this section. However, Wynn's $\eps$ method still has its limitations if the sequence converges to the desired value too slowly.

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

- Richardson Extrapolation. We take $1/h = 250, 500, 1000, 2000$ and calculate that the Richardson extrapolation three times would result in almost machine precision.
- Wynn's $\eps$ method. We take the sequence 

  $$
  s_k = \sum_{j=0}^k \frac{(-1)^j}{2 j + 1}
  $$
  
  as the truncated series at $z = 1$. With about 20 terms, we already reached machine precision.
```

## Newton-Cotes Quadrature

## Romberg Integration

## Gaussian Quadrature

## Monte Carlo Integration

## Extended Reading

```{bibliography}
:filter: docname in docnames
```
