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

For a variety of programming languages, recursive implementation can be made more efficient by using the so called **tail recusion optimization**. This is a compiler feature that allows the recursive program to be executed with the same efficiency as the iterative one. However, ``Python`` does **not** natively support tail recursion optimization.
```

For the first two cases, we can repeat the process with the new interval until certain stop criteria are met. Each iteration reduces the size of the interval by half (gaining one bit each iteration), the total number of iterations required to reduce the interval to a certain size is $\log_2\left(\frac{b - a}{\epsilon}\right)$, where $\epsilon$ is the desired tolerance.

Once the function $f$ has a sign change over the bracket $[a, b]$, the bisection method is guaranteed to converge to a root, but it is not very efficient. It is usually used to obtain a rough estimate of the root, which is then taken as an initial guess for a more efficient method.

### False Position Method

The bisection method only uses $\text{sgn}(f(a))$ and $\text{sgn}(f(b))$ instead of the function values. The **false position method** improves the bisection method by taking the function values into account. Instead of selecting the midpoint $c = \frac{a + b}{2}$, the false position method selects the point $c\in[a, b]$ that lies on the line connecting $(a, f(a))$ and $(b, f(b))$, that is 

$$c = \frac{a f(b) - b f(a)}{f(b) - f(a)}.$$ 

The false position method is also guaranteed to converge to a root if $f(a)f(b) < 0$, and it usually converges faster than the bisection method. The implementation is quite similar to the bisection method. 

```{prf:example}
Let $f(x) = x^3 - 2x^2 - 4$ in the interval $[1, 3]$. Using the bisection method, we obtain the sequence of midpoints $c_n$ and the corresponding function values $f(c_n)$.
```

```{margin}
The bracket methods need to first locate an interval $[a, b]$ such that $f(a)f(b) < 0$. A common approach is to sample a few equidistant points in a large interval and then use the sign of the function values to identify the bracket. This is a simple and robust approach, but it may require a large number of function evaluations.
```

## Iterative Methods

### Newton-Raphson Method

### Secant Method

## Applications in Optimization

## Exercises

### Theoretical Part

### Computational Part
