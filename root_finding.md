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

For the first two cases, we can repeat the process with the new interval until certain stop criteria are met. Each iteration reduces the size of the interval by half (gaining one bit each iteration), the total number of iterations required to reduce the interval to a certain size is $\log_2\left(\frac{b - a}{\epsilon}\right)$, where $\epsilon$ is the desired tolerance.

Once the function $f$ has a sign change over the bracket $[a, b]$, the bisection method is guaranteed to converge to a root, but it is not very efficient. It is usually used to obtain a rough estimate of the root, which is then taken as an initial guess for a more efficient method.

### False Position Method

The bisection method only uses $\text{sgn}(f(a))$ and $\text{sgn}(f(b))$ instead of the function values. The **false position method** (*Regula falsi* in Latin) improves the bisection method by taking the function values into account. Instead of selecting the midpoint $c = \frac{a + b}{2}$, the false position method selects the point $c\in[a, b]$ that lies on the line connecting $(a, f(a))$ and $(b, f(b))$, that is

$$c = \frac{a f(b) - b f(a)}{f(b) - f(a)}.$$

The false position method is also guaranteed to converge to a root if $f(a)f(b) < 0$ and the implementation is quite similar to the bisection method. It usually converges faster than the bisection method, but sometimes exceptions occur.

```{note}
Let's consider a convex function $f(x)$, 
```

``````{prf:example}
Let us try the aforementioned methods to find the root of $f(x) = x^3 - 2x^2 - 4$ on the interval $[1, 3]$. The root $x^{\ast}$ can be computed analytically through cubic root formula, which is roughly ``2.5943130163548496``.

Using the previous methods, we obtain the sequence of selection $c_n$ and the error $|c_n - x^{\ast}|$.

`````{tab-set}
````{tab-item} Bisection Method
```{code-block}
iter  1 | 2.0000000000000000000000000 | 4.00e+00
iter  2 | 2.5000000000000000000000000 | 8.75e-01
iter  3 | 2.7500000000000000000000000 | 1.67e+00
iter  4 | 2.6250000000000000000000000 | 3.07e-01
iter  5 | 2.5625000000000000000000000 | 3.06e-01
iter  6 | 2.5937500000000000000000000 | 5.52e-03
iter  7 | 2.6093750000000000000000000 | 1.49e-01
iter  8 | 2.6015625000000000000000000 | 7.15e-02
iter  9 | 2.5976562500000000000000000 | 3.29e-02
iter 10 | 2.5957031250000000000000000 | 1.37e-02
iter 11 | 2.5947265625000000000000000 | 4.06e-03
iter 12 | 2.5942382812500000000000000 | 7.33e-04
iter 13 | 2.5944824218750000000000000 | 1.66e-03
iter 14 | 2.5943603515625000000000000 | 4.65e-04
iter 15 | 2.5942993164062500000000000 | 1.34e-04
iter 16 | 2.5943298339843750000000000 | 1.65e-04
iter 17 | 2.5943145751953125000000000 | 1.53e-05
iter 18 | 2.5943069458007812500000000 | 5.96e-05
iter 19 | 2.5943107604980468750000000 | 2.21e-05
iter 20 | 2.5943126678466796875000000 | 3.42e-06
iter 21 | 2.5943136215209960937500000 | 5.94e-06
iter 22 | 2.5943131446838378906250000 | 1.26e-06
iter 23 | 2.5943129062652587890625000 | 1.08e-06
iter 24 | 2.5943130254745483398437500 | 8.95e-08
iter 25 | 2.5943129658699035644531250 | 4.95e-07
iter 26 | 2.5943129956722259521484375 | 2.03e-07
iter 27 | 2.5943130105733871459960938 | 5.67e-08
iter 28 | 2.5943130180239677429199219 | 1.64e-08
iter 29 | 2.5943130142986774444580078 | 2.02e-08
iter 30 | 2.5943130161613225936889648 | 1.90e-09
iter 31 | 2.5943130170926451683044434 | 7.24e-09
iter 32 | 2.5943130166269838809967041 | 2.67e-09
iter 33 | 2.5943130163941532373428345 | 3.86e-10
iter 34 | 2.5943130162777379155158997 | 7.57e-10
iter 35 | 2.5943130163359455764293671 | 1.86e-10
iter 36 | 2.5943130163650494068861008 | 1.00e-10
iter 37 | 2.5943130163504974916577339 | 4.27e-11
iter 38 | 2.5943130163577734492719173 | 2.87e-11
iter 39 | 2.5943130163541354704648256 | 7.00e-12
iter 40 | 2.5943130163559544598683715 | 1.09e-11
```
````

````{tab-item} False Position Method
```{code-block}
iter  1 | 2.0000000000000000000000000 | 4.00e+00
iter  2 | 2.4444444444444446418174266 | 1.34e+00
iter  3 | 2.5621621621621617492792211 | 3.10e-01
iter  4 | 2.5876913365185605364615640 | 6.47e-02
iter  5 | 2.5929610854818996301673906 | 1.33e-02
iter  6 | 2.5940374914642010395482430 | 2.70e-03
iter  7 | 2.5942568846837747997824408 | 5.51e-04
iter  8 | 2.5943015817106331866170876 | 1.12e-04
iter  9 | 2.5943106870264029950590157 | 2.29e-05
iter 10 | 2.5943125418534931370118102 | 4.66e-06
iter 11 | 2.5943129196954899384763849 | 9.49e-07
iter 12 | 2.5943129966646405470953596 | 1.93e-07
iter 13 | 2.5943130123438118417311671 | 3.94e-08
iter 14 | 2.5943130155377716050679737 | 8.02e-09
iter 15 | 2.5943130161884040418840414 | 1.63e-09
iter 16 | 2.5943130163209429106530024 | 3.33e-10
iter 17 | 2.5943130163479417582550468 | 6.78e-11
iter 18 | 2.5943130163534418031190398 | 1.38e-11
iter 19 | 2.5943130163545622401954915 | 2.81e-12
iter 20 | 2.5943130163547900579601446 | 5.77e-13
iter 21 | 2.5943130163548366873271789 | 1.19e-13
iter 22 | 2.5943130163548464572897956 | 2.31e-14
iter 23 | 2.5943130163548482336466350 | 5.33e-15
iter 24 | 2.5943130163548486777358448 | 0.00e+00
```
````
`````
``````

```{prf:remark}
```

```{prf:remark}
The bracket methods need to first locate an interval $[a, b]$ such that $f(a)f(b) < 0$. A common approach is to sample a few equidistant points in a large interval and then use the sign of the function values to identify the bracket. This is a simple and robust approach, but it may require a large number of function evaluations.
```

## Iterative Methods

### Newton-Raphson Method

### Secant Method

## Applications in Optimization

## Exercises

### Theoretical Part

### Computational Part
