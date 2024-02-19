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

````{note}
When $f''$ keeps the same sign over $[a, b]$, it is not hard to show that only one side of the bracket is updating. The bracket size will never decrease to zero, which is different from the bisection method. In the following, we use an example to illustrate this. The function $f(x) = x^2 - 1$ over the initial bracket $[0, 2]$. The left endpoint is updating to the root while the right endpoint is fixed at $2$.

```{image} images/doc/root_finding_img_0.png
:alt: false position method
:class: bg-primary 
:width: 600px
:align: center
```
````

It is actually easy to improve the false position method by forcing more weight towards the other endpoint. This is called the **Illinois method**. Once the same side is updated in two consecutive iterations, the Illinois method will use a slightly different formula for $c$.

$$
{c} = \frac{\lambda_b a f(b) - \lambda_a b f(a)}{\lambda_b f(b) - \lambda_a f(a)}.
$$

- Initially, $\lambda_a = \lambda_b = 1$.
- If the previous and current iteration update the left endpoint, then $\lambda_b = \frac{\lambda_b}{2}$, otherwise reset $\lambda_a = 1$.
- If the previous and current iteration update the right endpoint, then $\lambda_a = \frac{\lambda_a}{2}$, otherwise reset $\lambda_b = 1$.



````{note}
We use the previous example to illustrate the difference between the false position method and the Illinois method. At the second iteration, the Illinois method finds the updating is still on left side, so it modifies right endpoint $f(b)$ into $\frac{1}{2} f(b)$ to compute the new $c$, which makes the selected point ${c}$ closer to the right endpoint than the false position method.

```{image} images/doc/root_finding_img_1.png
:alt: Illinois method
:class: bg-primary 
:width: 600px
:align: center
```
````

``````{prf:example}
Let us try the aforementioned methods to find the root of $f(x) = x^3 - 2x^2 - 4$ on the interval $[1, 3]$. The root $x^{\ast}$ can be computed analytically through cubic root formula, which is roughly ``2.5943130163548496``.

Using the previous methods, we obtain the sequence of selection $c_n$ and the error $|c_n - x^{\ast}|$.

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
iter 25 | 2.5943129658699035644531250 | 5.05e-08 
iter 26 | 2.5943129956722259521484375 | 2.07e-08 
iter 27 | 2.5943130105733871459960938 | 5.78e-09 
iter 28 | 2.5943130180239677429199219 | 1.67e-09 
iter 29 | 2.5943130142986774444580078 | 2.06e-09 
iter 30 | 2.5943130161613225936889648 | 1.94e-10 
iter 31 | 2.5943130170926451683044434 | 7.38e-10 
iter 32 | 2.5943130166269838809967041 | 2.72e-10 
iter 33 | 2.5943130163941532373428345 | 3.93e-11 
iter 34 | 2.5943130162777379155158997 | 7.71e-11 
iter 35 | 2.5943130163359455764293671 | 1.89e-11 
iter 36 | 2.5943130163650494068861008 | 1.02e-11 
iter 37 | 2.5943130163504974916577339 | 4.35e-12 
iter 38 | 2.5943130163577734492719173 | 2.92e-12 
iter 39 | 2.5943130163541354704648256 | 7.14e-13 
iter 40 | 2.5943130163559544598683715 | 1.10e-12 
iter 41 | 2.5943130163550449651665986 | 1.95e-13 
iter 42 | 2.5943130163545902178157121 | 2.59e-13 
iter 43 | 2.5943130163548175914911553 | 3.20e-14 
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
iter 12 | 2.5943129966646405470953596 | 1.97e-08 
iter 13 | 2.5943130123438118417311671 | 4.01e-09 
iter 14 | 2.5943130155377716050679737 | 8.17e-10 
iter 15 | 2.5943130161884040418840414 | 1.66e-10 
iter 16 | 2.5943130163209429106530024 | 3.39e-11 
iter 17 | 2.5943130163479417582550468 | 6.91e-12 
iter 18 | 2.5943130163534418031190398 | 1.41e-12 
iter 19 | 2.5943130163545622401954915 | 2.87e-13 
iter 20 | 2.5943130163547900579601446 | 5.95e-14
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
iter  7 | 2.5943130163543197674869134 | 5.30e-13 
iter  8 | 2.5943130163553771438955664 | 5.28e-13 
iter  9 | 2.5943130163548486777358448 | 8.88e-16 
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
