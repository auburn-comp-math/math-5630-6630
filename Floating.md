# Floating Point Arithmetic

In this chapter, we will introduce some basics on the real number system for modern computers and discuss the arithmetic operations of the number system.

## Representation of Real Numbers

Any nonzero real number $x\in \mathbb{R}$ can be accurately represented with an infinite sequence of digits. This can be understood as the consequence that rational numbers are dense on any interval. ```{margin}
**Dense** means that between any two distinct real numbers, there is a rational number.
```Therefore, with the binary representation, we can write 

$$
x = \pm (0.d_1 d_2 d_3\dots d_{t-1} d_t d_{t+1}\dots) \times 2^e,
$$

where $e$ is an integer exponent and $d_1=1$, the other binary digits $d_i\in \{0, 1\}$. The mantissa part

$$
0.d_1 d_2 d_3\dots =  \frac{d_1}{2} + \frac{d_2}{2^2} + \frac{d_3}{2^3} + \cdots.
$$

