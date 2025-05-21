# Path Simulation Using Computational Methods

## About this project

This is an effort in modeling the path of a trajectory using multiple methods including

- [Epicycles and deferents (complex Fourier series)](https://en.wikipedia.org/wiki/Deferent_and_epicycle)
- [Forward kinematics transformation matrix (matrix manipulation, commonly used in robot manipulators)](https://en.wikipedia.org/wiki/Forward_kinematics)

## Background

### Fouries series in continuous time domain

Base form
$$f(t)=a_0+\sum_{n=1}^{\infty}a_ncos(\frac{2\pi nt}{T})+\sum_{n=1}^{\infty}b_nsin(\frac{2\pi nt}{T})$$
where:
$$
a_0=\frac{1}{T}\int_{0}^{T}f(t)dt
$$
$$
a_n=\frac{2}{T}\int_{0}^{T}f(t)cos(\frac{2\pi nt}{T})dt
$$
$$
b_n=\frac{2}{T}\int_{0}^{T}f(t)sin(\frac{2\pi nt}{T})dt
$$

Compact form
$$f(t)=c_0+\sum_{n=1}^{\infty}c_ncos(\frac{2\pi nt}{T}+\phi)$$
where:
$$
c_0=a_0=\frac{1}{T}\int_{0}^{T}f(t)dt
$$
$$
c_n=\sqrt{a_n^2+b_n^2}
$$

Exponential form
$$
f(t)=\sum_{n=-\infty}{\infty}e_0e^{-in\omega_0t}
$$
where:
$$
c_n=\begin{cases}
\frac{a_n-ib_n}{2}, n < 0 \\
a_0, n = 0 \\
\frac{a_n+ib_n}{2}, n > 0 \\
\end{cases}
$$

### Discrete Fourier transform in time domain

$$
X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j \frac{2\pi}{N}kn}
$$

Amplitude: $A_k = \frac{2}{N} |X[k]|$

Phase: $\phi_k = \arg(X[k])$

Angular frequency: $\omega_k = \frac{2\pi k}{T}$

## Carry out the project

### Modeling trajectory based on the Fourier series
