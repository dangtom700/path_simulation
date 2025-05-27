# Path Simulation Using Computational Methods

## About this project

This is an effort in modeling the path of a trajectory using multiple methods including

- [Epicycles and deferents (complex Fourier series)](https://en.wikipedia.org/wiki/Deferent_and_epicycle)
- [Forward kinematics transformation matrix (matrix manipulation, commonly used in robot manipulators)](https://en.wikipedia.org/wiki/Forward_kinematics)
- [Closed loop 4 bar linkage system](https://en.wikipedia.org/wiki/Four-bar_linkage)

## Modeling path with Fourier series

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

## Modeling path in close loop with 4 bar linkage system

Assuming a 4-side geometry that consist of 4 bar of different dimensions and they form 4 angles so that the total inner angles is 360 degree.

4 bar linkage problems often come in 2 types:

- close loop
- open loop

In a close loop system, there are a unique constraint of the sum of 4 inner angles must add up to 360 degrees

Formula used to compute the location of this system (also applied to open system):

- $\frac{BD}{sinBCD}=\frac{BC}{sinBDC}=\frac{CD}{sinCBD}$
- $\frac{BA}{sinBDA}=\frac{BD}{sinBAD}=\frac{AD}{sinABD}$
- $\frac{AB}{sinACB}=\frac{BC}{sinBAC}=\frac{AC}{sinABC}$
- $\frac{AD}{sinACD}=\frac{AC}{sinADC}=\frac{CD}{sinCAD}$
- $BD=\sqrt{AB^2+AD^2-2\times AB\times ADcosBAD}=\sqrt{BC^2+CD^2-2\times BC\times CDcosBCD}$
- $AC=\sqrt{AB^2+BC^2-2\times AB\times BCcosABC}=\sqrt{AD^2+CD^2-2\times AD\times CDcosADC}$

The known parameters are the previous positions of each linkage and their angles with one known angle change that drive the rest of the parameters.
