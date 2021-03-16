# heat

Implementation of n-dimensional solver for Laplace's equation with the following formulation:

![Formulation](img/eq.png)

using a second-order scheme in space and a Backward Euler scheme in time. To solve the resulting linear system,
conjugate gradient method is used.

## Examples:

## 1D solution with initial condition of y(t=0) = 10sin($\pi$ x) and diffusion coefficient of 0.1

![Formulation](img/sin.png)

## 2D solution with initial condition of y(t=0) = 10sin($\pi$ x)sin($\pi$ y) and diffusion coefficient of 0.1

![Alt Text](img/Laplace2D.gif)

