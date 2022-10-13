This directory contains the Matlab codes relative to the numerical experiments of the paper:

A. A. Aguiar, O. P. Ferreira and L. F. Prudente, Inexact gradient projection method with relative error tolerance, Computational Optimization and Applications, 2022.


The codes are designed to solve the problem:

Minimize 1/2 * |AX-B|_F^2
  s.t.   tr(X) = 1
         X >= 0,
         
where A and B are given n x m (m>n) matrices and X is the n x n matrix that we seek to find.

Implemented methods:

- Gradient inexact projection method employing constant step size;
- Gradient exact projection method employing constant step size;
- Gradient inexact projection method employing Armijo step size;
- Gradient exact projection method employing Armijo step size.

Instructions:

—————————————

1) File main.m contains the main program. Modify it to choose the problem data and the method to be used.

2) Go to the Matlab Command Window and type:

main

and see the output in the screen.
