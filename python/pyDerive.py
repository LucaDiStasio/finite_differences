# -*- coding: utf-8 -*-

'''
=====================================================================================

Copyright (c) 2017 Université de Lorraine & Luleå tekniska universitet
Author: Luca Di Stasio <luca.distasio@gmail.com>
                       <luca.distasio@ingpec.eu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=====================================================================================

DESCRIPTION

A class to implement explicit and implicit Finite Difference Method for derivatives of every order m over arbitrarily spaced grid

Tested with Python 2.7 Anaconda 2.4.1 (64-bit) distribution
       in Windows 7.

M                       Order of the highest derivative to approximate
N                       Number of grid points
z                       Location at which derivative is approximated
xs                      Vector of N grid points forming the stencil for the finite difference procedure
fs                      Vector of N functional evaluation at the corresponding grid points
weights                 Vector of finite difference weights
derivative              Approximated value of derivative at location z

Pseudo-code from
"Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", Bengt Fornberg.
Mathematics of Computation; Volume 51 (1988), Number 184, Pages 699-706.
see also
Calculation of Weights in Finite Difference Formulas, Bengt Fornberg.
SIAM Reviews; Volume 40 (1998), Number 3, Pages 685-691
for a FORTRAN version of the algorithm.

'''

from os.path import join
import sys
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction

def computeWeights(N,M,xs,x0):
    delta = []              # delta^m_n,n ==> delta[m][n][n] 0<=m<M+1, 0<=n<N+1, 0<=n<N+1
    for m in range(0,M+1):
        matrix = []
        for n in range(0,N+1):
            row = []
            for n in range(0,N+1):
                row.append(0.0)
            matrix.append(row)
        delta.append(matrix)
    delta[0][0][0] = 1.0
    c1 = 1.0
    for n in range(1,N+1):
        c2 = 1.0
        for nu in range(0,n+1):
            c3 = xs[n] - xs[nu]
            if c3!=0:
                c2 *= c3
                for m in range(0,min(n,M)+1):
                    delta[m][n][nu] = ((xs[n]-x0)*delta[m][n-1][nu]-m*delta[m-1][n-1][nu])/c3
        if c2!=0:
            for m in range(0,min(n,M)+1):
                delta[m][n][n] = c1*(m*delta[m-1][n-1][n-1]-(xs[n-1]-x0)*delta[m][n-1][n-1])/c2
        c1 = c2
    return delta


def computeDerivativeAtPoint(M,x0,stencil,fs):
    allWeights = computeWeights(len(stencil)-1,M,stencil,x0)
    weights = allWeights[-1][-1]
    derivative = 0.0
    for v,value in enumerate(fs):
        derivative += weights[v]*value
    return derivative

def ftest(xs):
    ys = []
    dys = []
    d2ys = []
    d3ys = []
    d4ys = []
    for x in xs:
        ys.append(3.0*x**4-7.0*x**3+x**2-8.0)
        dys.append(12.0*x**3-21.0*x**2+2.0*x)
        d2ys.append(36.0*x**2-42.0*x+2)
        d3ys.append(72.0*x-42.0)
        d4ys.append(72.0)
    return ys, dys, d2ys, d3ys, d4ys

def main(argv):

    alpha = [-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0]

    ys, dys, d2ys, d3ys, d4ys = ftest(alpha)

    y0, dy0, d2y0, d3y0, d4y0 = ftest([0.5])

    numDev = computeDerivativeAtPoint(0,0.5,alpha,ys)

    print("{:10.15f}".format(numDev))
    print("{:10.15f}".format(y0[0]))




if __name__ == "__main__":
    main(sys.argv[1:])
