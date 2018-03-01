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

def computeWeigths(N,M,xs,x0):
    delta = []              # delta^m_n,n ==> delta[m][n][n] 0<=m<M+1, 0<=n<N+1, 0<=n<N+1
    for m in range(0,M+1):
        matrix = []
        for n in range(0,N+1):
            row = []
            for n in range(0,N+1):
                row.append(0.0)
            matrix.append(row)
        delta.append(matrix)
    mn = 0 # !!!! check this !!!!
    delta[0][0][0] = 1.0
    c1 = 1.0
    #c4 = xs[0] - z;
    for n in range(1,N+1):
        c2 = 1.0
        #c5 = c4
        #c4 = xs[i] - z
        #mn = min(i,M)
        for nu in range(0,n+1):
            c3 = xs[n] - xs[nu]
            c2 *= c3
            for m in range(0,min(n,M)+1):
                delta[m][n][nu] = ((xs[n]-x0)*delta[m][n-1][nu]-m*delta[m-1][n-1][nu])/c3
        for m in range(0,min(n,M)+1):
            delta[m][n][n] = c1*(m*delta[m-1][n-1][n-1]-(xs[n-1]-x0)*delta[m][n-1][n-1])/c2
        c1 = c2


def computeDerivative(fs,weights):
    derivative = 0.0
    for i,f in enumerate(fs):
        derivative += weights[i]*f
    return derivative

def numDerive(N,M,xs,fs):
    dfs = []
    for n in range(0,N+1):
        df = []
        for x in xs:
            df.append(computeDerivative(fs,computeWeigths(N,M,xs,x)))
        dfs.append(df)
    return dfs

def main(argv):

    M = 4
    N = 8
    
    alphas = [-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0]

    computeWeigths(N,4,alphas,0.0)



if __name__ == "__main__":
    main(sys.argv[1:])
