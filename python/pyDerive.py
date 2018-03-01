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
        c5 = c4
        c4 = xs[i] - z
        mn = min(i,M)
        for j in range(0,i):
            c3 = xs[i] - xs[j]
            c2 = c2*c3
            if j==i-1:
                for k in range(mn,0,-1):
                    delta[i][k] = c1*(k*delta[i-1][k-1]-c5*delta[i-1][k])/c2
                delta[i][0] = -c1*c5*delta[i-1][0]/c2
            for k in range(mn,0,-1):
                delta[j][k] = (c4*delta[j][k]-k*delta[j][k-1])/c3
            delta[j][0] = c4*delta[j][0]/c3
        c1 = c2
    weights = []
    for i in range(0,N):
        weights.append(delta[i][M])
    return weights

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
    xtrain = np.arange(1,10,0.1)
    ytrain = xtrain**4 - xtrain**3 + 7

    stencilSize = 3

    devorderApprox = 1

    yapprox = numDerive(stencilSize,devorderApprox,xtrain,ytrain)

    plt.plot(xtrain,ytrain,'-k')
    plt.plot(xtrain,yapprox[0],'.r')
    plt.grid
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
