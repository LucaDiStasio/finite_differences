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
            if c3!=0:
                c2 *= c3
                for m in range(0,min(n,M)+1):
                    delta[m][n][nu] = ((xs[n]-x0)*delta[m][n-1][nu]-m*delta[m-1][n-1][nu])/c3
        if c2!=0:
            for m in range(0,min(n,M)+1):
                delta[m][n][n] = c1*(m*delta[m-1][n-1][n-1]-(xs[n-1]-x0)*delta[m][n-1][n-1])/c2
        c1 = c2
    return delta


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

#def main(argv):

M = 4

alpha = [0,1,2,3,4,5,6,7,8]
alpha1 = [-1.0,0.0,1.0]
alpha2 = [-2.0,-1.0,0.0,1.0,2.0]
alpha3 = [-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0]
alpha4 = [-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0]

weight = computeWeights(len(alpha)-1,4,alpha,0.0)
weight1 = computeWeights(len(alpha1)-1,4,alpha1,0.0)
weight2 = computeWeights(len(alpha2)-1,4,alpha2,0.0)
weight3 = computeWeights(len(alpha3)-1,4,alpha3,0.0)
weight4 = computeWeights(len(alpha4)-1,4,alpha4,0.0)

for m,matrix in enumerate(weight):
    for n,row in enumerate(matrix):
        for nu,value in enumerate(row):
            try:
                weight[m][n][nu] = str(Fraction(value).limit_denominator())
            except Exception:
                weight[m][n][nu] = str(value)
                sys.exc_clear()

for m,matrix in enumerate(weight1):
    for n,row in enumerate(matrix):
        for nu,value in enumerate(row):
            try:
                weight1[m][n][nu] = str(Fraction(value).limit_denominator())
            except Exception:
                weight1[m][n][nu] = str(value)
                sys.exc_clear()
for m,matrix in enumerate(weight2):
    for n,row in enumerate(matrix):
        for nu,value in enumerate(row):
            try:
                weight2[m][n][nu] = str(Fraction(value).limit_denominator())
            except Exception:
                weight2[m][n][nu] = str(value)
                sys.exc_clear()
for m,matrix in enumerate(weight3):
    for n,row in enumerate(matrix):
        for nu,value in enumerate(row):
            try:
                weight3[m][n][nu] = str(Fraction(value).limit_denominator())
            except Exception:
                weight3[m][n][nu] = str(value)
                sys.exc_clear()
for m,matrix in enumerate(weight4):
    for n,row in enumerate(matrix):
        for nu,value in enumerate(row):
            try:
                weight4[m][n][nu] = str(Fraction(value).limit_denominator())
            except Exception:
                weight4[m][n][nu] = str(value)
                sys.exc_clear()

for wset in weight:
    print(str(wset))

#print(str(weight1[1][-1])+'\n'+str(weight2[1][-1])+'\n'+str(weight3[1][-1])+'\n'+str(weight4[1][-1]))

#print(str(weight1[2][-1])+'\n'+str(weight2[2][-1])+'\n'+str(weight3[2][-1])+'\n'+str(weight4[2][-1]))

#print(str(weight1[3][-1])+'\n'+str(weight2[3][-1])+'\n'+str(weight3[3][-1])+'\n'+str(weight4[3][-1]))

#print(str(weight1[4][-1])+'\n'+str(weight2[4][-1])+'\n'+str(weight3[4][-1])+'\n'+str(weight4[4][-1]))

#if __name__ == "__main__":
#    main(sys.argv[1:])
