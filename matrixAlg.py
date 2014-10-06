#!/usr/bin/python
#####################################
# Written by Gavin Heverly-Coulson
# Email: gavin <at> quantumgeranium.com
#####################################
# A set of matrix algebra functions for performing
# basic matrix algebra operations. 
#
# Tested with Python 2.6/2.7
#
# This work is licensed under a Simplified BSD License
# Copyright (c) 2014, Gavin Heverly-Coulson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import math

# Print a nicely formatted matrix
def printMat(mat):
  newStr = ""
  for i in range(len(mat)):
    for j in range(len(mat[0])):
      newStr = newStr + str(mat[i][j]) + "  "
    newStr += "\n"
  print newStr


# Calculates the determinant of a 3x3 matrix, using the 2x2 sub-matrices method
def det3(mat):
  return ( ( mat[0][0]*det2([[mat[1][1], mat[1][2]], [mat[2][1], mat[2][2]]]) ) - ( mat[0][1]*det2([[mat[1][0], mat[1][2]], [mat[2][0], mat[2][2]]]) ) + (mat[0][2]*det2([[mat[1][0], mat[1][1]], [mat[2][0], mat[2][1]]])) )



# Calculates the determinant of a 2x2 matrix
def det2(mat):
  return ((mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0]))


# Calculates the transpose of a matrix
# Works for arbitrary NxM size
def transpose(mat):
  cols = len(mat) # number of rows in mat
  rows = len(mat[0])    # number of columns in mat
  transMat = [x[:] for x in [[None]*cols]*rows] # cols, rows
  
  for a in range(rows):
    for b in range(cols):
      transMat[a][b] = mat[b][a]

  return transMat


# Calculates the dot product of two vectors, A and B
def dotProduct(A, B):
  counter = 0
  product = 0
  while counter < len(A):
    product = product + (A[counter] * B[counter])
    counter += 1
  return product


# Calculates the length of a vector
def vectLength(A):
  sumSquares = 0
  for i in A:
    sumSquares = sumSquares + (i**2)
  return math.sqrt(sumSquares)


# Multiplies two matrices (A and B) and returns the result
def matMult(A, B):
  if len(A[0]) != len(B):
    print "Matrix dimensions don't match!\nA has {0} columns and B has {1} rows.".format(len(A[0]), len(B))
  else:
    newMat = [[0.0 for cols in range(len(B[0]))] for rows in range(len(A))]
    
    for i in range(len(A)):
      for j in range(len(B[0])):
        for k in range(len(B)):
          newMat[i][j] += A[i][k]*B[k][j]
    return newMat


# Converts a given matrix (not necessarily square) to
# reduced row echelon form
def toRedRowEchelon(mat):
  colPos = 0
  rows = len(mat)
  cols = len(mat[0])
  for r in range(rows):
    if colPos >= cols:
      return mat
    i = r
    while mat[i][colPos] == 0.0:
      i += 1
      if i == rows:
        i = r
        colPos += 1
        if colPos == cols:
          return mat
    mat[i], mat[r] = mat[r], mat[i] # swap rows i and r
    lv = mat[r][colPos]
    mat[r] = [mrx / lv for mrx in mat[r]]
    for i in range(rows):
      if i != r:
        lv = mat[i][colPos]
        mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
    colPos += 1
  return mat


# Finds the inverse of a given matrix
def invMat(mat):
  matDim = len(mat)
  idenMat = [[0.0 for col in range(matDim)] for row in range(matDim)]
  for i in range(matDim):
    idenMat[i][i] = 1.0
  
  newMat = [None] * matDim
  for i in range(matDim):
    newMat[i] = mat[i] + idenMat[i]
  
  solvedMat = toRedRowEchelon(newMat)
  invertMat = [None] * matDim
  for i in range(matDim):
    invertMat[i] = solvedMat[i][-1*matDim:]
  return invertMat
