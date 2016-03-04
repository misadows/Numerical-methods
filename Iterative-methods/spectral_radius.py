import random
import numpy as np
from math import sqrt
 
 
def fun_a(n, i, j):
    k = 11
    m = 2
    if i==j:
        return k
    return m/(n-i-j-2+0.5)

def fun_b(n, i, j):
    k = 11
    m = 2
    if i == j:
        return k
    return 1 / (abs(i - j) + m)

def fun_c(n,i,j):
    k = 11
    m = 2
    if i == j:
        return k
    if j>i:
        return ((-1)**(j+1))*m/(j+1)
    if j<i-1:
        return 0
    return m/(i+1)
    
 
 
def create_matrix_a(n, fun):
    A = np.zeros((n,n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            A[i][j] = fun(n,i,j)
    return A
 
 
def get_iteration_matrix(A):
    D = np.diag(A)
    R = A - np.diagflat(D)
    return R / D
 
 
def spectral_radius(A):
    return max(abs(np.linalg.eigvals(A)))
 
 
def get_statistics(N, fun):
    for n in range(2, N):
        A = create_matrix_a(n, fun)
        yield (n, spectral_radius(get_iteration_matrix(A)))
 
 
for x, y in get_statistics(501, fun_c):
    print("{x},{y}".format(x=x, y=y))