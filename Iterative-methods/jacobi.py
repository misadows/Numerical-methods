import numpy as np
from random import choice

def get_x(N):
    x = np.zeros(N)
    for i in range(N):
        x[i] = choice([-1,1])
    return x

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

def euclidean_distance(a, b):
    return np.sqrt(np.sum((a-b)**2))

def create_matrix_a(n, fun):
    a = np.zeros((n,n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            a[i][j] = fun(n,i,j)
    return a

def first_stop_condition(A, x1, x2, b, rho):
    norm = euclidean_distance(x1, x2)
    return norm<rho

def second_stop_condition(A, x1, x2, b, rho):
    p = np.dot(A, x1)
    norm = euclidean_distance(p, b)
    return norm<rho

def jacobi(A, b, N, rho, stop_condition, x_old=None):
    if x_old==None:
        x_old = np.ones(len(A[0]))
                                                                                                                                                             
    D = np.diag(A)
    R = A - np.diagflat(D)
    i=0

    while True:
        i+=1
        x = (b - np.dot(R, x_old)) / D
        if stop_condition(A, x, x_old, b, rho):
            break    
        x_old = x
    return i,x

def full_test(N):
    A = create_matrix_a(N, fun_a)
    x = get_x(N)
    b = np.dot(A, x)

    print('\n')
    print('First condition')
    num_it, res = jacobi(A, b, N, rho=0.00001, stop_condition=first_stop_condition)
    print('Finished in: ', num_it)
    print('Norm Ax-b: ', np.linalg.norm(np.dot(A, res)-b, np.inf))

    print('\n')
    print('Second condition')
    num_it, res = jacobi(A, b, N, rho=0.00001, stop_condition=second_stop_condition)
    print('Finished in: ', num_it)
    print('Norm Ax-b: ', np.linalg.norm(np.dot(A, res)-b, np.inf))
    print('\n') 

def small_test(N):
    A = create_matrix_a(N, fun_a)
    x = get_x(N)
    b = np.dot(A, x)

    num_it1, res1 = jacobi(A, b, N, rho=0.0001, stop_condition=first_stop_condition)

    num_it2, res2 = jacobi(A, b, N, rho=0.0001, stop_condition=second_stop_condition)
    print('{}, {}, {}, {}, {}'.format(N, num_it1, num_it2, euclidean_distance(np.dot(A, res1), b), euclidean_distance(np.dot(A, res2), b)))

def small_test_it(N):
    for i in range(2,N):
        small_test(i)

small_test_it(300)





