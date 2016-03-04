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

def sor(A, b, factor, rho, stop_condition, x=None):
    n = len(A[0])
    if x is None:
        x = np.zeros(n, dtype=np.float64)
 
    k = 0
    while True:
        y = np.copy(x)
        k += 1
        for i in range(n):
            p = 0
            for j in range(n):
                if i != j:
                    p += A[i][j]*x[j]
            x[i] = (1-factor)*x[i] + (factor/A[i][i])*(b[i] - p)
 
        if stop_condition(A, y, x, b, rho):
            break
 
    return x, k

def full_test(N, relaxation_factor = 1.2):
    A = create_matrix_a(N, fun_a)
    x = get_x(N)
    b = np.dot(A, x)

    print('\n')
    print('First condition')
    res, num_it = sor(A, b, relaxation_factor, rho=0.00001, stop_condition=first_stop_condition)
    print('Finished in: ', num_it)
    print('Norm Ax-b: ', np.linalg.norm(np.dot(A, res)-b, np.inf))

    print('\n')
    print('Second condition')
    res, num_it = sor(A, b, relaxation_factor, rho=0.00001, stop_condition=second_stop_condition)
    print('Finished in: ', num_it)
    print('Norm Ax-b: ', np.linalg.norm(np.dot(A, res)-b, np.inf))
    print('\n') 


def small_test(N, relaxation_factor):
    A = create_matrix_a(N, fun_a)
    x = get_x(N)
    b = np.dot(A, x)

    res1, num_it1= sor(A, b, relaxation_factor, rho=0.0001, stop_condition=first_stop_condition)

    res2, num_it2 = sor(A, b, relaxation_factor, rho=0.0001, stop_condition=second_stop_condition)
    print('{}, {}, {}, {}, {}'.format(relaxation_factor, num_it1, num_it2, euclidean_distance(np.dot(A, res1), b), euclidean_distance(np.dot(A, res2), b)))

def small_test_it(N):
    for i in range(40,200):
        small_test(150, i/100)

small_test_it(300)


