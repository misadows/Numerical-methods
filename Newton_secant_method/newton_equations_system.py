import numpy as np


apply_vectorized = np.vectorize(lambda f, x: f(x))
apply_vectorized_multi = np.vectorize(lambda f, x: f(*x), excluded=[1])


def first_stop_condition(F, x, d, rho=0.001):
    return np.linalg.norm(d) < rho

def second_stop_condition(F, x, d, rho=0.001):
    return np.linalg.norm(apply_vectorized_multi(F, x)) < rho

def solve_equation_system(F, dF, x_start, stop_condition):
    x=x_start
    i=0
    while True:
        i+=1
        try:
            d = np.linalg.solve(apply_vectorized(dF, x), -1*apply_vectorized_multi(F, x))
        except np.linalg.linalg.LinAlgError:
            return None, x_start, np.linalg.norm(d), i

        if stop_condition(F, x, d):
            break
        x = x+d

        if i > 100:
            return None, x_start, np.linalg.norm(d), i

    return x,x_start, np.linalg.norm(d), i