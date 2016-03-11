from math import exp
import random
import newton_method
import newton_system
import secant_method
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def first_stop_condition(func, x, next_x, rho):
    return abs(next_x - x) < rho


def second_stop_condition(func, x, next_x, rho):
    return abs(func(x)) < rho


def f_b(x):
    return (x-1)*exp(-16*x)+x**14


def df_b(x):
    return exp(-16*x)*(14*exp(16*x)*x**13-16*x+17)


def get_color(x):
    if np.linalg.norm(x-[0.5, 1, -0.5]) < 0.005: return 'r'
    if np.linalg.norm(x-[-1, 1, -1]) < 0.005: return 'b'
    if np.linalg.norm(x-[-1, 1, 1]) < 0.005: return 'g'
    if np.linalg.norm(x-[0.5, 1, 0.5]) < 0.005: return 'y'


condition = second_stop_condition


# print('Newton method')
# for start, x, value, iterations, in newton_method.test_ranges(0.2, 2.5, f_b, df_b, condition, 0.1, 0.0001):
#     print("{iterations},{x},{y},{start}".format(iterations=iterations, x=x, y=value, start=round(start,2)))

# print('\n\nSecant method')
# for start, end, iterations, x, value in secant_method.test_ranges(0.2, 2.5, f_b, condition, 0.1, 0.0001):
#     print("{iterations},{x},{y},{start},{end}".format(iterations=iterations, x=x, y=value, start=round(start,2), end=round(end,2)))


# for i_rho in range(0,10):
#     yn1, _, i_n = newton_method.solve_newton(f_b, df_b, 0.2, first_stop_condition, 10**(-i_rho))
#     i_s, ys1, _ = secant_method.solve_secant(f_b, 0.2, 2.5, first_stop_condition, 10**(-i_rho))
#     yn2, _, i_n2 = newton_method.solve_newton(f_b, df_b, 0.2, second_stop_condition, 10**(-i_rho))
#     i_s2, ys2, _ = secant_method.solve_secant(f_b, 0.2, 2.5, second_stop_condition, 10**(-i_rho))
#     print("{rho},{i_n1},{i_s1},{i_n2},{i_s2},{pre_n1},{pre_s1},{pre_n2},{pre_s2}".format(i_n1=i_n, rho=10**(-i_rho), i_s1=i_s, i_n2=i_n2, i_s2=i_s2, pre_n1=abs(yn1), pre_n2=abs(yn2), pre_s1=abs(ys1), pre_s2=abs(ys2)))


dF = np.array([[lambda x: 2*x, lambda x: 2*x, lambda x: -2*x],
              [lambda x: 1, lambda x: -6*x**2, lambda x: 4*x],
              [lambda x: 4*x, lambda x: 1, lambda x: -4*x]]
              )

F = np.array([lambda x1, x2, x3: x1**2+x2**2-x3**2-1,
              lambda x1, x2, x3: x1-2*x2**3+2*x3**2+1,
              lambda x1, x2, x3: 2*x1**2+x2-2*x3**2-1
              ])
xs=[]
xs_start=[]
for i in range(1,500):
    x,start_x,norm,i=newton_system.solve_equation_system(F, dF, np.array([random.choice([-1,1]),random.choice([-1,1]),random.choice([-1,1])])*np.random.rand(3)*100, newton_system.second_stop_condition)

    #if np.all(x) and np.linalg.norm(x-[0.5, 1, -0.5]) < 0.005:
    xs.append(x)
    xs_start.append(start_x)

xs = np.array(xs)
xs_start=np.array(xs_start)
xs_start_ok = xs_start[xs != np.array(None)]
xs_start_nok = xs_start[xs == np.array(None)]

xs_ok = xs[xs != np.array(None)]
xs_ok = xs_ok.tolist()
xs_ok = np.array(xs_ok)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('x3')
ax.set_xlim3d(-100, 100)
ax.set_ylim3d(-100,100)
ax.set_zlim3d(-100,100)

#colors = ['r' if x is None else 'b' for x in xs]
colors = [get_color(x) for x in xs_ok]
ax.scatter(xs_start_ok[:,0],xs_start_ok[:,1],xs_start_ok[:,2], s=40, c=colors)
plt.show()


