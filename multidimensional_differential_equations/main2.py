from __future__ import division
import math
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def euclidean_error(u1, u2, n, m):
    result = 0
    for i in range(0, m):
        for j in range(0, n):
            # result += math.sqrt((u1[j, i] - u2[j, i]) * (u1[j, i] - u2[j, i]))
            if abs((u1[j, i] - u2[j, i])) > result:
                result = abs((u1[j, i] - u2[j, i]))
    return result
def calculate_next_value(r, prev, n, u_0t, u_at):

    s = 1 - 2*r
    b = prev[1:-1]

    b[0] -= r * u_0t
    b[-1] -= r * u_at

    n -= 2

    a = [r for x in range(n-1)]
    d = [s for x in range(n)]
    c = [r for x in range(n-1)]

    for k in range(1, n-1):
        quot = a[k-1] / d[k-1]
        d[k] -= quot * c[k-1]
        b[k] -= quot * b[k-1]


    x0 = [0 for x in range(n)]
    x0[n-2] = b[n-2] / d[n-2]
    for k in range(n-3, -1, -1):
        x0[k] = (b[k] -c[k] * x0[k+1]) / d[k]

    x0.insert(0,u_0t)
    x0.append(u_at)
    return x0


def implicit_method(ax, bx, at, bt, n, m, phi, psi1, psi2,a):
    h = (bx - ax)/n
    k = (bt - at)/m

    xs = [x * h for x in range(0, n)]
    ts = [x * k for x in range(0, m)]

    r = -a*a*k/(h*h)

    u = np.zeros((n, m))

    for i in range(0, n):
        u[i, 0] = phi(xs[i])

    prev = [u[i, 0] for i in range(0, n)]

    for j in range(1, m):
        next = calculate_next_value(r, prev, n, psi1(ts[j]), psi2(ts[j]))
        for i in range(0, n):
            u[i, j] = next[i]
            prev = next
    return xs, ts, u, r



def explicit_method(ax, bx, at, bt, n, m, phi, psi1, psi2,a):
    h = (bx - ax)/n
    k = (bt - at)/m

    xs = [x * h for x in range(0, n)]
    ts = [x * k for x in range(0, m)]

    r = a*a*k/(h*h)
    s = 1 - 2*r

    u = np.zeros((n, m))

    for i in range(0, m):
        u[0, i] = psi1(ts[i])
        u[n-1, i] = psi2(ts[i])

    for i in range(1, n-1):
        u[i, 0] = phi(xs[i])

    for j in range(1, m):
        for i in range(1, n-1):
            u[i, j] = s * u[i, j-1] + r * (u[i-1, j-1] + u[i+1, j-1])


    return xs, ts, u, r


def plot(filename, title, x, t, u):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X, Y = np.meshgrid(t, x)
    surf = ax.plot_surface(X, Y, u, rstride=1, linewidth=0)
    #ax.set_zlim(-2, 10)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title(title)
    ax.zaxis.set_major_locator(LinearLocator(5))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.0e'))
    ax.view_init(azim=61, elev=0)
    #plt.show()
    #fig.colorbar(surf)
    plt.savefig(filename + '.png')
    plt.close()


def start(n,m):
    A = 3
    w = 2
    a = 0.5

    phi = lambda x: 0
    psi1 = lambda t:A*t;
    # psi1 = lambda t: A * t;
    # psi2 = lambda t: math.sin(w * t)
    psi2 = lambda t:  math.sin(w*t)
    # phi = lambda x: math.sin(x)
    # psi1= lambda t: 0
    # psi2 = lambda t: 0
    np.set_printoptions(threshold=np.nan)

    ax = 0
    bx = 10
    # bx = 3 * math.pi
    at = 0
    bt = 3 * math.pi/2.0
    # bt = 5
    # n = 3000
    # m = 3000
    x, t, u1, r = explicit_method(ax, bx, at, bt, n, m, phi, psi1, psi2,a)
    plot("Metoda jawna-{}-{}".format(n,m), "Jawna n = {}, m = {}, p = {:.5f}".format(m, n, abs(r)), x,t,u1)
    x, t, u2,r = implicit_method(ax, bx, at, bt, n, m, phi, psi1, psi2,a)
    plot("Metoda niejawna-{}-{}".format(n,m), "Niejawna, n = {}, m = {}, p = {:.5f}".format(m, n, abs(r)), x,t,u2)
    print("{},{},{},{}".format(n,m,euclidean_error(u1,u2,n,m), abs(r)))


for n in [5,10, 15, 20, 25, 30,35,40,45,50,55,60,65,70,75]:
    for m in [5,10, 15, 20, 25, 30,35,40,45,50,55,60,65,70,75]:
        start(n,m)