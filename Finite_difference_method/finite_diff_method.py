from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

font = {'family': 'serif',
        'color': 'darkred',
        'weight': 'normal',
        'size': 12,
        }

colors = cycle(['b','r'])


def add_to_plot(x, y,label):
    plt.plot(x, y, "--",color=next(colors), label=label, linewidth=2)


def draw_plot(title, file_name):
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 10
    fig_size[1] = 6
    plt.rcParams["figure.figsize"] = fig_size
    plt.legend(loc='upper left', prop={'size': 12})
    plt.title(title, fontdict = font)
    plt.grid(True)
    plt.show()
    plt.savefig(file_name + '.png')
    plt.close()


def calculate_mean_square_error(x, f_y, f_i):
    n = len(x)
    return sum((f_y[i] - f_i[i]) ** 2 for i in range(0, n)) / n

def calculate_max_error(x, f_y, f_i):
    return max(abs(f_y[i] - f_i[i]) for i in range(0, len(x)))

##########################

def calculate_f_x_1(x, k=2, m=1):
    return np.exp((-k) * np.cos(m * x)) - k * np.cos(m * x) + 1


def calculate_differential_equation(x, y, k=2, m=1):
    return k*k*m*np.sin(m*x)*np.cos(m*x) + k*m*y*np.sin(m*x)


@np.vectorize
def calculate_f_x_2b(x, k, m):
    return np.cos(m*x)-x*np.sin(m*x)


@np.vectorize
def calculate_differential_equation_2a(x, k, m):
    return -2*m*np.cos(m * x)
############################


def euler(k, m, x0, xk, n, h, t, w, x):
    for i in range(0, n):
        w.append(w[i] + h * calculate_differential_equation(t, w[i], k, m))
        t += h
        x.append(t)


def runge_kutta(k, m, x0, xk, n, h, t, w, x):
    for i in range(0, n):
        k1 = calculate_differential_equation(t, w[i], k, m)
        k2 = calculate_differential_equation(t + h/2, w[i] + k1 * h/2, k, m)
        k3 = calculate_differential_equation(t + h/2, w[i] + k2 * h/2, k, m)
        k4 = calculate_differential_equation(t + h, w[i] + h*k3, k, m)
        w.append(w[i] + h*(k1 + 2*k2 + 2*k3 + k4)/6)
        t += h
        x.append(t)


def solve_equation(k, m, x0, xk, n, h, t, w, x, solver, method_name, draw_flag):
    solver(k, m, x0, xk, n, h, t, w, x)

    xx = np.arange(x0, xk, 0.1)
    f_y = [calculate_f_x_1(i, k, m) for i in x]

    mean_square_error = calculate_mean_square_error(x, f_y, w)
    max_error = calculate_max_error(x, f_y, w)


    print(",".join((str(n), str(mean_square_error), str(max_error))))

    if draw_flag:
        add_to_plot(xx, calculate_f_x_1(xx, k, m), "Rozwiązanie dokładne")
        add_to_plot(x, w, "Metoda Eulera")
        file_name = method_name + "_" + str(n)
        title = method_name + "  n = " + str(n)
        draw_plot(title, file_name)

def solve_equation_mrs(p, q, r, x_a, x_b, y_a, y_b, n, method_name, draw_flag, m, k):
    j = k
    h = (x_b - x_a) / float(n)
    xs = list(map(lambda j: x_a + h * j, range(1, n)))
    # create matrix
    b = list(map(lambda x: -h * h * r(x), xs))
    b[0] += (1.0 + h / 2.0 * p(xs[0])) * y_a
    b[-1] += (1.0 - h / 2.0 * p(xs[-1])) * y_b

    a = list(map(lambda x: -1.0 + h * p(x) / 2.0, xs[:-1]))
    d = list(map(lambda x: 2 + h * h * q(x), xs))
    c = list(map(lambda x: -1.0 - h * p(x) / 2.0, xs[1:]))
    for k in range(1, n - 1):
        quot = a[k - 1] / d[k - 1]
        d[k] -= quot * c[k - 1]
        b[k] -= quot * b[k - 1]

    x0 = list(map(lambda x: 0, range(n - 1)))
    x0[n - 2] = b[n - 2] / d[n - 2]

    for k in range(n - 3, -1, -1):
        x0[k] = (b[k] - c[k] * x0[k + 1]) / d[k]

    xs.insert(0, x_a)
    xs.append(x_b)

    x0.insert(0, y_a)
    x0.append(y_b)

    f_y = []
    xx = np.linspace(x_a, x_b, 10000)

    for i in xs:
        f_y.append(calculate_f_x_2b(i, j, m))

    max_error = calculate_max_error(xs, f_y, x0)

    max_format = "{:.16f}"
    print(";".join((str(n), max_format.format(max_error))))
    if draw_flag:
        add_to_plot(xx, calculate_f_x_2b(xx, j, m), "Rozwiązanie dokładne")
        add_to_plot(xs, x0, "Rozwiązanie MRS")
        file_name = method_name + "_" + str(n)
        title = method_name + "  n = " + str(n)
        draw_plot(title, file_name)

def simple_test(n):
    k = 2
    m = 1
    x0 = np.pi / 2
    xk = 7 * np.pi / 2
    h = (xk - x0) / n
    t = x0
    a = calculate_f_x_1(x0, k, m)
    w = [float(a)]
    x = [float(t)]

    draw_flag = True

    solve_equation(k, m, x0, xk, n, h, t, w, x, euler, "Euler", draw_flag)

def simple_test_mrs(n):
    k = 4
    m = 2
    a = 0
    b = (2 * np.pi + 2) / m
    p = lambda x: 0.0
    q = lambda x: -m * m
    r = lambda x: -2 * m * np.cos(m * x)
    y_exact = lambda x: np.cos(m * x) - x * np.sin(m * x)
    y_a = 1

    # b = (2 * np.pi + k) / m
    # r = lambda x: m * m * k * x
    # y_exact = lambda x: -k * np.sin(m * x) + k * x
    # y_a = 0

    y_b = y_exact(b)
    solve_equation_mrs(p, q, r, a, b, y_a, y_b, n, "Metoda MRS", True, m, k)

if __name__ == '__main__':
    #for n in [5,8,10,13,15,18,20,30,40,50,60,70,80,100,125,150,175,200,300,400,600,1000,1500,2000,3000,4000,6000,8000,10000]:
    #  simple_test_mrs(n)
    simple_test_mrs(30)