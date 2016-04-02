import numpy as np
import matplotlib.pyplot as plt


def fun(x):
    return np.sin(2*x)*np.sin(2*x*x/np.pi)

def fun_der(x):
    return -6.*np.sin(2*x)*np.exp(3*np.cos(2*x))


def calculate_mean_square_error(n, function_y, interpolated_y):
    error = 0
    for i in range(0, n):
        error += (function_y[i] - interpolated_y[i])**2
    error /= n
    return error


def calculate_max_error(n, function_y, interpolated_y):
    max_error = 0
    for i in range(0, n):
        max_error = max(np.fabs(function_y[i] - interpolated_y[i]), max_error)
    return max_error


def print_error(N, x_len, y_values, y_interpolated, label):
    print('{0},{1},{2},{3}'.format(label, N,
                                   calculate_mean_square_error(x_len, y_values, y_interpolated),
                                   calculate_max_error(x_len, y_values, y_interpolated)
                                   ))


def get_s_matrix_pol(x, n):
    m = n + 1
    S = np.empty((m, m))

    for i in range(m):
        for j in range(m):
            S[i][j] = (sum([xp**(i+j) for xp in x]))

    return S


def get_t_vector_pol(points, n):
    m = n + 1
    t = []

    for i in range(m):
        t.append(sum([y*(x**i) for x, y in points]))

    return t


def polynomial_approximate(points, n):
    x, y = zip(*points)

    S = get_s_matrix_pol(x, n)
    t = get_t_vector_pol(points, n)

    A = np.linalg.solve(S, t)

    def function(x_point):
        result = 0
        for i in range(0, n+1):
            result += (x_point**i) * A[i]
        return result

    return function


def polynomial_fun(x, i):
    return np.cos(x * i / 2) if i % 2 == 0 else np.sin(x * (i + 1) / 2)


def get_s_matrix_trig(x, n):
    m = n + 1
    S = np.empty((m, m))
    indxs = range(m)

    for i in indxs:
        for j in indxs:
            S[i][j] = sum([polynomial_fun(x[k], j) * polynomial_fun(x[k], i) for k in range(len(x))])

    return S


def get_t_vector_trig(points, n):
    x, y = zip(*points)
    m = n + 1
    t = []

    for i in range(m):
        sum = 0
        for j in range(len(x)):
            sum += polynomial_fun(x[j], i) * y[j]
        t.append(sum)

    return t


def trigonometric_approximate(points, n):
    x, y = zip(*points)

    S = get_s_matrix_trig(x, n)
    t = get_t_vector_trig(points, n)

    A = np.linalg.solve(S, t)

    def function(x_point):
        result = 0
        for i in range(0, n+1):
            result += polynomial_fun(x_point, i) * A[i]
        return result

    return function

def get_regular_interpolation_nodes(n, fun, start, end):
    return [(x, fun(x)) for x in np.linspace(start, end, n)]


def get_czebyshev_interpolation_nodes(n, fun, start, end):
    points = []
    half_length = (end - start) / 2
    middle = start + half_length
    for i in range(0, n):
        x = middle + np.cos(((2 * i + 1) * np.pi) / (2 * n)) * half_length
        y = fun(x)
        points.append((x, y))
    return points

def simple_test(N, deg):
    start = -np.pi
    end = 1.*np.pi
    step = 0.05
    function = fun

    x = np.arange(start, end, step)
    y = list(map(function, x))

    x_len = len(x)

    points = get_regular_interpolation_nodes(N, function, start, end)
    #points = get_czebyshev_interpolation_nodes(N, function, start, end)

    x_nodes, y_nodes = zip(*points)

    polynomial = polynomial_approximate(points, deg)
    trigonometric = trigonometric_approximate(points, deg)

    polynomial_y = list(map(polynomial, x))
    trigonometric_y = list(map(trigonometric, x))

    #print('label,N,mean_sqaure_error,max_error')
    #print_error(deg, x_len, y, polynomial_y, 'wielomianowa')
    print_error(deg, x_len, y, trigonometric_y, 'trigonometric')
    return x,y,x_nodes,y_nodes,polynomial_y,trigonometric_y

def plot_test(N, deg):
    x,y,x_nodes,y_nodes,polynomial_y,trigonometric_y = simple_test(N, deg)
		
    plt.plot(x, y, "y", linewidth=5, label='Funkcja')
    #plt.plot(x, polynomial_y, "r--", linewidth=4, label='Wielomianowa')
    plt.plot(x, trigonometric_y, "g:", linewidth=4, label='Trygonometyczna')
    plt.plot(x_nodes, y_nodes, 'yo')
    plt.legend(loc='upper right', prop={'size': 10})
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    #for i in range(2, 31):
    #    simple_test(30, i)
    plot_test(10, 9)

