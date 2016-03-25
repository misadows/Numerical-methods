import numpy as np
import matplotlib.pyplot as plt


def calculate_cubic_value(t, y):
    n = len(t)
    h = []
    b = []
    u = [0]
    v = [0]

    for i in range(0, n - 1):
        h.append(t[i + 1] - t[i])
        b.append(6 * (y[i + 1] - y[i]) / h[i])

    u.append(2 * (h[0] + h[1]))
    v.append(b[1] - b[0])

    for i in range(2, n - 1):
        u.append(2 * (h[i - 1] + h[i]) - (h[i - 1] * h[i - 1]) / u[i - 1])
        v.append((b[i] - b[i - 1]) - (h[i - 1] * v[i - 1] / u[i - 1]))

    z = [0] * (n)

    for i in range(n - 2, 0, -1):
        z[i] = ((v[i] - h[i] * z[i + 1]) / u[i])

    z[0] = 0
    return z, h

def calculate_cubic_value_not_a_knot(t, y):
    n = len(t)
    h = []
    b = []
    v = [0]
    A = np.zeros(shape=(n, n))

    for i in range(0, n - 1):
        h.append(t[i + 1] - t[i])
        b.append(6 * (y[i + 1] - y[i]) / h[i])

    A[0][0] = h[1]
    A[0][1] = - h[1] - h[0]
    A[0][2] = h[0]

    for i in range(1, n - 1):
        A[i][i - 1] = h[i - 1]
        A[i][i] = 2 * (h[i - 1] + h[i])
        A[i][i + 1] = h[i]
        v.append(b[i] - b[i - 1])

    v.append(0)
    A[n - 1][n - 3] = h[n - 2]
    A[n - 1][n - 2] = - h[n - 2] - h[n - 3]
    A[n - 1][n - 1] = h[n - 3]

    z = np.linalg.solve(A, v)

    return z, h


def create_cubic_function(t, y, calc):
    z, h = calc(t, y)

    def spline_function(x):
        i = 0
        while float(x) - t[i] >= 0:
            i += 1
        i -= 1
        A = (z[i + 1] - z[i]) / (6 * h[i])
        B = z[i] / 2.
        C = -(h[i] * (z[i + 1] + 2 * z[i])) / 6 + (y[i + 1] - y[i]) / h[i]
        return y[i] + ((x - t[i]) * (C + (x - t[i]) * (B + (x - t[i]) * A)))

    return spline_function

def calculate_mean_square_error(n, function_y, interpolated_y):
    error = 0
    for i in range(0, n):
        error += (function_y[i] - interpolated_y[i])**2
    error /= n
    return error


def calculate_max_error(n, function_y, interpolated_y):
    max_error = 0
    for i in range(0, n):
        if np.fabs(function_y[i] - interpolated_y[i]) > max_error:
            max_error = np.fabs(function_y[i] - interpolated_y[i])
    return max_error


def print_error(N, x_len, y_values, y_interpolated, label):
    print('{0},{1},{2},{3}'.format(label, N,
                                   calculate_mean_square_error(x_len, y_values, y_interpolated),
                                   calculate_max_error(x_len, y_values, y_interpolated)
                                   ))


def fun(x):
    return np.sin(2*x)*np.sin(2*x*x/np.pi)


def get_points(n, start, end, function):
    points = []
    r = end - start
    step = r / n

    for i in range(n):
        points.append((start+i*step, function(start+i*step)))
    points.append((end, function(end)))

    return points


def calculate_quadratic_value(t, y):
    n = len(t)
    z = [0]

    for i in range(0, n - 1):
        z.append(-z[i] + 2 * ((y[i + 1] - y[i]) / (t[i + 1] - t[i])))
    return z


def create_quadratic_function(t, y, calc):
    z = calc(t, y)

    def spline_function(x):
        i = 0
        while x - t[i] >= 0:
            i += 1
        i -= 1
        return ((z[i + 1] - z[i]) * (x - t[i]) * (x - t[i])) / (2 * (t[i + 1] - t[i])) + z[i] * (x - t[i]) + y[i]

    return spline_function


def calculate_quadratic_value_not_a_knot(t, y):
    n = len(t)
    d = [(y[1] - y[0]) / (t[1] - t[0])]
    A = np.zeros(shape=(n, n))
    A[0][0] = 1

    for i in range(1, n):
        A[i][i-1] = 1
        A[i][i] = 1
        d.append(2 * ((y[i] - y[i-1]) / (t[i] - t[i-1])))

    return np.linalg.solve(A, d)


def create_quadratic_function(t, y, calc):
    z = calc(t, y)

    def spline_function(x):
        i = 0
        while x - t[i] >= 0:
            i += 1
        i -= 1
        return ((z[i + 1] - z[i]) * (x - t[i]) * (x - t[i])) / (2 * (t[i + 1] - t[i])) + z[i] * (x - t[i]) + y[i]

    return spline_function


def simple_test(N):
    start = -np.pi
    end = np.pi
    step = 0.05
    function = fun

    points = get_points(N-1, start, end, function)

    x_points, y_points = zip(*points)

    x = np.arange(start, end, step)
    y = [function(x) for x in x]
    x_len = len(x)

    quadratic_spline_natural = create_quadratic_function(x_points, y_points, calculate_quadratic_value)
    quadratic_spline_notaknot = create_quadratic_function(x_points, y_points, calculate_quadratic_value_not_a_knot)
    cubic_spline_natural = create_cubic_function(x_points, y_points, calculate_cubic_value)
    cubic_spline_not_a_knot = create_cubic_function(x_points, y_points, calculate_cubic_value_not_a_knot)

    quadratic_spline_natural_values = [quadratic_spline_natural(x) for x in x]
    quadratic_spline_notaknot_values = [quadratic_spline_notaknot(x) for x in x]
    cubic_spline_natural_values = [cubic_spline_natural(x) for x in x]
    cubic_spline_notaknot_values = [cubic_spline_not_a_knot(x) for x in x]

    print('{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(N,
                                   calculate_mean_square_error(x_len, y, quadratic_spline_natural_values),
                                   calculate_max_error(x_len, y, quadratic_spline_natural_values),
                                   calculate_mean_square_error(x_len, y, cubic_spline_natural_values),
                                   calculate_max_error(x_len, y, cubic_spline_natural_values),
                                   calculate_mean_square_error(x_len, y, quadratic_spline_notaknot_values),
                                   calculate_max_error(x_len, y, quadratic_spline_notaknot_values),
                                   calculate_mean_square_error(x_len, y, cubic_spline_notaknot_values),
                                   calculate_max_error(x_len, y, cubic_spline_notaknot_values)
                                   ))

    return x, function, cubic_spline_natural_values, quadratic_spline_natural_values, x_points, y_points, cubic_spline_notaknot_values, quadratic_spline_notaknot_values

def plot_test(x, function, cubic_spline_natural_values, quadratic_spline_natural_values, x_points, y_points, cubic_spline_notaknot_values, quadratic_spline_notaknot_values):
    plt.plot(x, function(x), 'k', linewidth=2.5, label="funkcja")
    plt.plot(x, cubic_spline_natural_values, 'r--', linewidth=4, label='3 stopnia - natural splines')
    plt.plot(x, quadratic_spline_natural_values, 'g--', linewidth=4, label='2 stopnia - natural splines')
    plt.plot(x, quadratic_spline_notaknot_values, 'g:', color='#808000', linewidth=4, label='2 stopnia - not a knot')
    plt.plot(x, cubic_spline_notaknot_values, 'r:', color='#880000', linewidth=4, label='3 stopnia - not a knot')
    plt.plot(x_points, y_points, 'yo')
    plt.legend(loc='lower right')
    plt.show()

def error_stat(N_max):
    print('Liczba wezlow,2 stopnia - natural splines,2 stopnia - natural splines,3 stopnia - natural splines,3 stopnia - natural splines, 2 stopnia - not a knot,2 stopnia - not a knot,3 stopnia - not a knot,3 stopnia - not a knot')
    for i in range(4, N_max+1):
        simple_test(i)

if __name__ == '__main__':
    plot_test(*simple_test(50))




# czy rozwiazuje iteracyjnie czy układem równań