import numpy as np
import matplotlib.pyplot as plt

def calculate_mean_square_error(x, function_y, interpolated_y):
    n = len(x)
    error = 0
    for i in range(0, n):
        error += (function_y[i] - interpolated_y[i])**2
    error /= n
    return error


def calculate_max_error(x, function_y, interpolated_y):
    n = len(x)
    max_error = 0
    for i in range(0, n):
        if np.fabs(function_y[i] - interpolated_y[i]) > max_error:
            max_error = np.fabs(function_y[i] - interpolated_y[i])
    return max_error


def fun(x):
    return np.sin(2*x)*np.sin(2*x*x/np.pi)


def d_fun(x):
    return 4*x*np.sin(2*x)*np.cos(2*x*x/np.pi)/np.pi + 2*np.sin(2*x*x/np.pi)*np.cos(2*x)


def lagrange_interpolate(points):
    def result_polynomial(x):
        sum = 0
        n = len(points)
        for i in range(n):
            xi, yi = points[i]
            product = 1
            for j in range(n):
                if i != j:
                    xj, yj = points[j]
                    product *= (x - xj) / float(xi - xj)
            sum += yi*product
        return sum
    return result_polynomial


def newton_interpolate(points):
    x1, y1 = zip(*points)
    n = len(points)
    dq = list(y1)
    x = list(x1)
    for i in range(1, n):
        for j in range(n - 1, i - 1, -1):
            dq[j] = (dq[j] - dq[j - 1]) / (x[j] - x[j - i])

    def result_polynomial(xpoint):
        val = dq[0]
        factor = 1.0
        for i in range(1, n):
            factor *= (xpoint - x[i - 1])
            val += (dq[i] * factor)
        return val

    return result_polynomial


def hermit_interpolate(input, df):
    n = len(input)
    points = np.zeros(shape=(2 * n + 1, 2 * n + 1))
    X, Y = zip(*input)
    X = list(X)
    Y = list(Y)

    for i in range(0, 2 * n, 2):
        points[i][0] = X[i // 2]
        points[i + 1][0] = X[i // 2]
        points[i][1] = Y[i // 2]
        points[i + 1][1] = Y[i // 2]

    for i in range(2, 2 * n + 1):
        for j in range(1 + (i - 2), 2 * n):
            if i == 2 and j % 2 == 1:
                points[j][i] = df(X[j // 2]);

            else:
                points[j][i] = (points[j][i - 1] - points[j - 1][i - 1]) / (
                    points[j][0] - points[(j - 1) - (i - 2)][0])

    def result_polynomial(xpoint):
        val = 0
        for i in range(0, 2 * n):
            factor = 1.
            j = 0
            while j < i:
                factor *= (xpoint - X[j // 2])
                if j + 1 != i:
                    factor *= (xpoint - X[j // 2])
                    j += 1
                j += 1
            val += factor * points[i][i + 1]
        return val

    return result_polynomial

def get_regular_interpolation_nodes(n, start, end):
    interval = end - start
    step = interval / n
    points = []
    for x in np.arange(start, end, step):
        points.append((x, fun(x)))
    points.append((end, fun(end)))
    return points


def get_czebyshev_interpolation_nodes(n, start, end):
    points = []
    half_length = (end - start) / 2
    middle = start + half_length
    for i in range(0, n):
        x = middle + np.cos(((2 * i + 1) * np.pi) / (2 * n)) * half_length
        y = fun(x)
        points.append((x, y))
    return points


def simple_test(N, interpolation_nodes):
    start = -np.pi
    end = np.pi
    x = np.arange(start, end, 0.1)
    x = np.append(x, end)


    points = interpolation_nodes(N, start, end)

    lagrange_polynomial = lagrange_interpolate(points)
    newton_polynomial = newton_interpolate(points)
    hermit_polynomial = hermit_interpolate(points, d_fun)

    function_y = [fun(x) for x in x]
    lagrange_y = [lagrange_polynomial(x) for x in x]
    newton_y = [newton_polynomial(x) for x in x]
    hermit_y = [hermit_polynomial(x) for x in x]

    x_points, y_points = zip(*points)

    print('{0},{1},{2},{3},{4},{5},{6}'.format(N,
                                       calculate_mean_square_error(x, function_y, lagrange_y),
                                       calculate_mean_square_error(x, function_y, newton_y),
                                       calculate_mean_square_error(x, function_y, hermit_y),
                                       calculate_max_error(x, function_y, lagrange_y),
                                       calculate_max_error(x, function_y, newton_y),
                                       calculate_max_error(x, function_y, hermit_y)))

    return x, function_y, lagrange_y, newton_y, hermit_y, x_points, y_points


def simple_test_plot(x, function_y, lagrange_y, newton_y, hermit_y, x_points, y_points):
    plt.plot(x, function_y, "k-", linewidth=1.5, label="Funkcja poczatkowa")
    plt.plot(x, lagrange_y, "y--", linewidth=1.5, label="Funkcja interpolujaca Lagrange")
    #plt.plot(x, newton_y, "b:", linewidth=1.5, label="Funkcja interpolujaca Newtona")
    plt.plot(x, hermit_y, "g:", linewidth=1.5, label="Funkcja interpolujaca Hermita")
    plt.plot(x_points, y_points, 'yo')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.grid(True)
    plt.show()


def error_stat(N_max, interpolation_nodes):
    print('Wezly,Lagrange_srednia_kwadratow,Newton_srednia_kwadratow,Hermit_srednia_kwadratow,Lagrange_maksymalna_roznica,Newton_maksymalna_roznica,Hermit_maksymalna_roznica')
    for i in range(1, N_max+1):
        simple_test(i, interpolation_nodes)

if __name__ == '__main__':
    simple_test_plot(*simple_test(9, get_czebyshev_interpolation_nodes))
    #error_stat(150, get_czebyshev_interpolation_nodes)
