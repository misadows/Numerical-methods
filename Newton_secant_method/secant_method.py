

def solve_secant(f, a, b, stop_condition, rho):
    val_a, val_b, iterations = f(a), f(b), 0

    while True:
        iterations += 1

        if abs(val_a) > abs(val_b):
            a, b = b, a
            val_a, val_b = val_b, val_a

        s = (b - a)/(val_b - val_a)
        b = a
        val_b = val_a
        a = a - val_a * s
        val_a = f(a)

        if stop_condition(f, b, a, rho):
            break

    return iterations, a, val_a


def test_ranges(start, end, f, stop_condition, step, rho):
    a, b = start, end

    while a <= b:
        yield (a, end) + solve_secant(f, a, end, stop_condition, rho)
        a += step

    a, b = start, end
    while a <= b:
        yield (start, b) + solve_secant(f, a, end, stop_condition, rho)
        b -= step
