def solve_newton(f, df, x, stop_condition, rho):
 
    y = f(x)
    i = 0

    while True:
        i += 1
        x_next = x - y/df(x)
        y = f(x_next)

        if stop_condition(f, x, x_next, rho):
            break

        x = x_next

    return x, y, i



def test_ranges(start, end, f, df, stop_condition, step, rho):

    while start <= end:
        yield (start,) + solve_newton(f, df, start, stop_condition, rho)
        #yield (end,) + solve_newton(f, df, end, stop_condition, rho)

        start += step
        #end -= step