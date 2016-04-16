import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

n = 30
m = 30

a = 0.5
A = 3
w = 2
l = 10
tc = 3 * np.pi / 2


def psi_1(t):
    return A * t


def psi_2(t):
    return np.sin(w * t)


def phi(x):
    return 0


h = l / n
k = tc / m

matrix = np.zeros((n + 1, m + 1))

for i in range(n+1):
    matrix[i][0] = phi(i * h)

for j in range(m + 1):
    matrix[0][j] = psi_1(j * k)
    matrix[n][j] = psi_2(j * k)

mlt = (a ** 2) * k / (h ** 2)
for j in range(0, m):
    for i in range(1, n):
        matrix[i][j + 1] = mlt * (matrix[i + 1][j] - 2 * matrix[i][j] + matrix[i - 1][j]) + matrix[i][j]

print(mlt)
x = np.linspace(0, l, n)
t = np.linspace(0, tc, m)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(range(0,n+1), range(0,m+1))

ax.plot_surface(X, Y, matrix)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
