import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sympy as sp
import numpy as np
import math
#------------------------------------------------------------------------------------------------
u_a = 1.5
u_b = 2.5
y_a = 3
y_b = -3
k_1 = 80
k_2 = 1
k_3 = 20
x_0 = 100
c = u_a + (u_b - u_a) / 2
h = (u_b - u_a) / 150
x, t = sp.symbols('x t')
k = sp.Piecewise( (k_1, x < u_a + (u_b - u_a) / 3), (k_3, u_a + 2 * (u_b - u_a) / 3 <= x), (k_2, True), )
source = [(100, u_a + (u_b - u_a) / 2)]
#------------------------------------------------------------------------------------------------
def z4_solution(a, b, k, T_t, delta_function_s_integral, g1, g2, f):
    x_h = (b - a) / 50
    x_h_step = int((b - a) / x_h) + 1
    t_h = 0.5 * x_h ** 2 / k
    t_h_step = int(T_t / t_h) + 1
    x_hs = np.linspace(a, b, x_h_step)
    t_hs = np.linspace(0, T_t, t_h_step)
    matrix = np.zeros(shape = (t_h_step, x_h_step))
    matrix[0, 1: -1] = np.array([delta_function_s_integral(x_hs[i], t_hs[0]) for i in range(1, x_h_step - 1)])
    matrix[ :, 0] = np.array([g1(x_hs[0], t_hs[i]) for i in range(t_h_step)])
    matrix[ :, -1] = np.array([g2(x_hs[-1], t_hs[i]) for i in range(t_h_step)])
    coefficient = np.array([k * t_h / x_h ** 2, 1 - 2 * k * t_h / x_h ** 2, k * t_h / x_h ** 2])
    for i in range(1, t_h_step):
        for j in range(1, x_h_step - 1):
            matrix[i][j] = matrix[i - 1, j - 1: j + 2].dot(coefficient) + t_h * f(x_hs[j], t_hs[i - 1])
    ax = plt.axes(projection='3d')
    ax.set_ylabel('$T:$ Временная ось')
    ax.set_xlabel('$X:$ Пространственная ось')
    ax.set_zlabel('$Y:$ Ось значений')
    for i in range(0, t_h_step, 10):
        ax.plot3D(x_hs, np.array([t_hs[i]] * x_h_step), matrix[i, :])
#------------------------------------------------------------------------------------------------
def z3_solution(x_h, t_h, a, b, k, T_t, g1, g2, delta_function_s_integral, f):
    x_h_step = 1 + int((b - a) / x_h)
    t_h_step = 1 + int(T_t / t_h)
    x_hs = np.linspace(a, b, x_h_step)
    t_hs = np.linspace(0, T_t, t_h_step)
    matrix = np.zeros(shape = (t_h_step, x_h_step))
    matrix[0, 1: -1] = np.array([delta_function_s_integral(x_hs[i]) for i in range(1, x_h_step - 1)])
    matrix[:, 0] = np.array([g1(x_hs[0], t_hs[i]) for i in range(t_h_step)])
    matrix[:, -1] = np.array([g2(x_hs[-1], t_hs[i]) for i in range(t_h_step)])
    for i in range(1, t_h_step):
        for j in range(1, x_h_step - 1):
            matrix[i,j] = sum([k(x_hs[j] - x_h / 2) * t_h / x_h ** 2 * matrix[i - 1, j - 1], (1 - (k(x_hs[j] - x_h / 2) + k(x_hs[j] + x_h / 2)) * t_h / x_h ** 2) * matrix[i - 1, j], k(x_hs[j] + x_h / 2) * t_h / x_h ** 2 * matrix[i - 1, j + 1], t_h * f(x_hs[j], t_hs[i]) * (1 - math.exp(-t_hs[i])) ] )
    ax = plt.axes(projection ='3d')
    ax.set_ylabel('$T:$ Временная ось')
    ax.set_xlabel('$X:$ Пространственная ось')
    ax.set_zlabel('$Y:$ Ось значений')
    for i in range(0, t_h_step, 100):
        ax.plot3D(x_hs, np.array([t_hs[i]] * x_h_step), matrix[i, :])
#------------------------------------------------------------------------------------------------
def delta_function_s_integral(x, x_0, c):
    if abs(x - x_0) - h / 2 < 1e-5:
        return c / 2
    elif 2 * x - h < 2 * x_0 < 2 * x + h:
        return c
    else:
        return 0
def z2_solution(a, b, y_a, y_b, h, delta_function_s_integral, k_expr, sources):
    l, r = sp.symbols('l r')
    a_b = sp.lambdify((l, r), h * (sp.integrate(1 / k_expr, (x, l, r))) ** (-1))
    n = int((b - a) / h) + 1
    matrix = np.zeros(shape = (n, n))
    t = np.zeros(shape = (n, 1))
    xs = np.linspace(a, b, n)
    matrix[0, 0] = matrix[-1, -1] = 1
    t[0] = y_a
    t[-1] = y_b
    for i in range(1, n - 1):
        matrix[i, i - 1] = a_b(xs[i - 1], xs[i])
        matrix[i, i] = -a_b(xs[i - 1], xs[i]) - a_b(xs[i], xs[i + 1])
        matrix[i, i + 1] =  a_b(xs[i], xs[i + 1])
        t[i] = -h * sum(delta_function_s_integral(xs[i], x_0_i, ci) for ci, x_0_i in sources)
    return xs, np.linalg.solve(matrix, t)

xs, ys = z2_solution(u_a, u_b, y_a, y_b, h, delta_function_s_integral, k, source)
plt.plot(xs, ys)
plt.grid()
plt.show()
#------------------------------------------------------------------------------------------------
u_a = 0
u_b = 1
y_a = 0
y_b = 0
x_0 = 100
c = u_a + (u_b - u_a) / 2
h = (u_b - u_a) / 150
sources = [ [(10, u_a + (u_b - u_a) / 2)], [(10, u_a + (u_b - u_a) / 4), (10, u_a + 3 * (u_b - u_a) / 4)], [(10, u_a + (u_b - u_a) / 4), (50, u_a + 3 * (u_b - u_a) / 4)], ]
for k_1, k_2 in [[1, 100], [100, 1]]:
    k = sp.Piecewise( (k_1, x < u_a + (u_b - u_a) / 2), (k_2, True), )
    for source in sources:
        xs, ys = z2_solution(u_a, u_b, y_a, y_b, h, delta_function_s_integral, k, source)
        plt.plot(xs, ys, label = '{} {} {}'.format(k_1, k_2, source))
        plt.legend()
        plt.grid()
        plt.show()
for k_1, k_2, k_3 in [ [1, 3, 9], [9, 3, 1], [1, 2, 1], [20, 1, 20] ]:
    k = sp.Piecewise( (k_1, x < u_a + (u_b - u_a) / 3), (k_3, u_a + 2 * (u_b - u_a) / 3 <= x), (k_2, True), )
    for source in sources:
        xs, ys = z2_solution(u_a, u_b, y_a, y_b, h, delta_function_s_integral, k, source)
        plt.plot(xs, ys, label = '{} {} {} {}'.format(k_1, k_2, k_3, source))
        plt.legend()
        plt.grid()
        plt.show()
#------------------------------------------------------------------------------------------------
a = 1.5
b = 5/2
g1 = sp.lambdify((x, t), 3)
g2 = sp.lambdify((x, t), 3)
f = sp.lambdify((x, t), x + x ** 0.5)
k = sp.lambdify(x, x ** (-1 / 3))
delta_function_s_integral = sp.lambdify(x, 12 * (x - 2) ** 2)
h_x = 0.05
h_t = 0.001
T_t = 500 * 0.001
z3_solution(h_x, h_t, a, b, k, T_t, g1, g2, delta_function_s_integral, f)
plt.show()
#------------------------------------------------------------------------------------------------
k = 1
T_t = 0.5
a = 0
b = 1
delta_function_s_integral = sp.lambdify((x, t), 0)
g1 = sp.lambdify((x, t), 0)
g2 = sp.lambdify((x, t), 0)
f = sp.lambdify((x, t), x)
z4_solution(a, b, k, T_t, delta_function_s_integral, g1, g2 , f)
plt.show()