import matplotlib.pyplot as plt
import numpy as np
import math

k = 2
n = 1
x_a = -1
x_b = 1

def solve(matrix, vector, n):
    result = np.zeros(n + 1)
    result[1: n] = np.linalg.solve(matrix, vector)
    return result

def solution_of_system(n, a_f, b_f):
    for i in range(1, 10):
        n = n + 2
        h = (x_b - x_a) / n
        matrix, vector, X1 = system_creation(n, h, a_f, b_f)
        Y1 = solve(matrix, vector, n)
        plt.plot(X1, Y1)
        #print(matrix)
        #print(vector)
    plt.show()

def system_creation(n, h, a_f, b_f):
    xi_list = [x_a + (h * i) for i in range(n + 1)]
    matrix = np.zeros((n - 1, n - 1))
    vector = np.linspace(-h * h, -h * h, n - 1)
    ##################
    for i in range(n - 1):
        line = np.linspace(0, 0, n - 1)
        a = a_f(xi_list[i + 1])
        b = b_f(xi_list[i + 1])
        if i == 0:
            line[i:i + 2] = np.array([-(2 * a - (1 + b * (xi_list[i + 1] ** 2)) * h ** 2), a])
        elif i == (n - 2):
            line[i - 1:i + 1] = np.array([a, -(2 * a - (1 + b * (xi_list[i + 1] ** 2)) * h ** 2)])
        else:
            line[i - 1:i + 2] = np.array([a, -(2 * a - (1 + b * (xi_list[i + 1] ** 2)) * h ** 2), a])
        matrix[i, :] = line

    return matrix, vector, xi_list

def main():
    solution_of_system(n, lambda x=1: 1, lambda x=1: 1)
    solution_of_system(n, lambda x: math.sin(k), lambda x: math.cos(k)) 
    solution_of_system(n, lambda x: math.sin(k * x), lambda x: math.cos(k * x))
    plt.show()

if __name__ == "__main__":
    main()
