import math
from prettytable import PrettyTable
import matplotlib.pyplot as plt
import random
import sympy as sp
import numpy as np 
#----------------------------------------------------------------------------------------------------------------------------------
def Overt_scheme_way1(N, t):
    H = (b - a) / N
    t_n = int(T / t) + 1
    t_k = N + 1
    system_M = np.zeros(shape = (t_n, t_k))
    system_M [: , 0] = np.array([G__1(ti) for ti in np.linspace(0, T, t_n)])
    system_M[0, : ] = np.array([delta_func(xi) for xi in np.linspace(a, b, t_k)])
    coefficient = np.array([k * t / H ** 2, 1 - 2 * k * t / H ** 2, k * t / H ** 2])
    for i in range(1, t_n):
        for j in range(1, t_k - 1):
            system_M[i][j] = system_M[i - 1, j - 1 : j + 2].dot(coefficient) + t * f(a + j * H, i * t)
        system_M[i][-1] = system_M[i][-2] + H * G__2(i * t)
    return system_M
#----------------------------------------------------------------------------------------------------------------------------------
def Overt_scheme_way2(N, t):
    H = (b - a) / N
    t_n = int(T / t) + 1
    t_k = N + 1
    system_M = np.zeros(shape = (t_n, t_k))
    system_M [: , 0] = np.array([G__1(ti) for ti in np.linspace(0, T, t_n)])
    system_M[0, : ] = np.array([delta_func(xi) for xi in np.linspace(a, b, t_k)])
    coefficient = np.array([k * t / H ** 2, 1 - 2 * k * t / H ** 2, k * t / H ** 2])
    for i in range(1, t_n):
        for j in range(1, t_k - 1):
            system_M[i, j] = system_M[i - 1, j-1:j + 2].dot(coefficient) + t * f(a + j * H, i * t)
        system_M[i, -1] = 2 * k * t / (H ** 2) * system_M[i - 1, -2] + (1 - 2 * k * t / (H ** 2)) * system_M[i - 1, -1] + t * f(b, i * t) + 2 * k * (t ** 3) / (H ** 2) * G__2((i - 1) * t)
    return system_M
#----------------------------------------------------------------------------------------------------------------------------------
def Unovert_scheme_way1(N, t):
    H = (b - a) / N
    t_n = int(T / t) + 1
    t_k = N + 1
    solution = np.zeros(shape =(t_n, t_k))
    solution [: , 0] = np.array([G__1(ti) for ti in np.linspace(0, T, t_n)])
    solution[0, : ] = np.array([delta_func(xi) for xi in np.linspace(a, b, t_k)])
    system_M = np.zeros(shape = (t_k - 1, t_k - 1))
    system_M[0, 0: 2] = [1 + 2 * k * t / H ** 2, -k * t / H ** 2]
    for j in range(1, t_k - 2):
            system_M[j, j - 1] = - k * t / H ** 2 
            system_M[j, j] = 1 + 2 * k * t / H ** 2
            system_M[j, j+1] = -k * t / H ** 2
    system_M[-1, -1], system_M[-1, -2] = 1, -1 
    ter = np.zeros(t_k - 1)
    for i in range(1, t_n):
        ter[0] = t * f(a + H, i * t) + solution[i - 1, 1] + k * t / H ** 2 * solution[i, 0]
        for j in range(1, t_k - 2):
            ter[j] = t * f(a + (j + 1) * H, i * t) + solution[i - 1, j + 1]
        ter[-1] = t * G__2(i * t)
        solution[i, 1: ] = np.linalg.solve(system_M, ter)
    return solution
#----------------------------------------------------------------------------------------------------------------------------------
def Unovert_scheme_way2(N, t):
    H = (b - a) / N
    t_n = int(T / t) + 1
    t_k = N + 1
    solution = np.zeros(shape = (t_n, t_k))
    solution [:, 0] = np.array([G__1(ti) for ti in np.linspace(0, T, t_n)])
    solution[0, : ] = np.array([delta_func(xi) for xi in np.linspace(a, b, t_k)])
    system_M = np.zeros(shape = (t_k - 1, t_k - 1))
    system_M[0, 0:2] = [1 + 2 * k * t / H ** 2, -k * t / H ** 2]
    for j in range(1, t_k - 2):
            system_M[j, j - 1] = - k * t / H ** 2 
            system_M[j, j] = 1 + 2 * k * t / H ** 2
            system_M[j, j + 1] = -k * t / H ** 2
    system_M[-1, -1], system_M[-1, -2] = -2 * k * t / (H ** 2), 1 + 2 * k * t / (H ** 2) 
    ter = np.zeros(t_k - 1)
    for i in range(1, t_n):
        ter[0] = t * f(a + H, i * t) + solution[i - 1, 1] + k * t / H ** 2 * solution[i, 0]
        for j in range(1, t_k - 2):
            ter[j] = t * f(a + (j + 1) * H, i * t) + solution[i - 1, j + 1]
        ter[-1] = t * f(a + (j + 1) * H, i * t) + solution[i - 1, j + 1] - 2 * k * t / H * G__2(i * t)
        solution[i, 1: ] = np.linalg.solve(system_M, ter)
    return solution

def T_Standart_deflection(m1, m2, Level1, Level2):
    if len(m1[0]) > len(m2[0]):
        m1, m2 = m2, m1
        Level1, Level2 = Level2, Level1
    return np.sqrt(sum((m1[Level1][i] - m2[Level2][2 * i]) ** 2 for i in range(len(m1[0]))) / len(m1[0]))

def H_Standart_deflection(m1, m2, Level1, Level2):
    return np.sqrt(sum((m1[Level1][i] - m2[Level2][i]) ** 2 for i in range(len(m1[0]))) / len(m1[0]))

def Find_deviation_Max_of_T(m1, m2, Level1, Level2):
    if len(m1[0]) > len(m2[0]):
        m1, m2 = m2, m1
        Level1, Level2 = Level2, Level1
    return max(np.abs([m1[Level1][i] - m2[Level2][2 * i] for i in range(len(m1[0]))]))

def Find_deviation_Max_of_H(m1, m2, Level1, Level2):
    return max(np.abs([m1[Level1][i] - m2[Level2][i] for i in range(len(m1[0]))]))

def Find_Fix_deviation_T(func):
    N_s = np.array([5 * 2 ** i for i in range(4, -1, -1)])
    table_t = PrettyTable()
    table_t.field_names = ['N', 'П„', 's(1)', 's(2)', 'max(1)', 'max(2)']
    Level1 = random.randint(1, N_s[-1] - 1)
    Level2 = random.randint(1, N_s[-1] - 1)
    prev_res = None
    Xs = np.linspace(a, b, N_s[0] + 1)
    П„ = 0.5 * (Xs[1] - Xs[0]) ** 2 / k
    for i, N in enumerate(N_s):
        h_x = (b - a) / N
        result = func(N, П„)
        if prev_res is None:
            table_t.add_row([Xs[1] - Xs[0], round(П„, man_len), * list('-'*4)])
        else:
            t11 = Level1 * N_s[i-1] // N_s[-1]
            t21 = Level2 * N_s[i-1] // N_s[-1]
            t22 = Level2 * N_s[i] // N_s[-1]
            t12 = Level1 * N_s[i] // N_s[-1]
            table_t.add_row(list(map(lambda x: round(x, man_len),[ h_x, П„, T_Standart_deflection(prev_res, result, t11, t12), T_Standart_deflection(prev_res, result, t21, t22), Find_deviation_Max_of_T(prev_res, result, t11, t12), Find_deviation_Max_of_T(prev_res, result, t21, t22), ])))
        prev_res = result 
    print(table_t)

def Find_Fix_deviation_H(func):
    T_s = np.array([100*2**i for i in range(5)])
    table_h = PrettyTable()
    table_h.field_names = ['N', 'П„', 's(1)', 's(2)', 'max(1)', 'max(2)']
    max_П„ = T / T_s[0]
    N = int((b-a) / np.sqrt(2 * k * max_П„))
    H = np.sqrt(2 * k * max_П„)
    Level1 = random.randint(1, T_s[0] - 1)
    Level2 = random.randint(1, T_s[0] - 1)
    prev_res = None
    for i, h_t in enumerate(T_s):
        П„ = T / h_t
        result = func(N, П„)
        if prev_res is None:
            table_h.add_row([H, round(П„, man_len), * list('-' * 4)])
        else:
            t21 = Level2 * T_s[i-1] // T_s[0]
            t11 = Level1 * T_s[i-1] // T_s[0]
            t22 = Level2 * T_s[i] // T_s[0]
            t12 = Level1 * T_s[i] // T_s[0]
            table_h.add_row(list(map(lambda x: round(x, man_len),[ H, П„, H_Standart_deflection(prev_res, result, t11, t12), H_Standart_deflection(prev_res, result, t21, t22), Find_deviation_Max_of_H(prev_res, result, t11, t12), Find_deviation_Max_of_H(prev_res, result, t21, t22), ])))
        prev_res = result
    print(table_h)
#----------------------------------------------------------------------------------------------------------------------------------
delta_func = lambda x: 0
G__1 = lambda t: 0
G__2 = lambda t: 0 
a = 0
b = 1
k = 1
T = 0.05
man_len = 10
f = lambda x, t: x
N = 100
t = (((b-a) / N)**2) / 6
#----------------------------------------------------------------------------------------------------------------------------------
xs = np.linspace(a, b, N + 1)
solution = Overt_scheme_way1(N, t)
for i in range(0, len(solution), 100):
    plt.plot(xs, solution[i], label = i)
plt.grid()
plt.show()
print("*******************************РЇРІРЅР°СЏ СЃС…РµРјР° (1 СЃРїРѕСЃРѕР±)*******************************")
Find_Fix_deviation_T(Overt_scheme_way1)
Find_Fix_deviation_H(Overt_scheme_way1)
xs = np.linspace(a, b, N + 1)
solution = Overt_scheme_way2(N, t)
for i in range(0, len(solution), 100):
    plt.plot(xs, solution[i], label = i)
plt.grid()
plt.show()
print("*******************************РЇРІРЅР°СЏ СЃС…РµРјР° (2 СЃРїРѕСЃРѕР±)*******************************")
Find_Fix_deviation_T(Overt_scheme_way2)
Find_Fix_deviation_H(Overt_scheme_way2)
xs = np.linspace(a, b, N + 1)
solution = Unovert_scheme_way1(N, t)
for i in range(0, len(solution), 100):
    plt.plot(xs, solution[i], label = i)
plt.grid()
plt.show()
print("*******************************РќРµСЏРІРЅР°СЏ СЃС…РµРјР° (1 СЃРїРѕСЃРѕР±)*******************************")
Find_Fix_deviation_T(Unovert_scheme_way1)
Find_Fix_deviation_H(Unovert_scheme_way1)
xs = np.linspace(a, b, N + 1)
solution = Unovert_scheme_way2(N, t)
for i in range(0, len(solution), 100):
    plt.plot(xs, solution[i], label = i)
plt.grid()
plt.show()
print("*******************************РќРµСЏРІРЅР°СЏ СЃС…РµРјР° (2 СЃРїРѕСЃРѕР±)*******************************")
Find_Fix_deviation_T(Unovert_scheme_way2)
Find_Fix_deviation_H(Unovert_scheme_way2)