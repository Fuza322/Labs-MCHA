import numpy as nump
from matplotlib import pyplot as graphic
import matplotlib
import seaborn as plate
#------------------------------------------------------------------------
n_x = 40
n_y = 50
u_0 = lambda x, y: nump.arctan(nump.cos(nump.pi * x / A))
u_0_part = lambda x, y: nump.sin(2 * nump.pi * x / A) * nump.sin(nump.pi * y / B)
B = 0.02
A = 0.01
n_t = 2000
T = 0.2
#------------------------------------------------------------
δ_x = A / n_x
δ_y = B / n_y
δ_t = T / n_t
Cond = δ_t / δ_x + δ_t / δ_y
if Cond > 1:
	print('---------------------------------------------------')
	print('ERROR!!! The condition is not satisfied')
	print('---------------------------------------------------')
#------------------------------------------------------------
x_value = nump.linspace(-A / 2, A / 2, n_x)
y_value = nump.linspace(-B / 2, B / 2, n_y)
mas = nump.zeros((n_t, n_x, n_y))
mas[1, 1: -1, 1] = mas[1, 1: -1, 0]
mas[1, 1: -1, -1] = mas[1, 1: -1, -2]
for i in range(1, n_x - 1):
	for j in range(1, n_y - 1):
		mas[0, i, j] = u_0(x_value[i], y_value[j])
for i in range(1, n_x - 1):
	for j in range(1, n_y - 1):
		mas[1, i, j] = (u_0(x_value[i], y_value[j]) + u_0_part(x_value[i], y_value[j]) * δ_t + δ_t ** 2 / (2 * δ_x ** 2) * (mas[0, i + 1, j] - 2 * mas[0, i, j] + mas[0, i - 1, j]) + δ_t ** 2 / (2 * δ_y ** 2) * (mas[0, i, j + 1] - 2 * mas[0, i, j] + mas[0, i, j - 1]))  
for v in range(1, n_t - 1):
	mas[v+1, 1: -1, 1: -1] = (2 * mas[v, 1: -1, 1: -1] - mas[v - 1, 1: -1, 1: -1] + δ_t ** 2 / δ_x ** 2 * (mas[v, : -2, 1: -1] - 2 * mas[v, 1: -1, 1: -1] + mas[v, 2:, 1: -1]) + δ_t ** 2 / δ_y ** 2 * (mas[v, 1: -1, : -2] - 2 * mas[v, 1: -1, 1: -1] + mas[v, 1: -1, 2:]))
#------------------------------------------------------------
for i in range(0, 200, 20):
    plate.heatmap(mas[i].T, cmap='BuPu', vmin=mas.min(), vmax=mas.max())
    graphic.show()
#------------------------------------------------------------