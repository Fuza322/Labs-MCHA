from matplotlib import pyplot as graphic
import numpy as nump
#------------------------------------------------------------------------
fun = lambda x: - nump.fabs((x - leng / 2) * 2 * Δ_u / leng) + Δ_u
ρ = 5900
n_x = 100
n_t = 10000
Elasticity = 120000000000
leng = 0.18
Δ_u = 0.002
T = 1e-3
#------------------------------------------------------------------------
dx = leng / n_x
dt = T / n_t
Cond = nump.sqrt(Elasticity / ρ) * dt / dx
if Cond > 1:
    print('---------------------------------------------------')
    print('ERROR!!! The convergence condition is not satisfied')
    print('---------------------------------------------------')
values_on_OX = nump.linspace(0, leng, n_x)
 #------------------------------------------------------------------------
mas = nump.zeros((n_t, n_x))
mas[0] = [fun(x) for x in values_on_OX]
mas[1, 1:-1] = [fun(values_on_OX[i]) * (1 - dt ** 2 / 2) for i in range(1, n_x - 1)]
#------------------------------------------------------------------------
for i in range(1, n_t - 1):
    mas[i + 1, 1: - 1] = Cond ** 2 * (mas[i, 2: ] - 2 * mas[i, 1: -1] + mas[i, : -2]) + 2 * mas[i, 1: -1] - mas[i-1, 1: -1]
values_on_OX = nump.linspace(0, leng, n_x)
#------------------------------------------------------------------------
for i in range(0, n_t, int(n_t / 20)):
	graphic.plot(values_on_OX, mas[i])
	graphic.xlabel('X', size = 20)
	graphic.ylabel('U(X, t)', size = 20)
graphic.grid()
graphic.show()
#------------------------------------------------------------