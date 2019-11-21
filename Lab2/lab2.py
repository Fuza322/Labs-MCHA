import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

k = 15
'''c = [lambda x: (1+math.cos(k)*(x ** 2)),
     lambda x: 0,
     lambda x: math.sin(k)]
'''
c = [lambda x: (1+(x ** 2)),
     lambda x: 0,
     lambda x: 1]

c0 = lambda x: -1
df_dx =[lambda f, x: f(x),
        lambda f, x: (f(x+h)-f(x-h)) / (2*h),
        lambda f, x: (f(x-h)-2*f(x)+f(x+h)) / (h**2)]
a0, b0 = -1, 1
n = 100
h = 0.0001
dh = (b0 - a0) / (n+2)

def integration(f, a, b, n = 50):
    h = (b - a) / n
    return h / 3 * sum([f(a + h * i) + 4 * f(a + h * (i + 1)) + f(a + h * (i + 2)) for i  in range(0, n - 1, 2)])

def gen_f(i):
    return lambda x: x**(i-1) * (1-x**2)

def collation(collation_point):
    matrix = []
    start = timer()
    b = []
    for x in collation_point:
        matrix.append([sum(c[j](x)*df_dx[j](f[i], x) for j in range(3)) for i in range(1, n)])
        b.append(c0(x)-sum(c[j](x)*df_dx[j](f[0], x) for j in range(3)))
    sol = np.linalg.solve(np.array(matrix), np.array(b))
    # print(sol)
    end = timer()
    result = end - start
    print("РњРµС‚РѕРґ РєРѕР»Р»РѕРєР°С†РёР№: {0} СЃРµРєСѓРЅРґ.".format(result))
    return lambda x: f[0](x) + sum(sol[i]*f[i+1](x) for i in range(n-1))

def part(q):
    return lambda x: sum(c[i](x)*df_dx[i](f[q], x) for i in range(3))

def mnk_int():
    matrix = [[integration(lambda x: part(i)(x)*part(j)(x), a0, b0)for j in range(1, n)]for i in range(1, n)]
    start = timer()
    b = [integration(lambda x: (c0(x)-part(0)(x))*part(i)(x), a0, b0) for i in range(1, n)]
	# print(matrix)
	# print(b)
    sol = np.linalg.solve(np.array(matrix), np.array(b))
    # print(sol)
    end = timer()
    result = end - start
    print("РРЅС‚РµРіСЂР°Р»СЊРЅС‹Р№ РјРµС‚РѕРґ РЅР°РёРјРµРЅСЊС€РёС… РєРІР°РґСЂР°С‚РѕРІ: {0} СЃРµРєСѓРЅРґ.".format(result))
    return lambda x: f[0](x) + sum(sol[i]*f[i+1](x) for i in range(n-1))

def mnk_discr(points):
    matrix = np.zeros((n-1, n-1))
    start = timer()
    b = np.zeros(n-1)
    for i in range(1, n):
        for j in range(1, n):
            der = lambda x: part(i)(x)*part(j)(x)
            matrix[i-1][j-1] = sum(der(x) for x in points)
        der = lambda x: (c0(x)-part(0)(x))*part(i)(x)
        b[i-1] = sum(der(x) for x in points)
    #print(matrix)
    #print(b)
    sol = np.linalg.solve(np.array(matrix), np.array(b))
    #print(sol)
    end = timer()
    result = end - start
    print("Р”РёСЃРєСЂРµС‚РЅС‹Р№ РјРµС‚РѕРґ РЅР°РёРјРµРЅСЊС€РёС… РєРІР°РґСЂР°С‚РѕРІ: {0} СЃРµРєСѓРЅРґ.".format(result))
    return lambda x: f[0](x) + sum(sol[i]*f[i+1](x) for i in range(n-1))

def galerkin():
    matrix = [[integration(lambda x: part(j)(x)*f[i](x), a0, b0) for j in range(1, n)] for i in range(1, n)]
    start = timer()
    #print(matrix)
    b = [integration(lambda x: (c0(x)-part(0)(x))*f[i](x), a0, b0) for i in range(1, n)]
    #print(b)
    sol = np.linalg.solve(np.array(matrix), np.array(b))
    #print(sol)
    end = timer()
    result = end - start
    print("РњРµС‚РѕРґ Р“Р°Р»С‘СЂРєРёРЅР°: {0} СЃРµРєСѓРЅРґ.".format(result))
    return lambda x: f[0](x) + sum(sol[i]*f[i+1](x) for i in range(n-1))

f = [lambda x: 0]
for i in range(1, n):
    f.append(gen_f(i))

x = np.arange(a0, b0+0.01, 0.01)
plt.plot(x, collation(np.linspace(a0, b0, n-1))(x), 'r')
plt.show()
plt.plot(x, galerkin()(x), 'b')
plt.show()
plt.plot(x, mnk_discr(np.linspace(a0, b0, n + 2))(x), 'y')
plt.show()
plt.plot(x, mnk_int()(x), 'g')
plt.show()
print()
print()
#plt.plot(x, collation(np.linspace(a0, b0, n-1))(x), 'r')
#plt.plot(x, mnk_int()(x), 'g')
#plt.plot(x, mnk_discr(np.linspace(a0, b0, n + 2))(x), 'y')
#plt.plot(x, galerkin()(x), 'b')
#plt.grid()
plt.show()