import matplotlib.pyplot as plt
import numpy as np
import requests as rqst
import re
import scipy.special as sp
import os


if not os.path.isdir("results"):
    os.mkdir("results")
url = 'http://www.jenyay.net/uploads/Student/Modelling/task_02.txt'
task = rqst.get(url)
z = re.search(r'^9\..+', task.text, flags=re.M)
z1 = (z.group().split(';'))
D = float(z1[0].split('=')[1])
fmin = float(z1[1].split('=')[1])
fmax = float(z1[2].split('=')[1])

# открытие файла с результатами и создание заголовка
fil = open('results/task_02_4o-506c_Shemshura_09.txt', 'w', encoding='utf-8')
print('f[ГГц]\t    σ[м^2]', file=fil)

# Вычислительная часть
c = 3e8
r = D / 2
f = np.linspace(fmin, fmax + 1, 500)
σ1 = []
for fx in f:
    λ = c / fx
    k = 2 * np.pi / λ
    h = []
    b = [0]
    a = []
    σN = []
    n = 0
    while n < 70:
        h.append(
            sp.spherical_jn(n, k * r) + 1j * sp.spherical_yn(n, k * r))
        a.append(sp.spherical_jn(n, k * r) / h[n])
        n += 1
    n = 1
    while n < 70:
        b.append((k * r * sp.spherical_jn(n - 1, k * r) - n *
                  sp.spherical_jn(n, k * r)) / (k * r * h[n - 1] - n * h[n]))
        n += 1
    n = 1
    while n < 70:
        σN.append(((-1)**n) * (n + 1 / 2) * (b[n] - a[n]))
        n += 1
    σ = (λ**2 / np.pi) * (abs(np.sum(σN)))**2
    print(str(fx * 1e-9) + '\t ' + str(σ), file=fil)
    σ1.append(σ / (np.pi * (r**2)))
fil.close()

# Построение графика
plt.figure()
plt.plot(2 * np.pi * f * r / c, σ1)
plt.grid()
plt.ylabel('ЭПР')
plt.xlabel('2πr/λ')
plt.savefig('img.png')
plt.show()
