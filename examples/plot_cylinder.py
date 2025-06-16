import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

# x y z n T Px Py Pz Mx My Mz type p_inf S_inf
data = np.loadtxt(sys.argv[1], skiprows=1)
print(data.shape)

Na = 6.02214129e+23
kB = 1.381e-23
Ru = 8.3144598

Mol = 40e-3
Rg = Ru / Mol
m = Mol / Na

g = 5.0 / 3.0
d = 3.418e-10

Pr = 2.0 / 3.0

C = 144.4
T_0 = 273.11
mu_0 = 2.125e-5

mu_suth = lambda T: mu_0 * ((T_0 + C) / (T + C)) * (pow(T / T_0, 3.0 / 2.0))
mu = lambda T: mu_suth(200.0) * (pow(T / 200.0, 0.734))

##################
l_s = 0.1524
n_s = 1.6967383084574999e+19
u_s = 2624.0
T_s = 200.0
##################
v_s = np.sqrt(2 * Rg * T_s)
rho_s = m * n_s
p_s = rho_s * Rg * T_s

f_s = n_s * v_s ** 3

mu_s = mu(T_s)

data = data[data[:, 11] == 3]

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

n = data[:, 3]
T = data[:, 4]

Px = data[:, 5]
Py = data[:, 6]
Pz = data[:, 7]

Mx = data[:, 8]# / (n_s * v_s ** 3)
My = data[:, 9]# / (n_s * v_s ** 3)
Mz = data[:, 10]# / (n_s * v_s ** 3)

p_inf = data[0, 12]
S_inf = data[0, 13]

r = np.sqrt(x**2 + y**2 + z**2)
R = np.mean(r)

normal = np.array([x/r, y/r, z/r])

Pn = normal[0, :] * Px + normal[1, :] * Py + normal[2, :] * Pz
Pt = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 - Pn ** 2)
En = normal[0, :] * Mx + normal[1, :] * My + normal[2, :] * Mz

cp = (Pn - p_inf) / (S_inf ** 2)
cf = (Pt) / (S_inf ** 2)
ch = 2 * (En) / (S_inf ** 3)

cp = abs(cp)
cf = abs(cf)
ch = abs(ch)

angle = -np.arctan(y / x)
angle = np.where(angle < 0, angle + np.pi, angle)
distance_surface = R * angle
angle *= 180.0 / np.pi

plt.figure(figsize = (10, 7))
# d = np.loadtxt('DONE/cp.txt', delimiter=', ', skiprows=1)
# plt.plot(d[:, 0], d[:, 1], 'k--', label='DSMC')
plt.plot(angle, cp, 'o', markersize=10)
plt.suptitle(r'$c_p$', fontsize=48)
plt.xlabel("Angle")
plt.savefig('cp.png')
plt.close()

plt.figure(figsize = (10, 7))
# d = np.loadtxt('DONE/cf.txt', delimiter=', ', skiprows=1)
# plt.plot(d[:, 0], d[:, 1], 'k--', label='DSMC')
plt.plot(angle, cf, 'o', markersize=10)
plt.suptitle(r'$c_f$', fontsize=48)
plt.xlabel("Angle")
plt.savefig('cf.png')
plt.close()

plt.figure(figsize = (10, 7))
# d = np.loadtxt('DONE/ch.txt', delimiter=', ', skiprows=1)
# plt.plot(d[:, 0], d[:, 1], 'k--', label='DSMC')
plt.plot(angle, ch, 'o', markersize=10)
plt.suptitle(r'$c_h$', fontsize=48)
plt.xlabel("Angle")
plt.savefig('ch.png')
plt.close()

data = np.loadtxt(sys.argv[1], skiprows=1)
data = data[data[:, 11] == 4]
data = data[data[:, 0] < 0]
data = data[data[:, 0] > -6]

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

n = data[:, 3]
T = data[:, 4]

plt.figure(figsize = (10, 7))
# d = np.loadtxt('DONE/ch.txt', delimiter=', ', skiprows=1)
# plt.plot(d[:, 0], d[:, 1], 'k--', label='DSMC')
plt.plot(x, T, 'o', markersize=10)
plt.suptitle(r'$T$', fontsize=48)
plt.xlabel("x")
plt.savefig('T.png')
plt.close()
