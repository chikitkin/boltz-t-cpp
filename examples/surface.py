import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

# x y z n T rho p Px Py Pz Mx My Mz 
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

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

n = data[:, 3]
T = data[:, 4]

rho = data[:, 5]# / rho_s
p   = data[:, 6]# / p_s

Px = data[:, 7]
Py = data[:, 8]
Pz = data[:, 9]

Mx = data[:, 10]# / (n_s * v_s ** 3)
My = data[:, 11]# / (n_s * v_s ** 3)
Mz = data[:, 12]# / (n_s * v_s ** 3)

r = np.sqrt(x**2 + y**2 + z**2)
R = np.mean(r)

normal = np.array([x/r, y/r, z/r])

Pn = normal[0, :] * Px + normal[1, :] * Py + normal[2, :] * Pz
Pt = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 - Pn ** 2)
En = normal[0, :] * Mx + normal[1, :] * My + normal[2, :] * Mz
Et = np.sqrt(Mx ** 2 + My ** 2 + Mz ** 2 - En ** 2)

cp = 2 * (Pn) / (rho_s * u_s ** 2)
cf = 2 * (Pt) / (rho_s * u_s ** 2)
ch = 2 * (En) / (rho_s * u_s ** 3)

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")

im = ax.scatter3D(x, y, z, c=cf)
ax.set_title(r'$c_f$', fontsize=48)

ax.view_init(35, 135, 0)
ax.set_aspect('equal')
fig.colorbar(im)

plt.savefig('surface.png')