import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

######################################
# Timings

time_full = 17+2+2+12+21+.72

time_2 = 6872.66 / 1000
time_3 = 39762.2 / 1000
time_4 = 99260.2 / 1000

plt.figure(figsize = (5, 3))
plt.semilogx([1e-2, 1e-3, 1e-4], [time_2, time_3, time_4], 'ko-', linewidth=1)
plt.hlines(time_full, 0, 1, color='k', linestyle="dashed", linewidth=1)
plt.suptitle("Average time of a step", fontsize=15)
plt.xlabel("Rounding error")
plt.xticks([1e-2, 1e-3, 1e-4], [r'$10^{-2}$', r'$10^{-3}$', r'$10^{-4}$'])
plt.ylabel("Time, s")
plt.grid()
plt.savefig('timings.png', dpi=300)
plt.close()

print("Time full:", time_full, "speed-up", time_full / time_full)
print("Time 1e-2:", time_2, "speed-up", time_full / time_2)
print("Time 1e-3:", time_3, "speed-up", time_full / time_3)
print("Time 1e-4:", time_4, "speed-up", time_full / time_4)

######################################
# Compression
filename_list = [""]

for filename in filename_list:
    data = np.loadtxt("./" + filename + "macro_restart.txt", skiprows=1)
    print(filename, "compression", 1/np.mean(data[:, -1]))

######################################

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

# x y z n T Px Py Pz Mx My Mz type p_inf S_inf
cp_list = []
cf_list = []
ch_list = []
n_list = {}
T_list = {}

for filename in filename_list:
    data = np.loadtxt("./" + filename + "boundary.txt", skiprows=1)
    print(data.shape)

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

    r = np.sqrt(x**2 + y**2)
    R = np.mean(r)
    print("mean radius", R)

    normal = np.array([x/r, y/r, 0*z])

    Pn = normal[0, :] * Px + normal[1, :] * Py + normal[2, :] * Pz
    Pt = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 - Pn ** 2)
    En = normal[0, :] * Mx + normal[1, :] * My + normal[2, :] * Mz

    cp = (Pn - p_inf) / (S_inf ** 2)
    cf = (Pt) / (S_inf ** 2)
    ch = 2 * (En) / (S_inf ** 3)

    cp = abs(cp)
    cf = abs(cf)
    ch = abs(ch)
    
    cp_list.append(cp)
    cf_list.append(cf)
    ch_list.append(ch)

    angle = -np.arctan(y / x)
    angle = np.where(angle < 0, angle + np.pi, angle)
    distance_surface = R * angle
    angle *= 180.0 / np.pi
    
    ###################################

    data = np.loadtxt("./" + filename + "boundary.txt", skiprows=1)
    data = data[data[:, 11] == 4]
    data = data[np.argsort(data[:, 0])]
    data = data[data[:, 0] < 0]
    data = data[data[:, 0] > -5.5]

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    n = data[:, 3]
    T = data[:, 4]
    
    n_list[filename] = n
    T_list[filename] = T
    
# data_ref = np.loadtxt("./plot-data-T.csv", skiprows=1, delimiter=",")
# print("data ref T shape", data_ref.shape)
# xT_ref = data_ref[:, 0]
# T_ref  = data_ref[:, 1]

# data_ref = np.loadtxt("./plot-data-cp.csv", skiprows=1, delimiter=",")
# print("data ref cp shape", data_ref.shape)
# xcp_ref = data_ref[:, 0]
# cp_ref  = data_ref[:, 1]

# data_ref = np.loadtxt("./plot-data-cf.csv", skiprows=1, delimiter=",")
# print("data ref cf shape", data_ref.shape)
# xcf_ref = data_ref[:, 0]
# cf_ref  = data_ref[:, 1]

# data_ref = np.loadtxt("./plot-data-ch.csv", skiprows=1, delimiter=",")
# print("data ref ch shape", data_ref.shape)
# xch_ref = data_ref[:, 0]
# ch_ref  = data_ref[:, 1]

labels = [r"$\varepsilon=10^{-2}$", r"$\varepsilon=10^{-3}$", r"$\varepsilon=10^{-4}$"]

custom_cycler = (matplotlib.cycler(color=['b', 'orange', 'g']) + matplotlib.cycler(lw=[1, 2, 3]))
#matplotlib.rcParams['axes.prop_cycle'] = custom_cycler
linewidth = 2
markersize = 1

plt.figure(figsize = (7, 4))

plt.plot(angle, cp_list[0], '--', linewidth=linewidth, markersize=markersize, label=labels[0])
# plt.plot(angle, cp_list[1], '-.', linewidth=linewidth, markersize=markersize, label=labels[1])
# plt.plot(angle, cp_list[2], ':',  linewidth=linewidth, markersize=markersize, label=labels[2])
    
# plt.scatter(xcp_ref, cp_ref, s=50, c="none", marker='o', label="DSMC", zorder=10, edgecolors='r')
plt.xlabel("Angle")
plt.ylabel(r"$c_p$", fontsize=18)
plt.xticks([0, 30, 60, 90, 120, 150, 180])
plt.legend()
#plt.grid()
plt.savefig('cp.png', dpi=300)
plt.close()

plt.figure(figsize = (7, 4))
#fig, axs = plt.subplot(2, 2, figsize = (14, 8))

plt.plot(angle, cf_list[0], '--', linewidth=linewidth, markersize=markersize, label=labels[0])
# plt.plot(angle, cf_list[1], '-.', linewidth=linewidth, markersize=markersize, label=labels[1])
# plt.plot(angle, cf_list[2], ':',  linewidth=linewidth, markersize=markersize, label=labels[2])
    
# plt.scatter(xcf_ref, cf_ref, s=50, c="none", marker='o', label="DSMC", zorder=10, edgecolors='r')
plt.xlabel("Angle")
plt.ylabel(r"$c_f$", fontsize=18)
plt.xticks([0, 30, 60, 90, 120, 150, 180])
plt.legend()
#plt.grid()
plt.savefig('cf.png', dpi=300)
plt.close()

plt.figure(figsize = (7, 4))

plt.plot(angle, ch_list[0], '--', linewidth=linewidth, markersize=markersize, label=labels[0])
# plt.plot(angle, ch_list[1], '-.', linewidth=linewidth, markersize=markersize, label=labels[1])
# plt.plot(angle, ch_list[2], ':',  linewidth=linewidth, markersize=markersize, label=labels[2])
    
# plt.scatter(xch_ref, ch_ref, s=50, c="none", marker='o', label="DSMC", zorder=10, edgecolors='r')
plt.xlabel("Angle")
plt.ylabel(r"$c_h$", fontsize=18)
plt.xticks([0, 30, 60, 90, 120, 150, 180])
plt.legend()
#plt.grid()
plt.savefig('ch.png', dpi=300)
plt.close()

plt.figure(figsize = (5, 4))
for name, n in n_list.items():
    plt.plot(x, n, '--', linewidth=linewidth, label=name)
plt.suptitle(r'$n$', fontsize=20)
plt.xlabel("x")
plt.legend()
plt.grid()
plt.savefig('n.png', dpi=300)
plt.close()

plt.figure(figsize = (7, 4))
for i, T in enumerate(T_list.values()):
    plt.plot(x, T, '--', linewidth=linewidth, markersize=markersize, label=labels[i])
# plt.scatter(xT_ref, T_ref, s=50, c="none", marker='o', label="DSMC", zorder=10, edgecolors='r')
plt.xlabel("x") #, fontsize=18)
plt.ylabel("T") #, fontsize=18)
#plt.tight_layout()
plt.legend()
#plt.grid()
plt.savefig('T.png', dpi=300)
plt.close()



