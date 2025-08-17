import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

data = np.loadtxt("./macro_restart.txt")

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

n = data[:, 3]
ux = data[:, 4]
uy = data[:, 5]
T = data[:, 7]

R = (ux**2 + uy**2) ** .5

plt.figure(figsize = (7, 7))
plt.scatter(x, y, c=n, s=30, cmap='coolwarm')
plt.colorbar()
plt.axis('equal')
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
plt.title("n", fontsize=20)
plt.savefig('n.png', dpi=300)
plt.close()

plt.figure(figsize = (7, 7))
plt.scatter(x, y, c=T, s=30, cmap='coolwarm')
plt.colorbar()
plt.axis('equal')
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
plt.title("T", fontsize=20)
plt.savefig('T.png', dpi=300)
plt.close()

plt.figure(figsize = (7, 7))
plt.quiver(x, y, ux, uy)
plt.axis('equal')
plt.xlabel("x", fontsize=18)
plt.ylabel("y", fontsize=18)
plt.title("u", fontsize=20)
plt.savefig('u.png', dpi=300)
plt.close()