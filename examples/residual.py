import sys
import numpy as np
from matplotlib import pyplot as plt

folder = sys.argv[1]

fig, ax = plt.subplots(figsize = (7, 5))

ax.plot(np.loadtxt(folder + '/res.txt'))

# data = np.loadtxt('../../PHD/timings.txt')
# print(data.shape)
#
# plt.plot(np.sum(data, axis=1))

ax.set_title('RHS')
plt.savefig('res.png')