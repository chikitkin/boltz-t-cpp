import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    data = np.loadtxt(sys.argv[1])
    
    if data.ndim != 1:
        fig, ax = plt.subplots(figsize = (10, 6))
        
        data = data[data[:, 0].argsort()]
        x = data[:, 0]
        n = data[:, 1]
        ux = data[:, 2]
        T = data[:, 3]
        
        n = (n - n[0]) / (n[-1] - n[0])
        ux = (ux - ux[-1]) / (ux[0] - ux[-1])
        T = (T - T[0]) / (T[-1] - T[0])
        
        ax.plot(x, n, label="n")
        ax.plot(x, ux, label="ux")
        ax.plot(x, T, label="T")
        
        ax.grid()
        plt.legend()
        plt.savefig("T.png", dpi=200)
        
    else:
        plt.figure(figsize=(12, 10))
        plt.semilogy(data)
        plt.title('res', fontsize=60)
        plt.savefig("res.png", dpi=200)
