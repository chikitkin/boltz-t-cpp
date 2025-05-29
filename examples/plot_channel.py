import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    coords = np.loadtxt(sys.argv[1])
    params = np.loadtxt(sys.argv[2])
    
    if params.ndim != 1:
        fig, ax = plt.subplots(figsize = (10, 6))
        
        dist = (coords[:, 0]**2 + coords[:, 1]**2)**.5
        coords = coords[dist == np.amin(dist)]
        params = params[dist == np.amin(dist)]
        
        params = params[coords[:, 2].argsort()]
        coords = coords[coords[:, 2].argsort()]
        
        z  = coords[:, 2]
        n  = params[:, 0]
        T  = params[:, 4]
        
        # n = (n - n[0]) / (n[-1] - n[0])
        # T = (T - T[0]) / (T[-1] - T[0])
        
        ax.plot(z, n, label="n")
        ax.plot(z, T, label="T")
        
        ax.grid()
        plt.legend()
        plt.savefig("channel.png", dpi=200)