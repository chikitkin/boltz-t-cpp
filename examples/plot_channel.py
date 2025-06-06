import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    params = np.loadtxt(sys.argv[1])

    fig, ax = plt.subplots(figsize = (10, 6))
    
    dist = (params[:, 0]**2 + params[:, 1]**2)**.5
    params = params[dist == np.amin(dist)]
    
    params = params[params[:, 2].argsort()]
    
    z  = params[:, 2]
    n  = params[:, 0+3]
    T  = params[:, 4+3]
    
    # n = (n - n[0]) / (n[-1] - n[0])
    # T = (T - T[0]) / (T[-1] - T[0])
    
    ax.plot(z, n, label="n")
    ax.plot(z, T, label="T")
    
    ax.grid()
    plt.legend()
    plt.savefig("channel.png", dpi=200)