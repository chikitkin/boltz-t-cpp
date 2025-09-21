import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    data = np.loadtxt("T.txt")

    fig, ax = plt.subplots(figsize = (10, 6))
    
    data = data[data[:, 2].argsort()]
    z  = data[:, 2]
    n  = data[:, 3]
    u  = data[:, 6]
    T  = data[:, 7]
    
    n  = (n - n[0] ) / (n[-1] - n[0] )
    u  = (u - u[-1]) / (u[0]  - u[-1])
    T  = (T - T[0] ) / (T[-1] - T[0] )
    
    lw = 3
    
    ax.plot(z, n,  'b-', linewidth=lw, label="Density")
    ax.plot(z, u,  'g-', linewidth=lw, label="Velocity")
    ax.plot(z, T,  'r-', linewidth=lw, label="Temperature")
    
    ax.set_xlabel('x', fontsize=16)
    ax.grid()
    plt.legend()
    plt.savefig("T.png", dpi=200)
    
    fig, ax = plt.subplots(figsize = (10, 6))
    
    ax.plot(np.diff(n), 'b-', linewidth=lw, label=r"$\frac{\partial n}{\partial x}$")
    
    ax.set_xlabel('x', fontsize=16)
    ax.grid()
    plt.legend(fontsize=30)
    plt.savefig("stress.png", dpi=200)