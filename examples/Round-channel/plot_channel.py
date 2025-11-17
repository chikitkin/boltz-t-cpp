import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    params = np.loadtxt("./macro_restart.txt")

    fig, ax = plt.subplots(figsize = (10, 6))
    
    dist = (params[:, 0]**2 + params[:, 1]**2)**.5
    params = params[np.isclose(dist, np.amin(dist), rtol=1e-2)]
    
    params = params[params[:, 2].argsort()]
    
    z   = params[:, 2]
    n   = params[:, 3]
    uz  = params[:, 6]
    T   = params[:, 7]
    
    # n  = (n  -  n[0]) / (n[0])
    # uz = (uz - uz[0]) / (uz[0])
    # T  = (T  -  T[0]) / (T[0])
    
    ax.plot(z,  n, 'o-', label="n")
    ax.plot(z, uz, 'o-', label="uz")
    ax.plot(z,  T, 'o-', label="T")
    
    ax.grid()
    plt.legend()
    plt.savefig("channel.png", dpi=200)
    
    
    params = np.loadtxt("./macro_restart_1e-3.txt")
    
    Mol = 40e-3
    Na = 6.02214129e+23
    m = Mol / Na
    Ru = 8.3144598
    Rg = Ru / Mol
    
    n_s = 1.5217740541258367e+20
    v_s = (2 * Rg * 200) ** .5
    
    x = params[:, 0]
    y = params[:, 1]
    z = params[:, 2]
    n = params[:, 3]
    
    ux = params[:, 4]
    uy = params[:, 5]
    uz = params[:, 6]
    
    T = params[:, 7]
    V = params[:, 8]
    
    compression = params[:, 9]
    print("compression", np.mean(compression))
    
    rho = n
    uz = uz
    z = z.round(decimals=4)
    z_unique = np.unique(z)
    
    Q_0 = .5 * np.pi ** -.5
    
    Q = []
    
    for z_ in z_unique:
        rho_ = rho[z == z_]
        uz_  = uz [z == z_]
        V_   = V  [z == z_]
        
        Q_ = np.sum(rho_ * uz_ * V_) / np.sum(V_)
    
        Q.append(Q_)

    Q = np.array(Q)
    Q_norm = Q / Q_0
    
    N = len(Q)
    
    print((Q_norm[N // 2] + Q_norm[(N + 1) // 2]) / 2)
    print(np.mean(Q_norm))

    fig, ax = plt.subplots(figsize = (10, 6))
    ax.plot(z_unique, Q_norm, 'o-')
    ax.grid()
    
    plt.xlabel('z')
    plt.ylabel('Q')
    plt.savefig("rate.png", dpi=200)
    
    