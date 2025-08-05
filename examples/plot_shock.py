import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    data = np.loadtxt(sys.argv[1])

    fig, ax = plt.subplots(figsize = (10, 6))
    
    data = data[data[:, 0].argsort()]
    x  = data[:, 0]
    n  = data[:, 1]
    u  = data[:, 2]
    T  = data[:, 3]
    
    n  = (n - n[0]) / (n[-1] - n[0])
    u  = (u - u[-1]) / (u[0] - u[-1])
    T  = (T - T[0]) / (T[-1] - T[0])
    
    delta = (n[1:] - n[:-1]) / (x[1:] - x[:-1])
    print("delta =", np.max(delta))
    
    ax.plot(x, n,  'b-', label="Density")
    ax.plot(x, u, 'g-', label="Velocity")
    ax.plot(x, T,  'r-', label="Temperature")
    
    
    
    data = [-10.0, 1.000, 1.000, 1.000,
    -9.2, 1.000, 1.000, 1.000,
    -8.4, 1.000, 1.000, 1.000,
    -7.6, 1.000, 1.000, 1.020,
    -6.8, 1.000, 1.020, 1.040,
    -6.2, 1.000, 1.040, 1.079,
    -5.6, 1.000, 1.079, 1.139,
    -5.0, 1.003, 1.159, 1.298,
    -4.4, 1.006, 1.318, 1.596,
    -3.8, 1.011, 1.676, 2.252,
    -3.2, 1.025, 2.371, 3.563,
    -2.8, 1.042, 3.226, 5.153,
    -2.4, 1.071, 4.557, 7.598,
    -2.0, 1.121, 6.564, 11.174,
    -1.6, 1.206, 9.386, 15.924,
    -1.2, 1.350, 12.824, 21.210,
    -0.8, 1.587, 16.301, 25.443,
    -0.4, 1.945, 18.964, 27.112,
    0.0, 2.410, 20.415, 26.277,
    0.4, 2.896, 20.932, 24.449,
    0.8, 3.288, 21.051, 22.879,
    1.2, 3.542, 21.031, 21.905,
    1.6, 3.683, 20.971, 21.369,
    2.0, 3.753, 20.932, 21.110,
    2.4, 3.787, 20.912, 20.991,
    2.8, 3.804, 20.892, 20.932,
    3.2, 3.813, 20.892, 20.892,
    3.8, 3.818, 20.872, 20.872,
    4.4, 3.818, 20.872, 20.872,
    5.0, 3.821, 20.872, 20.872,
    5.6, 3.821, 20.872, 20.872,
    6.2, 3.821, 20.872, 20.872,
    6.8, 3.821, 20.872, 20.872,
    7.6, 3.821, 20.872, 20.872,
    8.4, 3.821, 20.872, 20.872,
    9.2, 3.821, 20.872, 20.872,
    10.0, 3.821, 20.872, 20.872]
    
    data = np.array(data).reshape((-1, 4))
    
    # x_ref = data[:, 0]
    # n_ref = data[:, 1]
    # # u_ref = data[:, 2]
    # T_ref = data[:, 2]
    
    # n_ref  = (n_ref - n_ref[0])  / (n_ref[-1] - n_ref[0])
    # # u_ref  = (u_ref - u_ref[-1]) / (u_ref[0] - u_ref[-1])
    # T_ref  = (T_ref - T_ref[0])  / (T_ref[-1] - T_ref[0])

    # ax.plot(x_ref, n_ref,  'ko', label="Reference Density")
    # # ax.plot(x_ref, u_ref,  'kD', label="Reference Velocity")
    # ax.plot(x_ref, T_ref,  'kx', label="Reference Temperature")
    
    ax.set_xlabel('x', fontsize=16)
    ax.grid()
    plt.legend()
    plt.savefig("T.pdf", dpi=200)
    
    
    fig, ax = plt.subplots(figsize = (10, 6))
    ax.plot((x[1:] + x[:-1]), delta, 'b-', label=r"$\frac{\partial n}{\partial x}$")
    ax.set_xlabel('x', fontsize=16)
    ax.grid()
    plt.legend()
    plt.savefig("stress.pdf", dpi=200)