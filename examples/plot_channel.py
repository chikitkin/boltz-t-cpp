import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    params = np.loadtxt(sys.argv[1])

    fig, ax = plt.subplots(figsize = (10, 6))
    
    dist = (params[:, 0]**2 + params[:, 1]**2)**.5
    params = params[np.isclose(dist, np.amin(dist), rtol=1e-2)]
    
    params = params[params[:, 2].argsort()]
    
    z   = params[:, 2]
    n   = params[:, 3]
    uz  = params[:, 6]
    T   = params[:, 7]
    
    n  = (n  -  n[0]) / (n[0])
    uz = (uz - uz[0]) / (uz[0])
    T  = (T  -  T[0]) / (T[0])
    
    ax.plot(z,  n, 'o-', label="n")
    ax.plot(z, uz, 'o-', label="uz")
    ax.plot(z,  T, 'o-', label="T")
    
    ax.grid()
    plt.legend()
    plt.savefig("channel.png", dpi=200)
    
    
    params = np.loadtxt(sys.argv[1])
    
    
    
    
    # params = np.loadtxt(sys.argv[2], skiprows=1)
    # print(params.shape)
    # params = params[np.isclose(abs(params[:, -1]), 1, rtol=1e-2)]
    # print(params.shape)
    # params = params[np.argsort(params[:, 2])]
    # print("n faces", params.shape[0] / 75)
    
    # z_list = []
    # M_0 = (2. / np.pi) * (np.pi /4)
    # M_list = []
    
    # for i in range(33):
    #     z_tmp = 0.0
    #     M_tmp = 0.0
    #     for j in range(75):
    #         rho = params[i*75+j,  3]
    #         un  = params[i*75+j,  6]
    #         S   = params[i*75+j, -4]
    #         z_tmp += params[i*75+j, 2]
    #         M_tmp += np.sum(rho * un * S)
    #     z_list.append(z_tmp / 75)
    #     M_list.append(M_tmp / M_0)
    

    # fig, ax = plt.subplots(figsize = (10, 6))
    # ax.plot(z_list, M_list, 'o-')
    # ax.grid()
    
    # plt.xlabel('z')
    # plt.ylabel('Q')
    # plt.savefig("rate.png", dpi=200)