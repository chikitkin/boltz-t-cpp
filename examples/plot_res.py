import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == "__main__":

    data = np.loadtxt("res.txt")
    
    plt.figure(figsize=(12, 10))
    plt.semilogy(data)
    plt.title('res', fontsize=60)
    plt.savefig("res.png", dpi=200)
    
    print("Time:", np.sum(np.loadtxt("timings.txt")), "seconds")
