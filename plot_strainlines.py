import numpy as np
import matplotlib.pyplot as plt
import glob

for file in glob.glob('strainlines/*'):
    data = np.genfromtxt(file, dtype=float)
    num_samples = data.shape[0] // 3
    data = data.reshape((num_samples, 3))
    # print(data)
    plt.plot(data[:, 0], data[:,1], marker='*')

plt.show()