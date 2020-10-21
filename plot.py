import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Generate new colormap

# N = 256
# vals = np.ones((N, 4))
# vals[:, 0] = np.linspace(1, 0/256, N)
# vals[:, 1] = np.linspace(1, 92/256, N)
# vals[:, 2] = np.linspace(1, 132/256, N)
#vals[247:-1,0] = 0.
#vals[247:-1,1] = 92/256.
#vals[247:-1,2] = 132./256.
# newcmp = ListedColormap(vals)

ftle = np.genfromtxt('helicity.data', skip_header=0).reshape((500,500))
# ftleMask = np.genfromtxt('helicityMask.data', skip_header=0, dtype=int).reshape((500, 500))
# ftle = np.ma.array(ftle, mask=ftleMask == 1)
x = np.linspace(0, 2*np.pi,num=500)
y = np.linspace(0, 2*np.pi,num=500)
X, Y = np.meshgrid(x, y, indexing='ij')
ax = plt.subplot(111)
CS = ax.contourf(X, Y, ftle, 100, cmap=plt.cm.jet)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.colorbar(CS)
# Plot seed points
seeds = np.genfromtxt('bin/seedPoints.data', skip_header=0)
num_seeds = seeds.size//3
seeds = seeds.reshape( (num_seeds, 3) )

ax.plot(seeds[:, 0], seeds[:,1], 'ro')
plt.title('Helicity field; ABC flow, z = 0')
plt.show()
