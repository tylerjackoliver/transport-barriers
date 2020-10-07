import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Generate new colormap

N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(1, 0/256, N)
vals[:, 1] = np.linspace(1, 92/256, N)
vals[:, 2] = np.linspace(1, 132/256, N)
#vals[247:-1,0] = 0.
#vals[247:-1,1] = 92/256.
#vals[247:-1,2] = 132./256.
newcmp = ListedColormap(vals)

ftle = np.genfromtxt('helicity.data', skip_header=0).reshape((400, 400))
ftleMask = np.genfromtxt('helicityMask.data', skip_header=0, dtype=int).reshape((400, 400))
ftle = np.ma.array(ftle, mask=ftleMask == 1)
x = np.linspace(0, 2*np.pi,num=400)
y = np.linspace(0, 2*np.pi,num=400)
X, Y = np.meshgrid(x, y, indexing='ij')
ax = plt.subplot(111)
ax.contourf(X, Y, ftle, 100, cmap=newcmp)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
