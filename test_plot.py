import numpy as np
import matplotlib.pyplot as plt 

data = np.genfromtxt('helicity.data')
x = np.linspace(0, 2. * np.pi, data.size) 

plt.plot(x, data)
plt.ylim((0., .32))
plt.show()

print("The minimum point is " + str(np.min(data)))