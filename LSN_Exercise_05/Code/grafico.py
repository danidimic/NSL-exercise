import numpy as np
import matplotlib.pyplot as plt

n,x = np.loadtxt('data.out', delimiter='  ', unpack=True)

plt.hist(x, bins=80)
#plt.plot(n,x)
plt.show()
