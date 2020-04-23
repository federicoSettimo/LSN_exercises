import matplotlib
import matplotlib.pyplot as plt
import numpy as np

filein = open("output_temp.dat")
N = int(int(filein.readline())/10)
t = np.zeros(N)
for i in range(N):
    t[i] = float(filein.readline()) - 1.4

x = range(N)
plt.plot(x,t)
plt.ylim(-1,1)
plt.show()
