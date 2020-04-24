import matplotlib
import matplotlib.pyplot as plt
import numpy as np


filein = open("output_temp.dat")
#N = int(int(filein.readline())/10)
N = 1000
y = np.zeros(N)
x = np.arange(N)
for i in range(N):
    y[i] = float(filein.readline())
filein.close()
plt.plot(x, y, label=r'$T$')

filein = open("output_ekin.dat")
filein.readline()
for i in range(N):
    y[i] = float(filein.readline())
filein.close()
plt.plot(x, y, label=r'$K$')

filein = open("output_epot.dat")
filein.readline()
for i in range(N):
    y[i] = float(filein.readline())
filein.close()
plt.plot(x, y, label=r'$U$')

filein = open("output_etot.dat")
filein.readline()
for i in range(N):
    y[i] = float(filein.readline())
filein.close()
plt.plot(x, y, label=r'$E$')



plt.xlabel('step')
plt.grid(True)
plt.legend(loc='lower right')
plt.show()
