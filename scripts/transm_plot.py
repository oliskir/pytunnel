import matplotlib.pyplot as plt
import numpy as np

ylog = True


x, y = np.loadtxt('out_5fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -5$ fm')

x, y = np.loadtxt('out_10fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -10$ fm')

x, y = np.loadtxt('out_20fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -20$ fm')

x, y = np.loadtxt('out_30fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -30$ fm')

x, y = np.loadtxt('out_20fm_4fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -20$ fm, $\sigma=4$ $\mu$m')

x, y = np.loadtxt('out_30fm_6fm.csv', delimiter=',', unpack=True)
plt.plot(x,y,label=r'$x_0 = -30$ fm, $\sigma=6$ $\mu$m')


dy = np.max(y) - np.min(y)
if ylog:
    ymin = 0.5 * np.min(y)
    ymax = 10 * np.max(y)
else:
    ymin = np.min(y) - 0.2 * dy
    ymax = np.max(y) + 0.2 * dy

plt.xlabel('Classical position when barrier height is increased (fm)')

plt.ylabel('Transmission probability')
plt.title('Transmission of Gaussian wave-packet')

if ylog: plt.yscale('log')
plt.ylim(ymin, ymax)
plt.grid(True)

plt.legend(loc='lower right')

plt.show()
