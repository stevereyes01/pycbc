import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math

n_bins = 20
data = np.genfromtxt('test1.txt', delimiter =', ',names=['m1','m2','s1a','s1b','s1c','s2a','s2b','s2c'])

c = data['s1a']
d = data['s1b']
e = data['s1c']
f = data['s2a']
g = data['s2b']
h = data['s2c']

s1=[]
s2=[]
for i in range(0,len(data['s1a'])):
    s1.append(math.sqrt(c[i]**2 + d[i]**2 + e[i]**2))
    s2.append(math.sqrt(f[i]**2 + g[i]**2 + h[i]**2))

fig, axes = plt.subplots(nrows=2, ncols=4)
ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7 = axes.flat

ax2.hist(c, bins = n_bins)
ax3.hist(d, bins = n_bins)
ax4.hist(e, bins = n_bins)
ax5.hist(f, bins = n_bins)
ax6.hist(g, bins = n_bins)
ax7.hist(h, bins = n_bins)

ax0.hist(s1, bins = n_bins)
ax1.hist(s2, bins = n_bins)

ax2.set_title('Spin 1 x')
ax3.set_title('Spin 1 y')
ax4.set_title('Spin 1 z')
ax5.set_title('Spin 2 x')
ax6.set_title('Spin 2 y')
ax7.set_title('Spin 2 z')

ax0.set_title('Spin 1 mag')
ax1.set_title('Spin 2 mag')

spin = 1
spinD = 1 + 0.2
tick = 1 

ax0.set_xticks(np.arange(-1*spin, spinD, tick))
ax0.set_xlim(-0.2, spinD)

ax1.set_xticks(np.arange(-1*spin, spinD, tick))
ax1.set_xlim(-0.2,spinD)

ax2.set_xticks(np.arange(-1*spin, spinD, tick))
ax2.set_xlim(-1*spinD,spinD)
ax3.set_xticks(np.arange(-1*spin, spinD, tick))
ax3.set_xlim(-1*spinD,spinD)
ax4.set_xticks(np.arange(-1*spin, spinD, tick))
ax4.set_xlim(-1*spinD,spinD)
ax5.set_xticks(np.arange(-1*spin, spinD, tick))
ax5.set_xlim(-1*spinD,spinD)
ax6.set_xticks(np.arange(-1*spin, spinD, tick))
ax6.set_xlim(-1*spinD,spinD)
ax7.set_xticks(np.arange(-1*spin, spinD, tick))
ax7.set_xlim(-1*spinD,spinD)

plt.tight_layout()
plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/spinsDistr.png")
