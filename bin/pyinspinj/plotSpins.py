import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math

n_bins = 25
data = np.genfromtxt('test1.csv', delimiter =',',names=True)

c = data['spin1x']
d = data['spin1y']
e = data['spin1z']
f = data['spin2x']
g = data['spin2y']
h = data['spin2z']

s1=[]
s2=[]

file1 = open('spin1.txt', 'w')
file2 = open('spin2.txt', 'w')

for i in range(0,len(data['spin1x'])):
    string1 = str(c[i]) + ' ' + str(d[i]) + ' ' + str(e[i]) + '\n'
    string2 = str(f[i]) + ' ' + str(g[i]) + ' ' + str(h[i]) + '\n'

    file1.write(string1)
    file2.write(string2)

    s1.append(math.sqrt(c[i]**2 + d[i]**2 + e[i]**2))
    s2.append(math.sqrt(f[i]**2 + g[i]**2 + h[i]**2))

file1.close()
file2.close()

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
#plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/spinsDistr.png")
plt.savefig("/home/steven.reyes/WWW/LSC/ER10/spinsDistr.png")
