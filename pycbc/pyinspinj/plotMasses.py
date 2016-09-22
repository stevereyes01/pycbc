import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#plt.figure()
n_bins = 20
data = np.genfromtxt('test1.csv', delimiter =',',names=True)

a = data['mass1']
b = data['mass2']
c = data['mchirp']
d = data['eta']
e = data['q']

fig, axes = plt.subplots(nrows=3, ncols=2)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flat

ax0.hist(a, bins = n_bins)
ax1.hist(b, bins = n_bins)
ax2.hist(c, bins = n_bins)
ax3.hist(d, bins = n_bins)
ax4.hist(e, bins = n_bins)
ax5.hist(e, bins = n_bins)

ax0.set_title('Distribution of mass 1 masses')
ax1.set_title('Distribution of mass 2 masses')
ax2.set_title('Distribution of chirp masses')
ax3.set_title('Distribution of eta')
ax4.set_title('Distribution of q')
ax5.set_title('Distribution of q')

plt.tight_layout()
plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/massDistr.png")
