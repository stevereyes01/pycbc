import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#plt.figure()
n_bins = 20
data = np.genfromtxt('test1.csv', delimiter =',',names=True)

a = data['distance']

fig, axes = plt.subplots(nrows=1, ncols=2)
ax0, ax1 = axes.flat

ax0.hist(a, bins = n_bins)
ax1.hist(a, bins = n_bins)

ax0.set_title('Distribution of distances')
ax1.set_title('Distribution of distances')

plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/distDistr.png")
