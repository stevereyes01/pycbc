import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.figure()
n_bins = 20
data = np.genfromtxt('test1.csv', delimiter =',',names=True)

a = data['time']

plt.hist(a, bins = n_bins)

plt.title('Distribution of times')
plt.xlabel('GPS times')

plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/timeDistr.png")
