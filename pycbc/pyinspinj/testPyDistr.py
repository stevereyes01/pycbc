import pycbc.inference.distributions as distr
import injection_params as inj

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict

n_bins = 20
numberOfInj = 10000
mBounds1 = distr.Uniform(mass1=(1.0,3.0)) 
mBounds2 = distr.Uniform(mass2=(4.0,5.0))
#mDist = mBounds.rvs(size=numberOfInj)

#for i in range(0,numberOfInj):
#    print mDist[i]

#myInjection = inj.InjParams()
#a=[]
#b=[]
#c=[]
#for i in range(0,numberOfInj):
#    mDist = mBounds.rvs(size=1)
#    myInjection.mass1 = (mDist[0])[0]
#    mDist = mBounds.rvs(size=1)
#    myInjection.mass2 = (mDist[0])[0]
#    a.append(myInjection)
#    print a[i].mass1+1.0, a[i].mass2+1.0
#    b.append(a[i].mass1+1.)
#    c.append(a[i].mass2+1.)

injDict = defaultdict(list)

mDist1 = mBounds1.rvs(size=numberOfInj)
mDist2 = mBounds2.rvs(size=numberOfInj)

injDict['mass1'] = mDist1
injDict['mass2'] = mDist2

print type(mDist1)
print mDist1
print mDist1[0][0]
print mDist2[5][0]

mass1=[]
mass2=[]
for i in range(0,numberOfInj):
    mass1.append( injDict['mass1'][i][0])
    mass2.append(injDict['mass2'][i][0])

fig, axes = plt.subplots(nrows=1, ncols=2)
ax0, ax1 = axes.flat

ax0.hist(mass1, bins = n_bins)
ax1.hist(mass2, bins = n_bins)

ax0.set_title('Distribution of mass 1 masses')
ax1.set_title('Distribution of mass 2 masses')

plt.savefig("/home/steven.reyes/public_html/ER10/py_inspinj/massDistrTest1.png")

