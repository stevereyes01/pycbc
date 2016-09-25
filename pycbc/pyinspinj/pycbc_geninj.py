import sys
import argparse
import logging
import time
import numpy as np
import math
from pycbc.inference import distributions as distr
from collections import defaultdict
import csv

'''Functions for py_inspinj'''

'''GPS and Time functions'''

# Draw an injection time between GPS START and GPS END 


# Draw an injection time at specified step interval between GPS START & END

'''Source Distribution functions'''

# Draw a distance via uniform in distance (Mpc)
def uniformDist(lowerBound, upperBound, num_injections):
    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i] = np.random.uniform(lowerBound, upperBound)   
    return array1

# Draw a distance via uniform in distance squared (Mpc**2)
def uniformArea(lowerBound, upperBound, num_injections):
    shellMin = lowerBound**2
    shellMax = upperBound**2

    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i] = math.sqrt(np.random.uniform(shellMin, shellMax))
    return array1
 
# Draw a distance via uniform in volume (distance cubed) (Mpc**3)
def uniformVolume(lowerBound, upperBound, num_injections):
    sphereMin = lowerBound**3
    sphereMax = upperBound**3

    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i] = pow(np.random.uniform(sphereMin, sphereMax),(1./3))
    return array1

# Draw a distance via uniform in log 10 distance (Mpc)
def log10Dist(lowerBound, upperBound, num_injections):
# Logarithms return NaN if you're at 0!
    if lowerBound <= 0:
       lowerBound = 0.000000000001

    logDMin = math.log10(lowerBound)
    logDMax = math.log10(upperBound)

    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i] = pow(10, np.random.uniform(logDMin, logDMax))
    return array1

# Draw a distance uniform in Chirp Distance

# Draw a sky location for the binary
#def uniformTheta(lowerBound, upperBound):

#def uniformPhi(lowerBound, upperBound):
# Draw a polarization for the binary

# Draw an inclination for the binary

# Draw a coalescence phase for the binary

'''Mass Distribution functions'''

# Parse a Uniform Distribution from distributions.py to simple numpy array
def parseUniformDistr(numpyDumpy, num_injections):
    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i]=numpyDumpy[i][0]
    return array1

# Parse a Gaussian Distribution from distributions.py to a simple numpy array
def parseGaussDistr(numpyDumpy, num_injections):
    array1 = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        array1[i]=numpyDumpy[i][0]
    return array1

'''Spin Distribution functions'''

# Draw a random spin from a uniform spin distribution
def uniformSpin(lowerBound, upperBound):
    spinSample = np.random.uniform(lowerBound, upperBound)
    return spinSample

# Check that the spin magnitude is less than 1
def checkSpinMagBig(spin1, spin2, spin3):
    spinMag = math.sqrt(spin1**2 + spin2**2 + spin3**2)

    if spinMag < 1 :
       return False
   
    return True

'''Tapering Injected Waveform functions'''


parser = argparse.ArgumentParser(description='A mock up of lalapps_inspinj' \
                                             'for Python.')
# GPS Inputs and Time Options
parser.add_argument('--gps-start-time', type=float, required=False, default=0,
                     help='Optional: A beginning GPS time where the injection' \
                          ' can be bounded by. (default: 0)')
parser.add_argument('--gps-end-time', type=float, required=False, default=100,
                     help='Optional: An end GPS time where the injection can' \
                          ' be bounded by. (default: 100)')

parser.add_argument('--time-step', type=float, required=False,
                     default=2630/math.pi, help='Optional: Sets the interval' \
                    ' at which injections will uniformly be distributed' \
                    ' across (default: 2630/pi)')

# Detector PSD Inputs
#parser.add_arguments('--L1-PSD', type, required=True,
#                      help='Required: A Power Spectral Density Ascii text' \
#                           ' file for the Livingston detector.')

#parser.add_arguments('--H1-PSD', type, required=True,
#                      help='Required: A Power Spectral Density Ascii txt' \
#                           ' file for the Hanford detector.')

# Source Distribution Options
parser.add_argument('--num-injections', type=int, required=False, default=10000,
                     help='Optional: The number of injections to create' \
                          ' (default: 100)')

parser.add_argument('--dist-distr', type=str, required=False,
                    default='uniform', help='Choose distance distribution as ' \
                                            '<uniform>, <uniform_area>, ' \
                                            '<uniformVolume>, or ' \
                                            '<uniformLog10>. For uniform in ' \
                                            'distance, uniform in area, ' \
                                            'uniform in volume, or uniform ' \
                                            'uniform in log10 distance.')

parser.add_argument('--dist-max', type=float, required=False, default=1000,
                     help='Optional: Choose the maximum distance of ' \
                          'injections in Mpc. (default: 1000)')

parser.add_argument('--dist-min', type=float, required=False, default=0,
                    help='Optional: Choose the minimum distance of ' \
                         'injections in Mpc. (default: 0)')

#parser.add_argument(RIGHT ASCENSION BOUNDS)

#parser.add_argument(DECLINATION BOUNDS)

# Mass Distribution Options
parser.add_argument('--min-mass1', type=float, required=False, default=1.0,
                     help='Optional: Minimum mass for mass 1 distribution' \
                          ' (default:1.0)')

parser.add_argument('--max-mass1', type=float, required=False, default=3.0,
                     help='Optional: Maximum mass for mass 1 distribution' \
                          ' (default:3.0)')

parser.add_argument('--mean-mass1', type=float, required=False, default=1.35,
                     help='Optional: Mean mass for mass 1 distribution' \
                          ' (default:1.35)')

parser.add_argument('--stdev-mass1', type=float, required=False, default=0.13,
                     help='Optional: Standard deviation for mass 1 distribution' \
                          ' (default:0.13)')

parser.add_argument('--mass1-uniform', type=int, required=False, default=0,
                     help='Optional: Set mass 1 distribution to uniform' \
                          ' (default:0, False)')

parser.add_argument('--min-mass2', type=float, required=False, default=1.0,
                     help='Optional: Minimum mass for mass 2 distribution' \
                          ' (default:1.0)')

parser.add_argument('--max-mass2', type=float, required=False, default=3.0,
                     help='Optional: Maximum mass for mass 2 distribution' \
                          ' (default:3.0)')

parser.add_argument('--mean-mass2', type=float, required=False, default=1.35,
                     help='Optional: Mean mass for mass 2 distribution' \
                          ' (default:1.35)')

parser.add_argument('--stdev-mass2', type=float, required=False, default=0.13,
                     help='Optional: Standard deviation for mass 2 distribution' \
                          ' (default:0.13)')

parser.add_argument('--mass2-uniform', type=int, required=False, default=0,
                     help='Optional: Set mass 2 distribution to uniform' \
                          ' (default:0, False)')

# Spin Distribution Options
parser.add_argument('--min-spin1', type=float, required=False, default=-1.0,
                     help='Optional: Minimum spin for spin1 (default:-1.0)')

parser.add_argument('--max-spin1', type=float, required=False, default=1.0,
                     help='Optional: Maximum spin for spin1 (default:1.0)')

parser.add_argument('--min-spin2', type=float, required=False, default=-1.0,
                     help='Optional: Minimum spin for spin2 (default:-1.0)')

parser.add_argument('--max-spin2', type=float, required=False, default=1.0,
                     help='Optional: Maximum spin for spin2 (default:1.0)')

# Tapering the Injected Waveform
# ?????

# Output & Misc Options
parser.add_argument('--seed', required=False, default=1,
                     help='Optional: The seed to feed to the random number' \
                          ' generator (default:1')

parser.add_argument('--output', required=True,
                     help='Required: The file to output the injections to.')

parser.add_argument('--verbose', action='store_true', help='Optional: Give a' \
                    ' verbose output of operations.')

opts = parser.parse_args()

# Logging information.
log_fmt = '%(asctime)s %(message)s'
log_date_fmt = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(level=logging.INFO, format=log_fmt, datefmt=log_date_fmt)

logging.info('Generating Injections...')

# Create a dictionary for all of the parameters of the injections

injDict = defaultdict(list)

##############################################################################
#                      Write injections to a dictionary                      #
##############################################################################

# Mass 1
if opts.mass1_uniform == True:
   mBounds1 = distr.Uniform(mass1=(opts.min_mass1,opts.max_mass1))
   injDict['mass1'] = parseUniformDistr(mBounds1.rvs(size=opts.num_injections),
                                        opts.num_injections)
else :
#   mBounds1 = distr.Gaussian('mass1',opts.min_mass1,opts.max_mass1,
#                             opts.mean_mass1, opts.stdev_mass1**2)
   mBounds1 = distr.Gaussian(['mass1'],[opts.min_mass1],[opts.max_mass1],
                            [opts.mean_mass1],[opts.stdev_mass1**2])
   injDict['mass1'] = parseGaussDistr(mBounds1.rvs(size=opts.num_injections),
                                      opts.num_injections)

logging.info('Mass 1 parameters written.')

# Mass 2
if opts.mass2_uniform == True:
   mBounds2 = distr.Uniform(mass2=(opts.min_mass2,opts.max_mass2))
   injDict['mass2'] = parseUniformDistr(mBounds2.rvs(size=opts.num_injections),
                                        opts.num_injections)
else :
#   mBounds2 = distr.Gaussian(mass1=(opts.min_mass2,opts.max_mass2,
#                             opts.mean_mass2, opts.stdev_mass2**2))
   mBounds2 = distr.Gaussian(['mass2'],[opts.min_mass2],[opts.max_mass2],
                            [opts.mean_mass2],[opts.stdev_mass2**2])
   injDict['mass2'] = parseGaussDistr(mBounds2.rvs(size=opts.num_injections),
                                      opts.num_injections)

logging.info('Mass 2 parameters written.')

# MTotal = m1 + m2
injDict['mtotal'] = injDict['mass1'] + injDict['mass2']
logging.info('Total mass written.')

# MChirp = (m1*m2)**(3./5.) / (m1 + m2)**(1./5.)
injDict['mchirp'] = ((injDict['mass1']*injDict['mass2'])**(3./5.)) / \
                     ((injDict['mass1'] + injDict['mass2'])**(1./5.))

logging.info('Chirp mass written.')

# Eta
injDict['eta'] = (injDict['mass1']*injDict['mass2']) / \
                 ((injDict['mass1'] + injDict['mass2'])**2)

logging.info('Eta written.')

# q mass ratio (Convention: write it as bigger mass divided by smaller mass!)

for i in range(0,opts.num_injections):
    if injDict['mass1'][i] > injDict['mass2'][i]:
       injDict['q'].append(injDict['mass1'][i]/injDict['mass2'][i])
    else:
       injDict['q'].append(injDict['mass2'][i]/injDict['mass1'][i])
    if injDict['q'][i] < 1.0:
       print 'Q is too small!'

# Spin 1
spin1x = np.ndarray(shape=(opts.num_injections), dtype=float)
spin1y = np.ndarray(shape=(opts.num_injections), dtype=float)
spin1z = np.ndarray(shape=(opts.num_injections), dtype=float)

for i in range(0, opts.num_injections):
    spinDistr1 = distr.Uniform(spin1=(opts.min_spin1, opts.max_spin1))
    s1x = spinDistr1.rvs(size=1)[0][0]
    s1y = spinDistr1.rvs(size=1)[0][0]
    s1z = spinDistr1.rvs(size=1)[0][0]

    while checkSpinMagBig(s1x,s1y,s1z) == True:
          s1x = spinDistr1.rvs(size=1)[0][0]
          s1y = spinDistr1.rvs(size=1)[0][0]
          s1z = spinDistr1.rvs(size=1)[0][0] 
    
    spin1x[i] = s1x
    spin1y[i] = s1y
    spin1z[i] = s1z

injDict['spin1x'] = spin1x
injDict['spin1y'] = spin1y
injDict['spin1z'] = spin1z

logging.info('Spin 1 parameters written.')

# Spin 2
spin2x = np.ndarray(shape=(opts.num_injections), dtype=float)
spin2y = np.ndarray(shape=(opts.num_injections), dtype=float)
spin2z = np.ndarray(shape=(opts.num_injections), dtype=float)

for i in range(0, opts.num_injections):
    spinDistr2 = distr.Uniform(spin2=(opts.min_spin2, opts.max_spin2))
    s2x = spinDistr2.rvs(size=1)[0][0]
    s2y = spinDistr2.rvs(size=1)[0][0]
    s2z = spinDistr2.rvs(size=1)[0][0]

    while checkSpinMagBig(s2x,s2y,s2z) == True:
          s2x = spinDistr2.rvs(size=1)[0][0]
          s2y = spinDistr2.rvs(size=1)[0][0]
          s2z = spinDistr2.rvs(size=1)[0][0]                 

    spin2x[i] = s2x
    spin2y[i] = s2y
    spin2z[i] = s2z

injDict['spin2x'] = spin2x
injDict['spin2y'] = spin2y
injDict['spin2z'] = spin2z

logging.info('Spin 2 parameters written.')

# Write the projection of the spin on the angular momentum axis
# x1 = (c*vec{S1}/(G m1**2) dot vec{L}

# Write effective spin to the dictionary (x1m1 + x2m2)/(m1+m2)

#injDict['chi_eff'] = 

# Draw the distances to the binaries
if opts.dist_distr == 'uniform':
    injDict['distance'] = uniformDist(opts.dist_min, opts.dist_max,
                                          opts.num_injections)
   
elif opts.dist_distr == 'uniformArea':
    injDict['distance'] = uniformArea(opts.dist_min, opts.dist_max,
                                      opts.num_injections)

elif opts.dist_distr == 'uniformVolume':
    injDict['distance'] = uniformVolume(opts.dist_min, opts.dist_max,
                                        opts.num_injections)
elif opts.dist_distr == 'uniformLog10':
    injDict['distance'] = log10Dist(opts.dist_min, opts.dist_max,
                                    opts.num_injections)

# Throw an error if the input doesn't look good
else:
    raise ValueError('Unrecognized distance distribution! Could not find %s.' \
                     'Please try uniform, uniformArea, uniformVolume, or ' \
                     'uniformLog10 distance distributions.' %(opts.dist_distr))

logging.info('Distance parameters written.')

# Draw the sky location of a binary (Right Ascension first)
# (Declination second)
# In radians!

ra = np.ndarray(shape=(opts.num_injections), dtype=float)
dec = np.ndarray(shape=(opts.num_injections), dtype=float)

# Sampling dec in -pi/2 to pi/2 gives a uniform
# distribution over S2 surface of sphere. Sampling
# from 0 to pi gives bunching at the poles.
for i in range(0,opts.num_injections):
    ra[i] = 2*np.pi*np.random.uniform(0, 1)
    dec[i] = np.arccos(2*np.random.uniform(0, 1)-1)

injDict['ra'] = ra
injDict['dec'] = dec

# Draw the GPS times for each injection

deltaT = opts.gps_end_time - opts.gps_start_time

time = np.ndarray(shape=(opts.num_injections), dtype=int)
for i in range(0, opts.num_injections):
    time [i] = np.random.uniform(opts.gps_start_time, opts.gps_end_time)
injDict['time'] = time

logging.info('Time stamps for injections written.')

with open(opts.output, "wb") as outfile:
   writer = csv.writer(outfile)
   writer.writerow(injDict.keys())
   writer.writerows(zip(*injDict.values()))
    
logging.info('100% done')
