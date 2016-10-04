import sys
import argparse
import logging
import time
import numpy as np
import math
from pycbc.inference import distributions as distr
from collections import defaultdict
import pycbc.pnutils as pnu
import csv
import pycbc.coordinates as coord

'''GPS and Time functions'''

# Draw an injection time between GPS START and GPS END 
def uniformTime(lowerBound, upperBound, num_injections):
    deltaT = opts.gps_end_time - opts.gps_start_time

    time = np.ndarray(shape=(opts.num_injections), dtype=int)
    for i in range(0, opts.num_injections):
        time[i] = np.random.uniform(opts.gps_start_time, opts.gps_end_time)
    return time
# Draw an injection time at specified step interval between GPS START & END

'''Source Distribution functions'''

# Draw a distance via uniform in distance (Mpc)
def uniformDist(lowerBound, upperBound, num_injections):
    distance = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        distance[i] = np.random.uniform(lowerBound, upperBound)   
    return distance

# Draw a distance via uniform in distance squared (Mpc**2)
def uniformArea(lowerBound, upperBound, num_injections):
    shellMin = lowerBound**2
    shellMax = upperBound**2

    distance = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        distance[i] = math.sqrt(np.random.uniform(shellMin, shellMax))
    return distance
 
# Draw a distance via uniform in volume (distance cubed) (Mpc**3)
def uniformVolume(lowerBound, upperBound, num_injections):
    sphereMin = lowerBound**3
    sphereMax = upperBound**3

    distance = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        distance[i] = pow(np.random.uniform(sphereMin, sphereMax),(1./3))
    return distance

# Draw a distance via uniform in log 10 distance (Mpc)
def log10Dist(lowerBound, upperBound, num_injections):
    # Logarithms return NaN if you're at 0!
    if lowerBound <= 0:
       lowerBound = 0.000000000001

    logDMin = math.log10(lowerBound)
    logDMax = math.log10(upperBound)

    distance = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        distance[i] = pow(10, np.random.uniform(logDMin, logDMax))
    return distance

# Draw a distance uniform in Chirp Distance

# Draw a sky location for the binary
# Sampling dec in cos 0 to cos pi gives a uniform
# distribution over S2 surface of sphere. Sampling
# from 0 to 2pi in ra.
# See: http://mathworld.wolfram.com/SpherePointPicking.html 

def uniformDec(num_injections):
    dec = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        dec[i] = np.arccos(2*np.random.uniform(0,1)-1)
    return dec

def uniformRA(num_injections):
    ra = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        ra[i] = 2*np.pi*np.random.uniform(0,1)
    return ra

# Draw a polarization for the binary from 0 to pi
def uniformPolariAngle(num_injections):
    psi = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        psi[i] = 2.*np.pi*np.random.uniform(0.,1.)
    return psi
 
# Draw an inclination angle relative to the line of sight
# 0 to pi
def uniformIncAngle(num_injections):
    incAng = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        incAng[i] = np.arccos(2*np.random.uniform(0,1)-1)
    return incAng

# Draw a coalescence phase for the binary

'''Mass Distribution functions'''

# Parse a Uniform Distribution from distributions.py to simple numpy array
def parseUniformDistr(collinArray, num_injections):
    mass = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        mass[i] = collinArray[i][0]
    return mass

# Parse a Gaussian Distribution from distributions.py to a simple numpy array
def parseGaussDistr(collinArray, num_injections):
    mass = np.ndarray(shape=(num_injections), dtype=float)
    for i in range(0,num_injections):
        mass[i] = collinArray[i][0]
    return mass

'''Spin Distribution functions'''

# Check that the spin magnitude is less than 1
def drawSpinComponents(spin_mag_min, spin_mag_max, theta_min, theta_max,
                       phi_min, phi_max, num_injections):
    # Distributions.py uses angles in units of pi so convert or units
    theta_low = theta_min / np.pi
    theta_high = theta_max / np.pi

    phi_low = phi_min / np.pi
    phi_high = phi_max / np.pi

    spin_mag = np.ndarray(shape=(opts.num_injections), dtype=float)
    theta = np.ndarray(shape=(opts.num_injections), dtype=float)
    phi = np.ndarray(shape=(opts.num_injections), dtype=float)

    # Grab a distribution of polar and azimuthal angles,
    # polar between 0 and pi and azimuthal between 0 and 2pi
    anglePair = distr.UniformSolidAngle(polar_bounds=(theta_low,theta_high),
                                       azimuthal_bounds=(phi_low,phi_high))
    angles = anglePair.rvs(size=num_injections)

    # For some reason Collin's uniformSolidAngle inputs theta and phi but
    # outputs tuples of phi and theta  
    for i in range(0, num_injections):
        spin_mag[i] = np.random.uniform(spin_mag_min, spin_mag_max)
        phi[i] = angles[i][0]
        theta[i] = angles[i][1]

    return coord.spherical_to_cartesian(spin_mag, phi, theta)        

parser = argparse.ArgumentParser(description='A Python code for generating ' \
                                             'a astrophysical population of ' \
                                             'compact binaries.')
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
parser.add_argument('--min-spin1', type=float, required=False, default=0.0,
                     help='Optional: Minimum spin magnitude for spin1' \
                          ' (default:0.0)')

parser.add_argument('--max-spin1', type=float, required=False, default=1.0,
                     help='Optional: Maximum spin magnitude for spin1' \
                          ' (default:1.0)')

parser.add_argument('--min-spin1-theta', type=float, required=False,
                     default=0, help='Optional: Minimum angle to' \
                                     ' place the spin vector angle.' \
                                     ' Angle is relative to the' \
                                     ' azimuth (0 to pi).')
                                           
parser.add_argument('--max-spin1-theta', type=float, required=False,
                     default=np.pi,help='Optional: Maximum angle to' \
                                        ' place the spin vector angle.' \
                                        ' Angle is relative to the' \
                                        ' azimuth (0 to pi).')

parser.add_argument('--min-spin1-phi', type=float, required=False, 
                     default=0.,help='Optional: Minimum angle to place the' \
                                     ' spin vector angle. Angle is in the' \
                                     ' x-y plane (0 to 2pi) (default:0)') 
                                            
parser.add_argument('--max-spin1-phi', type=float, required=False,  
                     default=2.*np.pi, help='Optional: Maximum angle to place' \
                                            ' the spin vector angle. Angle is' \
                                            ' x-y plane (0 to 2pi)' \
                                            ' (default:2pi)')

parser.add_argument('--min-spin2', type=float, required=False, default=0.0,
                     help='Optional: Minimum spin magnitude for spin2' \
                          ' (default:0.0)')

parser.add_argument('--max-spin2', type=float, required=False, default=1.0,
                     help='Optional: Maximum spin magnitude for spin2' \
                          ' (default:1.0)')

parser.add_argument('--min-spin2-theta', type=float, required=False, 
                     default=0,help='Optional: Minimum angle to' \
                                            ' place the spin vector angle.' \
                                            ' Angle is relative to the' \
                                            ' azimuth (-pi/2 to pi/2).')
                                            
parser.add_argument('--max-spin2-theta', type=float, required=False,
                     default=np.pi,help='Optional: Maximum angle to' \
                                           ' place the spin vector angle.' \
                                           ' Angle is relative to the' \
                                           ' azimuth (-pi/2 to pi/2).')

parser.add_argument('--min-spin2-phi', type=float, required=False,
                     default=0.,help='Optional: Minimum angle to place the' \
                                     ' spin vector angle. Angle is in the' \
                                     ' x-y plane (0 to 2pi) (default:0)')
                                            
parser.add_argument('--max-spin2-phi', type=float, required=False,
                     default=2.*np.pi, help='Optional: Maximum angle to place' \
                                            ' the spin vector angle. Angle is' \
                                            ' x-y plane (0 to 2pi)' \
                                            ' (default:2pi)')

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

if opts.verbose:
    logging_level = logging.DEBUG
else:
    logging_level = logging.WARN

log_fmt = '%(asctime)s %(message)s'
log_date_fmt = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(level=logging_level, format=log_fmt, datefmt=log_date_fmt)

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
   mBounds1 = distr.Gaussian(mass1=(opts.min_mass1,opts.max_mass1,
                             opts.mean_mass1,opts.stdev_mass1**2))
   injDict['mass1'] = parseGaussDistr(mBounds1.rvs(size=opts.num_injections),
                                      opts.num_injections)

logging.info('Mass 1 parameters written.')

# Mass 2
if opts.mass2_uniform == True:
   mBounds2 = distr.Uniform(mass2=(opts.min_mass2,opts.max_mass2))
   injDict['mass2'] = parseUniformDistr(mBounds2.rvs(size=opts.num_injections),
                                        opts.num_injections)
else :
   mBounds2 = distr.Gaussian(mass2=(opts.min_mass2,opts.max_mass2,
                             opts.mean_mass2, opts.stdev_mass2**2))
   injDict['mass2'] = parseGaussDistr(mBounds2.rvs(size=opts.num_injections),
                                      opts.num_injections)

logging.info('Mass 2 parameters written.')

# MTotal = m1 + m2
injDict['mtotal'] = injDict['mass1'] + injDict['mass2']

logging.info('Total mass written.')

# MChirp = (m1*m2)**(3./5.) / (m1 + m2)**(1./5.)
# Eta = (m1*m2)/(m1+m2)**2
injDict['mchirp'], injDict['eta'] = pnu.mass1_mass2_to_mchirp_eta(injDict['mass1'],
                                                     injDict['mass2'])
                                   
logging.info('Chirp mass and Eta written.')

# q mass ratio (Convention: write it as bigger mass divided by smaller mass!)
for i in range(0,opts.num_injections):
    if injDict['mass1'][i] > injDict['mass2'][i]:
       injDict['q'].append(injDict['mass1'][i]/injDict['mass2'][i])
    else:
       injDict['q'].append(injDict['mass2'][i]/injDict['mass1'][i])

logging.info('Mass ratio written')

# Spin 1
injDict['spin1x'], injDict['spin1y'], injDict['spin1z'] = drawSpinComponents(
                                                          opts.min_spin1,
                                                          opts.max_spin1,
                                                          opts.min_spin1_theta,
                                                          opts.max_spin1_theta,
                                                          opts.min_spin1_phi,
                                                          opts.max_spin1_phi,
                                                          opts.num_injections)

logging.info('Spin 1 parameters written.')

# Spin 2
injDict['spin2x'], injDict['spin2y'], injDict['spin2z'] = drawSpinComponents(
                                                          opts.min_spin2,
                                                          opts.max_spin2,
                                                          opts.min_spin2_theta,
                                                          opts.max_spin2_theta,
                                                          opts.min_spin2_phi,
                                                          opts.max_spin2_phi,
                                                          opts.num_injections)

logging.info('Spin 2 parameters written.')

# Write the projection of the spin on the angular momentum axis
# x1 = (c*vec{S1}/(G m1**2) dot vec-hat{L}

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

else:
    raise ValueError('Unrecognized distance distribution! Could not find %s.' \
                     'Please try uniform, uniformArea, uniformVolume, or ' \
                     'uniformLog10 distance distributions.' %(opts.dist_distr))

logging.info('Distance parameters written.')

# Draw the sky location of a binary (Right Ascension first)
# (Declination second)
# In radians!
injDict['ra'] = uniformRA(opts.num_injections)
injDict['dec'] = uniformDec(opts.num_injections)

logging.info('RA and Dec written.')

# Draw the GPS times for each injection
injDict['time'] = uniformTime(opts.gps_start_time, opts.gps_end_time,
                              opts.num_injections)

logging.info('Time stamps for injections written.')

# Draw a random inclination of the binary to the line of sight
# between -pi/2 and pi/2
injDict['inclination'] = uniformIncAngle(opts.num_injections)

logging.info('Inclination angles written.')

# Draw a random polarization angle of the binary to the line of sight
# between 0 and 2 pi
injDict['polarization_angle'] = uniformPolariAngle(opts.num_injections)

logging.info('Polarization angles written.')

# Write the injection dictionary to file
logging.info('Writing to file.')

with open(opts.output, "wb") as outfile:
   writer = csv.writer(outfile)
   writer.writerow(injDict.keys())
   writer.writerows(zip(*injDict.values()))

logging.info('100% done')
