# Copyright (C) 2016  Steven Reyes
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module provides a class for injection parameters
"""

class InjParams(object):
   """
   The injection parameters class that stores all of the necessary parameters
   for a software injection:

                coa_phase              :
                mass1                  :
                mass2                  :
                eta                    :
                mchirp                 :
                latitude               :
                geocent_end_time       :
                spin1x                 :
                spin1y                 :
                spin1z                 :
                spin2x                 :
                spin2y                 :
                spin2z                 :
                h_end_time             :
                l_end_time             :
                distance               :
                eff_dist_l             :
                eff_dist_h             :
                lattitude              :
                longitude              :
                polarization           :
                waveform_approximant   :
                inclination            :
                simulation_id          :
                f_lower                :
                f_final                :
   """
   def __init__(self):
       self.coa_phase = 0
       self.mass1 = 0
       self.mass2 = 0
       self.eta = 0
       self.mchirp = 0
       self.lattitude = 0
       self.geocent_end_time = 0
       self.spin1x = 0
       self.spin1y = 0
       self.spin1z = 0
       self.spin2x = 0
       self.spin2y = 0
       self.spin2z = 0
       self.h_end_time = 0
       self.l_end_time = 0
       self.distance = 0
       self.eff_dist_l = 0
       self.eff_dist_h = 0
       self.lattitude = 0
       self.longitude = 0
       self.polarization = 0
       self.waveform_approximant = ""
       self.inclination = 0
       self.simulation_id = 0
       self.f_lower = 0
       self.f_final = 0

