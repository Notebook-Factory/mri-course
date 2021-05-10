# Utilities for design of RF pulses 
#
# Dr. John Pauly's RF Toolbox written in MATLAB
#
# written by John Pauly, 1992
# (c) Board of Trustees, Leland Stanford Junior University
#
# Ported from MATLAB to Python by Mahesh Keerthivasan, 
# University of Arizona, 2018

import numpy as np
import scipy as sp

def rfscaleg(rf,t):
    # scale RF to Gaus
    #  rf  -- rf waveform, scaled so sum(rf) = flip angle
    #  t   -- duration of RF pulse in ms
    #  rfs -- rf waveform scaled to Gauss

    gamma = 2*np.pi*4.257 #kHz/G
    dt = t/np.size(rf)
    rfs = np.divide(rf,(gamma*dt))
    return rfs