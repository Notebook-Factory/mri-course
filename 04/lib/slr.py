# SLR Transform for RF Pulse Design
#
# Dr. John Pauly's RF Toolbox written in MATLAB
#
# written by John Pauly
# (c) Board of Trustees, Leland Stanford Junior University
#
# Ported from MATLAB to Python by Mahesh Keerthivasan, 
# University of Arizona, 2018

import numpy as np
import scipy as sp

def ab2ex(a,b):
    # mxy = ab2ex(a,b)    -- or --    mxy = ab2ex(ab)
    # Computes the excitation profile 2*conj(a).b
    #  written by John Pauly, 1992
    #  (c) Board of Trustees, Leland Stanford Junior University
    #
    # Ported from MATLAB to Python by Mahesh Keerthivasan, 
    # University of Arizona, 2018

    mxy = np.multiply(2*np.conj(a),b)
    return mxy

def ab2inv(a,b):
    # mz = ab2inv(a,b)    -- or --    mz = ab2inv(ab)
    # Computes the inversion profile 1-2*b.*conj(b)
    #  written by John Pauly, 1992
    #  (c) Board of Trustees, Leland Stanford Junior University
    #
    # Ported from MATLAB to Python by Mahesh Keerthivasan, 
    # University of Arizona, 2018
    
	mz = 1 - 2 * np.multiply(np.conj(b),b)
	return mz

def abrm(rf,g,*posVec):
    #  [a b] = abrm(rf,[g],[x [,y])
    #  simulate an rf pulse, returning the Cayley-Klein parameters, alpha and
    #  beta. This is the .m file version of abr(), which is compiled and faster.
    #  Inputs: 
    #   rf -- rf scaled so that sum(rf) = flip angle
    #   g -- optional gradient waveform, scaled so that (gamma/2*pi)*sum(g) = k 
    #           in cycles/cm
    #   x -- position vector
    #   y -- optional position vector for 2D pulses (assumes imag(g)  = gy)
    #
    #  Written by John Pauly, Dec 22, 2000
    #  (c) Boaard of Trustees, Leland Stanford Jr. University
    #
    #  Translated from octave, and modified to scale gradient by 
    #  2pi, so k = cumsum(g) is cycles/cm
    #  Sept 27, 2004
    #
    # Ported from MATLAB to Python by Mahesh Keerthivasan, 
    # University of Arizona, 2018
    # Python port: all inputs are numpy arrays
    
    if len(posVec) == 1:
        x = posVec[0]
        y = np.zeros(np.shape(x))
    else:
        x = posVec[0]
        y = posVec[1]

    a = np.zeros((np.size(x),np.size(y)),dtype=np.complex_)
    b = np.zeros((np.size(x),np.size(y)),dtype=np.complex_)

    for jj in range(0,np.size(y)):
        for kk in range(0,np.size(x)):
            # make sure om isn't exactly zero, so n doesn't blow up
            om = np.multiply(x[kk],np.real(g)) + np.multiply(y[jj], np.imag(g))
            om = om + (np.absolute(om)< np.finfo(float).eps)*np.finfo(float).eps            
            phi = np.sqrt(np.multiply(rf,np.conj(rf))+ np.square(om))
            n1 = np.divide(np.real(rf), phi)
            n2 = np.divide(np.imag(rf),phi)
            n3 = np.divide(om,phi)
            av = np.cos(phi/2) - 1j *np.multiply(n3, np.sin(phi/2))
            bv = np.multiply(-1j*(n1+1j*n2),np.sin(phi/2))
            abt = np.array([[1], [0]])
            for m in range(0,np.size(phi)):
                tmp = np.array([[av[m], -np.conj(bv[m])], [bv[m], np.conj(av[m])]])
                abt = np.matmul(tmp, abt)
            
            a[kk,jj] = abt[0,0]
            b[kk,jj] = abt[1,0]

    return a,b

def gt2cm(x,g,t):
    #  Takes dimensionless x used by abr, and scales it to cm based
    #  on a gradient strength g (G/cm) and pulse duration t (ms)
    #
    #   xs = gt2cm(x,g,t)
    #
    #     x  -- normalized frequency x vector used by abr
    #     g  -- Gradient strength, G/cm
    #     t  -- pulse duration in ms
    #
    #     xs -- scaled spatial axis, in cm
    #
    #  written by John Pauly, 1992
    #  (c) Board of Trustees, Leland Stanford Junior University
    #
    # Ported from MATLAB to Python by Mahesh Keerthivasan, 
    # University of Arizona, 2018

    xs = np.divide(x,(4.257*g*t))
    return xs





