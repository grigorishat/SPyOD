# -*- coding: utf-8 -*-
"""
This finds the coupled mode pairs.

Created on Mon Mar 28 13:54:15 2022

@author: Hatzissawidis

The code is based on https://www.researchgate.net/publication/280933681_SPOD_matlab_example.

"""

import numpy as np

def findpairs(a):
    """
    This function finds mode pairs from a Spectral Proper Orthogonal Decomposition.

    Function input:
        a:          SPOD mode coefficients -> output of spod()

    Function output:
    mode: Structure where single element (for one mode pair k) consists of:
        mode['ind']: Indices [i,j] of the coressponding SPOD modes
        mode['c']:   Magnitude of the harmonic corellation
        mode['a']:   Analytic mode coefficient a[:,i] + 1j*a[:,i]
        mode['at']:  Derivative of analytic mode coefficient
        mode['K']:   Relative energy content of a single mode pair
        mode['f']:   Estimated frequency of the mode pair [1/sample]
    """
    Nsnap, Npod = np.shape(a)
    # relative energy content, without the part which is cut off 
    lmbd = np.diag(np.dot(a.transpose(),a))
    lmbd = lmbd/np.sum(lmbd)

    ag = np.gradient(a) # gradient
    ag = ag[0]
    # find combined modes from DMD

    T = np.linalg.lstsq(a[:Nsnap-1],a[1:Nsnap],rcond=-1)[0].conj().T


    mu, V = np.linalg.eig(T)
    
    # compute frequency

    omega = np.imag(np.log(mu))

    # linked modes are only phase shifted -> correlate real and imaginary part
    # of DMD eigenvectors
    C = np.imag(V @ np.diag(np.sign(omega)) @ V.conj().T)/2

    # problem in relative mode energy: last mode is the first, so change
    Nmode = np.trunc(Npod/2)
    Nmode = Nmode.astype(int)
    
    mode = {'c' : np.zeros(Nmode),'ind' : np.zeros((Nmode,2)),
            'a' : np.zeros((Nmode, Nsnap), dtype=complex),
            'at': np.zeros((Nmode, Nsnap), dtype=complex),
            'K' : np.zeros(Nmode),'f' : np.zeros(Nmode)}

    for i in range(Nmode):
        #pick maximum
        tmp = np.max(C)
        mode['c'][i] = tmp
        indij = np.unravel_index(np.argmax(C, axis=None),C.shape)
        # delete rows and column
        C[indij,:] = 0
        C[:,tuple(reversed(indij))] = 0
        # save indices to output
        mode['ind'][i] = indij

        # compose analytical mode coefficient
        mode_sign = np.sign(np.sum(a[:,indij[0]]*ag[:,indij[1]]))
        mode['a'][i] = a[:,indij[0]] + mode_sign*1j*a[:,indij[1]]
        mode['at'][i] = ag[:,indij[0]] + mode_sign*1j*ag[:,indij[1]]

        # calculate common energy content

        mode['K'][i] = lmbd[indij[0]] + lmbd[indij[1]]

        # estimate mode frequency

        V_comb = np.abs(V[indij[0],:])**2 + np.abs(V[indij[1],:])**2
        ind = np.argsort(-V_comb)
        tmp = -np.sort(-V_comb)
        # sum the first 3 frequencies (each frequency appears twice)
        mode['f'][i] = np.sum(np.abs(omega[ind[:6]]) * tmp[:6])/np.sum(tmp[:6])/(2*np.pi)

    return mode