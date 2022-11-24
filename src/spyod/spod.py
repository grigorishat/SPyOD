# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 20:57:12 2022

@author: Hatzissawidis

The code is based on https://www.researchgate.net/publication/280933681_SPOD_matlab_example.
"""

__author__ = "Grigorios Hatzissawidis"
__authors__ = ["Grigorios Hatzissawidis", "Moritz Sieber"]
__contact__ = "grigorios.hatzissawidis@fst.tu-darmstadt.de"
__date__ = "2022/10/07"
__email__ =  "grigorios.hatzissawidis@fst.tu-darmstadt.de"
__license__ = "MIT License"
__maintainer__ = "Grigorios Hatzissawidis"
__status__ = "Development"
__version__ = "0.1.0"


import time
import warnings

import numpy as np
from scipy.linalg import toeplitz
from scipy.signal import convolve2d, fftconvolve


def spod(*args, **kwargs):
    """
    This function calculates the Spectral Proper Orthogonal Decomposition according
    to Sieber et al., 2016, JFM.
    SPOD provides basis vectors (modes) and temporal coefficients of the input data.

    The input data can be recomposed from the output data
        If Npod < Nsnap the expression is only approximately equal

    Function inputs:

    Unspecific number of arguments (one argument required at least):
     U,V,W,,...:    Arbitrary number of input components e.g. velocity components.
                    Each numpy array may have any number of dimensions but all must have
                    the same number of dimensions. 
                    For example: The first axis og a 3D numpy array spans the temporal 
                    direction of the data, e.g. U[i,:,:] is the i-th snapshot of a 2D data set.

    Optional keyword arguments:
    *args are the input experimental or simulation data and *kwargs are the settings,
    if the input setting parameter is not defined, default values are used. SPOD(*args,*kwargs)


    Nfilt:      The length of the SPOD filter. Scalar integer greater than or equal
                to 0. Default is 0 (classical POD).
    Npod:       Number of modes to return. Scalar integer greater than 0. Default
                is number of snapshots.
    Wxyz:       Spatial weighting for the inner product. Double array with the same
                dimensions as the spatial dimensions of U,V,W... . The weighted
                inner product reads <u,v>_w = transpose(u)*diag(w)*v. Default is
                uniform spatial weighting e.g. Wxyz = np.ones(U[i,:,:])
    boundary:   Handling of missing data at start and end of the time series.
                conditions. String with value:
                "periodic" boundary conditions
                "zeros": pad with zeros
                Default is periodic.
    corr_type:  Type of correlation used to calculate the SPOD. String with value:
                "spatial" correlation (use when Ngrid*Nfilt < Nsnap),
                "temporal" correlation (use when Ngrid*Nfilt > Nsnap).
                Default is "temporal" (snapshot POD).

    Function outputs:

    Ui:             POD modes (spatial modes). Number of output modes depends
                    on the number of input velocity components. The dimension
                    of the output data is the same as the input data, whereas
                    the last dimension contains individual modes instead of
                    snapshots. The output is a list.
    a:              POD coefficients (temporal modes). 2D array where the first
                    dimension represents temporal snapshots and the second the
                    single modes, e.g. a[:,1] is the time-series of the first mode
                    coefficient.
    lambda:         Energy of individual modes. 1D array of length Npod.
                    np.diag(transpose(a)*a)/Nsnap = lambda.
    mode_norm:      Norm of the single spatial modes. 1D array of length Npod. For
                    classical POD (Nfilt=0) the norm is 1 for all modes.
                    sqrt(<Ui,Ui>_w + <Vi,Vi>_w + <Wi,Wi>_w + ...) = mode_norm.
    """
    t = time.time()
    # parse input arguments
    Nsnap, Ncomp, Ngrid, Nfilt, Npod, input_size, Wxyz, boundary, corr_type, f = parse_input(args,kwargs)

    print("There are {} components.".format(Ncomp)) #pylint: disable=consider-using-f-string

    if corr_type == 'temporal':

    # calculate temporal correlation matrix
        R = np.zeros((Nsnap,Nsnap))

        for ii in range(Ncomp):
            U = np.transpose(args[ii]).reshape((Ngrid,Nsnap))
            R = R + U.transpose() * Wxyz @ U
        R = R/Nsnap/np.sum(Wxyz)

        # calculate filtered correlation matrix S

        if boundary == "periodic":
            # convolve periodically extended matrix with filtered diagonal matrix

            ind = np.array(range(1-Nfilt, Nsnap+Nfilt+1))
            ind = np.mod(ind-1,Nsnap)
            #S = convolve2d(R[np.ix_(ind, ind)], np.diag(f), mode = "valid")
            S = fftconvolve(R[np.ix_(ind, ind)], np.diag(f), mode = "valid")


        # dft case
        elif boundary == "DFTcase":
            rr = np.zeros(Nsnap)
            for ii in range(Nsnap):
                rr[ii] = np.sum(np.diag(R,ii))*f[0]

            # periodic boundary condition -> circulant matrix
            rr[1:] = rr[1:] + rr[:0:-1]
            S = toeplitz(rr)

        # convolve zero padded matrix with filter diagonal matrix
        # ("center" differs from matlab if number of rows/ columns is odd)

        elif boundary == "zeros":
            #S = convolve2d(R,np.diag(f), mode = "same")
            S = fftconvolve(R, np.diag(f), mode = "same") # much faster than convolve2d

        # calculate eigenvalues and eigenvectors of S and sort in descending order
        # distinguish between DFTcase and periodic/ zeros, since eigh is used for symmetric matrices

        if boundary == "DFTcase":

            w, v = np.linalg.eigh(S, UPLO='U')
        else:
            w, v = np.linalg.eig(S)

        idx = w.argsort()[::-1]
        lmbd = w[idx] # mode variance or energy in case of velocity data
        v = v[:,idx] # normalised temporal coefficients

        # cut off POD modes according to the rank of matrix S
        Nrank = np.linalg.matrix_rank(S)

        if not Npod:
            Npod = Nrank
        elif Npod > Nrank:
            Npod = Nrank
            print('SPODrankDeficit: Correlation matrix rank is less than number of request POD modes. Npod is set to Nrank. Npod = '+ str(Nrank))

# compute scaled temporal coefficients
        a = np.multiply(v[:,:Npod],np.sqrt(Nsnap*lmbd[:Npod]))

# calculate spatial modes for all components
        a_proj = np.divide(v[:,:Npod],np.sqrt(Nsnap*lmbd[:Npod]))
        output_size = list(input_size)
        output_size[0] = Npod
        mode_norm = np.zeros(Npod)
        Ui = list()

        for i in range(Ncomp):

            U = np.transpose(args[i]).reshape((Ngrid,Nsnap))
            #Ui.append(np.dreshape(np.dot(U,a_proj),output_size))
            #Ui_vec = np.dot(U,a_proj)
            mode_norm = mode_norm + np.sum(np.dot(U,a_proj).T**2 * Wxyz,axis=1)/np.sum(Wxyz)
            Ui.append(np.dot(U,a_proj).T.reshape(output_size, order = 'F'))

    if corr_type == "spatial":
        Nfilt2 = f.size
        Ncorr = Ngrid*Nfilt2*Ncomp
        R = np.zeros((Ngrid,Nfilt2,Ncomp,Ngrid,Nfilt2,Ncomp))

        # first loop over input component
        for k in range(Ncomp):
            Uk = np.transpose(args[k]).reshape((Ngrid,Nsnap))* np.sqrt(Wxyz)[:,None]

            # second loop over input component
            for l in range(Ncomp):
                Ul = np.transpose(args[l]).reshape((Ngrid,Nsnap))* np.sqrt(Wxyz)[:,None]

                # first loop over filter coefficient
                for i in range(-Nfilt, Nfilt + 1, 1):
                    for j in range(-Nfilt, Nfilt + 1, 1):
                        indi = i + Nfilt
                        indj = j + Nfilt
                        ij_shift = i-j
                        f_ij = np.sqrt(f[indi]*f[indj])

                        if boundary == "periodic":
                            # correlation with periodic boundary conditions
                            R[:,indi,k,:,indj,l] = np.multiply(f_ij, Uk @ np.roll(Ul,[0, ij_shift], axis=(0, 1)).transpose())

                        else:
                            # correlation with finite signal (zero padded)
                            u_lim = min([Nsnap-i, Nsnap-j, Nsnap])
                            l_lim = max([1-i, 1-j, 1])
                            subi = np.array(range(l_lim, u_lim,1)) + i
                            subj = np.array(range(l_lim, u_lim,1)) + j
                            R[:,indi,k,:,indj,l] = np.multiply(f_ij, Uk[:,subi] @ Ul[:,subj].transpose())

        del Uk, Ul
        R = np.reshape(R,(Ncorr,Ncorr),order='F') / (Nsnap*sum(Wxyz))

# calculate eigenvalues and eigenvector of R
        lmbd, v = np.linalg.eig(R) # energy of modes and normalized spatial modes

# cut POD modes according to rank of matrix R if necessary
        Nrank = np.linalg.matrix_rank(R)
        if not Npod:
            Npod = Nrank
        elif Npod > Nrank:
            Npod = Nrank
            print('SPODrankDeficit: Correlation matrix rank is less than number of request POD modes. Npod is set to Nrank.')

        # compute scaled temporal coefficients
        output_size = np.array(input_size)
        output_size[0] = Npod
        a = np.zeros((Nsnap,Npod),dtype = 'complex_')
        mode_norm = np.zeros((1,Npod))
        Nconvfilt = Ngrid * Nfilt2
        Ui = list()
        for j in range(Ncomp):
            U = np.transpose(args[j]).reshape((Ngrid,Nsnap))
            if boundary == "periodic":
                # create periodically extended time series
                Uext = np.concatenate((U[:,-Nfilt:], U, U[:,0:Nfilt]),1)
            else: 
                # create zero padded time series
                Uext = np.concatenate((np.zeros((Ngrid,Nfilt)),U,np.zeros((Ngrid,Nfilt))),1)

            idx = j * Nconvfilt + range(Nconvfilt) 

            for i in range(Npod):
                # pick apropriate spatial mode
                f_proj = v[idx,i].reshape(Ngrid,Nfilt2, order = 'F')
                # apply spatial weihting
                f_proj = f_proj * np.sqrt(Wxyz/np.sum(Wxyz))[:,None]
                # apply temporal weighting (filter)
                f_proj = f_proj * np.sqrt(f)
                # flip dimensions before convolution to obtain a correlation
                f_proj = np.rot90(f_proj,2)
                # calculate temporal coefficient
                a[:,i][:,None] = a[:,i][:,None] + convolve2d(Uext, f_proj, mode = "valid").T

            # pick central spatial mode for the output and scale it
            idx = j * Nconvfilt + Nfilt * Ngrid + range(Ngrid)
            scale = np.sqrt(np.sum(Wxyz)/ f[Nfilt+1]/ Wxyz)
            scale[Wxyz==0] = 1 # avoid singularities
            Ui_vec = v[idx,:Npod] * scale[:,None]
            mode_norm = mode_norm + np.sum(Ui_vec**2 * Wxyz[:,None])/np.sum(Wxyz)
            Ui.append(Ui_vec.T.reshape(output_size, order = 'F'))

    print("Elapsed: %.2f min" %(np.double((time.time() - t))/60)) #pylint: disable=consider-using-f-string

    return lmbd, a, mode_norm, Ui

def parse_input(in_data,in_set):
    """
    this function checks the data for validity

    inputs: in_data is a non-keyworded variable-length argument list of the data to be analysed
            in_set is a keyworded, variable-length argument list of the SPOD settings

    output: Nsnap, Ncomp, Ngrid, Nfilt, Npod, input_size, Wxyz, boundary, corr_type, f
    """
    Ncomp = len(in_data)
    assert not Ncomp < 1, 'Input components U,V,W,.. are missing.'

    # determine size of input data
    for ii in range(Ncomp):
        U = in_data[ii]

        if ii == 0:
            input_size = U.shape
            Nsnap = input_size[0]
        else:
            if U.shape[0] != Nsnap:
                raise Exception('Number of snapshots in additional component is not consistent with first component.')
    Ngrid = np.prod(input_size[1:])



    if 'Nfilt' in in_set:
        if np.isscalar(in_set.get('Nfilt')):
            if in_set.get('Nfilt') <= Nsnap:
                if in_set.get('Nfilt') >= 0:
                    Nfilt = np.round(in_set.get('Nfilt'))
                else:
                    raise Exception('filter size must be a positive number.')

            else:
                warnings.warn('filter size is larger than number of snapshots')
                Nfilt = np.round(in_set.get('Nfilt'))
        else:
            raise Exception('filter size must be scalar number.')

    else:
        Nfilt = 0



    if 'Npod' in in_set:
        if np.isscalar(in_set.get('Npod')):
            if in_set.get('Npod') >= 1:
                Npod = np.round(in_set.get('Npod'))
            else:
                raise Exception('number of modes must be a positive number.')

        else:
            raise Exception('number of modes must be scalar number.')

    else:
        Npod = []



    if 'Wxyz' in in_set:
        if in_set.get('Wxyz').size == Ngrid:
            Wxyz = in_set.get('Wxyz')
            Wxyz = np.concatenate(Wxyz.transpose())
        else:
            raise Exception('dimension of weighting does not match the input data.')
    else:
        Wxyz = np.ones(Ngrid)


    if 'boundary' in in_set:
        if in_set.get('boundary') == "zeros" or in_set.get('boundary') == "periodic":
            boundary = in_set.get('boundary')

        else:
            raise Exception('boundary parameter must be "zeros" or "periodic".')
    else:
        boundary = "periodic"

    # temporal or spatial correlation?
    if 'corr_type' in in_set:
        if in_set.get('corr_type') == 'temporal' or in_set.get('corr_type') == 'spatial':
            corr_type = in_set.get('corr_type')
        else:
            raise Exception('correlation type must be temporal or spatial.')  
    else:
        corr_type = 'temporal'

    # define filter
    if Nfilt == Nsnap and corr_type == 'temporal':
    # choose box filter to obtain DFT
        f = filter_coefficient(Nfilt, 'box')
        boundary = "DFTcase"
    else:
        f = filter_coefficient(Nfilt, 'gauss')



    # check problem size and display appropriate warning
    if corr_type == 'temporal':
        Ncorr = Nsnap
    else: # spatial
        Ncorr = Ngrid*len(f)
    if Ncorr > 10000:
        warnings.warn('SPOD:LargeProblem computation takes long time >1 hour')
    elif Ncorr > 2000:
        warnings.warn('SPOD:LargeProblem computation takes some time >1 min')
    elif Ncorr > 1000:
        warnings.warn('SPOD:LargeProblem computation may take some time')

    #print(Nsnap, Ncomp, Ngrid, Nfilt, Npod, input_size, Wxyz, boundary, corr_type, f)
    print('Input parameter checked.')
    return Nsnap, Ncomp, Ngrid, Nfilt, Npod, input_size, Wxyz, boundary, corr_type, f


def filter_coefficient(N_filt, filt_type):
    """
    calculates the filter coefficient for a box or a gaussian filter

    input: filter size N_filt and filt_type
    output: filter coefficient f
    """
    if filt_type == 'box':
        f = np.ones(N_filt)
    elif filt_type == 'gauss':
        f = np.exp(-np.linspace(-2.285,2.285,2*N_filt+1)**2)
    else:
        raise Exception('unknown filter type.')
    f = f/np.sum(f)
    return f