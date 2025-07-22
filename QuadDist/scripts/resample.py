#!/usr/bin/env python3

import numpy as np

# For a description of the jackknife method as implemented here, see
# "Everything you wanted to know about Data Analysis and Fitting but were
# afraid to ask" by Peter Young, ArXiv:1210.3781


def jackknife_resample(samples, k=1):
    '''
    Computes the averages of resampled sets of samples for k-removed jackknife
    method. The length of the list of samples must be divisible by k.

    Parameters:
      samples    list or NumPy array of samples
      k          block size (default: 1)

    Returns a NumPy array of the partial averages.
    '''

    # Assert that the samples can be evenly split to blocks
    nr_samples = len(samples)
    if nr_samples % k != 0:
        raise ValueError("Number of samples not divisible by the jackknife "
            "block size!")

    nr_blocks = nr_samples // k
    sum_all = np.sum(samples)
    sum_block = np.array([np.sum(samples[i*k:(i+1)*k]) for i in range(nr_blocks)])
    jackknife_averages = (sum_all - sum_block)/(nr_samples - k)

    return jackknife_averages


def jackknife_observable(samples, signs, k=1):
    '''
    Computes the jackknifed averages of samples for k-removed jackknife
    from Monte Carlo samples with signs.
    '''
    # Assert that the samples can be evenly split to blocks
    nr_samples = len(samples)
    if nr_samples % k != 0:
        raise ValueError("Number of samples not divisible by the jackknife "
            "block size!")

    if nr_samples != len(signs):
        raise ValueError("Number of samples is not consistent with the number "
            "of signs!")
    
    nr_blocks = nr_samples // k
    sum_all = np.sum(samples*signs)
    sum_all_signs = np.sum(signs)
    sum_block = np.array([np.sum(samples[i*k:(i+1)*k]*signs[i*k:(i+1)*k]) for i in range(nr_blocks)])
    sum_block_signs = np.array([np.sum(signs[i*k:(i+1)*k]) for i in range(nr_blocks)])

    jackknife_averages = (sum_all - sum_block)/(sum_all_signs - sum_block_signs)

    if np.any(np.abs(jackknife_averages.imag) > 1e-6):
        raise ValueError("Expected imaginary parts to vanish -- check this!")
    
    return jackknife_averages.real


def jackknife_results(function_values):
    '''
    Computes the average and the standard error from function values computed
    from the jackknife-resampled inputs.

    Parameters:
      function_values   Numpy array of function values computed with

    Returns (average, standard_error)
    '''
    nr_blocks = len(function_values)
    average = sum(function_values)/nr_blocks

    # Note: the resampled standard error has a prefactor that takes into
    # account the fact that the jackknifed samples are strongly correlated
    # (i.e. this is not meant to be the "usual" standard error formula)
    standard_error = np.sqrt((nr_blocks - 1.0)/nr_blocks
        * sum( (average - function_values)**2 ))

    return average, standard_error

