""" The module contains functions for calculating the
Kullback-Leibler divergence and other information theory
quantities.
"""

import numpy
from scipy import stats


def k_l_div(samples_1, samples_2, pdf_1 = False, pdf_2 = False,
            x_min=None, x_max=None, bandwidth_1="scott",
            bandwidth_2="scott", base=numpy.e)
    """ Computes the Kullback-Leibler divergence for a single parameter
    from two distributions. The first sample set is the reference set,
    while the second sample set is the comparison set. This is an
    asymmetric comparison of the informational content or divergence
    between the first sample set and the second sample set. This can
    also be thought of as the relative entropy between the two sample
    sets.
    (https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)

    This module uses a kernel density estimator (KDE) to estimate the
    probability distribution across all samples. Note, this function is
    not defined if the pdf of the two distributions is 0 within the
    domain of the two distributions.

    Parameters
    ----------
    samples1 : numpy.array
        Samples from the reference set (use as pdf if pdf_1 = True).
    samples2 : numpy.array
        Samples from the comparison set (use as pdf if pdf_2 = True).
    pdf_1 : bool
        Use True if inputting pdf values and not samples.
        [Default = False].
    pdf_2 : bool
        Use True if inputting pdf values and not samples.
        [Default = False].
    x_min : float
        Use x_min bounds for the minimum bound when pdf_1(2) is given.
        [Default = None]
    x_max : float
        Use x_max bounds for the maximum bound when pdf_1(2) is given.
        [Default = None]
    bandwidth : str, float
        The KDE used to estimate the probability distribution relies
        on a bandwidth or standard deviation for the estimation process.
        Default is [Scott] method.
    base : float
        The base to calculate the Kullback-Liebler divergence in.
        The default is in the base of the Euler constant (~2.718)

    Returns
    -------
    numpy.float64
        The Kullback-Leibler divergence value.
    """

    if pdf_1 is False:
        kernel_1 = stats.gaussian_kde(samples_1)
        kernel_1.set_bandwidth(bandwidth_1)

        if pdf_2 is True:
            x_vals_idx = numpy.logical_and(samples_1 > x_min,
                                           samples_2 < x_max)

            x_vals = samples_1[x_vals]

        else :
            x_vals = numpy.concatenate([samples_1, samples_2])
            x_vals = numpy.unique(x_vals)

        pdf_samples_1 = kernel_1(x_vals)

    else :
        pdf_samples_1 = samples_1

    if pdf_2 is False:
        kernel_2 = stats.gaussian_kde(samples_2)
        kernel_2.set_bandwidth(bandwidth_2)

        if pdf_1 is True:
            x_vals_idx = numpy.logical_and(samples_1 > x_min,
                                           samples_2 < x_max)
            x_vals = samples_2[x_vals_idx]

        else :
            x_vals = numpy.concatenate([samples_1, samples_2])
            x_vals = numpy.unique_1(x_vals) 

        pdf_samples_2 = kernel_2(x_vals)

    else :
        pdf_samples_2 = samples_2

    return stats.entropy(pdf_samples_1, qk=pdf_samples_2, base=base)
