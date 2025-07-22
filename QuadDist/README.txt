Supplementary material for "Statistical theory of deformation distributions in nuclear spectra"
===============================================================================================

This README file contains the description for the files included in the supplementary material.
For details on how to use the included scripts, see INSTRUCTIONS.txt.


Contents of the input/ directory
--------------------------------

The files resampled_sm*_db64.txt contain the delete-k jackknife averages (k=40 corresponding to a single
independent MC walk in this work) of the moments <Q_20^n> (n=1,2,...,5), the (unprojected) thermal energy, and the (unprojected) heat capacity. The first column is the inverse temperature 1/T and the subsequent columns list the values of the observables in the above order. The sets with a different MC walk omitted are separated by two empty lines.

For instructions how to reproduce the results presented in this work from the data, see the file INSTRUCTIONS.txt.


Contents of the data/ directory
-------------------------------

The files out_abc_sm*_db64.txt contain six data sets each, separated from each other by two empty
lines and numbered 0-5 (in correspondence with the indexing used gnuplot):

(0) the expansion parameters a, b, and c, as well as tau=a*c/b^2, computed from the AFMC data;
(1) the temperature derivatives of a, b, and c computed from the AFMC data, calculated from the
    three-point formula for the derivative;
(2) the second temperature derivatives of a, b, and c using the AFMC data;
(3) the spline fit to a, b, and c evaluated at discrete values of the inverse temperature;
(4) the first derivatives of the spline fits evaluated at the discrete points in (2); and
(5) the second derivatives of the spline fits evaluated at the discrete points.

All errors are estimated using the jackknife method as described in the manuscript. The order of
the columns is explained in the comment lines in the beginning of each dataset.


The files out_prob_sm*_db64.txt contain two datasets each:

(0) the integral of the probability distribution over the (beta, gamma) plane, and
(1) the integrals of the probability distribution over the different shape regions.


The files pq20.sm*.db64.txt contain three datasets each: these are the P_{q_20} distributions
obtained from the AFMC code at three distinct temperatures (specified in the initial comments of
each dataset). The values of q_{20} in this file are to be scaled by a factor 2 to obtain
the actual values of q_20 (to take into account core polarization) and, correspondingly, the values
of P_{q_20} should be divided by 2).


The files qmoments.sm*.db64.txt contain the moments <Q_20^k>/2^k computed from the SMMC data, with
their uncertainties. That is, the actual moments (and their uncertainties) are obtained by multiplying
the values in this file by 2^k.


The files sm*.db64.lvl2.txt contain the total state densities computed directly from the AFMC thermal energy and heat capacity.


The files zvb.sm*.db64.txt contain the P_{q20} distributions computed from the a, b, c fits at
three different temperatures. The values in these files do not need scaling.


The files state_density_ratios_sm*_db64.txt contain the ratios of shape-dependent state
densities (integrated over the three shape regions) to the total state density computed directly
from the AFMC energy and heat capacity.


Contents of the plots/ directory
--------------------------------

This directory contains the scripts that generate the figures in the manuscript. The .plt files are
gnuplot scripts tested to work with gnuplot version 5.2.2. The .py files for generating the more
complicated (beta, gamma) plots are compatible with Python 3.6, matplotlib 2.1, and NumPy 1.13.


Contents of the scripts/ directory
----------------------------------

This directory contains the main scripts used for post-processing of the AFMC results. The Python files are tested with Python 3.6, NumPy 1.13, and SciPy 1.0. The Julia files are compatible with Julia 0.6.

qlvl_abcfit.jl: the main routine for determining the Landau expansion parameters a, b, c from the moments 
                <Q_20^k> (k=2,3,4).

q20dist.jl:    the main routine for computing the marginal distribution P(q_{20}) from the Landau-like
               distribution with the given parameters a, b, and c.

qdist.jl:      the functions used by q20dist.jl and qlvl_abcfit.jl.

qlvl_main.py: the main program for computing the various projected observables and their integrals
              in the (beta, gamma) plane given the Landau expansion parameters a, b, c.

qlevden.py:   functions required by qlvl_main.py, doing the actual computational work.

resample.py:  auxiliary functions for handling the jackknife resampling.


