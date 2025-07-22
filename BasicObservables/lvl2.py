#!/usr/bin/env python
# Encoding: UTF-8
# Simple code to compute the level density from E vs. β with statistical
# errors. If the heat capacity is available, will use that too.
#
# TO USE THIS SCRIPT,
#  1. Create an input file as described below
#  2. Run
#       lvl2.py <inputfile>
#
#  The input file should contain:
#
#     np = <np>     # Number of protons
#     nn = <nn>     # Number of neutrons
#     nsp = <nsp>   # Number of proton single-particle states
#     nsn = <nsn>   # Number of neutron single-particle states
#     E0 = <E0>     # ground-state energy
#     energies = [...]  # Excitation energies to interpolate to ([] for none)
#     filename = <filename>  # Data file
#     zxz = <zxz>   # (Optional) Z_X/Z
#
#  Where Z_X/Z (optional) is the ratio of the projected to the unprojected
#  partition function at β=0, and <filename> is for the data file.
#
#  The data file should have the format
#
#        <beta> <E> <dE> [<C> <dC>]
#        <beta> <E> <dE> [<C> <dC>]
#          :     :   :
#
#   in order of increasing beta. The heat capacity <C> <dC> is optional.
import sys, os
from math import *
from numpy import *
from numpy import savetxt


################################################################################
# Main code                                                                    #
################################################################################

if len(sys.argv) != 2:
    sys.stderr.write("Error: Must specify input file, no other parameters.\n")
    sys.exit(1)

# "Read" input file
#execfile(sys.argv[1])
exec(open(sys.argv[1]).read())

if not 'zxz' in globals():
    zxz = 1.

################################################################################
# Functions                                                                    #
################################################################################

def lnbinomial(n,k):
    """Compute log (n choose k)"""
    if k > n:
        raise Exception("Error: k > n in lnbinomial")
    lnb = 0.
    for m in range(1,min(k,n-k)+1):
        lnb = lnb + log(float(n-m+1)/float(m))
    return lnb


def lfit(x,y,dy):
    """Linear least-squares fit y = a x + b to data with errors.

    Input:
      x:   List of x values
      y:   List of y values
      dy:  List of errors in y values

    Returns:
      a:  Slope of linear fit
      da: Error in slope
      b:  Intercept
      db: Error in the intercept

    Note that we don't compute the covariance between a and b.
    The equations below can be derived by minimizing
       chi^2 = \sum (a x[i] + b - y[i])**2 / dy[i]**2
    with respect to the parameters a, b.
    """
    A11 = sum([ x[i]**2 / dy[i]**2 for i in range(0,len(x)) ])
    A12 = sum([    x[i] / dy[i]**2 for i in range(0,len(x)) ])
    A21 = A12
    A22 = sum([     1   / dy[i]**2 for i in range(0,len(x)) ])
    detA = A11*A22 - A12*A21

    a = sum([ y[i] * (A22 * x[i] - A12) / dy[i]**2 for i in range(0,len(x)) ])
    a = a / detA
    b = sum([ y[i] * (A11 - A21 * x[i]) / dy[i]**2 for i in range(0,len(x)) ])
    b = b / detA

    da = sum([ (A22 * x[i] - A12)**2 / dy[i]**2 for i in range(0,len(x)) ])
    da = sqrt(da)/abs(detA)
    db = sum([ (A11 - A21 * x[i])**2 / dy[i]**2 for i in range(0,len(x)) ])
    db = sqrt(db)/abs(detA)

    return a, da, b, db


def interpolate(Ei, betas, E, dE, lnrho, dlnrho):
    """Interpolate ln(ρ) to a particular energy E.

    Input:
      Ei:  Energy to interpolate to
      betas, E, dE, lnrho, dlnrho: Arrays of length 2

    Output:
      lnrhoi, dlnrhoi:  Interpolated ln(ρ) and the error"""

    if not all(map(lambda x: len(x) == 2, [E,dE,lnrho,dlnrho,betas])):
        raise Exception("Error in interpolate()")

    a = (E[1] - Ei)/(E[1] - E[0])
    b = (Ei - E[0])/(E[1] - E[0])
    R = (lnrho[0] - lnrho[1])/(E[1] - E[0])

    lnr = a*lnrho[0] + b*lnrho[1]
    dlnr = a**2 * R**2 * dE[0]**2 + b**2 * R**2 * dE[1]**2 + a**2 * \
        dlnrho[0]**2 + b**2 * dlnrho[1]**2 + a**2 * R * betas[0] * dE[0]**2 + \
        b**2 * R * betas[1] * dE[1]**2
    dlnr = sqrt(dlnr)

    return lnr, dlnr


################################################################################
# Read in data                                                                 #
################################################################################

# Read in E vs. β data
# Must be in order of increasing β
E = []
dE = []
C = []
dC = []
betas = []
for line in open(filename):
    s = line.strip()  # Get rid of whitespace
    # Ignore comments and blank lines
    if len(s) > 0 and s[0] != "#":
        vals = s.split()  # beta  E  dE [C dC]
        betas.append(float(vals[0]))
        E.append(float(vals[1]))
        dE.append(float(vals[2]))
        if len(vals) > 3:
            C.append(float(vals[3]))
            dC.append(float(vals[4]))
        else:
            C.append(float('NaN'))
            dC.append(float('NaN'))

if len(betas) < 2:
    raise Exception("Error: need more than one beta point.")

################################################################################
# Level density                                                                #
################################################################################

lnrho = zeros(len(betas))
dlnrho = zeros(len(betas))
lnrho[:] = float('NaN')
dlnrho[:] = float('NaN')
S = zeros(len(betas))
dS = zeros(len(betas))
lnz0 = lnbinomial(nsp,np) + lnbinomial(nsn,nn) + log(zxz) # ln(Z) at beta = 0
lnz = lnz0
vlnz = 0.
# Running sum for the error in ln(Z)
errsum = dE[0]**2 * (betas[1]-betas[0])**2/4

for i in range(1,len(betas)):
    lnz = lnz - (E[i-1] + E[i]) * (betas[i] - betas[i-1]) / 2.
    # Variance of ln(Z)
    vlnz = errsum + (dE[i] * (betas[i]-betas[i-1])/2)**2
    if i < len(betas)-1:
        errsum = errsum + (dE[i]*(betas[i+1]-betas[i-1])/2)**2
    S[i] = betas[i] * E[i] + lnz  # Entropy
    dS[i] = sqrt(betas[i]**2 * dE[i]**2 + vlnz)

    # Compute heat capacity if necessary
    if isnan(C[i]):
        # Derivative: 3-point rule
        #deriv = (E[i+1] - E[i-1])/(betas[i+1] - betas[i-1])
        # Derivative: Linear fit w/ error bars
        # (We add a small constant to dE to avoid problems when dE = 0)
        if i < len(betas)-1:
            deriv, dderiv, tmp1, tmp2 = \
                lfit(betas[i-1:i+2],E[i-1:i+2], [dE[j] + 1E-10 for j in range(i-1,i+2)])
            C[i] = -deriv*betas[i]**2
            dC[i] = dderiv*betas[i]**2

    if C[i] > 0.:
        #lnrho[i] = S[i] - 0.5 * log(-2. * pi * deriv)
        # TODO: Include covariance term?
        #dlnrho[i] = sqrt(dS[i]**2 + 0.5**2 / deriv**2  * dderiv**2)
        lnrho[i] = S[i] - 0.5 * log(2. * pi * C[i]/betas[i]**2)
        dlnrho[i] = sqrt(dS[i]**2 + 0.25 / C[i]**2  * dC[i]**2)


################################################################################
# Interpolation                                                                #
################################################################################

lnrhoi = zeros(len(energies))
dlnrhoi = zeros(len(energies))
lnrhoi[:] = float('NaN')
dlnrhoi[:] = float('NaN')
i1s = [0]*len(energies)
i2s = [0]*len(energies)

for i in range(0,len(energies)):
    Ei = energies[i] + E0

    i2 = 1
    while i2 < len(betas) and E[i2] > Ei:
        # i1: previous index s.t. E[i1] > E
        # i2: next valid index
        i1 = i2
        i2 = i2 + 1
        # Skip invalid points
        while i2 < len(betas) and isnan(lnrho[i2]):
            i2 = i2 + 1
    if i2 == 1 or i2 >= len(betas): continue # Couldn't bracket the energy
    if i2-i1 > 4: continue # Avoid extreme interpolations

    lnrhoi[i], dlnrhoi[i] = interpolate(Ei, [betas[i1],betas[i2]],
                                        [E[i1],E[i2]], [dE[i1],dE[i2]],
                                        [lnrho[i1],lnrho[i2]],
                                        [dlnrho[i1],dlnrho[i2]])
    i1s[i], i2s[i] = i1,i2

################################################################################
# Output                                                                       #
################################################################################

sys.stdout.write("# Level density calculation from lvl2.py\n")
sys.stdout.write("# Input file: %s\n" % sys.argv[1])
sys.stdout.write("# Z_X/Z = %14.7E at beta=0\n" % zxz)
sys.stdout.write("#%11s %16s %14s %14s %14s %14s %14s %14s\n" %
                 ("beta", "E", "ln(rho)", "dln(rho)", "S", "dS", "C", "dC"))
for i in range(1,len(betas)-1):
    sys.stdout.write("%12.6F %16.8E %14.6E %14.6E"
                     " %14.6E %14.6E %14.6E %14.6E\n" %
                     (betas[i], E[i]-E0, lnrho[i], dlnrho[i],
                      S[i], dS[i], C[i], dC[i]))
i = len(betas)-1
sys.stdout.write("%12.6F %16.8E %14s %14s"
                 " %14.6E %14.6E %14s %14s\n" %
                 (betas[i], E[i]-E0, "--", "--",
                  S[i], dS[i], "--", "--"))

if len(energies) > 0:
    sys.stdout.write("\n\n")
    sys.stdout.write("# Interpolated level densities\n")
    sys.stdout.write("%16s %16s %16s %12s %12s\n" % ("Ex", "ln(rho)", "+/-",
                                                     "beta1", "beta2"))
    for i in range(0,len(energies)):
        sys.stdout.write("%16.8E %16.8E %16.8E %12.6F %12.6F\n" %
                         (energies[i], lnrhoi[i], dlnrhoi[i], betas[i1s[i]],
                          betas[i2s[i]]))



os.chdir("../Data")
savetxt("leveldensity.dat", stack([betas,array(E)-E0,lnrho,dlnrho,S,dS,C,dC],axis=1), delimiter="\t", header="Beta\tEx\tlog(rho),dlog(rho),S,dS,C,dC")
