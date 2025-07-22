#!/usr/bin/env python3
# Encoding: UTF-8
#
# A script for computing q-projected level densities with uncertainties
# estimated by the block jackknife method. Designed to work together with
# the scripts qlvl_resample.py and qlvl_abcfit.jl.
#
# USAGE:
#  python qlvl_final.py qlvl_final.in
#
# See the example input file for details of input parameters.


import sys
import numpy as np
import qlevden
import resample


# Execute the input parameter file to "read" the parameters into the global
# namespace
if len(sys.argv) != 2:
    sys.stderr.write("Error: Must specify input file, no other parameters.\n")
    sys.exit(1)

exec(open(sys.argv[1]).read())

# Translate lists of energies and inverse temperatures into NumPy arrays
# (if they weren't already).
Evals = np.array(energy_values)
bevals = np.array(inv_temperature_values)


def deriv(y, x):
    '''
    Returns derivatives dy/dx approximated with the three-point formula.
    '''
    d = np.full(len(x), np.NaN)
    for i in range(1,len(x)-1):
      d[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    return d


def jackknife(resampled):
    '''
    Does resample.jackknife_results for a whole range of inverse temperatures.
    '''
    a_aerr = [resample.jackknife_results([resampled[iblock][ibeta]
        for iblock in range(len(resampled))]) for ibeta in range(len(resampled[0]))]
    return [x[0] for x in a_aerr], [x[1] for x in a_aerr]


# Initialize lists
rhos = []
lnbeta = np.log(bevals)
rhos_sph, rhos_pro, rhos_obl = [], [], []
as_orig, bs_orig, cs_orig, taus_orig = [], [], [], []
das_orig, dbs_orig, dcs_orig = [], [], []
d2as_orig, d2bs_orig, d2cs_orig = [], [], []
as_spline, bs_spline, cs_spline = [], [], []
das_spline, dbs_spline, dcs_spline = [], [], []
d2as_spline, d2bs_spline, d2cs_spline = [], [], []
as_splinetest, bs_splinetest, cs_splinetest = [], [], []
Es, Cs, lnZs = [], [], []
Ps = []
Ps_sph, Ps_pro, Ps_obl = [], [], []

if output_projected:
    Eqs = [[] for x in range(len(deformations))]
    Cqs = [[] for x in range(len(deformations))]
    lnrhoqs = [[] for x in range(len(deformations))]
    Fqs = [[] for x in range(len(deformations))]
    dFqs = [[] for x in range(len(deformations))]
    d2Fqs = [[] for x in range(len(deformations))]
    Pqs = [[] for x in range(len(deformations))]
    Sqs = [[] for x in range(len(deformations))]


def read_blocks(filename):
    '''
    Yields blocks of the input data, separated by empty lines. Lines starting
    with # are stripped away as comment lines.
    '''
    with open(filename) as f:
        block = []
        for line in f:
            ln = line.strip()
            if ln:
                # Append non-comment lines to the current block
                if not ln.startswith("#"):
                    block.append(ln)
            else:
                # If the current block is empty, this is still one of the empty
                # lines separating blocks
                if block:
                    yield block
                    block = []


def read_data_blocks(filename):
    '''
    Returns the values of a, b, c, tau, energy, and heat capacity for every
    inverse temperature beta in each block in the input file.
    '''
    for block in read_blocks(filename):
        data = list(zip(*[ln.split() for ln in block]))
        beta, a, b, c, E, C, tau = [np.array([float(x) for x in column]) for column in data]
        yield beta, a, b, c, E, C, tau


# The first pass: Determine the jackknife errors of a, b, and c obtained from
# the moments. These will be needed for properly weighting the points when
# doing the spline fit.
for beta, a, b, c, E, C, tau in read_data_blocks(inputfile):
    as_orig.append(a)
    bs_orig.append(b)
    cs_orig.append(c)
    taus_orig.append(tau)

# Compute the jackknifed values and errors
a, aerr = jackknife(as_orig)
b, berr = jackknife(bs_orig)
c, cerr = jackknife(cs_orig)
tau, tauerr = jackknife(taus_orig)

abc_spline_weights = [np.array(aerr), np.array(berr), np.array(cerr)]
spline_dof = len(beta) - nr_spline_knots - 1

# Determine the optimal knot placement
abc_spline_knots = qlevden.place_spline_knots(nr_spline_knots, beta)
print("Number of spline knot points: {}".format(nr_spline_knots))
print("Knot points placed at temperatures:")
print(" ".join(str(x) for x in abc_spline_knots))

# The second pass
# f = open("samples16.txt", "w")
for beta, a, b, c, E, C, tau in read_data_blocks(inputfile):

    # initialize the level density computation with these inputs
    pq_fit = qlevden.PqFitData(beta=beta, E=E, a=a, b=b, c=c, C=C)
    ldsolver = qlevden.QProjectedLevelDensitySolver(pq_fit, A, nr_states_n,
        nrn, nr_states_p, nrp, 3, abc_spline_knots, abc_spline_weights)

    # a, b, c and their temperature derivatives
    if output_abc:
        
        # Naive, inaccurate derivatives with the three-point rule
        da = -beta**2*deriv(a, beta)
        db = -beta**2*deriv(b, beta)
        dc = -beta**2*deriv(c, beta)
        d2a = -beta**2*deriv(da, beta)
        d2b = -beta**2*deriv(db, beta)
        d2c = -beta**2*deriv(dc, beta)

        das_orig.append(da)
        dbs_orig.append(db)
        dcs_orig.append(dc)
        d2as_orig.append(d2a)
        d2bs_orig.append(d2b)
        d2cs_orig.append(d2c)

        # splines for a, b, c, and derivatives
        as_spline.append( ldsolver.a_spline(bevals) )
        bs_spline.append( ldsolver.b_spline(bevals) )
        cs_spline.append( ldsolver.c_spline(bevals) )
        das_spline.append( -bevals**2*ldsolver.da_spline(bevals) )
        dbs_spline.append( -bevals**2*ldsolver.db_spline(bevals) )
        dcs_spline.append( -bevals**2*ldsolver.dc_spline(bevals) )
        d2as_spline.append( bevals**3 * (2*ldsolver.da_spline(bevals)
            + bevals*ldsolver.d2a_spline(bevals)) )
        d2bs_spline.append( bevals**3 * (2*ldsolver.db_spline(bevals)
            + bevals*ldsolver.d2b_spline(bevals)) )
        d2cs_spline.append( bevals**3 * (2*ldsolver.dc_spline(bevals)
            + bevals*ldsolver.d2c_spline(bevals)) )

        # splines evaluated at the SMMC points
        # as beta=0 gives divisions by zero, suppress that warning
        with np.errstate(divide = 'ignore'):
            as_splinetest.append( ldsolver.a_spline(beta) )
            bs_splinetest.append( ldsolver.b_spline(beta) )
            cs_splinetest.append( ldsolver.c_spline(beta) )

    # <E>, Cv, ln Z
    if output_unprojected:
        Es.append( E )
        Cs.append( C )
        lnZs.append( ldsolver.lnZ )

    # q-projected E, F, dF/dT, d^2F/dT^2, C, P, S, rho
    # when varying beta or gamma
    if output_projected:
        for i, q in enumerate(deformations):
            be, ga = q
            rs = [ldsolver.evaluate_at_temperature(be, ga, 1.0/x) for x in bevals]
            Eqs[i].append( [x.E for x in rs] )
            Pqs[i].append( [x.P for x in rs] )
            Cqs[i].append( [x.C for x in rs] )
            Sqs[i].append( [x.S for x in rs] )
            Fqs[i].append( [x.F for x in rs] )
            dFqs[i].append( [x.dFdT for x in rs] )
            d2Fqs[i].append( [x.d2FdT2 for x in rs] )
            with np.errstate(divide = 'ignore'):
                lnrhoqs[i].append( [np.log(x.rho) for x in rs] )
            # DEBUG: save resampled values for computing the correlation matrix
            # if i == 0: # selected deformation
            #     x = rs[16] # selected temp
            #     f.write(f"{x.T} {x.E} {x.P} {x.C} {x.S} {x.F} {x.dFdT} {x.d2FdT2} {x.rho} {x.P} {x.lnZq} {x.P} {x.lnZ} {x.dlnP} {x.d2lnP}\n")
    
    # normalization probability and probabilities of shapes
    if output_probability:
        results = [ldsolver.integrate_at_temperature(1.0/x, (0.0, BETA_MAX), (0.0, np.pi/3)) for x in bevals]
        Ps.append([x.P for x in results])
        results = [ldsolver.integrate_at_temperature(1.0/x, (0.0, BETA_LIMIT), (0.0, np.pi/3)) for x in bevals]
        Ps_sph.append([x.P for x in results])
        results = [ldsolver.integrate_at_temperature(1.0/x, (BETA_LIMIT, BETA_MAX), (0.0, np.pi/6)) for x in bevals]
        Ps_pro.append([x.P for x in results])
        results = [ldsolver.integrate_at_temperature(1.0/x, (BETA_LIMIT, BETA_MAX), (np.pi/6, np.pi/3)) for x in bevals]
        Ps_obl.append([x.P for x in results])
    
    # Integrated level densities over shapes
    if output_integrated:
        # the total level density
        results = [ldsolver.integrate_at_energy(x + E_gs, (0.0, BETA_MAX), (0.0, np.pi/3)) for x in Evals]
        rhos.append([x.rho for x in results])
    
        # spherical, prolate, oblate level densities
        results = [ldsolver.integrate_at_energy(x + E_gs, (0.0, BETA_LIMIT), (0.0, np.pi/3)) for x in Evals]
        rhos_sph.append([x.rho for x in results])
        results = [ldsolver.integrate_at_energy(x + E_gs, (BETA_LIMIT, BETA_MAX), (0.0, np.pi/6)) for x in Evals]
        rhos_pro.append([x.rho for x in results])
        results = [ldsolver.integrate_at_energy(x + E_gs, (BETA_LIMIT, BETA_MAX), (np.pi/6, np.pi/3)) for x in Evals]
        rhos_obl.append([x.rho for x in results])

# f.close()

# At this point, the computations are complete for resampled sets; for each
# type of observable computed, take their average and error and print them
# out to the file they belong to.

# a, b, c and derivatives
if output_abc:
    with open(output_abc, 'w') as f:
        a, aerr = jackknife(as_orig)
        b, berr = jackknife(bs_orig)
        c, cerr = jackknife(cs_orig)
        tau, tauerr = jackknife(taus_orig)
        da, daerr = jackknife(das_orig)
        db, dberr = jackknife(dbs_orig)
        dc, dcerr = jackknife(dcs_orig)
        d2a, d2aerr = jackknife(d2as_orig)
        d2b, d2berr = jackknife(d2bs_orig) 
        d2c, d2cerr = jackknife(d2cs_orig)

        f.write("# (0) original a, b, c, tau\n")
        f.write(("# {:^12}"+8*"  {:^12}"+"\n").format("beta", "a", "+/-", "b", "+/-", "c", "+/-", "tau", "+/-"))
        for ib, be in enumerate(beta):
            f.write(("{:12}"+8*"  {:12.4e}"+"\n").format(be, a[ib], aerr[ib],
                b[ib], berr[ib], c[ib], cerr[ib], tau[ib], tauerr[ib]))
        f.write("\n\n")

        f.write("# (1) original da/dT, db/dT, dc/dT\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "da/dT", "+/-", "db/dT", "+/-", "dc/dT", "+/-"))
        for ib, be in enumerate(beta):
            f.write(("{:12}"+6*"  {:12.4e}"+"\n").format(be, da[ib], daerr[ib], db[ib], dberr[ib], dc[ib], dcerr[ib]))
        f.write("\n\n")

        f.write("# (2) original d^2a/dT^2, d^2b/dT^2, d^2c/dT^2\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "d^2a/dT^2", "+/-", "d^2b/dT^2", "+/-", "d^2c/dT^2", "+/-"))
        for ib, be in enumerate(beta):
            f.write(("{:12}"+6*"  {:12.4e}"+"\n").format(be, d2a[ib], d2aerr[ib], d2b[ib], d2berr[ib], d2c[ib], d2cerr[ib]))
        f.write("\n\n")

        a, aerr = jackknife(as_spline)
        b, berr = jackknife(bs_spline)
        c, cerr = jackknife(cs_spline)
        da, daerr = jackknife(das_spline)
        db, dberr = jackknife(dbs_spline)
        dc, dcerr = jackknife(dcs_spline)
        d2a, d2aerr = jackknife(d2as_spline)
        d2b, d2berr = jackknife(d2bs_spline)
        d2c, d2cerr = jackknife(d2cs_spline)

        f.write("# (3) spline-fit a, b, c\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "a", "+/-", "b", "+/-", "c", "+/-"))
        for ib, be in enumerate(bevals):
            f.write(("{:12.4f}"+6*"  {:12.4e}"+"\n").format(be, a[ib], aerr[ib], b[ib], berr[ib], c[ib], cerr[ib]))
        f.write("\n\n")

        f.write("# (4) spline-fit da/dT, db/dT, dc/dT\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "da/dT", "+/-", "db/dT", "+/-", "dc/dT", "+/-"))
        for ib, be in enumerate(bevals):
            f.write(("{:12.4f}"+6*"  {:12.4e}"+"\n").format(be, da[ib], daerr[ib], db[ib], dberr[ib], dc[ib], dcerr[ib]))
        f.write("\n\n")

        f.write("# (5) spline-fit d^2a/dT^2, d^2b/dT^2, d^2c/dT^2\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "d^2a/dT^2", "+/-", "d^2b/dT^2", "+/-", "d^2c/dT^2", "+/-"))
        for ib, be in enumerate(bevals):
            f.write(("{:12.4f}"+6*"  {:12.4e}"+"\n").format(be, d2a[ib], d2aerr[ib], d2b[ib], d2berr[ib], d2c[ib], d2cerr[ib]))
        
        # Check how well the spline fit works
        a, aerr = jackknife(as_orig)
        b, berr = jackknife(bs_orig)
        c, cerr = jackknife(cs_orig)
        
        atest, atesterr = jackknife(as_splinetest)
        btest, btesterr = jackknife(bs_splinetest)
        ctest, ctesterr = jackknife(cs_splinetest)
        
        print("Number of temperature values: {}".format(len(beta)))
        print("\nchi^2 of the spline fits:")
        nrvals = len(beta)
        sigma = [(x1-x2)/np.sqrt(dx**2 + dx2**2 + 1e-10) for x1, x2, dx, dx2 in zip(a, atest, aerr, atesterr)]
        print("a:  {}".format(sum([x**2 for x in sigma[1:]])/spline_dof))
        sigma = [(x1-x2)/np.sqrt(dx**2 + dx2**2 + 1e-10) for x1, x2, dx, dx2 in zip(b, btest, berr, btesterr)]
        print("b:  {}".format(sum([x**2 for x in sigma[1:]])/spline_dof))
        sigma = [(x1-x2)/np.sqrt(dx**2 + dx2**2 + 1e-10) for x1, x2, dx, dx2 in zip(c, ctest, cerr, ctesterr)]
        print("c:  {}".format(sum([x**2 for x in sigma[1:]])/spline_dof))
        

# Total energy, heat capacity, and log of the partition function
if output_unprojected:
    with open(output_unprojected, 'w') as f:
        E, dE = jackknife(Es)
        C, dC = jackknife(Cs)
        lnZ, dlnZ = jackknife(lnZs)

        f.write("# non-projected <E>, Cv, ln Z\n")
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "<E>", "+/-", "Cv", "+/-", "ln Z", "+/-"))
        for ib, be in enumerate(beta):
            f.write(("{:12.4f}"+6*"  {:12.4e}"+"\n").format(be, E[ib], dE[ib], C[ib], dC[ib], lnZ[ib], dlnZ[ib]))


# q-projected quantities at given deformations (beta, gamma)
if output_projected:
    with open(output_projected, 'w') as f:
        for i, q in enumerate(deformations):
            f.write("# ({}) beta, gamma = {}, {}\n".format(i, *q))
            f.write(("# {:^12}"+16*"  {:^12}"+"\n").format("beta", "ln rho", "+/-", "Eq", "+/-", "Pq", "+/-",
                "Sq", "+/-", "Cq", "+/-", "Fq", "+/-", "dFq", "+/-", "d2Fq", "+/-"))
            E, dE = jackknife(Eqs[i])
            C, dC = jackknife(Cqs[i])
            S, dS = jackknife(Sqs[i])
            P, dP = jackknife(Pqs[i])
            F, dF = jackknife(Fqs[i])
            derF, dderF = jackknife(dFqs[i])
            der2F, dder2F = jackknife(d2Fqs[i])
            with np.errstate(invalid='ignore'): # ignore NaNs
                lnrho, dlnrho = jackknife(lnrhoqs[i])
            for ib, be in enumerate(bevals):
                f.write(("{:12.4f}"+16*"  {:12.4e}"+"\n").format(be, lnrho[ib], dlnrho[ib], 
                    E[ib], dE[ib], P[ib], dP[ib], S[ib], dS[ib], C[ib], dC[ib],
                    F[ib], dF[ib], derF[ib], dderF[ib], der2F[ib], dder2F[ib]))
            f.write("\n\n")


# Probabilities of the different shapes at given temperatures
if output_probability:
    with open(output_probability, 'w') as f:
        P, dP = jackknife(Ps)

        f.write("# (0) normalization (integrated Pq)\n")
        f.write(("# {:^12}"+2*"  {:^12}"+"\n").format("beta", "P", "+/-"))
        for ib, be in enumerate(bevals):
            f.write(("{:12.4f}"+2*"  {:12.4e}"+"\n").format(be, P[ib], dP[ib]))
        f.write("\n\n")

        Psph, dPsph = jackknife(Ps_sph)
        Ppro, dPpro = jackknife(Ps_pro)
        Pobl, dPobl = jackknife(Ps_obl)

        f.write("# (1) probabilities of shapes (integrated Pq) (beta0 = {})\n".format(BETA_LIMIT))
        f.write(("# {:^12}"+6*"  {:^12}"+"\n").format("beta", "P (sph.)", "+/-", "P (prol.)", "+/-", "P (obl.)", "+/-"))
        for ib, be in enumerate(bevals):
            f.write(("{:12.4f}"+6*"  {:12.4e}"+"\n").format(be, Psph[ib], dPsph[ib], Ppro[ib], dPpro[ib], Pobl[ib], dPobl[ib]))


# Integrated level density for the different shapes (and all deformations)
if output_integrated:
    with open(output_integrated, 'w') as f:
        f.write("# (0) integrated level density\n")
        f.write(("# {:^12}"+2*"  {:^12}"+"\n").format("E", "rho", "+/-"))
        for i, E in enumerate(Evals):
            lnrho_blocks = np.array([np.log(r[i]) for r in rhos])
            lnrho, dlnrho = resample.jackknife_results(lnrho_blocks)
            f.write("{:12.4f}  {:12.4e}  {:12.4e}\n".format(E, lnrho, dlnrho))
        f.write("\n\n")

        f.write("# (1) integrated level density over shapes\n")
        for i, E in enumerate(Evals):
            lnrho_blocks = np.array([np.log(r[i]) for r in rhos_sph])
            lnrho_sph, dlnrho_sph = resample.jackknife_results(lnrho_blocks)
            lnrho_blocks = np.array([np.log(r[i]) for r in rhos_pro])
            lnrho_pro, dlnrho_pro = resample.jackknife_results(lnrho_blocks)
            lnrho_blocks = np.array([np.log(r[i]) for r in rhos_obl])
            lnrho_obl, dlnrho_obl = resample.jackknife_results(lnrho_blocks)
            f.write("{:12.4f}  {:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}\n".format(
                E, lnrho_sph, dlnrho_sph, lnrho_pro, dlnrho_pro, lnrho_obl, dlnrho_obl))
