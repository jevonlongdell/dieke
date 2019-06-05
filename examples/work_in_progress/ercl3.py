import numpy as np
import matplotlib.pyplot as plt
import dieke


# Calculates the energy levels of Pr:LaF3 as the strength of the
# crystal field terms are varied from zero to 150% or their actual
# values.

nf = 11  # 11 f-electrons means we're dealing with Er
Er = dieke.RareEarthIon(nf)


# Read in a a set of crystal field parameters from Er:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)

#Add in crystal field parameters for ErCl3.6H2O
#Values are the Er-I-Fit-R set from
#M. Karbowiak, P. Gnutek and C. Rudowicz, Physica B, 405, (2010), 1927–1940
#doi:10.1016/j.physb.2010.01.086
cfparams['nf'] =         11.0
cfparams['F2'] =      96869.00
cfparams['F4'] =      67842.00
cfparams['F6'] =      55715.00
cfparams['ZETA'] =     2380.00
cfparams['ALPHA'] =      19.17
cfparams['BETA'] =     -611.00
cfparams['GAMMA'] =    1571.00
cfparams['T2'] =        518.00
cfparams['T3'] =         34.00
cfparams['T4'] =         76.00
cfparams['T6'] =       -340.00
cfparams['T7'] =        317.00
cfparams['T8'] =        393.00
cfparams['M0'] =          4.58
cfparams['P2'] =        740.00
cfparams['B20'] =       170
cfparams['B40'] =      -360
cfparams['B60'] =      -201
cfparams['B22'] =      -349
cfparams['B42'] =       451 + 32j
cfparams['B44'] =      -519 - 158j
cfparams['B62'] =       182 + 87j
cfparams['B64'] =      -370 + 162j
cfparams['B66'] =        86 + 233j
cfparams['M2'] =          2.56
cfparams['M4'] =          1.74
cfparams['P4'] =        555.00
cfparams['P6'] =        370.00

# # Number of levels and states
numLSJ = Er.numlevels()
numLSJmJ = Er.numstates()
# # Variable multiplier for the crystal field strength
# # this will be the x-axis of the graph we will draw
# cfstrvals = np.linspace(0, 1.5, 15)


# Make a free-ion Hamiltonian
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Er.FreeIonMatrix:
        print("Adding free ion parameter \'%s\' = %g" % (k, cfparams[k]))
        H0 = H0+cfparams[k]*Er.FreeIonMatrix[k]

# Add in the crystal field terms and diagonalise the result
H = H0
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Er.Cmatrix(k, q)
            else:
                Bkq = cfparams['B%d%d' % (k, q)]
                Bkmq = (-1)**q*np.conj(Bkq)
                Ckq = Er.Cmatrix(k, q)
                Ckmq = Er.Cmatrix(k, -q)
                #See page 44, eq 3.1 of the crystal field handbook
                H = H + Bkq*Ckq + Bkmq*Ckmq

(evals, evects) = np.linalg.eig(H)
E0 = np.min(evals)
calc_nrg_levels = np.sort(evals-E0)
calc_nrg_levels = calc_nrg_levels[::2]  # ignore every second element   

print(calc_nrg_levels)

# # Make a matrix for the energy levels vs crystal field strength
# # (this is what we will plot)
# numLSJmJ = len(LSJmJstates)
# nrglevels = np.zeros([len(cfstrvals), numLSJmJ])


# # Work out the energy levels as the crystal field strengths are varied.
# level_index = 0
# cfparams2 = cfparams
# for cfstr in cfstrvals:
#     H = H0
#     for k in [2, 4, 6]:
#         for q in range(0, k+1):
#             if 'B%d%d' % (k, q) in cfparams:
#                 # print("adding in B %d %d"%(k, q))
#                 if q == 0:
#                     H = H+cfstr*cfparams['B%d%d' % (k, q)] * Ckq[(k, q)]
#                 else:
#                     H = H+cfstr*cfparams['B%d%d' % (k, q)] * Ckq[(k, q)]
#                     H = H+cfstr*((-1)**q) \
#                         * np.conj(cfparams['B%d%d' % (k, q)]) * Ckq[(k, -q)]
#     (evals, evects) = np.linalg.eig(H)
#     nrglevels[level_index, :] = np.sort(np.real(evals))-E0
#     level_index = level_index+1


# # Plot the results
# plt.ion()
# plt.axis([0, max(cfstrvals), -2000, 25000])
# plt.plot(cfstrvals, nrglevels, 'x-', [1, 1], [-10000, +50000])
# plt.draw()


# # Print out the 20 lowest enery levels
# for k in range(20):
#     print(calc_nrg_levels[k])
