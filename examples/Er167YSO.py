import numpy as np
import matplotlib.pyplot as plt
import dieke
import gzip
import pickle

# Calculates the energy levels of Er:YSO
nf = 11  # 11 f-electrons means we're dealing with Er
I = 7/2
# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"

# try:
#     fname = 'Er167matricies.dat.gz'
#     f = gzip.open(fname, 'rb')
#     print("Loading matricies")
#     Er = pickle.load(f)
# except FileNotFoundError:
#     print("Making matricies")
#     Er = dieke.RareEarthIon(nf,I)
#     pickle.dump(Er, gzip.open(fname, 'wb'))
Er = dieke.RareEarthIon(11,7/2)

# Read in a a set of crystal field parameters from Pr:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Er.numlevels()
numLSJmJ = Er.numstates()

# Variable multiplier for the crystal field strength
# this will be the x-axis of the graph we will draw
# cfstrvals = np.linspace(0, 1.5, 15)


cfparams['E0'] = 35503.5 #this is ignored
cfparams['ZETA'] = 2362.9
cfparams['F2'] = 96029.6
cfparams['F4'] = 67670.6
cfparams['F6'] = 53167.1

cfparams['B20'] = -149.8
cfparams['B21'] = 420.6+396.0j
cfparams['B22'] = -228.5+27.6j
cfparams['B40'] = 1131.2
cfparams['B41'] = 985.7+34.2j
cfparams['B42'] = 296.8+145.0j
cfparams['B43'] = -402.3-381.7j
cfparams['B44'] = -282.3+1114.3j
cfparams['B60'] = -263.2
cfparams['B61'] = 111.9+222.9j
cfparams['B62'] = 124.7+195.9j
cfparams['B63'] = -97.9+139.7j
cfparams['B64'] = -93.7-145.0j
cfparams['B65'] = 13.9+109.5j
cfparams['B66'] = 3.0-108.6j
cfparams['A'] = 0.005466 #this is ignored
cfparams['Q'] = 0.0716 #this is ignored


cfparams['HF'] = cfparams['A']

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

#energy levels from Sebastians paper 
seb_levels = [15, 47, 75, 130, 199, 388, 462, 508,
              6522, 6560, 6583, 6640, 6777, 6833, 6867,
              10206, 10236, 10267, 10339, 10381, 10398]
seb_levels = np.array(seb_levels)
E0seb = np.min(seb_levels)
seb_levels = seb_levels-E0seb


print('\n    Jevon  Sebastian Difference')
for k in range(len(seb_levels)):
    print("%9.1f %9.1f %6.1f"%(
        np.real(calc_nrg_levels[k]),
        seb_levels[k],
        np.real(calc_nrg_levels[k]) - seb_levels[k]))
