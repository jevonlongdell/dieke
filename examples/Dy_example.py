import numpy as np
import dieke
import pickle



# Calculates the energy levels of Pr:LaF3 as the strength of the
# crystal field terms are varied from zero to 150% or their actual
# values.


nf = 9  # 9 f-electrons means we're dealing with Dy

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"

try:
    f = open('Dymatricies.dat', 'rb')
    print("Loading matricies")
    Dy = pickle.load(f)
except FileNotFoundError:
    print("Making matricies")
    Dy = dieke.RareEarthIon(nf)
    pickle.dump(Dy,open('Dymatricies.dat', 'wb'))


# Read in a a set of crystal field parameters from Dy:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Dy.numlevels()
numLSJmJ = Dy.numstates()

# Variable multiplier for the crystal field strength
# this will be the x-axis of the graph we will draw
cfstrvals = np.linspace(0, 1.5, 15)


# Make a free-ion Hamiltonian
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Dy.FreeIonMatrix:
        print("Adding free ion parameter %s\n" % (k))
        H0 = H0+cfparams[k]*Dy.FreeIonMatrix[k]

# Add in the crystal field terms and diagonalise the result
H = H0
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Dy.Cmatrix(k, q)
            else:
                H = H+cfparams['B%d%d' % (k, q)]*Dy.Cmatrix(k, q)
                H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
                    * Dy.Cmatrix(k, -q)


# zeroify all cf params that are incompatible with YPO4s symmetry
H = H0
for k, q in [(2, 0), (4, 0), (6, 0), (4, 4), (6, 4)]:
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Dy.Cmatrix(k, q)
            else:
                H = H+cfparams['B%d%d' % (k, q)]*Dy.Cmatrix(k, q)
                H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
                    * Dy.Cmatrix(k, -q)

(evals, vec) = np.linalg.eig(H)

idx = np.argsort(evals)

evals = evals[idx]
evals = evals - evals[0]
vec = vec[:, idx]

for st in range(8):
    for i in range(numLSJmJ):
        if abs(vec[i, st]) > 1e-2:
            print("%5.5f %+5.5f  %s" % (vec[i, st].real,
                                        vec[i, st].imag,
                                        Dy.LSJmJstateLabels[i]))
    print('\n\n')


#evals = np.matrix(np.diag(evals))
#evects = np.matrix(evects)





# evals = evals[index]
# evects = evects[:,index]


# E0 = np.min(evals)
# calc_nrg_levels = np.sort(evals-E0)

